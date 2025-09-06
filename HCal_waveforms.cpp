// hcal_allblock_waves_pdf_leaf.C
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TLatex.h"

static void print_progress(Long64_t done, Long64_t total){
  double f = (total>0) ? double(done)/double(total) : 1.0;
  int width = 40, filled = int(f*width);
  std::cout << "\r[";
  for(int i=0;i<width;i++) std::cout << (i<filled? '=' : (i==filled? '>' : ' '));
  std::cout << "] " << int(f*100.0) << "%  (" << done << "/" << total << ")" << std::flush;
  if(done == total) std::cout << std::endl;
}

// Draw all blocks for one event; tile across pages if needed.
// Uses TLeaf to read dynamic leaf-list Double_t arrays.
static int draw_event_pages_all_blocks_leaf(
    TTree* T,
    TLeaf* l_row, TLeaf* l_col, TLeaf* l_ns, TLeaf* l_sidx, TLeaf* l_samps,
    Long64_t ievt, TCanvas* c, const char* outpdf,
    int max_per_page, bool fixY, double ymin, double ymax, bool sort_by_rc)
{
  T->GetEntry(ievt);

  // Per-block metadata (length = #blocks in this event)
  int nb_blocks = std::min({ l_row->GetNdata(), l_col->GetNdata(),
                             l_ns->GetNdata(),  l_sidx->GetNdata() });
  if(nb_blocks <= 0){
    c->Clear(); c->cd();
    TLatex t; t.SetNDC(true);
    t.DrawLatex(0.12,0.5,Form("Entry %lld: no HCAL blocks", ievt));
    c->Print(outpdf);
    return 1;
  }

  // Pointers to the internal buffers (Double_t*)
  auto rowv  = (const Double_t*) l_row->GetValuePointer();
  auto colv  = (const Double_t*) l_col->GetValuePointer();
  auto nsv   = (const Double_t*) l_ns->GetValuePointer();     // stored as double; cast to int on use
  auto sidxv = (const Double_t*) l_sidx->GetValuePointer();   // stored as double; cast to int on use

  // Sample array for the whole event (flattened across all blocks)
  auto samps = (const Double_t*) l_samps->GetValuePointer();
  int  nsamp_total = l_samps->GetNdata();

  // Optional order: sort blocks by (row,col)
  std::vector<int> idx(nb_blocks);
  for(int i=0;i<nb_blocks;i++) idx[i]=i;
  if(sort_by_rc){
    std::sort(idx.begin(), idx.end(), [&](int a, int b){
      if(int(rowv[a]) != int(rowv[b])) return int(rowv[a]) < int(rowv[b]);
      return int(colv[a]) < int(colv[b]);
    });
  }

  int pages = 0;
  for(int pos=0; pos<nb_blocks; pos += max_per_page){
    int end   = std::min(pos + max_per_page, nb_blocks);
    int nplot = end - pos;

    // near-square tiling
    int nx = int(std::ceil(std::sqrt(double(nplot))));
    int ny = (nplot + nx - 1)/nx;

    c->Clear();
    c->Divide(nx, ny, 0.002, 0.002);

    c->cd();
    TLatex header; header.SetNDC(true); header.SetTextSize(0.03);
    header.DrawLatex(0.02, 0.98,
      Form("HCAL waveforms  |  entry=%lld  |  blocks on this page=%d  |  total blocks in event=%d",
           ievt, nplot, nb_blocks));

    for(int p=0; p<nplot; ++p){
      int i = idx[pos + p];

      int r   = int(rowv[i]);
      int ccol= int(colv[i]);
      int n   = int(nsv[i]);
      int off = int(sidxv[i]);

      c->cd(p+1);
      if(n<=0 || off<0 || (off+n) > nsamp_total){
        TLatex w; w.SetNDC(true); w.SetTextSize(0.05);
        w.DrawLatex(0.1,0.5,Form("r=%d c=%d  (invalid waveform)", r, ccol));
        continue;
      }

      TGraph *g = new TGraph(n);
      double lo=1e300, hi=-1e300;
      for(int k=0;k<n;k++){
        double y = samps[off+k];
        g->SetPoint(k, k, y);
        if(y<lo) lo=y; if(y>hi) hi=y;
      }
      g->SetLineWidth(1);
      g->SetTitle(Form("r=%d c=%d  off=%d  n=%d;sample;ADC", r, ccol, off, n));
      g->Draw("AL");

      if(fixY){
        g->GetYaxis()->SetRangeUser(ymin, ymax);
      }else{
        double pad = 0.06*(hi - lo + (hi==lo));
        g->GetYaxis()->SetRangeUser(lo - pad, hi + pad);
      }
    }

    c->Update();
    c->Print(outpdf);
    pages++;
  }

  return pages;
}

// Entry point
// root -l -q 'hcal_allblock_waves_pdf_leaf.C("file.root",0,-1,"Tout","HCAL_all_blocks.pdf")'
void hcal_allblock_waves_pdf_leaf(const char* rootfile,
                                  Long64_t first_entry = 0,
                                  Long64_t nentries    = -1,          // -1 => to end
                                  const char* treename = "Tout",
                                  const char* outpdf   = "hcal_all_blocks.pdf",
                                  int max_per_page     = 25,          // waveforms per page
                                  bool fixY            = false,       // true => fix y-range
                                  double ymin          = 0,
                                  double ymax          = 0,
                                  bool sort_by_rc      = true)        // sort plots by (row,col)
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(rootfile, "READ");
  if(!f || f->IsZombie()){ std::cerr << "ERROR: cannot open " << rootfile << "\n"; return; }
  TTree* T = dynamic_cast<TTree*>(f->Get(treename));
  if(!T){ std::cerr << "ERROR: tree \"" << treename << "\" not found.\n"; return; }

  // Grab leaves (leaf-list Double_t)
  TLeaf *l_row   = T->GetLeaf("sbs.hcal.adcrow");
  TLeaf *l_col   = T->GetLeaf("sbs.hcal.adccol");
  TLeaf *l_ns    = T->GetLeaf("sbs.hcal.nsamps");
  TLeaf *l_sidx  = T->GetLeaf("sbs.hcal.samps_idx");
  TLeaf *l_samps = T->GetLeaf("sbs.hcal.samps");
  if(!l_row || !l_col || !l_ns || !l_sidx || !l_samps){
    std::cerr << "ERROR: missing one or more leaves:"
              << " adcrow="<<(bool)l_row<<" adccol="<<(bool)l_col
              << " nsamps="<<(bool)l_ns<<" samps_idx="<<(bool)l_sidx
              << " samps="<<(bool)l_samps << "\n";
    return;
  }

  Long64_t N = T->GetEntries();
  if(first_entry < 0) first_entry = 0;
  Long64_t last = (nentries < 0) ? N : std::min(first_entry + nentries, N);
  if(first_entry >= last){ std::cerr << "Nothing to do (empty range).\n"; return; }

  TCanvas* c = new TCanvas("c_hcal_all","HCAL all block waveforms", 1200, 800);
  c->Print(TString::Format("%s[", outpdf)); // open multi-page

  Long64_t total = last - first_entry;
  Long64_t pages_written = 0;

  for(Long64_t ie = first_entry; ie < last; ++ie){
    pages_written += draw_event_pages_all_blocks_leaf(
        T, l_row, l_col, l_ns, l_sidx, l_samps,
        ie, c, outpdf, max_per_page, fixY, ymin, ymax, sort_by_rc);

    Long64_t done = (ie - first_entry + 1);
    if(done == total || (done % std::max<Long64_t>(1, total/100)) == 0)
      print_progress(done, total);
  }

  c->Print(TString::Format("%s]", outpdf)); // close PDF
  std::cout << "Wrote " << pages_written << " page(s) to " << outpdf
            << " covering entries [" << first_entry << ", " << (last-1) << "].\n";

  delete c;
  f->Close(); delete f;
}

