// hcal_primarycluster_allwaves_pdf_leaf.C
//
// Goal:
//  - Plot ALL waveforms belonging to the PRIMARY cluster using clus_blk.* arrays.
//  - NEW: Use ONLY sbs.hcal.clus_blk.id to locate waveform blocks.
//         Treat clus_blk.id as 1-based index into event-level block arrays.
//         => evt_index = clus_blk.id[k] - 1
//  - Use clus_blk.e for per-block energy annotation.
//  - Extract waveforms via event-level arrays:
//        sbs.hcal.adcrow, sbs.hcal.adccol,
//        sbs.hcal.nsamps, sbs.hcal.samps_idx, sbs.hcal.samps
//  - Dedicated header pad to prevent overlap.
//  - Keep cout statements for debugging.
//  - Optional global cut via TTreeFormula.
//
// Usage (no cut):
//   root -l -q 'hcal_primarycluster_allwaves_pdf_leaf.C("file.root",0,-1,"Tout","HCAL_primary_cluster_blocks.pdf")'
//
// Usage (with cut):
//   root -l -q 'hcal_primarycluster_allwaves_pdf_leaf.C("file.root",0,-1,"Tout","HCAL_primary_cluster_blocks.pdf",16,false,0,0,true,"sbs.hcal.clus.e[0]>0.05")'
//

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TTreeFormula.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TString.h"
#include "TLatex.h"

using std::cout;
using std::endl;

// Progress bar
static void print_progress(Long64_t done, Long64_t total){
  double f = (total>0) ? double(done)/double(total) : 1.0;
  int width = 40, filled = int(f*width);
  std::cout << "\r[";
  for(int i=0;i<width;i++) std::cout << (i<filled? '=' : (i==filled? '>' : ' '));
  std::cout << "] " << int(f*100.0) << "%  (" << done << "/" << total << ")" << std::flush;
  if(done == total) std::cout << std::endl;
}

// Generic leaf fetch with fallbacks
static TLeaf* getLeafAny(TTree* T, const std::vector<std::string>& names){
  for(const auto& n : names){
    if(!T) continue;
    TLeaf* l = T->GetLeaf(n.c_str());
    if(l) return l;
  }
  return nullptr;
}

// Draw primary-cluster block waveforms for one event using clus_blk.id -> evt_index.
static int draw_event_pages_primarycluster_blocks_leaf(
    TTree* T,
    // event-level per-block waveform leaves:
    TLeaf* l_row, TLeaf* l_col, TLeaf* l_ns, TLeaf* l_sidx, TLeaf* l_samps,
    // primary cluster block-list leaves:
    TLeaf* l_cblk_e, TLeaf* l_cblk_id,
    // optional cluster-level primary energy (nice for header):
    TLeaf* l_clusE,
    // optional global cut formula:
    TTreeFormula* fcut,
    Long64_t ievt, TCanvas* canv, const char* outpdf,
    int max_per_page, bool fixY, double ymin, double ymax, bool sort_by_rc)
{
  // Load event
  T->GetEntry(ievt);

  // Apply global cut if requested
  if(fcut){
    double pass = fcut->EvalInstance();
    if(pass == 0) return 0; // no pages written
  }

  // Event-level block counts
  int nb_blocks = std::min({ l_row->GetNdata(), l_col->GetNdata(),
                             l_ns->GetNdata(),  l_sidx->GetNdata() });

  if(nb_blocks <= 0){
    canv->Clear(); canv->cd();
    TLatex t; t.SetNDC(true);
    t.DrawLatex(0.12,0.5,Form("Entry %lld: no HCAL blocks", ievt));
    canv->Print(outpdf);
    return 1;
  }

  // Event-level buffers
  auto rowv  = (const Double_t*) l_row->GetValuePointer();
  auto colv  = (const Double_t*) l_col->GetValuePointer();
  auto nsv   = (const Double_t*) l_ns->GetValuePointer();
  auto sidxv = (const Double_t*) l_sidx->GetValuePointer();
  auto samps = (const Double_t*) l_samps->GetValuePointer();
  int  nsamp_total = l_samps->GetNdata();

  // Primary-cluster block-list buffers (ID + energy)
  const Double_t* cblkEv  = l_cblk_e  ? (const Double_t*) l_cblk_e->GetValuePointer()  : nullptr;
  const Double_t* cblkIDv = l_cblk_id ? (const Double_t*) l_cblk_id->GetValuePointer() : nullptr;

  int ncblk_e  = l_cblk_e  ? l_cblk_e->GetNdata()  : 0;
  int ncblk_id = l_cblk_id ? l_cblk_id->GetNdata() : 0;

  // Cluster primary energy (if available)
  double Eclus0 = std::numeric_limits<double>::quiet_NaN();
  if(l_clusE){
    const Double_t* clusEv = (const Double_t*) l_clusE->GetValuePointer();
    int nclusE = l_clusE->GetNdata();
    if(clusEv && nclusE > 0) Eclus0 = clusEv[0];
  }

  // Debug prints
  cout << "==============================" << endl;
  cout << "Entry = " << ievt << endl;
  if(fcut) cout << "Global cut PASS" << endl;
  cout << "nb_blocks(event) = " << nb_blocks << endl;
  cout << "clus_blk id N = " << ncblk_id << endl;
  cout << "clus_blk e  N = " << ncblk_e  << endl;
  cout << "Eclus[0] (if present) = " << Eclus0 << endl;

  if(!cblkIDv || ncblk_id <= 0){
    canv->Clear(); canv->cd();
    TLatex t; t.SetNDC(true); t.SetTextSize(0.04);
    t.DrawLatex(0.12,0.55,Form("Entry %lld: sbs.hcal.clus_blk.id not found or empty.", ievt));
    canv->Print(outpdf);
    return 1;
  }

  // Number of primary-cluster blocks we will attempt
  int ncblk = ncblk_id;

  // Build list of event-level block indices that belong to primary cluster
  struct BlockPick {
    int evt_index;   // index in event-level arrays
    int row;
    int col;
    double eblk;     // energy from clus_blk.e (matched by same k)
    int clusid;      // id from clus_blk.id (1-based)
  };

  std::vector<BlockPick> picks;
  picks.reserve(ncblk);

  for(int k=0; k<ncblk; ++k){
    int id1 = int(cblkIDv[k]);      // assumed 1-based
    int i   = id1 - 1;             // convert to 0-based event index

    double e = std::numeric_limits<double>::quiet_NaN();
    if(cblkEv && k < ncblk_e) e = cblkEv[k];

    cout << "clus_blk k=" << k
         << "  id(1-based)=" << id1
         << "  -> evt_index=" << i
         << "  Eblk=" << e
         << endl;

    if(i < 0 /*|| i >= nb_blocks*/){
      cout << "   WARNING: evt_index out of range for this event. Skipping." << endl;
      continue;
    }

    int r = int(rowv[i]);
    int c = int(colv[i]);

    picks.push_back({i, r, c, e, id1});
  }

  if(picks.empty()){
    canv->Clear(); canv->cd();
    TLatex t; t.SetNDC(true); t.SetTextSize(0.04);
    t.DrawLatex(0.12,0.55,Form("Entry %lld: no valid clus_blk.id mapped to event blocks.", ievt));
    canv->Print(outpdf);
    return 1;
  }

  // Optional sort by (row,col) derived from event-level arrays
  if(sort_by_rc){
    std::sort(picks.begin(), picks.end(), [](const BlockPick& a, const BlockPick& b){
      if(a.row != b.row) return a.row < b.row;
      return a.col < b.col;
    });
  }

  int nb_sel = (int)picks.size();
  int pages = 0;

  for(int pos=0; pos<nb_sel; pos += max_per_page){
    int end   = std::min(pos + max_per_page, nb_sel);
    int nplot = end - pos;

    int nx = int(std::ceil(std::sqrt(double(nplot))));
    int ny = (nplot + nx - 1)/nx;

    // ---------------- Layout with header pad (fix overlap) ----------------
    canv->Clear();

    double headerFrac = 0.08; // 8% reserved for header text

    TPad *phead = new TPad("phead","phead", 0.0, 1.0-headerFrac, 1.0, 1.0);
    TPad *pgrid = new TPad("pgrid","pgrid", 0.0, 0.0,           1.0, 1.0-headerFrac);

    phead->SetFillStyle(0);
    pgrid->SetFillStyle(0);

    phead->SetTopMargin(0.10);
    phead->SetBottomMargin(0.10);
    phead->SetLeftMargin(0.02);
    phead->SetRightMargin(0.02);

    pgrid->SetTopMargin(0.02);
    pgrid->SetBottomMargin(0.08);
    pgrid->SetLeftMargin(0.06);
    pgrid->SetRightMargin(0.02);

    phead->Draw();
    pgrid->Draw();

    // Header text
    phead->cd();
    TLatex header;
    header.SetNDC(true);
    header.SetTextSize(0.45); // relative to header pad

    TString htxt;
    if(std::isfinite(Eclus0)){
      htxt = Form("HCAL primary cluster block waveforms (ID-mapped) | entry=%lld | E_{clus}[0]=%.3f GeV | blocks in primary=%d | blocks this page=%d",
                  ievt, Eclus0, nb_sel, nplot);
    } else {
      htxt = Form("HCAL primary cluster block waveforms (ID-mapped) | entry=%lld | E_{clus}[0]=N/A | blocks in primary=%d | blocks this page=%d",
                  ievt, nb_sel, nplot);
    }
    header.DrawLatex(0.01, 0.5, htxt.Data());

    // Divide grid area
    pgrid->cd();
    pgrid->Divide(nx, ny, 0.002, 0.002);

    // Draw each block waveform
    for(int p=0; p<nplot; ++p){
      const auto& bp = picks[pos + p];
      int i = bp.evt_index;

      int r   = int(rowv[p]);
      int ccol= int(colv[p]);
      int n   = int(nsv[i]);
      int off = int(sidxv[i]);

      cout << "Plot p=" << p
           << "  evt_index=" << i
           //<< "  row=" << r
           //<< "  col=" << ccol
           << "  nsamps=" << n
           << "  sidx=" << off
           << "  Eblk(clus_blk)=" << bp.eblk
           << "  id(clus_blk,1-based)=" << bp.clusid
           << endl;

      pgrid->cd(p+1);

      if(n<=0 || off<0 || (off+n) > nsamp_total){
        TLatex w; w.SetNDC(true); w.SetTextSize(0.08);
        w.DrawLatex(0.1,0.5,Form("id=%d (invalid waveform)", bp.clusid));
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

      // Short title for small pads
      if(std::isfinite(bp.eblk)){
        g->SetTitle(Form("id=%d E_{blk}=%.3f GeV;sample;ADC", bp.clusid, bp.eblk));
      } else {
        g->SetTitle(Form("id=%d;sample;ADC", bp.clusid));
      }

      g->Draw("AL");

      if(fixY){
        g->GetYaxis()->SetRangeUser(ymin, ymax);
      } else {
        double pad = 0.06*(hi - lo + (hi==lo));
        g->GetYaxis()->SetRangeUser(lo - pad, hi + pad);
      }
    }

    canv->Update();
    canv->Print(outpdf);
    pages++;
  }

  return pages;
}

// Entry point
void hcal_primarycluster_allwaves_pdf_leaf(const char* rootfile,
                                           Long64_t first_entry = 0,
                                           Long64_t nentries    = -1,
                                           const char* treename = "Tout",
                                           const char* outpdf   = "HCAL_primary_cluster_blocks.pdf",
                                           int max_per_page     = 16,
                                           bool fixY            = false,
                                           double ymin          = 0,
                                           double ymax          = 0,
                                           bool sort_by_rc      = true,
                                           const char* cutexpr  = "")
{
  gROOT->SetBatch(kTRUE);
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(rootfile, "READ");
  if(!f || f->IsZombie()){
    std::cerr << "ERROR: cannot open " << rootfile << "\n";
    return;
  }

  TTree* T = dynamic_cast<TTree*>(f->Get(treename));
  if(!T){
    std::cerr << "ERROR: tree \"" << treename << "\" not found.\n";
    return;
  }

  // -------- Required event-level waveform leaves ----------
  TLeaf *l_row   = T->GetLeaf("sbs.hcal.adcrow");
  TLeaf *l_col   = T->GetLeaf("sbs.hcal.adccol");
  TLeaf *l_ns    = T->GetLeaf("sbs.hcal.nsamps");
  TLeaf *l_sidx  = T->GetLeaf("sbs.hcal.samps_idx");
  TLeaf *l_samps = T->GetLeaf("sbs.hcal.samps");

  if(!l_row || !l_col || !l_ns || !l_sidx || !l_samps){
    std::cerr << "ERROR: missing one or more required leaves:"
              << " adcrow="<<(bool)l_row<<" adccol="<<(bool)l_col
              << " nsamps="<<(bool)l_ns<<" samps_idx="<<(bool)l_sidx
              << " samps="<<(bool)l_samps << "\n";
    return;
  }

  // -------- Primary cluster block-list leaves ----------
  // ID is REQUIRED in this version.
  TLeaf *l_cblk_id = getLeafAny(T, {
    "sbs.hcal.clus_blk.id",
    "sbs.hcal.clus_blk_id"
  });

  TLeaf *l_cblk_e = getLeafAny(T, {
    "sbs.hcal.clus_blk.e",
    "sbs.hcal.clus_blk_e"
  });

  // Optional cluster energy for header
  TLeaf *l_clusE = getLeafAny(T, {
    "sbs.hcal.clus.e",
    "sbs.hcal.clus_e",
    "sbs.hcal.eclus"
  });

  cout << "==== Leaf check ====" << endl;
  cout << "l_cblk_id  = " << l_cblk_id  << endl;
  cout << "l_cblk_e   = " << l_cblk_e   << endl;
  cout << "l_clusE    = " << l_clusE    << endl;

  if(!l_cblk_id){
    std::cerr << "ERROR: This version requires sbs.hcal.clus_blk.id.\n";
    return;
  }

  // -------- Global cut formula ----------
  TTreeFormula* fcut = nullptr;
  if(cutexpr && std::strlen(cutexpr) > 0){
    cout << "Global cut expression: " << cutexpr << endl;
    fcut = new TTreeFormula("globalcut", cutexpr, T);
  }

  Long64_t N = T->GetEntries();
  if(first_entry < 0) first_entry = 0;
  Long64_t last = (nentries < 0) ? N : std::min(first_entry + nentries, N);
  if(first_entry >= last){
    std::cerr << "Nothing to do (empty range).\n";
    if(fcut) delete fcut;
    return;
  }

  TCanvas* c = new TCanvas("c_hcal_primary_all",
                           "HCAL primary cluster block waveforms (ID-mapped)",
                           1200, 800);

  c->Print(TString::Format("%s[", outpdf)); // open multi-page

  Long64_t total = last - first_entry;
  Long64_t pages_written = 0;

  for(Long64_t ie = first_entry; ie < last; ++ie){
    pages_written += draw_event_pages_primarycluster_blocks_leaf(
        T,
        l_row, l_col, l_ns, l_sidx, l_samps,
        l_cblk_e, l_cblk_id,
        l_clusE,
        fcut,
        ie, c, outpdf,
        max_per_page, fixY, ymin, ymax, sort_by_rc);

    Long64_t done = (ie - first_entry + 1);
    if(done == total || (done % std::max<Long64_t>(1, total/100)) == 0)
      print_progress(done, total);
  }

  c->Print(TString::Format("%s]", outpdf)); // close PDF

  std::cout << "Wrote " << pages_written << " page(s) to " << outpdf
            << " covering entries [" << first_entry << ", " << (last-1) << "].\n";

  delete c;
  f->Close(); delete f;
  if(fcut) delete fcut;
}
