// plot_hcal_clus.C
// Usage:
//   root -l -q 'plot_hcal_clus.C("QE_data_GEN3_sbs100p_nucleon_np_model2_sbstrackingon.root")'
//   root -l -q 'plot_hcal_clus.C("file.root","Tout","hcal_plots.pdf")'

#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TTreeReaderValue.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TString.h"
#include "TSystem.h"
#include <iostream>

static void Style1D(TH1 *h){
  h->SetLineWidth(2);
  h->SetStats(1);
}

void plot_hcal_clus(const char* infile,
                    const char* treename="Tout",
                    const char* outpdf="hcal_clus_plots.pdf")
{
  gStyle->SetOptStat(1110);

  // ---- Input (supports wildcards and multiple files if you want later)
  TChain chain(treename);
  int nadded = chain.Add(infile);
  if(nadded <= 0){
    std::cerr << "[plot_hcal_clus] ERROR: no files matched: " << infile << "\n";
    return;
  }
  std::cout << "[plot_hcal_clus] Added " << nadded << " file(s)\n";
  std::cout << "[plot_hcal_clus] Entries: " << chain.GetEntries() << "\n";

  // ---- Reader (safe for variable-length arrays like hcal_clus_*[Ndata_clus_id])
  TTreeReader R(&chain);

  // Cluster arrays
  TTreeReaderArray<Double_t> hcal_clus_e      (R, "hcal_clus_e");
  TTreeReaderArray<Double_t> hcal_clus_x      (R, "hcal_clus_x");
  TTreeReaderArray<Double_t> hcal_clus_y      (R, "hcal_clus_y");
  TTreeReaderArray<Double_t> hcal_clus_atime  (R, "hcal_clus_atime");
  TTreeReaderArray<Double_t> hcal_clus_tdctime(R, "hcal_clus_tdctime");
  TTreeReaderArray<Double_t> hcal_clus_nblk   (R, "hcal_clus_nblk");

  // Scalar
  TTreeReaderValue<Double_t> hcal_time(R, "hcal_time");
  TTreeReaderValue<Double_t> bbcal_time(R,"bbcal_time");
  TTreeReaderValue<Double_t> coin_time(R,"coin_time");
  TTreeReaderValue<Double_t> W2(R,"W2");
  TTreeReaderValue<Double_t> dx(R,"dx");
  TTreeReaderValue<Double_t> dy(R,"dy");
  TTreeReaderValue<Double_t> vz(R,"vz");
  TTreeReaderValue<Double_t> ePS(R,"ePS");
  TTreeReaderValue<Double_t> eSH(R,"eSH");
  TTreeReaderValue<Double_t> trP(R,"trP");


  // ---- Histograms (set broad ranges now; you can tighten after a first look)
  TH1D *hE   = new TH1D("hE",   "HCAL cluster E;hcal_clus_e;Counts", 400, 0.0, 1.0);
  TH1D *hX   = new TH1D("hX",   "HCAL cluster x;hcal_clus_x;Counts", 400, -4.0, 4.0);
  TH1D *hY   = new TH1D("hY",   "HCAL cluster y;hcal_clus_y;Counts", 400, -4.0, 4.0);
  TH1D *hAt  = new TH1D("hAt",  "HCAL cluster ADC time;hcal_clus_atime;Counts", 600, -300.0, 300.0);
  TH1D *hTdc = new TH1D("hTdc", "HCAL cluster TDC time;hcal_clus_tdctime;Counts", 600, -300.0, 300.0);
  TH1D *hNb  = new TH1D("hNb",  "HCAL cluster nblk;hcal_clus_nblk;Counts", 50, 0.0, 50.0);

  TH1D *hT   = new TH1D("hT",   "HCAL time (scalar);hcal_time;Counts", 600, -300.0, 300.0);

  TH2D *hXY  = new TH2D("hXY",  "HCAL cluster x vs y;hcal_clus_x;hcal_clus_y",
                        240, -4.0, 4.0, 240, -4.0, 4.0);

  TH2D *hEHCal_vs_coin_time = new TH2D("hEHCal_vs_coin_time","HCal cluster enegry vs (hcal_clus_atime-bbcal_time)",200,170,210,300,0,1.5);

  Style1D(hE); Style1D(hX); Style1D(hY); Style1D(hAt); Style1D(hTdc); Style1D(hNb); Style1D(hT);

  // ---- Event loop
  Long64_t iev = 0;
  while(R.Next()){
    // Placeholder for event-level cuts you add later:
    // bool passEvent = (some condition);
    bool passEvent = abs(*vz)<0.27&&*ePS>0.2&&abs((*ePS+*eSH)/(*trP)-1)<0.2&&abs(*coin_time-185)<5&&*W2<1.6;
    if(!passEvent) continue;

    // fill scalar
    hT->Fill(*hcal_time);

    // number of clusters in this event == size of any of the arrays
    const auto nclus = hcal_clus_e.GetSize();

    for(Long64_t ic = 1; ic < (Long64_t)nclus; ic++){
      // Placeholder for per-cluster cuts you add later:
      // if(hcal_clus_e[ic] < 0.02) continue;

      hE->Fill(hcal_clus_e[ic]);
      hX->Fill(hcal_clus_x[ic]);
      hY->Fill(hcal_clus_y[ic]);
      hAt->Fill(hcal_clus_atime[ic]);
      hTdc->Fill(hcal_clus_tdctime[ic]);
      hNb->Fill(hcal_clus_nblk[ic]);
      hXY->Fill(hcal_clus_x[ic], hcal_clus_y[ic]);
      hEHCal_vs_coin_time->Fill((hcal_clus_atime[ic]-*bbcal_time),hcal_clus_e[ic]);
    }

    iev++;
    if(iev % 2000000 == 0) std::cout << "  processed " << iev << " events\n";
  }

  // ---- Draw to a multipage PDF
  TCanvas *c = new TCanvas("c","c",1000,800);

  TString pdf(outpdf);
  c->Print(pdf + "[");

  c->Clear(); hE->Draw("hist"); c->Print(pdf);
  c->Clear(); hX->Draw("hist"); c->Print(pdf);
  c->Clear(); hY->Draw("hist"); c->Print(pdf);
  c->Clear(); hNb->Draw("hist"); c->Print(pdf);

  c->Clear(); hAt->Draw("hist"); c->Print(pdf);
  c->Clear(); hTdc->Draw("hist"); c->Print(pdf);
  c->Clear(); hT->Draw("hist"); c->Print(pdf);
  c->Clear(); hEHCal_vs_coin_time->Draw("COLZ"); c->Print(pdf);
  //hEHCal_vs_coin_time->SetStats(0);
  hEHCal_vs_coin_time->GetXaxis()->SetTitle("(hcal_clus_atime - bbcal_time) [ns]");
  hEHCal_vs_coin_time->GetYaxis()->SetTitle("HCAL cluster energy [GeV]");

  c->Clear();
  gStyle->SetOptStat(0);
  hXY->Draw("colz");
  c->Print(pdf);
  gStyle->SetOptStat(1110);

  c->Print(pdf + "]");

  std::cout << "[plot_hcal_clus] Wrote: " << outpdf << "\n";
}

