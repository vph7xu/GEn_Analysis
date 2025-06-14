// -----------------------------------------------------------------------------
// ClusterDistance.C
// -----------------------------------------------------------------------------
// Compile with  .L ClusterDistance.C+   (the “+” makes ACLiC compile it)
//
// Call        ComputeClusterDistances("myfile.root","Tout");
//
// Produces a histogram "hDist" with the distances (cm, mm …
// whatever units x,y are in your tree) of every secondary cluster’s
// centroid to the primary (cid==0) centroid, event-by-event.
// -----------------------------------------------------------------------------

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>

#include <unordered_map>
#include <cmath>        // std::hypot
#include <iostream>

#include "parse.h"

void ComputeClusterDistances(const char* fname          = "input.root",
                             const char* kin            = "GEN4b_He3",
                             const char* treename       = "Tout",
                             int         maxHitsPerEvt  = 288)
{
  //---------------------------------------------------------------------------
  // 1. Open file and get tree
  //---------------------------------------------------------------------------
  TFile *f   = TFile::Open(fname,"READ");
  if(!f || f->IsZombie()){
    std::cerr << "Error opening file " << fname << '\n';
    return;
  }

  TTree *T = nullptr;
  f->GetObject(treename, T);
  if(!T){
    std::cerr << "Tree " << treename << " not found in " << fname << '\n';
    return;
  }


  std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin));

  double coin_time_L = getDoubleValue(config,"coin_time_L");
  double coin_time_H = getDoubleValue(config,"coin_time_H");

  double W2_L = getDoubleValue(config,"W2_L");
  double W2_H = getDoubleValue(config,"W2_H");

  double dy_L = getDoubleValue(config,"dy_L");
  double dy_H = getDoubleValue(config,"dy_H");

  double dx_L = getDoubleValue(config,"dx_L");
  double dx_H = getDoubleValue(config,"dx_H");

  double eHCAL_L = getDoubleValue(config,"eHCAL_L");

  double dy_ac_L = getDoubleValue(config,"dy_ac_L");
  double dy_ac_H = getDoubleValue(config,"dy_ac_H");

  double run_num_L = getDoubleValue(config,"run_num_L");
  double run_num_H = getDoubleValue(config,"run_num_H");

  double IHWP_flip = getDoubleValue(config,"IHWP_flip");


  //---------------------------------------------------------------------------
  // 2. Set up branch addresses
  //---------------------------------------------------------------------------
  int    Nhit = 0;                                    // Ndata_goodblock_id
  double cid  [512];                                  // goodblock_cidHCAL
  double eHit [512];                                  // goodblock_eHCAL
  double xHit [512];                                  // goodblock_xHCAL
  double yHit [512];                                  // goodblock_yHCAL
  double dx;
  double dy;
  double W2;
  double ePS;
  double vz;
  double eHCAL;
  double coin_time;


  // NOTE: 512 > your max 288, so we’re safe; change if you have more hits.
  T->SetBranchAddress("Ndata_goodblock_id", &Nhit);
  T->SetBranchAddress("goodblock_cidHCAL",  cid );
  T->SetBranchAddress("goodblock_eHCAL",    eHit);
  T->SetBranchAddress("goodblock_xHCAL",    xHit);
  T->SetBranchAddress("goodblock_yHCAL",    yHit);
  T->SetBranchAddress("dx", &dx);
  T->SetBranchAddress("dy", &dy);
  T->SetBranchAddress("W2", &W2);
  T->SetBranchAddress("ePS", &ePS);
  T->SetBranchAddress("vz", &vz);
  T->SetBranchAddress("eHCAL", &eHCAL);
  T->SetBranchAddress("coin_time", &coin_time);


  //---------------------------------------------------------------------------
  // 3. Histogram for distances
  //---------------------------------------------------------------------------
  TH1D *hDist = new TH1D("hDist","Distance of secondary cluster centroids to primary (cid==0); d = sqrt((x-x_{0})^{2}+(y-y_{0})^{2});Counts",200, 0, 5);      // tune bins/range to your detector
  TH1D *hx = new TH1D("hx","hx;",200, -2, 2); 
  TH1D *hy = new TH1D("hy","hy;",200, -2, 2); 
  TH1D *hdx = new TH1D("hdx","hdx",200,-4,3);
  TH1D *hdx_1 = new TH1D("hdx_1","hdx_1",200,-4,3);
  TH2D *hdxvsdist = new TH2D("hdxvsdist","hdxvsdist ; dist ; dx",200,0.05,5,200,-4,4);
  TH2D *henergyvsdist = new TH2D("henergyvsdist","henergyvsdist ; dist (m) ; energy (GeV)",200,0.05,5,200,0,2);
  TH2D *henergyvsdist1 = new TH2D("henergyvsdist1","henergyvsdist1 ; dist (m) ; energy (GeV)",200,0.05,5,200,0,0.5);
  TH2D *henergy0vsdist = new TH2D("henergy0vsdist","henergy0vsdist ; dist (m) ; energy (GeV)",200,0.05,5,200,0,0.7);
  TH1D *hnclus = new TH1D("hnclus","hnclus;number of clusters", 100, 0, 7);

  //---------------------------------------------------------------------------
  // 4. Event loop
  //---------------------------------------------------------------------------
  const Long64_t nEntries = T->GetEntries();
  std::cout << "Processing " << nEntries << " events …\n";

  for(Long64_t ev=0; ev<nEntries; ++ev){
    T->GetEntry(ev);

    //bool goodHelicity = (lookupValue(HelicityCheck, runnum) == 1);
    //bool goodMoller   = (lookupValue(MollerQuality, runnum) == 1);
    bool goodVz       = abs(vz) < 0.27;
    bool goodPS       = (ePS > 0.2);
    //bool goodRunRange = (run_num_L < runnum && runnum < run_num_H);
    bool goodEHCAL    = (eHCAL > eHCAL_L); 
    bool goodW2       = (W2>W2_L && W2<W2_H);
    bool goodCointime = (coin_time>coin_time_L && coin_time<coin_time_H);
    bool goodDx       = (dx>dx_L && dx<dx_H);
    bool goodDy        = (dy>dy_L && dy<dy_H);
    bool goodQE       = goodW2 && goodCointime /*&& goodDx*/ && goodDy;
    bool goodGlobal   = (goodVz && goodPS && goodEHCAL);
    //bool validHel     = (helicity == -1 || helicity == 1);
    //bool goodGrinch = (grinch_track == 0) && (grinch_clus_size>2);
    //bool goodSbs_track = ntrack_sbs>0 && abs(vz_sbs)<0.27;

    if (!goodGlobal || !goodQE) continue;

    if(Nhit <= 0) continue;

    // Accumulators keyed by cid
    std::unordered_map<int,double> sumE, sumEx, sumEy;

    for(int i=0; i<Nhit && i<maxHitsPerEvt; ++i){
      int    ic = static_cast<int>(cid[i]);
      double e  = eHit[i];
      sumE [ic] += e;
      sumEx[ic] += e * xHit[i];
      sumEy[ic] += e * yHit[i];
    }

    const double Ecut = 0.10;              // choose your threshold

    int nStrongClusters = 0;
    for (const auto& kv : sumE) {
      if (kv.second > Ecut) ++nStrongClusters;   // kv.second == ΣE of that cid
    }

    hnclus->Fill(nStrongClusters);

    // --------- identify primary cluster (largest ΣE) ---------------------
    int    primaryCid     = -1;
    double primaryEnergy  = -1.0;

    int    secondaryCid    = -1;   
    double secondaryEnergy = -1.0;    

    for (const auto &kv : sumE) {

      double eTot = kv.second;

      if (eTot > primaryEnergy) {
        // shift the current primary down to secondary slot
        secondaryCid    = primaryCid;
        secondaryEnergy = primaryEnergy;

        primaryCid      = kv.first;
        primaryEnergy   = eTot;

      } else if (eTot > secondaryEnergy) {
        // new runner-up
        secondaryCid    = kv.first;
        secondaryEnergy = eTot;
      }
    }
    if (primaryCid < 0) continue;                       // should not happen

    const double x0 = sumEx[primaryCid] / sumE[primaryCid];
    const double y0 = sumEy[primaryCid] / sumE[primaryCid];

    // --------- loop over all SECONDARY clusters --------------------------
    for (const auto &kv : sumE) {
      int cidNow = kv.first;
      if (cidNow == primaryCid) continue;               // skip primary itself

      if (cidNow != secondaryCid) continue;
      //if (sumE[cidNow]<Ecut) continue;

      double xc = sumEx[cidNow] / sumE[cidNow];
      double yc = sumEy[cidNow] / sumE[cidNow];
      double d  = std::hypot(xc - x0, yc - y0);

      hDist          ->Fill(d);
      hx             ->Fill(xc);
      hy             ->Fill(yc);
      hdxvsdist      ->Fill(d, dx);
      henergyvsdist  ->Fill(d, sumE[cidNow]);
      henergyvsdist1 ->Fill(d, sumE[cidNow]);
      henergy0vsdist ->Fill(d, primaryEnergy);          // primary ΣE vs distance
    }

    int nCluster = static_cast<int>(sumE.size());
    // --------- optional: keep your special dx cut for cid==1 -------------
    
    if (!sumE.count(secondaryCid) || sumE[secondaryCid] > Ecut) {
      hdx_1->Fill(dx);
      //hnclus->Fill(nCluster);
    }

    hdx->Fill(dx);

    if (ev % 1000 == 0) {
      std::cout << (ev * 100.0 / nEntries) << "% \r";
      std::cout.flush();
    }

  }

  //---------------------------------------------------------------------------
  // 5. Draw result
  //---------------------------------------------------------------------------
  TCanvas *c = new TCanvas("cDist","Secondary-to-primary distance",3000,3600);
  TCanvas *c1 = new TCanvas("c1","c1",3000,3600);
  TCanvas *c2 = new TCanvas("c2","c2",3000,3600);
  //c->SetLogy();
  c->Divide(2,2);

  c->cd(1);
  hDist->Draw();
  c->cd(2);
  hx->Draw();
  c->cd(3);
  hy->Draw();


  gPad->SetLogz();
  c1->Divide(2,2);
  c1->cd(1);
  hdxvsdist->Draw("COLZ");
  c1->cd(2);
  henergyvsdist->Draw("COLZ");
  c1->cd(3);
  henergyvsdist1->Draw("COLZ");
  c1->cd(4);
  henergy0vsdist->Draw("COLZ");

  c2->Divide(2,2);
  c2->cd(1);
  hdx->Draw();
  c2->cd(2);
  hdx_1->Draw();
  c2->cd(3);
  hdx->SetLineColor(kBlue);
  hdx_1->SetLineColor(kRed);
  hdx->Draw();
  hdx_1->Draw("same");
  c2->cd(4);
  hnclus->Draw();

  c->SaveAs("cluster_distance.pdf(");      // optional
  c1->SaveAs("cluster_distance.pdf");
  c2->SaveAs("cluster_distance.pdf)");

  std::cout << "Finished.  Histogram saved as 'hDist' (and pdf).\n";

}

