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

#include "models.h"
#include "parse.h"
#include "plotdxdy.h"

void secondary_cluster_analysis(const char* fname          = "input.root",
                             const char* printfilename     = "GEN4b_He3",
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
  std::map<int, int> HelicityCheck = readCSVToMap("DB/Helicity_quality.csv");
  std::map<int, int> MollerQuality = readCSVToMap("DB/Moller_quality.csv");

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
  int    nclus = 0;                                    // Ndata_goodblock_id
  double dx;
  double dy;
  double W2;
  double ePS;
  double vz;
  double eHCAL;
  double hcal_atime;
  double coin_time;
  int runnum;
  double hcal_clus_e[100];
  double hcal_clus_x[100];
  double hcal_clus_y[100];
  double hcal_clus_atime[100];
  double hcal_clus_nblk[100];
  double dx_all_clus[100];
  double dy_all_clus[100];
  int helicity;
  double ntrack_sbs;
  double vz_sbs;


  // NOTE: 512 > your max 288, so we’re safe; change if you have more hits.
  T->SetBranchAddress("Ndata_clus_id", &nclus);
  T->SetBranchAddress("dx", &dx);
  T->SetBranchAddress("dy", &dy);
  T->SetBranchAddress("W2", &W2);
  T->SetBranchAddress("ePS", &ePS);
  T->SetBranchAddress("vz", &vz);
  T->SetBranchAddress("eHCAL", &eHCAL);
  T->SetBranchAddress("runnum", &runnum);
  T->SetBranchAddress("helicity", &helicity);
  T->SetBranchAddress("hcal_atime",&hcal_atime);
  T->SetBranchAddress("coin_time", &coin_time);
  T->SetBranchAddress("hcal_clus_e",hcal_clus_e);
  T->SetBranchAddress("hcal_clus_x",hcal_clus_x);
  T->SetBranchAddress("hcal_clus_y",hcal_clus_y);
  T->SetBranchAddress("hcal_clus_atime",hcal_clus_atime);
  T->SetBranchAddress("hcal_clus_nblk",hcal_clus_nblk);
  T->SetBranchAddress("dx_all_clus",&dx_all_clus);
  T->SetBranchAddress("dy_all_clus",&dy_all_clus);
  T->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
  T->SetBranchAddress("vz_sbs",&vz_sbs);

  TH1D* h_eratio = new TH1D("h_eratio","eratio distribution ; eratio ; number of secondary clusters",100,0,2);
  TH1D* h_eratio_test = new TH1D("h_eratio_test","h_eratio_test",100,0,2);
  TH1D* h_eratio_test1 = new TH1D("h_eratio_test1","h_eratio_test1",100,0,2);
  TH1D* h_eratio_wide_bins = new TH1D("h_eratio_wide_bins","eratio distribution ; eratio ; number of clusters",10,0.0,1.0);
  TH1D* h_eratio_wide_bins_tdiff_cut = new TH1D("h_eratio_wide_bins_tdiff_cut","eratio distribution with a tdiff cut ; eratio ; number of clusters",10,0.0,1.0);
  TH1D* h_tdiff = new TH1D("h_tdiff","tdiff distribution; tdiff (ns);number of clusters",100,-100,100); 
  TH1D* h_nclusters = new TH1D("h_nclusters", "h_nclusters", 100, 0, 10);  
  TH1D* h_nclusters_1 = new TH1D("h_nclusters_1", "h_nclusters_1", 100, 0, 10);
  TH1D* h_nclusters_2 = new TH1D("h_nclusters_2", "h_nclusters_2", 100, 0, 10);
  TH1D* h_nstrong_clusters = new TH1D("h_nstrong_clusters", "number of clusters with in tdiff per event;number of clusters ; number of events", 20, 0, 5);
  TH1D* h_nstrong_clusters_2 = new TH1D("h_nstrong_clusters_2", "number of clusters with in tdiff with eratio>0.2 per event;number of clusters ; number of events", 20, 0, 5);
  TH1D* h_nstrong_clusters_4 = new TH1D("h_nstrong_clusters_4", "number of clusters with in tdiff with eratio>0.4 per event;number of clusters ; number of events", 20, 0, 5);
  TH1D* h_nstrong_clusters_6 = new TH1D("h_nstrong_clusters_6", "number of clusters with in tdiff with eratio>0.6 per event;number of clusters ; number of events", 20, 0, 5);
  TH1D* h_nstrong_clusters_8 = new TH1D("h_nstrong_clusters_8", "number of clusters with in tdiff with eratio>0.8 per event;number of clusters ; number of events", 20, 0, 5);
  TH1D* hdist = new TH1D("hdist","distance from the primary cluster to the secondaries (QE + tdiff cut); dist (m)", 100, 0 , 5);
  TH1D* hdist2 = new TH1D("hdist2","distance from the primary cluster to the secondaries (QE + tdiff cut + eration>0.2); dist (m)", 100, 0 , 5); 
  TH1D* hdist4 = new TH1D("hdist4","distance from the primary cluster to the secondaries (QE + tdiff cut + eration>0.4); dist (m)", 100, 0 , 5);
  TH1D* hdist6 = new TH1D("hdist6","distance from the primary cluster to the secondaries (QE + tdiff cut + eration>0.6); dist (m)", 100, 0 , 5);
  TH1D* hdist8 = new TH1D("hdist8","distance from the primary cluster to the secondaries (QE + tdiff cut + eration>0.8); dist (m)", 100, 0 , 5);

  TH2D* h_energyvsdist = new TH2D("h_energyvsdist","secondary cluster energy vs distance from the primary cluster (QE + tdiff cut); dist(m); energy(GeV)",100,0,5,100,0,0.4);
  TH2D* h_energyvsdist2 = new TH2D("h_energyvsdist2","secondary cluster energy vs distance from the primary cluster (QE + tdiff cut + eration>0.2); dist(m); energy(GeV)",100,0,5,100,0,0.4);
  TH2D* h_energyvsdist4 = new TH2D("h_energyvsdist4","secondary cluster energy vs distance from the primary cluster (QE + tdiff cut + eration>0.4); dist(m); energy(GeV)",100,0,5,100,0,0.4);
  TH2D* h_energyvsdist6 = new TH2D("h_energyvsdist6","secondary cluster energy vs distance from the primary cluster (QE + tdiff cut + eration>0.6); dist(m); energy(GeV)",100,0,5,100,0,0.4);
  TH2D* h_energyvsdist8 = new TH2D("h_energyvsdist8","secondary cluster energy vs distance from the primary cluster (QE + tdiff cut + eration>0.8); dist(m); energy(GeV)",100,0,5,100,0,0.4);

  TH1D* h_dx = new TH1D("h_dx","delta-x distribution with QE cuts ; delta-x (m)",100,-4,3);
  TH1D* h_dx_er2 = new TH1D("h_dx_er2","delta-x distribution for eratio<0.2 with other QE cuts ; delta-x (m)",100,-4,3);
  TH1D* h_dx_er4 = new TH1D("h_dx_er4","delta-x distribution for eratio<0.4 with other QE cuts ; delta-x (m)",100,-4,3);
  TH1D* h_dx_er2_antisbs = new TH1D("h_dx_er2_antisbs","delta-x distribution for eratio<0.2 with other QE cuts and ntrack_sbs==0 ; delta-x (m)",100,-4,3);
  TH1D* h_dx_er4_antisbs = new TH1D("h_dx_er4_antisbs","delta-x distribution for eratio<0.4 with other QE cuts and ntrack_sbs==0 ; delta-x (m)",100,-4,3);

  TH1D* h_W2 = new TH1D("h_W2","W^{2} ;W^{2} (GeV^{2}) ",100,-3,5);
  TH1D* h_coin_time = new TH1D("h_coin_time","coincidence time; coincidence time (ns)",100, coin_time_L-30,coin_time_H+30);
  TH2D* h_dxdy = new TH2D("h_dxdy","delta-x delta-y distribution; delta-y (m); delta-x (m)",100,-4,3,100,-4,3);

  TH2D* h_dxeratio = new TH2D("h_dxeratio","eratio vs delta-x; delta-x (m); eratio",100,-4,3,100,0,1);

  TH2D* h_secondary_dx_dist = new TH2D("h_secondary_dx_dist","distance vs delta-x (QE cuts); dx(m); distance(m)",100,-4,3,100,0,5);
  TH2D* h_secondary_dy_dist = new TH2D("h_secondary_dy_dist","distance vs delta-y (QE cuts); dy(m); distance(m)",100,-4,3,100,0,5);

  TH2D* h_secondary_dx_dist2 = new TH2D("h_secondary_dx_dist2","distance vs delta-x (QE + tdiff + eratio>0.2); dx(m); distance(m)",100,-4,3,100,0,5);
  TH2D* h_secondary_dy_dist2 = new TH2D("h_secondary_dy_dist2","distance vs delta-y (QE + tdiff + eratio>0.2); dy(m); distance(m)",100,-4,3,100,0,5);

  TH2D* h_secondary_dx_dist4 = new TH2D("h_secondary_dx_dist4","distance vs delta-x (QE + tdiff + eratio>0.4); dx(m); distance(m)",100,-4,3,100,0,5);
  TH2D* h_secondary_dy_dist4 = new TH2D("h_secondary_dy_dist4","distance vs delta-y (QE + tdiff + eratio>0.4); dy(m); distance(m)",100,-4,3,100,0,5);

  TH2D* h_secondary_dxdy = new TH2D("h_secondary_dxdy","delta-x vs delta-y for secondary clusters (QE + tdiff) ; dy (m); dx(m)",100,-4,3,100,-4,3);
  TH2D* h_secondary_dxdy2 = new TH2D("h_secondary_dxdy2","delta-x vs delta-y for secondary clusters (QE + tdiff + eratio>0.2) ; dy (m); dx(m)",100,-4,3,100,-4,3);
  TH2D* h_secondary_dxdy4 = new TH2D("h_secondary_dxdy4","delta-x vs delta-y for secondary clusters (QE + tdiff + eratio>0.4) ; dy (m); dx(m)",100,-4,3,100,-4,3);
  TH2D* h_secondary_dxdy6 = new TH2D("h_secondary_dxdy6","delta-x vs delta-y for secondary clusters (QE + tdiff + eratio>0.6) ; dy (m); dx(m)",100,-4,3,100,-4,3);
  TH2D* h_secondary_dxdy8 = new TH2D("h_secondary_dxdy8","delta-x vs delta-y for secondary clusters (QE + tdiff + eratio>0.8) ; dy (m); dx(m)",100,-4,3,100,-4,3);

  int nentries = T->GetEntries();

  for (int i = 0; i<nentries; i++){

    T->GetEntry(i);

    bool goodHelicity = (lookupValue(HelicityCheck, runnum) == 1);
    bool goodMoller   = (lookupValue(MollerQuality, runnum) == 1);
    bool goodRunRange = (run_num_L < runnum && runnum < run_num_H);
    bool goodVz       = abs(vz) < 0.27;
    bool goodPS       = (ePS > 0.2);
    bool goodEHCAL    = (eHCAL > eHCAL_L); 
    bool goodW2       = (W2>W2_L && W2<W2_H);
    bool goodCointime = (coin_time>coin_time_L && coin_time<coin_time_H);
    bool validHelicity = (helicity==-1 || helicity==1);
    bool goodDx       = (dx>dx_L && dx<dx_H);
    bool goodDy       = (dy>dy_L && dy<dy_H);
    bool goodQE       = goodW2 && goodCointime /*&& goodDx*/ && goodDy;
    bool goodGlobal   = (goodVz && goodPS && goodEHCAL);
    bool antisbs = (ntrack_sbs == 0 || abs(vz-vz_sbs)>0.05 || abs(vz)>0.27);

    if(!goodHelicity || !goodMoller || !goodRunRange || !goodGlobal || !validHelicity) continue;

    if (goodW2 and goodDx and goodDy){
      h_coin_time->Fill(coin_time);
    }

    if (goodCointime and goodDx and goodDy){
      h_W2->Fill(W2);
    }

    if(goodCointime and goodW2){
      h_dxdy->Fill(dy,dx);
    }

    if(!goodHelicity || !goodMoller || !goodRunRange || !goodGlobal || !goodQE || !validHelicity) continue;

    h_dx->Fill(dx);

    h_nclusters_2->Fill(nclus);

    if(nclus<2){
      h_nclusters->Fill(nclus);
      h_nclusters_1->Fill(nclus);
      h_nstrong_clusters->Fill(0);
      h_nstrong_clusters_2->Fill(0);
      h_nstrong_clusters_4->Fill(0);
      h_nstrong_clusters_6->Fill(0);
      h_nstrong_clusters_8->Fill(0);
      //h_nclusters_2->Fill(nclus);

    }
    else{
      std::vector<int> order(nclus);
      std::iota(order.begin(), order.end(), 0);          // 0,1,2,…

      std::sort(order.begin(), order.end(),
          [&](int a, int b){ return hcal_clus_e[a] > hcal_clus_e[b]; });

      const int i0 = order[0];          // highest-E (“primary”)
      const int i1 = order[1];          // next-highest (“secondary”)
      double  e0   = hcal_clus_e[i0];                     // cache primary energy
      double  t0   = hcal_clus_atime[i0];
      double  x0   = hcal_clus_x[i0];
      double  y0   = hcal_clus_y[i0];
      int nstrong = 0;
      int nstrong2 = 0;
      int nstrong4 = 0;
      int nstrong6 = 0;
      int nstrong8 = 0;

      for (int k = 1; k < nclus; ++k) {                   // start at 1 ⇒ “secondary”
          int    ik   = order[k];                         // kth-highest-E cluster
          double er   = hcal_clus_e[ik]   / e0;           // E_secondary / E_primary
          double dt   = t0 - hcal_clus_atime[ik];         // Δt = t_primary – t_sec
          double er0  = hcal_clus_e[ik]   / eHCAL;        // optional cross-check
          double dist = std::hypot(hcal_clus_x[ik] - x0, hcal_clus_y[ik] - y0);

          double dxi = dx_all_clus[ik];
          double dyi = dy_all_clus[ik];

          h_eratio     ->Fill(er);                        // all secondary ratios
          h_eratio_wide_bins->Fill(er);
          h_tdiff      ->Fill(dt);                        // all secondary Δt
          h_eratio_test->Fill(er0);                       // sec / total-cluster-E
          //hdist->Fill(dist);
          h_secondary_dx_dist->Fill(dxi,dist);
          h_secondary_dy_dist->Fill(dyi,dist);

          bool gooddt = abs(dt+3)<2;

          if(!gooddt){
            h_dxeratio->Fill(dxi,er);

          }

          if(gooddt){
            hdist->Fill(dist);
            h_energyvsdist->Fill(dist,hcal_clus_e[ik]);
            h_eratio_wide_bins_tdiff_cut->Fill(er);
            h_secondary_dxdy->Fill(dyi,dxi);
            nstrong++;
          }

          if(er > 0.2 and gooddt ){
            hdist2->Fill(dist);
            h_energyvsdist2->Fill(dist,hcal_clus_e[ik]);
            h_secondary_dxdy2->Fill(dyi,dxi);
            h_secondary_dx_dist2->Fill(dxi,dist);
            h_secondary_dy_dist2->Fill(dyi,dist);
            nstrong2++;
          }
          if (er > 0.4 and gooddt ){
            hdist4->Fill(dist);
            h_energyvsdist4->Fill(dist,hcal_clus_e[ik]);
            h_secondary_dxdy4->Fill(dyi,dxi);
            h_secondary_dx_dist4->Fill(dxi,dist);
            h_secondary_dy_dist4->Fill(dyi,dist);
            nstrong4++;
          }
          if (er > 0.6 and gooddt){
            hdist6->Fill(dist);
            h_energyvsdist6->Fill(dist,hcal_clus_e[ik]);
            h_secondary_dxdy6->Fill(dyi,dxi);
            nstrong6++;
          }
          if (er > 0.8 and gooddt){
            hdist8->Fill(dist);
            h_energyvsdist8->Fill(dist,hcal_clus_e[ik]);
            h_secondary_dxdy8->Fill(dyi,dxi);
            nstrong8++;
          }

          // If you still want the “eratio_test1” histogram to hold the *same*
          // information but for clarity keep it, just duplicate the fill:
          //h_eratio_test1->Fill(er0);
      }
      h_nstrong_clusters->Fill(nstrong);
      h_nstrong_clusters_2->Fill(nstrong2);
      h_nstrong_clusters_4->Fill(nstrong4);
      h_nstrong_clusters_6->Fill(nstrong6);
      h_nstrong_clusters_8->Fill(nstrong8);

      if (nstrong2 == 0){
        h_dx_er2->Fill(dx);
        if (antisbs){
          h_dx_er2_antisbs->Fill(dx);
        }

      }
      if (nstrong4 == 0){
        h_dx_er4->Fill(dx);
        if (antisbs){
          h_dx_er4_antisbs->Fill(dx);
        }
      }


      double eratio_secondary = hcal_clus_e[i1] / hcal_clus_e[i0];
      double tdiff_secondary  = hcal_clus_atime[i0] - hcal_clus_atime[i1];

      bool goodtdiff_secondary = abs(tdiff_secondary+3)<2;

      //double eratio_test = hcal_clus_e[i0]/eHCAL;

      double eratio_test1 = hcal_clus_e[i1]/eHCAL;
      h_eratio_test1->Fill(eratio_test1);

      //h_eratio->Fill(eratio);
      //h_tdiff ->Fill(tdiff);
      //h_eratio_test->Fill(eratio_test);

      if(eratio_secondary > 0.2 and  goodtdiff_secondary){
        h_nclusters->Fill(nclus);
      }

      if(hcal_clus_e[i1] > 0.1 and goodtdiff_secondary ){
        h_nclusters_1->Fill(nclus);
      }

      //h_nclusters_2->Fill(nclus);

    }

    if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
    std::cout.flush();

  }

  // ---- Put this near the top of the macro, before you create canvases ----
  //gStyle->SetPadLeftMargin(0.18);   // 18 % of the pad width – tweak to taste
  //gStyle->SetPadRightMargin(0.03);  // keep the right side slim
  //gStyle->SetTitleOffset(1.3,"Y");  // pushes the Y–axis title in a bit

  TCanvas* c = new TCanvas("c","c",3000,3600);
  TCanvas* c1 = new TCanvas("c1","c1",3000,3600);
  TCanvas* c2 = new TCanvas("c2","c2",3000,3600);
  TCanvas* c3 = new TCanvas("c3","c3",3000,3600);
  TCanvas* c4 = new TCanvas("c4","c4",3000,3600);
  TCanvas* c5 = new TCanvas("c5","c5",3000,3600);
  TCanvas* c6 = new TCanvas("c6","c6",3000,3600);
  TCanvas* c7 = new TCanvas("c7","c7",3000,3600);
  TCanvas* c8 = new TCanvas("c8","c8",3000,3600);
  TCanvas* c9 = new TCanvas("c9","c9",3000,3600);
  TCanvas* c10 = new TCanvas("c10","c10",3000,3600);
  TCanvas* c11 = new TCanvas("c11", "c11",3000,3600);


  TLine *line1 = new TLine(-5,0.0,-5,h_tdiff->GetMaximum());
  line1->SetLineColor(kRed);
  line1->SetLineWidth(2);
  
  TLine *line2 = new TLine(-1,0.0,-1,h_tdiff->GetMaximum());
  line2->SetLineColor(kRed);
  line2->SetLineWidth(2);

  c->Divide(2,2);
  c->cd(1);
  h_eratio->Draw();
  c->cd(2);
  h_tdiff->Draw();
  line1->Draw("same");
  line2->Draw("same");
  c->cd(3);
  h_eratio_test->Draw();
  c->cd(4);
  h_eratio_test1->Draw();

  c1->Divide(2,2);
  c1->cd(1);
  h_nclusters->Draw();
  c1->cd(2);
  hdist->Draw();
  c1->cd(3);
  h_energyvsdist->Draw("COLZ");
  c1->cd(4);
  h_nclusters_1->Draw();

  c2->Divide(2,2);
  c2->cd(1);
  h_eratio_wide_bins->Draw();
  c2->cd(2);
  h_eratio_wide_bins_tdiff_cut->Draw();
  c2->cd(3);
  h_nclusters_2->Draw();
  c2->cd(4);
  h_nstrong_clusters->Draw();

  TLine *line3 = new TLine(coin_time_L,0.0,coin_time_L,h_coin_time->GetMaximum());
  line3->SetLineColor(kRed);
  line3->SetLineWidth(2);
  
  TLine *line4 = new TLine(coin_time_H,0.0,coin_time_H,h_coin_time->GetMaximum());
  line4->SetLineColor(kRed);
  line4->SetLineWidth(2);

  TLine *line5 = new TLine(W2_L,0.0,W2_L,h_W2->GetMaximum());
  line5->SetLineColor(kRed);
  line5->SetLineWidth(2);
  
  TLine *line6 = new TLine(W2_H,0.0,W2_H,h_W2->GetMaximum());
  line6->SetLineColor(kRed);
  line6->SetLineWidth(2);

  TBox *box = new TBox(dy_L,dx_L,dy_H,dy_H);
  box->SetLineColor(kRed);
  box->SetFillStyle(0);
  box->SetLineWidth(3);

  c3->Divide(2,2);
  c3->cd(1);
  h_coin_time->Draw();
  line3->Draw("same");
  line4->Draw("same");
  c3->cd(2);
  h_W2->Draw();
  line5->Draw("same");
  line6->Draw("same");
  c3->cd(3);
  h_dxdy->Draw("COLZ");
  box->Draw("same");

  c4->Divide(2,2);
  c4->cd(1);
  h_nstrong_clusters_2->Draw();
  c4->cd(2);
  h_nstrong_clusters_4->Draw();
  c4->cd(3);
  h_nstrong_clusters_6->Draw();
  c4->cd(4);
  h_nstrong_clusters_8->Draw();

  c5->Divide(2,2);
  c5->cd(1);
  h_energyvsdist2->Draw("COLZ");
  c5->cd(2);
  h_energyvsdist4->Draw("COLZ");
  c5->cd(3);
  h_energyvsdist6->Draw("COLZ");
  c5->cd(4);
  h_energyvsdist8->Draw("COLZ");

  c6->Divide(2,2);
  c6->cd(1);
  h_dx_er2->Draw();
  c6->cd(2);
  h_dx_er4->Draw();
  c6->cd(3);
  h_dx->SetLineColor(kAzure);
  h_dx_er2->SetLineColor(kOrange);
  h_dx_er4->SetLineColor(kViolet);
  h_dx_er2_antisbs->SetLineColor(kBlack);
  h_dx_er4_antisbs->SetLineColor(kCyan);
  h_dx->Draw();
  h_dx_er2->Draw("same");
  h_dx_er4->Draw("same");
  h_dx_er2_antisbs->Draw("same");
  h_dx_er4_antisbs->Draw("same");

  c6->cd(4);
  // -------------------------------------------------------------
  TLegend *leg = new TLegend(0.3, 0.3,   // x1, y1  (lower-left corner)
                             0.8, 0.8);  // x2, y2  (upper-right corner)
  leg->SetBorderSize(0);      // no frame
  leg->SetFillStyle(0);       // no grey background
  leg->SetTextFont(42);       // same font as axis labels

  leg->AddEntry(h_dx,      "primary clusters",                "l");
  leg->AddEntry(h_dx_er2,  "E_{sec}/E_{prim} > 0.2",          "l");
  leg->AddEntry(h_dx_er4,  "E_{sec}/E_{prim} > 0.4",          "l");

  leg->Draw();
  gPad->Update();    

  c7->Divide(2,2);
  c7->cd(1);
  hdist2->Draw();
  c7->cd(2);
  hdist4->Draw();
  c7->cd(3);
  hdist6->Draw();
  c7->cd(4);
  hdist8->Draw();

  c8->Divide(2,2);
  c8->cd(1);
  h_dxeratio->Draw("COLZ");
  c8->cd(2);
  h_secondary_dxdy->Draw("COLZ");
  c8->cd(3);
  h_secondary_dx_dist->Draw("COLZ");
  c8->cd(4);
  h_secondary_dy_dist->Draw("COLZ");

  c9->Divide(2,2);
  c9->cd(1);
  h_secondary_dxdy2->Draw("COLZ");
  c9->cd(2);
  h_secondary_dxdy4->Draw("COLZ");
  c9->cd(3);
  h_secondary_dxdy6->Draw("COLZ");
  c9->cd(4);
  h_secondary_dxdy8->Draw("COLZ");

  c10->Divide(2,2);
  c10->cd(1);
  h_secondary_dx_dist2->Draw("COLZ");
  c10->cd(2);
  h_secondary_dy_dist2->Draw("COLZ");
  c10->cd(3);
  h_secondary_dx_dist4->Draw("COLZ");
  c10->cd(4);
  h_secondary_dy_dist4->Draw("COLZ");


  c11->Divide(2,2);
  c11->cd(1); 
  gPad->SetLogy(1);
  h_dx->SetLineColor(kAzure);
  h_dx_er2->SetLineColor(kOrange);
  h_dx_er4->SetLineColor(kViolet);
  h_dx_er2_antisbs->SetLineColor(kBlack);
  h_dx_er4_antisbs->SetLineColor(kCyan);
  h_dx->Draw();
  h_dx_er2->Draw("same");
  h_dx_er4->Draw("same");
  //h_dx_er2_antisbs->Draw("same");
  //h_dx_er4_antisbs->Draw("same");

  c11->cd(2);
  // -------------------------------------------------------------
  TLegend *leg1 = new TLegend(0.3, 0.3,   // x1, y1  (lower-left corner)
                             0.8, 0.8);  // x2, y2  (upper-right corner)
  leg1->SetBorderSize(0);      // no frame
  leg1->SetFillStyle(0);       // no grey background
  leg1->SetTextFont(42);       // same font as axis labels

  leg1->AddEntry(h_dx,      "primary clusters",                "l");
  leg1->AddEntry(h_dx_er2,  "E_{sec}/E_{prim} < 0.2",          "l");
  leg1->AddEntry(h_dx_er4,  "E_{sec}/E_{prim} < 0.4",          "l");
  leg1->AddEntry(h_dx_er2_antisbs, "E_{sec}/E_{prim} < 0.2 && antisbs", "l");
  leg1->AddEntry(h_dx_er4_antisbs, "E_{sec}/E_{prim} < 0.4 && antisbs", "l");

  leg1->Draw();
  gPad->Update();    

  auto fixMargins = [](TCanvas* c){
    for (int ipad = 1; ipad <= c->GetListOfPrimitives()->GetSize(); ++ipad){
      c->cd(ipad);
      gPad->SetLeftMargin(0.18);   // or any value that keeps the labels visible
      gPad->SetBottomMargin(0.12); // optional: room for X–axis labels
    }
  };

  fixMargins(c);
  fixMargins(c1);
  fixMargins(c2);
  fixMargins(c3);
  fixMargins(c4);
  fixMargins(c5);
  fixMargins(c6);
  fixMargins(c7);
  fixMargins(c8);
  fixMargins(c9);
  fixMargins(c10);
  //fixMargins(c11);

  c->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf(",printfilename));
  c1->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c2->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c3->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c4->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c5->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c6->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c7->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c8->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c9->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c10->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf",printfilename));
  c11->SaveAs(Form("plots/%s_secondary_cluster_analysis.pdf)",printfilename));

}

