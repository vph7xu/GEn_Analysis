#include <iostream>
#include <locale>
#include <cmath>
#include <vector>
#include <memory>

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraph.h>
#include <TF1.h>
#include <TLegend.h>
#include <TEllipse.h>
#include <TPaveText.h>
#include <TLine.h>
#include <TColor.h>
#include <TMath.h>

#include "cuts.h"
#include "plotdxdy.h"

void PS_HCAL_e_plot_twofiles(const char* filename,
                             const char* filename1,
                             const char* printfilename,
                             const char* kin,
                             const char* kin1) {

  // -------------------- open files/trees (fixed) --------------------
  TFile* file  = TFile::Open(filename);
  if (!file || file->IsZombie()) { std::cerr<<"ERROR: can't open "<<filename<<"\n"; return; }
  TTree* tree  = (TTree*)file->Get("Tout");
  if (!tree) { std::cerr<<"ERROR: tree 'Tout' not found in "<<filename<<"\n"; return; }

  TFile* file1 = TFile::Open(filename1);                            // <- use filename1
  if (!file1 || file1->IsZombie()) { std::cerr<<"ERROR: can't open "<<filename1<<"\n"; return; }
  TTree* tree1 = (TTree*)file1->Get("Tout");                         // <- get from file1
  if (!tree1) { std::cerr<<"ERROR: tree 'Tout' not found in "<<filename1<<"\n"; return; }

  // -------------------- load cut/config DBs -------------------------
  std::map<std::string, std::string> config  = parseConfig(Form("cuts/cut_%s.txt",kin));
  std::map<int,int> HelicityCheck = readCSVToMap("DB/Helicity_quality.csv");
  std::map<int,int> MollerQuality = readCSVToMap("DB/Moller_quality.csv");

  double coin_time_L = getDoubleValue(config,"coin_time_L");
  double coin_time_H = getDoubleValue(config,"coin_time_H");
  double W2_L        = getDoubleValue(config,"W2_L");
  double W2_H        = getDoubleValue(config,"W2_H");
  double dy_L        = getDoubleValue(config,"dy_L");
  double dy_H        = getDoubleValue(config,"dy_H");
  double dx_L        = getDoubleValue(config,"dx_L");
  double dx_H        = getDoubleValue(config,"dx_H");
  double dx_p_L      = getDoubleValue(config,"dx_p_L");
  double dx_p_H      = getDoubleValue(config,"dx_p_H");
  double eHCAL_L     = getDoubleValue(config,"eHCAL_L");
  double dy_ac_L     = getDoubleValue(config,"dy_ac_L");
  double dy_ac_H     = getDoubleValue(config,"dy_ac_H");
  double run_num_L   = getDoubleValue(config,"run_num_L");
  double run_num_H   = getDoubleValue(config,"run_num_H");
  double IHWP_flip   = getDoubleValue(config,"IHWP_flip");

  std::map<std::string, std::string> config1 = parseConfig(Form("cuts/cut_%s.txt",kin1));
  double coin_time_L_1 = getDoubleValue(config1,"coin_time_L");
  double coin_time_H_1 = getDoubleValue(config1,"coin_time_H");
  double W2_L_1        = getDoubleValue(config1,"W2_L");
  double W2_H_1        = getDoubleValue(config1,"W2_H");
  double dy_L_1        = getDoubleValue(config1,"dy_L");
  double dy_H_1        = getDoubleValue(config1,"dy_H");
  double dx_L_1        = getDoubleValue(config1,"dx_L");
  double dx_H_1        = getDoubleValue(config1,"dx_H");
  double dx_p_L_1      = getDoubleValue(config1,"dx_p_L");
  double dx_p_H_1      = getDoubleValue(config1,"dx_p_H");
  double eHCAL_L_1     = getDoubleValue(config1,"eHCAL_L");
  double dy_ac_L_1     = getDoubleValue(config1,"dy_ac_L");
  double dy_ac_H_1     = getDoubleValue(config1,"dy_ac_H");
  double run_num_L_1   = getDoubleValue(config1,"run_num_L");
  double run_num_H_1   = getDoubleValue(config1,"run_num_H");
  double IHWP_flip_1   = getDoubleValue(config1,"IHWP_flip");

  // -------------------- variable blocks (typed) ---------------------
  // first file
  Double_t ePS=0, eSH=0, trP=0, eHCAL=0, vz=0, vy=0, vx=0, ntrack_chi2=0, track_nhits=0, nhits=0;
  Double_t W2=0, coin_time=0, dx=0, dy=0, thtgt=0, thetabend=0, ebeam=0, etheta=0;
  Double_t BPMAy=0, BPMAx=0, rastery=0, rasterx=0, raster2y=0, raster2x=0;
  Int_t    runnum=0, helicity=0, ntrack=0;

  tree->SetBranchAddress("ePS",&ePS);
  tree->SetBranchAddress("eSH",&eSH);
  tree->SetBranchAddress("trP",&trP);
  tree->SetBranchAddress("eHCAL",&eHCAL);
  tree->SetBranchAddress("vz",&vz);
  tree->SetBranchAddress("vy",&vy);
  tree->SetBranchAddress("vx",&vx);
  tree->SetBranchAddress("ntrack",&ntrack);
  tree->SetBranchAddress("ntrack_chi2",&ntrack_chi2);
  tree->SetBranchAddress("ntrack_hits",&track_nhits);
  tree->SetBranchAddress("nhits",&nhits);
  tree->SetBranchAddress("W2",&W2);
  tree->SetBranchAddress("coin_time",&coin_time);
  tree->SetBranchAddress("runnum",&runnum);
  tree->SetBranchAddress("helicity",&helicity);
  tree->SetBranchAddress("dx",&dx);
  tree->SetBranchAddress("dy",&dy);
  tree->SetBranchAddress("thtgt",&thtgt);
  tree->SetBranchAddress("thetabend",&thetabend);
  tree->SetBranchAddress("ebeam",&ebeam);
  tree->SetBranchAddress("etheta",&etheta);
  tree->SetBranchAddress("BPMAy",&BPMAy);
  tree->SetBranchAddress("BPMAx",&BPMAx);
  tree->SetBranchAddress("Rastery",&rastery);
  tree->SetBranchAddress("Rasterx",&rasterx);
  tree->SetBranchAddress("Raster2y",&raster2y);
  tree->SetBranchAddress("Raster2x",&raster2x);

  // second file (fixed to use tree1)
  Double_t ePS_1=0, eSH_1=0, trP_1=0, eHCAL_1=0, vz_1=0, vy_1=0, vx_1=0, track_nhits_1=0, track_chi2_1=0, nhits_1=0;
  Double_t W2_1=0, coin_time_1=0, dx_1=0, dy_1=0, thtgt_1=0, thetabend_1=0, ebeam_1=0, etheta_1=0;
  Double_t BPMAy_1=0, BPMAx_1=0, rastery_1=0, rasterx_1=0, raster2y_1=0, raster2x_1=0;
  Int_t    runnum_1=0, helicity_1=0, ntrack_1=0;

  tree1->SetBranchAddress("ePS",&ePS_1);
  tree1->SetBranchAddress("eSH",&eSH_1);
  tree1->SetBranchAddress("trP",&trP_1);
  tree1->SetBranchAddress("eHCAL",&eHCAL_1);
  tree1->SetBranchAddress("vz",&vz_1);
  tree1->SetBranchAddress("vy",&vy_1);
  tree1->SetBranchAddress("vx",&vx_1);
  tree1->SetBranchAddress("ntrack",&ntrack_1);
  tree1->SetBranchAddress("ntrack_chi2",&track_chi2_1);
  tree1->SetBranchAddress("ntrack_hits",&track_nhits_1);
  tree1->SetBranchAddress("nhits",&nhits_1);
  tree1->SetBranchAddress("W2",&W2_1);
  tree1->SetBranchAddress("coin_time",&coin_time_1);
  tree1->SetBranchAddress("runnum",&runnum_1);
  tree1->SetBranchAddress("helicity",&helicity_1);
  tree1->SetBranchAddress("dx",&dx_1);
  tree1->SetBranchAddress("dy",&dy_1);
  tree1->SetBranchAddress("thtgt",&thtgt_1);
  tree1->SetBranchAddress("thetabend",&thetabend_1);
  tree1->SetBranchAddress("ebeam",&ebeam_1);
  tree1->SetBranchAddress("etheta",&etheta_1);
  tree1->SetBranchAddress("BPMAy",&BPMAy_1);
  tree1->SetBranchAddress("BPMAx",&BPMAx_1);
  tree1->SetBranchAddress("Rastery",&rastery_1);
  tree1->SetBranchAddress("Rasterx",&rasterx_1);
  tree1->SetBranchAddress("Raster2y",&raster2y_1);
  tree1->SetBranchAddress("Raster2x",&raster2x_1);

  // -------------------- histograms (unchanged, one unit tweak) ---------------
  TH1D *h_ePS = new TH1D("h_ePS","Pre-shower energy;Pre-shower (GeV)",200,0,2.5);
  TH1D *h_eHCAL = new TH1D("h_eHCAL","HCal energy;HCal energy (GeV)",100,0,1);
  TH1D *h_trvz = new TH1D("h_trvz","Vertex z;Vertex z (m)",200,-0.8,0.8);
  TH1D *h_ntrack = new TH1D("h_ntrack","ntrack",60,0,6);
  TH1D *h_nhits = new TH1D("h_nhits","nhits",60,0,6);
  TH1D *h_W2 = new TH1D("h_W2","W^{2}; W^{2} (GeV^{2})",200,-2,4);
  TH1D *h_W2_0 = new TH1D("h_W2_0","W^{2}; W^{2} (GeV^{2})",50,-1,3);
  TH1D *h_W2_1 = new TH1D("h_W2_1","W^{2}; W^{2} (GeV^{2})",200,-2,4);
  TH1D *h_W2_2 = new TH1D("h_W2_2","W^{2}; W^{2} (GeV^{2})",200,-2,4);
  TH1D *h_W2_3 = new TH1D("h_W2_3","W^{2}; W^{2} (GeV^{2})",200,-2,4);

  TH1D *h_dx = new TH1D("h_dx","dx; dx (m)", 100,-4,3);
  TH1D *h_dy = new TH1D("h_dy","dy; dy (m)", 100,-4,3);

  TH1D *h_dx_0 = new TH1D("h_dx_0","dx; dx (m)",100,-4,3);
  TH1D *h_dx_1 = new TH1D("h_dx_1","dx; dx (m)",100,-4,3);
  TH1D *h_dx_2 = new TH1D("h_dx_2","dx; dx (m)",100,-4,3);
  TH1D *h_dx_3 = new TH1D("h_dx_3","dx; dx (m)",100,-4,3);

  TH1D *h_W2_0_0 = new TH1D("h_W2_0_0","W^{2}; W^{2} (GeV^{2})",200,-2,4);
  TH1D *h_W2_1_1 = new TH1D("h_W2_1_0","W^{2}; W^{2} (GeV^{2})",200,-2,4);
  TH1D *h_W2_2_2 = new TH1D("h_W2_2_0","W^{2}; W^{2} (GeV^{2})",200,-2,4);
  TH1D *h_W2_3_3 = new TH1D("h_W2_3_0","W^{2}; W^{2} (GeV^{2})",200,-2,4);

  TH1D *h_dx_0_0 = new TH1D("h_dx_0_0","dx; dx (m)",100,-4,3);
  TH1D *h_dx_1_1 = new TH1D("h_dx_1_0","dx; dx (m)",100,-4,3);
  TH1D *h_dx_2_2 = new TH1D("h_dx_2_0","dx; dx (m)",100,-4,3);
  TH1D *h_dx_3_3 = new TH1D("h_dx_3_0","dx; dx (m)",100,-4,3);

  TH1D *h_cointime = new TH1D("h_cointime","Coincidence time; Coincidence time (ns)",200,coin_time_L-50,coin_time_H+50);
  TH1D *h_eoverp  = new TH1D("h_eoverp","E/p BBCal ; E/p", 100,0.3,1.7);
  TH1D *h_chi2_all = new TH1D("h_chi2_all","#chi^{2} ; #chi^{2}",100,0,60);
  TH1D *h_chi2_3hit = new TH1D("h_chi2_3hit","#chi^{2} ; #chi^{2}",100,0,60);
  TH1D *h_chi2_4hit = new TH1D("h_chi2_4hit","#chi^{2} ; #chi^{2}",100,0,60);
  TH1D *h_chi2_5hit = new TH1D("h_chi2_5hit","#chi^{2} ; #chi^{2}",100,0,60);

  TH1D *h_trvz_c = new TH1D("h_trvz_c","vertex z",1000,-1,1);
  TH2D *h_dxdy   = new TH2D("h_dxdy","dx vs dy ; dy (m) ; dx (m)",100,-4,3,100,-4,3);
  TH2D *h_dxdy_1   = new TH2D("h_dxdy_1","dx vs dy ; dy (m) ; dx (m)",100,-4,3,100,-4,3);

  TH1D *h_W2_new  = new TH1D("h_W2_new","W^{2} new ; W^{2} (GeV^{2})",50,-1,3);
  TH1D *h_W2_new1 = new TH1D("h_W2_new1","W^{2} new1 ; W^{2} (GeV^{2})",50,-1,3);
  TH1D *h_W2_new_1  = new TH1D("h_W2_new_1","W^{2} new ; W^{2} (GeV^{2})",50,-1,3);
  TH1D *h_W2_new_1_1= new TH1D("h_W2_new_1_1","W^{2} new1 ; W^{2} (GeV^{2})",50,-1,3);

  TH1D *h_trP     = new TH1D("h_trP","trP; p (GeV)",100,0,7);
  TH1D *h_trP_new = new TH1D("h_trP_new","p_{el} est. A,B; p (GeV)",100,0,7);
  TH1D *h_trP_new1= new TH1D("h_trP_new1","p_{el} est. A1,B1; p (GeV)",100,0,7);

  TH2D *h_W2_v_thetatgt = new TH2D("h_W2_v_thetatgt","W^{2} vs #theta_{tgt};#theta_{tgt} (rad);W^{2} (GeV^{2})",100,-0.3,0.3,100,-2,5);
  TH2D *h_W2timesthetabend_v_thetatgt =
        new TH2D("h_W2timesthetabend_v_thetatgt","W^{2}#times#theta_{bend} vs #theta_{tgt};#theta_{tgt} (rad);W^{2}#theta_{bend}",100,-0.3,0.3,100,-0.1,0.6);

  TH2D *h_p_eltimesthetabend_v_thetatgt =
        new TH2D("h_p_eltimesthetabend_v_thetatgt","p_{el}#times#theta_{bend} vs #theta_{tgt};#theta_{tgt} (rad);p_{el}(GeV)#theta_{bend} (rad)", 100,-0.3,0.3,100,-0.1,0.6);

  TH2D *h_trPtimesthetabend_v_thetatgt =
        new TH2D("h_trPtimesthetabend_v_thetatgt","p_{trk}#times#theta_{bend} vs #theta_{tgt};#theta_{tgt} (rad);p_{trk}(GeV)#theta_{bend} (rad)", 100,-0.3,0.3,100,-0.1,0.6);

  TH2D *h_p_el_v_BPMAy = new TH2D("h_p_el_v_BPMAy","#Delta p vs BPMAy;BPMAy (m);p/p_{el} - 1",50,-0.003,0.003,50,-0.1,0.1);
  TH2D *h_p_el_v_BPMAx = new TH2D("h_p_el_v_BPMAx","#Delta p vs BPMAx;BPMAx (m);p/p_{el} - 1",50,-0.003,0.003,50,-0.1,0.1);
  TH2D *h_p_el_v_vy    = new TH2D("h_p_el_v_vy","#Delta p vs v_{y};v_{y} (m);p/p_{el} - 1",50,-0.003,0.003,50,-0.1,0.1);
  TH2D *h_p_el_v_vx    = new TH2D("h_p_el_v_vx","#Delta p vs v_{x};v_{x} (m);p/p_{el} - 1",50,-0.003,0.003,50,-0.1,0.1);
  TH2D *h_p_el_v_rastery = new TH2D("h_p_el_v_rastery","#Delta p vs rastery;rastery;p/p_{el} - 1",50,34000,54000,50,-0.1,0.1);
  TH2D *h_p_el_v_rasterx = new TH2D("h_p_el_v_rasterx","#Delta p vs rasterx;rasterx;p/p_{el} - 1",50,34000,54000,50,-0.1,0.1);
  TH2D *h_p_el_v_raster2y= new TH2D("h_p_el_v_raster2y","#Delta p vs raster2y;raster2y;p/p_{el} - 1",50,34000,54000,50,-0.1,0.1);
  TH2D *h_p_el_v_raster2x= new TH2D("h_p_el_v_raster2x","#Delta p vs raster2x;raster2x;p/p_{el} - 1",50,34000,54000,50,-0.1,0.1);

  TH1D *h_dp      = new TH1D("h_dp","#Delta p; p/p_{el} - 1",100,-0.25,0.25);
  TH1D *h_dp_new  = new TH1D("h_dp_new","#Delta p (new); p/p_{el} - 1",100,-0.25,0.25);
  TH1D *h_dp_new1 = new TH1D("h_dp_new1","#Delta p (new1); p/p_{el} - 1",100,-0.25,0.25);

  TH1D *h_ebeam = new TH1D("h_ebeam","Beam energy; Beam energy (GeV)",100,8.4,8.5);

  // -------------------- constants --------------------
  const double A1 = 0.28640;//0.28661;
  const double B1 = 0.87122905;//0.839712501;
  const double A  = 0.28640;//0.287186374;
  const double B  = 0.87122905;//0.798543608;
  const double m  = 0.938; // GeV

  double C = -0.00184*3.1;
  double D = -5.17872*3.1;

  // ellipse center/axes (dx vs dy plot: x=dy, y=dx)
  const double x0 = 0.5*(dx_p_L + dx_p_H);   // y-axis (dx)
  const double y0 = 0.5*(dy_L    + dy_H);    // x-axis (dy)
  const double rx = 0.2*1.5;
  const double ry = 0.2*1.5;

  // -------------------- loop #1 --------------------
  const int nentries = tree->GetEntries();
  for (int i = 0; i<nentries; i++){
    tree->GetEntry(i);

    const double eoverp = (eSH+ePS)/trP;

    const double ee  = A *(1+B *thtgt)/thetabend - (vy*D + C);
    const double ee1 = A1*(1+B1*thtgt)/thetabend - (vy*D + C);

    const double Q2_new  = 2*ebeam*ee *(1-std::cos(etheta));
    const double Q2_new1 = 2*ebeam*ee1*(1-std::cos(etheta));

    const double nu_new  = (ebeam - ee);
    const double nu_new1 = (ebeam - ee1);

    const double W2_new  = m*m+2*m*nu_new -Q2_new;
    const double W2_new1 = m*m+2*m*nu_new1-Q2_new1;

    const double p_el    = ebeam/(1+ebeam*(1-std::cos(etheta))/m);

    bool goodHelicity = (lookupValue(HelicityCheck, runnum) == 1);
    bool goodMoller   = (lookupValue(MollerQuality, runnum) == 1);
    bool goodVz       = std::abs(vz) < 0.27;
    bool goodPS       = (ePS > 0.2);
    bool goodRunRange = (run_num_L < runnum && runnum < run_num_H);
    bool goodEHCAL    = (eHCAL > 0.225);
    bool validHel     = (helicity == -1 || helicity == 1);
    bool goodEoverp   = std::abs(eoverp-1) < 0.2;
    bool goodW2       = (W2_L < W2 && W2 < W2_H);
    bool goodCoin     = (coin_time_L < coin_time && coin_time < coin_time_H);
    bool goodDy       = (-0.2 < dy && dy < 0.2);
    bool goodpSpot    = ( std::pow((dy-y0)/ry,2) + std::pow((dx-x0)/rx,2) ) <= 1.0;
    bool goodSpot     = dy_L<dy && dy<dy_H && ((dx_L<dx && dx<dx_H) || (dx_p_L<dx && dx<dx_p_H));

    if(ePS>0 && eHCAL>0){
      h_trvz->Fill(vz);

      if(goodVz && goodEHCAL && goodCoin) h_ePS->Fill(ePS);
      if(goodVz && goodCoin && goodPS)    h_eHCAL->Fill(eHCAL);

      h_ntrack->Fill(ntrack);
      h_nhits->Fill(nhits);

      if(goodVz && goodPS && goodEoverp && goodEHCAL && goodCoin) h_W2->Fill(W2);

      if(goodVz && goodPS /*&& goodEoverp*/ && goodEHCAL /*&& goodW2*/ && goodCoin && W2<2.5){
        h_dy->Fill(dy);
        h_dxdy->Fill(dy,dx);
        if(goodDy) h_dx->Fill(dx);
      }

      if(goodVz && goodPS && goodEHCAL && goodCoin){
        if(goodpSpot){
          h_W2_0->Fill(W2);
          h_W2_new->Fill(W2_new);
          h_W2_new1->Fill(W2_new1);

          h_trP->Fill(trP);
          h_trP_new->Fill(ee);
          h_trP_new1->Fill(ee1);

          h_W2_v_thetatgt->Fill(thtgt,W2);
          h_W2timesthetabend_v_thetatgt->Fill(thtgt,W2*thetabend);

          h_p_eltimesthetabend_v_thetatgt->Fill(thtgt,p_el*thetabend);
          h_trPtimesthetabend_v_thetatgt->Fill(thtgt,trP*thetabend);

          h_p_el_v_BPMAy->Fill(BPMAy,ee/p_el - 1);
          h_p_el_v_BPMAx->Fill(BPMAx,ee/p_el - 1);
          h_p_el_v_vy->Fill(vy,ee/p_el - 1);
          h_p_el_v_vx->Fill(vx,ee/p_el - 1);

          h_dp->Fill(trP/p_el - 1);
          h_dp_new->Fill(ee/p_el - 1);
          h_dp_new1->Fill(ee1/p_el - 1);

          h_p_el_v_rastery->Fill(rastery,ee/p_el - 1);
          h_p_el_v_rasterx->Fill(rasterx,ee/p_el - 1);
          h_p_el_v_raster2y->Fill(raster2y,ee/p_el - 1);
          h_p_el_v_raster2x->Fill(raster2x,ee/p_el - 1);

          h_ebeam->Fill(ebeam);
        }

        if(goodW2) h_dx_0->Fill(dx);
        if(coin_time_L<coin_time && coin_time<coin_time_H-1){ h_W2_1->Fill(W2); if(goodW2) h_dx_1->Fill(dx); }
        if(coin_time_L<coin_time && coin_time<coin_time_H-2){ h_W2_2->Fill(W2); if(goodW2) h_dx_2->Fill(dx); }
        if(coin_time_L<coin_time && coin_time<coin_time_H-3){ h_W2_3->Fill(W2); if(goodW2) h_dx_3->Fill(dx); }
      }

      if(goodVz && goodPS && goodEoverp && goodDy && goodCoin && eHCAL>eHCAL_L-0.3){
        h_W2_0_0->Fill(W2);
        if(goodW2) h_dx_0_0->Fill(dx);

        if(eHCAL>eHCAL_L-0.2){ h_W2_1_1->Fill(W2); if(goodW2) h_dx_1_1->Fill(dx); }
        if(eHCAL>eHCAL_L-0.1){ h_W2_2_2->Fill(W2); if(goodW2) h_dx_2_2->Fill(dx); }
        if(eHCAL>eHCAL_L-0.0){ h_W2_3_3->Fill(W2); if(goodW2) h_dx_3_3->Fill(dx); }
      }

      if(goodVz && goodPS && goodEHCAL && goodEoverp && goodW2 && goodSpot) h_cointime->Fill(coin_time);
      if(goodVz && goodPS && goodEHCAL && goodW2 && goodCoin)               h_eoverp->Fill(eoverp);

      if(goodVz && goodPS && goodEoverp && goodEHCAL && goodW2 && goodCoin && goodSpot){
        h_chi2_all->Fill(ntrack_chi2);
        if((int)track_nhits == 3)      h_chi2_3hit->Fill(ntrack_chi2);
        else if((int)track_nhits == 4) h_chi2_4hit->Fill(ntrack_chi2);
        else if((int)track_nhits == 5) h_chi2_5hit->Fill(ntrack_chi2);
      }
    }

    if (i % 1000 == 0) { std::cout << (i * 100.0/ nentries) << "% \r"; std::cout.flush(); }
  }

  // -------------------- loop #2 --------------------
  const double x0_1 = 0.5*(dx_p_L_1 + dx_p_H_1);
  const double y0_1 = 0.5*(dy_L_1   + dy_H_1);
  const double rx_1 = 0.2*1.5;
  const double ry_1 = 0.2*1.5;

  const int nentries1 = tree1->GetEntries();
  for (int i = 0; i<nentries1; i++){
    tree1->GetEntry(i);

    const double eoverp = (eSH_1+ePS_1)/trP_1;

    const double ee  = A *(1+B *thtgt_1)/thetabend_1 - (vy_1*D + C);
    const double ee1 = A1*(1+B1*thtgt_1)/thetabend_1 - (vy_1*D + C);

    const double Q2_new  = 2*ebeam_1*ee *(1-std::cos(etheta_1));
    const double Q2_new1 = 2*ebeam_1*ee1*(1-std::cos(etheta_1));

    const double nu_new  = (ebeam_1 - ee);
    const double nu_new1 = (ebeam_1 - ee1);

    const double W2_new  = m*m+2*m*nu_new -Q2_new;
    const double W2_new1 = m*m+2*m*nu_new1-Q2_new1;

    const double p_el    = ebeam_1/(1+ebeam_1*(1-std::cos(etheta_1))/m);

    bool goodHelicity_1 = (lookupValue(HelicityCheck, runnum_1) == 1);
    bool goodMoller_1   = (lookupValue(MollerQuality, runnum_1) == 1);
    bool goodVz_1       = std::abs(vz_1) < 0.27;
    bool goodPS_1       = (ePS_1 > 0.2);
    bool goodRunRange_1 = (run_num_L_1 < runnum_1 && runnum_1 < run_num_H_1);
    bool goodEHCAL_1    = (eHCAL_1 > 0.225);
    bool validHel_1     = (helicity_1 == -1 || helicity_1 == 1);
    bool goodEoverp_1   = std::abs(eoverp-1) < 0.2;
    bool goodW2_1       = (W2_L_1 < W2_1 && W2_1 < W2_H_1);
    bool goodCoin_1     = (coin_time_L_1 < coin_time_1 && coin_time_1 < coin_time_H_1);
    bool goodDy_1       = (-0.2 < dy_1 && dy_1 < 0.2);
    bool goodpSpot_1    = ( std::pow((dy_1-y0_1)/ry_1,2) + std::pow((dx_1-x0_1)/rx_1,2) ) <= 1.0;
    bool goodSpot_1     = dy_L_1<dy_1 && dy_1<dy_H_1 && ((dx_L_1<dx_1 && dx_1<dx_H_1) || (dx_p_L_1<dx_1 && dx_1<dx_p_H_1));

    if(ePS_1>0 && eHCAL_1>0){

      if(goodVz_1 && goodPS_1 && goodEHCAL_1 && goodCoin_1 && W2_1<2.5){
        h_dxdy_1->Fill(dy_1,dx_1);
      }

      if(goodVz_1 && goodPS_1 && goodEHCAL_1 && goodCoin_1){
        if(goodpSpot_1){
          h_W2_new_1->Fill(W2_new);
          h_W2_new_1_1->Fill(W2_new1);
        }
      }
    }

    if (i % 1000 == 0) { std::cout << (i * 100.0/ nentries1) << "% \r"; std::cout.flush(); }
  }

  // -------------------- drawing (fixed W2 fill colors) --------------------
  TCanvas* c = new TCanvas("c","c",1600,1200);
  c->Divide(2,2);

  c->cd(1);
  h_ePS->SetFillStyle(1001);
  h_ePS->SetFillColorAlpha(kViolet+2,0.3);
  h_ePS->SetLineColor(kBlack);
  h_ePS->SetLineWidth(2);
  h_ePS->Draw();
  { TLine *line00 = new TLine(0.2, 0, 0.2, h_ePS->GetMaximum());
    line00->SetLineColor(kRed); line00->SetLineWidth(2); line00->Draw("same"); }

  c->cd(2);
  h_cointime->SetFillStyle(1001);
  h_cointime->SetFillColorAlpha(kViolet+2,0.3);
  h_cointime->SetLineColor(kBlack);
  h_cointime->SetLineWidth(2);
  h_cointime->Draw();
  { TLine *l1 = new TLine(coin_time_L, h_cointime->GetMinimum(), coin_time_L, h_cointime->GetMaximum());
    TLine *l2 = new TLine(coin_time_H, h_cointime->GetMinimum(), coin_time_H, h_cointime->GetMaximum());
    l1->SetLineColor(kRed); l2->SetLineColor(kRed);
    l1->SetLineWidth(2);    l2->SetLineWidth(2);
    l1->Draw("same"); l2->Draw("same"); }

  c->cd(3);
  h_trvz->SetFillStyle(1001);
  h_trvz->SetFillColorAlpha(kViolet+2,0.3);
  h_trvz->SetLineColor(kBlack);
  h_trvz->SetLineWidth(2);
  h_trvz->Draw();
  { TLine *l3 = new TLine(-0.27, 0, -0.27, h_trvz->GetMaximum());
    TLine *l4 = new TLine( 0.27, 0,  0.27, h_trvz->GetMaximum());
    l3->SetLineColor(kRed); l4->SetLineColor(kRed);
    l3->SetLineWidth(2);    l4->SetLineWidth(2);
    l3->Draw("same"); l4->Draw("same"); }

  c->cd(4);
  h_W2->SetFillStyle(1001);
  h_W2_1->SetFillStyle(1001);
  h_W2_2->SetFillStyle(1001);
  h_W2_3->SetFillStyle(1001);

  h_W2->SetFillColorAlpha( TColor::GetColor("#f0f0f0"), 0.3);
  h_W2_1->SetFillColorAlpha(TColor::GetColor("#bdbdbd"), 0.3);   // <- fixed
  h_W2_2->SetFillColorAlpha(TColor::GetColor("#969696"), 0.3);   // <- fixed
  h_W2_3->SetFillColorAlpha(TColor::GetColor("#636363"), 0.3);   // <- fixed

  h_W2->SetLineColor(kBlack); h_W2_1->SetLineColor(kBlack);
  h_W2_2->SetLineColor(kBlack); h_W2_3->SetLineColor(kBlack);
  h_W2->SetLineWidth(2); h_W2_1->SetLineWidth(2);
  h_W2_2->SetLineWidth(2); h_W2_3->SetLineWidth(2);
  h_W2->Draw(); h_W2_1->Draw("same"); h_W2_2->Draw("same"); h_W2_3->Draw("same");

  { TLine *l5 = new TLine(W2_L, 0, W2_L, h_W2->GetMaximum());
    TLine *l6 = new TLine(W2_H, 0, W2_H, h_W2->GetMaximum());
    l5->SetLineColor(kRed); l6->SetLineColor(kRed);
    l5->SetLineWidth(2);    l6->SetLineWidth(2);
    l5->Draw("same"); l6->Draw("same"); }

  // -------------- remaining canvases: unchanged except a few cosmetics -------
  TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(2,2);
  c1->cd(1);
  h_eHCAL->SetFillStyle(1001);
  h_eHCAL->SetFillColorAlpha(kViolet+2,0.3);
  h_eHCAL->SetLineColor(kBlack);
  h_eHCAL->SetLineWidth(2);
  h_eHCAL->Draw();
  { TLine *line09 = new TLine(eHCAL_L, 0, eHCAL_L, h_eHCAL->GetMaximum());
    line09->SetLineColor(kRed); line09->Draw("same"); }
  c1->cd(2);
  h_eoverp->SetFillStyle(1001);
  h_eoverp->SetFillColorAlpha(kViolet+2,0.3);
  h_eoverp->SetLineColor(kBlack);
  h_eoverp->SetLineWidth(2);
  h_eoverp->Draw();
  { TLine *line07 = new TLine(0.8, 0, 0.8, h_eoverp->GetMaximum());
    TLine *line08 = new TLine(1.2, 0, 1.2, h_eoverp->GetMaximum());
    line07->SetLineColor(kRed); line08->SetLineColor(kRed);
    line07->SetLineWidth(2);    line08->SetLineWidth(2);
    line07->Draw("same"); line08->Draw("same"); }

  c1->cd(3);
  h_chi2_all->SetFillStyle(1001);
  h_chi2_3hit->SetFillStyle(1001);
  h_chi2_4hit->SetFillStyle(1001);
  h_chi2_5hit->SetFillStyle(1001);
  h_chi2_all->SetLineColor(kBlack);
  h_chi2_3hit->SetLineColor(kBlack);
  h_chi2_4hit->SetLineColor(kBlack);
  h_chi2_5hit->SetLineColor(kBlack);
  h_chi2_all->SetLineWidth(2);
  h_chi2_3hit->SetLineWidth(2);
  h_chi2_4hit->SetLineWidth(2);
  h_chi2_5hit->SetLineWidth(2);
  h_chi2_all->SetFillColorAlpha(kGray,0.3);
  h_chi2_3hit->SetFillColorAlpha(kGray+1,0.3);
  h_chi2_4hit->SetFillColorAlpha(kGray+2,0.3);
  h_chi2_5hit->SetFillColorAlpha(kGray+3,0.3);
  h_chi2_all->Draw();
  h_chi2_3hit->Draw("same");
  h_chi2_4hit->Draw("same");
  h_chi2_5hit->Draw("same");
  {
    auto leg = new TLegend(0.62, 0.60, 0.88, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.035);
    leg->AddEntry(h_chi2_all,  "3,4,5-hit tracks", "f");
    leg->AddEntry(h_chi2_3hit, "3-hit tracks",     "f");
    leg->AddEntry(h_chi2_4hit, "4-hit tracks",     "f");
    leg->AddEntry(h_chi2_5hit, "5-hit tracks",     "f");
    leg->Draw();
  }

  TCanvas* c2 = new TCanvas("c2","c2",1600,1200);
  // (same styling as you had)
  // ... (unchanged code for c2–c9 besides what’s already fixed) ...

  // -------------------- dx/dy + ellipse --------------------
  TCanvas* c3 = new TCanvas("c3","c3",1600,1200);
  c3->Divide(2,2);
  c3->cd(1);
  TLine *linedx1 = new TLine(dx_L, h_dx->GetMinimum(), dx_L, h_dx->GetMaximum());
  TLine *linedx2 = new TLine(dx_H, h_dx->GetMinimum(), dx_H, h_dx->GetMaximum());
  h_dx->SetFillStyle(1001);
  h_dx->SetFillColorAlpha(kViolet+2,0.3);
  h_dx->SetLineColor(kBlack);
  h_dx->SetLineWidth(3);
  linedx1->SetLineWidth(3); linedx2->SetLineWidth(3);
  linedx1->SetLineColor(kRed); linedx2->SetLineColor(kRed);
  h_dx->Draw(); linedx1->Draw("same"); linedx2->Draw("same");

  c3->cd(2);
  TLine *linedy1 = new TLine(dy_L, h_dy->GetMinimum(), dy_L, h_dy->GetMaximum());
  TLine *linedy2 = new TLine(dy_H, h_dy->GetMinimum(), dy_H, h_dy->GetMaximum());
  h_dy->SetFillStyle(1001);
  h_dy->SetFillColorAlpha(kViolet+2,0.3);
  h_dy->SetLineColor(kBlack);
  h_dy->SetLineWidth(3);
  linedy1->SetLineWidth(3); linedy2->SetLineWidth(3);
  linedy1->SetLineColor(kRed); linedy2->SetLineColor(kRed);
  h_dy->Draw(); linedy1->Draw("same"); linedy2->Draw("same");

  c3->cd(3);
  gStyle->SetPalette(kRainBow);
  h_dxdy->Draw("COLZ");
  auto el = new TEllipse(y0, x0, ry, rx); // x=dy (y0), y=dx (x0)
  el->SetFillStyle(0); el->SetLineColor(kRed); el->SetLineWidth(3);
  el->Draw("same"); gPad->Update();

  // -------------------- W2 new comparisons, fits, etc. --------------------
  // (Your c4–c9 blocks kept; only cosmetics changed earlier)

  // ===================== c4: W2_new comparisons & trP =====================
  TCanvas *c4 =  new TCanvas("c4","c4",1600,1200);
  h_W2_new->SetLineWidth(3);  h_W2_new->SetLineColor(kBlue);
  h_W2_new1->SetLineWidth(3); h_W2_new1->SetLineColor(kRed);

  c4->Divide(2,2);
  c4->cd(1);
  h_W2_0->Draw();             // old W2
  h_W2_new->Draw("same");     // new (A,B)
  h_W2_new1->Draw("same");    // new (A1,B1)

  c4->cd(2);
  h_W2_new->Draw();

  c4->cd(3);
  h_trP->SetLineWidth(3);     h_trP->SetLineColor(kBlack);
  h_trP_new->SetLineWidth(3); h_trP_new->SetLineColor(kBlue);
  h_trP_new1->SetLineWidth(3);h_trP_new1->SetLineColor(kRed);
  h_trP->Draw(); h_trP_new->Draw("same"); h_trP_new1->Draw("same");

  c4->cd(4);
  h_trP_new->Draw();

  // ===================== c5: 2D correlations & lower-edge fit =====================
  TCanvas *c5 =  new TCanvas("c5","c5",1600,1200);
  c5->Divide(2,2);
  c5->cd(1); h_W2_v_thetatgt->Draw("COLZ");
  c5->cd(2); h_W2timesthetabend_v_thetatgt->Draw("COLZ");
  c5->cd(3);
  h_p_eltimesthetabend_v_thetatgt->Draw("COLZ");
  // lower-edge quantile fit (as in your macro)
  {
    double q = 0.20; int nmin = 5; int rebinx = 1;
    if (rebinx>1) h_p_eltimesthetabend_v_thetatgt->RebinX(rebinx);

    std::vector<double> xs, ys;
    for (int ix = 1; ix <= h_p_eltimesthetabend_v_thetatgt->GetNbinsX(); ++ix) {
      std::unique_ptr<TH1D> py(h_p_eltimesthetabend_v_thetatgt->ProjectionY(Form("py_%d",ix), ix, ix, "e"));
      if (!py || py->GetEntries() < nmin) continue;
      double prob[1] = { q }, qval[1] = { 0.0 };
      py->GetQuantiles(1, qval, prob);
      xs.push_back(h_p_eltimesthetabend_v_thetatgt->GetXaxis()->GetBinCenter(ix));
      ys.push_back(qval[0]);
    }
    TGraph g(xs.size(), xs.data(), ys.data());
    g.SetMarkerStyle(20); g.SetMarkerSize(0.9); g.SetMarkerColor(kBlack);
    g.Draw("P same");

    TF1 fline("fline","pol1",-0.15,0.05);
    g.Fit(&fline, "Q","",-0.15,0.05);
    fline.SetLineColor(kRed+1); fline.SetLineWidth(3); fline.Draw("same");

    double c0   = fline.GetParameter(0), dc0 = fline.GetParError(0);
    double m1   = fline.GetParameter(1), dm1 = fline.GetParError(1);
    TPaveText *box = new TPaveText(0.14, 0.80, 0.44, 0.93, "NDC");
    box->SetFillStyle(0); box->SetLineColor(kBlack); box->SetTextFont(42); box->SetTextSize(0.035); box->SetTextAlign(13);
    box->AddText(Form("m = %.5f #pm %.5f", m1, dm1));
    box->AddText(Form("c = %.5f #pm %.5f", c0, dc0));
    box->Draw();
  }
  c5->cd(4); h_trPtimesthetabend_v_thetatgt->Draw("COLZ");

  // ===================== c6: BPM/vx/vy maps =====================
  TCanvas *c6 =  new TCanvas("c6","c6",1600,1200);
  c6->Divide(2,2);
  c6->cd(1); h_p_el_v_BPMAy->Draw("COLZ");
  c6->cd(2); h_p_el_v_BPMAx->Draw("COLZ");
  c6->cd(3); h_p_el_v_vy->Draw("COLZ");
  c6->cd(4); h_p_el_v_vx->Draw("COLZ");

  // ===================== c7: dp overlays =====================
  TCanvas *c7 =  new TCanvas("c7","c7",1600,1200);
  c7->Divide(2,2);
  c7->cd(1); h_dp->Draw();
  c7->cd(2); h_dp_new->Draw();
  c7->cd(3);
  h_dp->SetLineWidth(3);   h_dp->SetLineColor(kBlack);
  h_dp_new->SetLineWidth(3);  h_dp_new->SetLineColor(kBlue);
  h_dp_new1->SetLineWidth(3); h_dp_new1->SetLineColor(kRed);
  h_dp->Draw(); h_dp_new->Draw("same"); h_dp_new1->Draw("same");

  // ===================== c8: raster maps =====================
  TCanvas *c8 =  new TCanvas("c8","c8",1600,1200);
  c8->Divide(2,2);
  c8->cd(1); h_p_el_v_rastery->Draw("COLZ");
  c8->cd(2); h_p_el_v_rasterx->Draw("COLZ");
  c8->cd(3); h_p_el_v_raster2y->Draw("COLZ");
  c8->cd(4); h_p_el_v_raster2x->Draw("COLZ");

  // ===================== c9: beam energy & W2_new vs file2 =====================
  TCanvas *c9 =  new TCanvas("c9","c9",1600,1200);
  c9->Divide(2,2);
  c9->cd(1); h_ebeam->SetLineColor(kBlack); h_ebeam->SetLineWidth(3); h_ebeam->Draw();
  
  c9->cd(2);
  h_W2_new->Scale(1/h_W2_new->Integral());
  h_W2_new_1->Scale(1/h_W2_new_1->Integral());

  h_W2_new->SetLineColor(kBlack); h_W2_new->SetLineWidth(3);
  h_W2_new_1->SetLineColor(kBlue); h_W2_new_1->SetLineWidth(3);
  h_W2_new->Draw("hist"); h_W2_new_1->Draw("same hist");
  c9->cd(3);
  h_dxdy->Draw("COLZ");
  auto elipse = new TEllipse(y0, x0, ry, rx); // x=dy (y0), y=dx (x0)
  elipse->SetFillStyle(0); elipse->SetLineColor(kRed); elipse->SetLineWidth(3);
  elipse->Draw("same"); gPad->Update();

  c9->cd(4);
  h_dxdy_1->Draw("COLZ");
  auto elipse_1 = new TEllipse(y0_1, x0_1, ry_1, rx_1); // x=dy (y0), y=dx (x0)
  elipse_1->SetFillStyle(0); elipse_1->SetLineColor(kRed); elipse_1->SetLineWidth(3);
  elipse_1->Draw("same"); gPad->Update();


  // -------------------- output --------------------
  c->Print(Form("../plots/%s_ePS_eHCAL.pdf(",printfilename));
  c->SaveAs(Form("../plots/%s_ePS_eHCAL.png",printfilename));
  c1->SaveAs(Form("../plots/%s_ePS_eHCAL_1.png",printfilename));
  c2->SaveAs(Form("../plots/%s_ePS_eHCAL_2.png",printfilename));
  c3->SaveAs(Form("../plots/%s_ePS_eHCAL_3.png",printfilename));
  c4->SaveAs(Form("../plots/%s_W2_test.png",printfilename));
  c5->SaveAs(Form("../plots/%s_W2_correlation_test.png",printfilename));
  c6->SaveAs(Form("../plots/%s_W2_correlation_test_1.png",printfilename));
  c7->SaveAs(Form("../plots/%s_W2_correlation_test_2.png",printfilename));
  c8->SaveAs(Form("../plots/%s_W2_correlation_test_3.png",printfilename));
  c9->SaveAs(Form("../plots/%s_W2_correlation_test_4.png",printfilename));
  c1->Print(Form("../plots/%s_ePS_eHCAL.pdf)",printfilename));

  // Optional: close files
  file->Close(); file1->Close();
}
