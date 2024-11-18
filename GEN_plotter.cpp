#include "models.h"
#include "parse.h"
#include "plotdxdy.h"
#include "FF_param.cpp"
#include "Constants.h"
//#include "ExpConstants.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <TDatime.h>

vector<World_Data> data_points;

//const int npoints_new = 3;
const int npoints_new = 1;

//double Q2_new[npoints_new] = {2.90,6.62,9.48};
double Q2_new[npoints_new] = {6.6};

//double GEGM_new[npoints_new] = {-0.1952,-0.3527,-0.5127};
//double GEGM_stat_err_new[npoints_new] = {0.0194,0.1026,0.2629};
//double GEGM_sys_err_new[npoints_new] = {0.0130,0.0262,0.0552};
double GEGM_new[npoints_new] = {-0.34469};
double GEGM_stat_err_new[npoints_new] = {0.090659};
double GEGM_sys_err_new[npoints_new] = {0.0214585};

double GEGM_err_new[npoints_new];

//double GE_new[npoints_new] = {0.0157,0.0067,0.0046};
//double GE_stat_err_new[npoints_new] = {0.0016,0.0019,0.0023};
//double GE_sys_err_new[npoints_new] = {0.0011,0.0005,0.0005};
double GE_new[npoints_new] = {0.0065146788};
double GE_stat_err_new[npoints_new] = {0.001714851292};
double GE_sys_err_new[npoints_new] = {0.0004055680032};

double GE_err_new[npoints_new];

double kappa_d = -2.03;
double kappa_u = 1.67;


Double_t PSM_theory(Double_t *x,Double_t *par){
  return (0.4133*x[0] - 0.1038*pow(x[0],2) + 0.01252*pow(x[0],3) + 0.007748*pow(x[0],4))/(1 +  0.9771*x[0] - 0.3506*pow(x[0],2) + 0.09863*pow(x[0],3) + 0.009927*pow(x[0],4) + 0.000002173*pow(x[0],5));
}

double F1F2_to_GEGM(double F1, double F2, double Q2, bool is_GEGM){

  double tau = Q2 / (4 * constant::Mn * constant::Mn);
  double GE = F1 - tau*F2;
  double GM = F1 + F2;
  double GD = pow( 1 + Q2 / 0.71,-2);

  double val;
  if(is_GEGM) val = constant::mun*GE / GM;
  else val = GE * GD;

  return val;
}


TGraph *make_graph(TString var_name){
  
  TString DB_dir = "DB/";
  fstream file; file.open(DB_dir + var_name);

  vector<double> Q2, FF;

  string line;
  int iline = 0;
      
  while (getline(file, line)) {
    // Create a stringstream to parse the line
    stringstream ss(line);
    string cell;
    iline++;
	
    vector<string> val;

    // Split the line into cells using a comma as a delimiter
    while (getline(ss, cell, ',')) {
      val.push_back(cell);  // Put one line into vectros
    }
    if(val[0].substr(0, 1) == "#") continue;
    
    Q2.push_back(stod(val[0]));
    FF.push_back(stod(val[1]));
   
  }

  TGraph *g = new TGraph(Q2.size(),&Q2[0],&FF[0]);
  
  return g;
}


TGraph *RCQM_theory(bool is_GEGM){

  vector<double> Q2_plot;
  double Q2_max = 12;
  int nQ2 = 30;
  double Q2_bin_size = Q2_max / nQ2;


  Q2_plot.push_back(0);
  for(int iQ2 = 0; iQ2 <= nQ2; iQ2++)
    Q2_plot.push_back(iQ2*Q2_bin_size + 0.15);
  
  TGraph *g_F1u = make_graph("F1u.csv");
  TGraph *g_F1d = make_graph("F1d.csv");
  TGraph *g_F2u = make_graph("F2u.csv");
  TGraph *g_F2d = make_graph("F2d.csv");

  vector<double> GEGM_RCQM;

  for(int iQ2 = 0; iQ2 < Q2_plot.size(); iQ2++){
    double Q2 = Q2_plot[iQ2];

    double F1u = g_F1u->Eval(Q2) / (pow(Q2,4));
    double F1d = g_F1d->Eval(Q2) / (pow(Q2,4));
    double F2u = g_F2u->Eval(Q2) * kappa_u / (pow(Q2,4));
    double F2d = g_F2d->Eval(Q2) * kappa_d / (pow(Q2,4));
    
    double F1 = 2.0 / 3.0 * F1d - 1.0 / 3.0 * F1u;
    double F2 = 2.0 / 3.0 * F2d - 1.0 / 3.0 * F2u;
    
    if(Q2 > 0.0001) GEGM_RCQM.push_back(F1F2_to_GEGM(F1, F2, Q2, is_GEGM));
    else GEGM_RCQM.push_back(0);
    
  }

  TGraph *g = new TGraph(Q2_plot.size(),&Q2_plot[0],&GEGM_RCQM[0]);

  return g;
}


TGraph *Global_fit(bool is_GEGM){

  vector<double> Q2_plot, plot_val;
  int nQ2 = 100;
  double Q2max = 12;
  double Q2min = 1e-3;
  double Q2step = (Q2max - Q2min) / nQ2;

  for(int i=0; i < nQ2;i++){
    Q2_plot.push_back(Q2min + Q2step*i);
   
    double GE_temp, GM_temp,G_err_temp;
    GetGENFit(Q2_plot[i],GE_temp,G_err_temp);
    GetGMNFit(Q2_plot[i],GM_temp,G_err_temp);

    double val;
    if(is_GEGM) val = constant::mun*GE_temp / GM_temp;
    else val = GE_temp;
    
    plot_val.push_back(val);
  }

  TGraph *g = new TGraph(Q2_plot.size(),&Q2_plot[0],&plot_val[0]);
  return g;
}


void LoadWorldData(){
  
  ifstream file("DB/World_Data.dat");
  if (!file.is_open()) {
    cout << "Error opening file!" << endl;
  }

  string name_old = "first";
  vector<double> Q2_data;
  vector<double> GEGM_data;
  vector<double> GEGM_err_data;
  vector<double> GEn_data;
  vector<double> GEn_err_data;
  
  // Read the file line by line
  string line;
  int iline = 0;
  while (getline(file, line)) {
    iline++;
    if(iline == 1) continue; // skip the first line
    
    // Use stringstream to parse each line
    stringstream ss(line);
    string name;
    double Q2_temp;
    double GEGM_temp;
    double GEGM_err_temp;
    double GEn_temp;
    double GEn_err_temp;
    
    // Extract data from each column
    ss >> Q2_temp >> name >> GEGM_temp >> GEGM_err_temp >> GEn_temp >> GEn_err_temp;

    if((name != name_old && name_old != "first") || (file.eof())){
      if(name_old == "first") name_old = name;
      
      World_Data one_exp;
      one_exp.name = name_old;
      one_exp.Q2 = Q2_data;
      one_exp.GEGM = GEGM_data;
      one_exp.GEGM_err = GEGM_err_data;
      one_exp.GEn = GEn_data;
      one_exp.GEn_err = GEn_err_data;

      data_points.push_back(one_exp);
      
      Q2_data.clear();
      GEGM_data.clear();
      GEGM_err_data.clear();
      GEn_data.clear();
      GEn_err_data.clear();
    }

    Q2_data.push_back(Q2_temp);
    GEGM_data.push_back(GEGM_temp);
    GEGM_err_data.push_back(GEGM_err_temp);
    GEn_data.push_back(GEn_temp);
    GEn_err_data.push_back(GEn_err_temp);
    
    name_old = name;

  }

  // Close the file
  file.close();

}


void GEN_plotter(){

  gStyle->SetCanvasPreferGL();
  
  TF1 *PSM_theory_curve = new TF1("PSM_theory_curve", PSM_theory, 0, 12,1);

  TGraph *gGEGM_RCQM = RCQM_theory(1);
  TGraph *gGE_RCQM = RCQM_theory(0);
  TGraph *gGEGM_global = Global_fit(1);
  TGraph *gGE_global = Global_fit(0);
  LoadWorldData();

  
  for(int i=0; i < npoints_new; i++){
    GEGM_new[i] *= constant::mun;
    GEGM_err_new[i] = sqrt( GEGM_stat_err_new[i]*GEGM_stat_err_new[i] + GEGM_sys_err_new[i]*GEGM_sys_err_new[i] );
    GEGM_err_new[i] *= abs(constant::mun);
    GE_err_new[i] = sqrt( GE_stat_err_new[i]*GE_stat_err_new[i] + GE_sys_err_new[i]*GE_sys_err_new[i] );
  }

  
  TGraphErrors *g_new = new TGraphErrors(npoints_new,Q2_new,GEGM_new,0,GEGM_err_new);
  g_new->SetMarkerStyle(20);
  g_new->SetMarkerSize(1.5);
  g_new->SetMarkerColor(kRed);  
  g_new->SetLineColor(kRed);  
  
  TCanvas *c = new TCanvas("c","",800,600);
  g_new->SetTitle("G_{E}/G_{M} Neutron Results;Q^{2} (GeV^{2});#mu_{n}G_{E}^{n}/G_{M}^{n}");
  g_new->Draw("AP");


  TLegend *legend1 = new TLegend(0.11,0.64,0.41,0.89);
  legend1->SetLineColor(0);
  
  TLegend *legend3 = new TLegend(0.292,0.641,0.593,0.891);
  legend3->SetLineColor(0);
    
  int icolor = 1;
  int istyle = 21;
  for(int idata=0;idata < data_points.size(); idata++){
    World_Data exp = data_points[idata];

    if(icolor == 2 || icolor == 5) icolor++;  //Skip red and yellow
    
    TGraphErrors *g = new TGraphErrors(exp.Q2.size(),&exp.Q2[0],&exp.GEGM[0],0,&exp.GEGM_err[0]);
    g->SetMarkerColor(icolor);
    g->SetLineColor(icolor);
    g->SetMarkerStyle(istyle);
    g->SetMarkerSize(1.5);
    g->Draw("P");
    
    legend1->AddEntry(g,exp.name,"p");
    legend3->AddEntry(g,exp.name,"p");
    
    icolor++;
    istyle++;
  }
  
  PSM_theory_curve->Draw("same");
  gGEGM_RCQM->Draw("same");
  gGEGM_global->Draw("same");

  gGEGM_RCQM->SetLineColor(kBlue);
  gGEGM_RCQM->SetLineWidth(2);
  gGEGM_RCQM->SetLineStyle(2);

  PSM_theory_curve->SetLineStyle(9);
  
  gGEGM_global->SetLineWidth(2);

  g_new->GetXaxis()->SetLimits(0,12);
  g_new->GetYaxis()->SetRangeUser(-0.5,1.5);

  legend3->AddEntry(g_new,"This Work","p");
  legend1->AddEntry(g_new,"This Work","p");
  legend1->Draw("same");

  TLegend *legend2 = new TLegend(0.395,0.736,0.709,0.891);
  legend2->AddEntry(gGEGM_RCQM,"RCQM - Miller","l");
  legend2->AddEntry("PSM_theory_curve","DSE - Roberts","l");
  legend2->AddEntry(gGEGM_global,"Global Fit - Ye","l");
  legend2->SetLineColor(0);
  legend2->Draw("same");

  TLegend *legend4 = new TLegend(0.585,0.734,0.898,0.891);
  legend4->AddEntry(gGEGM_global,"Global Fit - Ye","l");
  legend4->SetLineColor(0);

  TPaveText *pt = new TPaveText(.139,.260,.860,.627,"ndc");
  pt->AddText("Exploratory");
  pt->SetFillStyle(0);
  pt->SetLineColorAlpha(kWhite,0);
  pt->SetTextColorAlpha(kRed,0.3);
  pt->Draw("same");

  g_new->GetXaxis()->SetLimits(0,12);
  g_new->GetYaxis()->SetRangeUser(0,1.5);

  g_new->GetXaxis()->SetTitleSize(0.05);
  g_new->GetXaxis()->SetTitleOffset(0.8);

  g_new->GetYaxis()->SetTitleSize(0.05);
  g_new->GetYaxis()->SetTitleOffset(0.8);

  TGraphErrors *gGE_new = new TGraphErrors(npoints_new,Q2_new,GE_new,0,GE_err_new);
  gGE_new->SetMarkerStyle(20);
  gGE_new->SetMarkerSize(1.5);
  gGE_new->SetMarkerColor(kRed);  
  gGE_new->SetLineColor(kRed); 


  TCanvas *c2 = new TCanvas("c2","",800,600);
  gGE_new->SetTitle("G_{E} Neutron Results;Q^{2} (GeV^{2});G_{E}^{n}");
  gGE_new->Draw("AP");

  icolor = 1;
  istyle = 21;
  for(int idata=0;idata < data_points.size(); idata++){
    World_Data exp = data_points[idata];

    if(icolor == 2 || icolor == 5) icolor++;  //Skip red and yellow
    
    TGraphErrors *g = new TGraphErrors(exp.Q2.size(),&exp.Q2[0],&exp.GEn[0],0,&exp.GEn_err[0]);
    g->SetMarkerColor(icolor);
    g->SetLineColor(icolor);
    g->SetMarkerStyle(istyle);
    g->SetMarkerSize(1.5);
    g->Draw("P");

    icolor++;
    istyle++;
  }

  gGE_new->GetXaxis()->SetTitleSize(0.05);
  gGE_new->GetXaxis()->SetTitleOffset(0.8);

  gGE_new->GetYaxis()->SetTitleSize(0.05);
  gGE_new->GetYaxis()->SetTitleOffset(0.95);
  
  gGE_global->Draw("same");
  legend3->Draw("same");
  legend4->Draw("same");
  pt->Draw("same");

  gGE_global->SetLineWidth(2);
  
  gGE_new->GetXaxis()->SetLimits(0,12);
  gGE_new->GetYaxis()->SetRangeUser(0,0.08);
}