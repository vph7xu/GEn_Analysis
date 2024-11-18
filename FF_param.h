#ifndef FF_param_H
#define FF_param_H

#include "Constants.h"

struct World_Data{
  TString name;
  vector<double> Q2;
  vector<double> GEGM;
  vector<double> GEGM_err;
  vector<double> GEn;
  vector<double> GEn_err;

  void Print(){
    cout<<"Name: "<<name<<endl;
    cout<<"Q2 = ";
    for(int i=0; i < Q2.size(); i++)
      cout<<Q2[i]<<",";
    cout<<endl;
  }
};


void GetFitFF(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err);

double GetGD(double Q2){
  return pow(1./(1. + Q2/0.71), 2);
}

void GetGEPFit(double Q2,double &G,double &Gerr){
  GetFitFF(1, Q2, &G, &Gerr); // Result is GEp/GD
  G *= GetGD(Q2);  // Remove GD
  Gerr *= GetGD(Q2);  // Remove GD
}

void GetGMPFit(double Q2,double &G,double &Gerr){
  GetFitFF(2, Q2, &G, &Gerr); // Result is GMp/u_pGD
  G *= constant::mup * GetGD(Q2);  // Remove u_pGD
  Gerr *= constant::mup * GetGD(Q2);  // Remove u_pGD
}

void GetGENFit(double Q2,double &G,double &Gerr){
  GetFitFF(3, Q2, &G, &Gerr); // Result is GEn/GD
  G *= GetGD(Q2);  // Remove GD
  Gerr *= GetGD(Q2);  // Remove GD
}

void GetGMNFit(double Q2,double &G,double &Gerr){
  GetFitFF(2, Q2, &G, &Gerr); // Result is GMn/u_nGD
  G *= constant::mun * GetGD(Q2);  // Remove u_nGD
  Gerr *= constant::mun * GetGD(Q2);  // Remove u_nGD
}

#endif
