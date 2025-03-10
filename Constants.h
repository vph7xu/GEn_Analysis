#ifndef CONST_H
#define CONST_H

#include "TMath.h"

namespace constant {

  // math
  static const double pi = TMath::Pi();
  
  // electron
  static const double qe = 1.602176634E-19; // C
  static const double Me = 0.5109989461E-03; // +/- 31E-13 GeV

  // muon
  static const double Mmu = 0.1056583745; // +/- 24E-10 GeV

  // pi-
  static const double Mpim = 0.13957039; // +/- 24E-10 GeV
  
  // proton
  static const double Mp = 0.938272081; // +/- 6E-9 GeV
  static const double mup = 2.7928473446; // +/- 7E-11 muN

  // neutron
  static const double Mn = 0.939565413; // +/- 6E-9 GeV
  static const double mun = -1.9130427; // +/- 5E-7 muN

  // atomic 
  static const double N_A = 6.02214076E23; // 1/mol
  static const double H2_Amass = 1.00784; // u
  static const double D2_Amass = 2.013553212745; // u
  static const double He3_Amass = 3.01603; // u

  // Ideal Gas
  static const double R_ideal = 0.08314;  //Divided by 100 because we are using cm instead of m

}

#endif