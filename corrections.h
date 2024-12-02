#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>

#pragma once

#include "parse.h"

//cut variables
class cuts{
public:
  double Aacc = 0.0;
  double Api = 0.0;
  double Ain = 0.0;
  double Ap = 0.0;
  double Afsi = 0.0;

  double facc = 0.0;
  double fN2 = 0.0;
  double fpi = 0.0;
  double fin = 0.0;
  double fp = 0.0;
  double ffsi = 0.0;

//function to get the assign cuts
  void parsecuts(const std::map<std::string, std::string> &config){

	Aacc = getDoubleValue(config,"Aacc");
	Api = getDoubleValue(config,"Api");	
	Ain = getDoubleValue(config,"Ain");
	Ap = getDoubleValue(config,"Ap");
	Afsi = getDoubleValue(config,"Afsi");
	
	facc = getDoubleValue(config,"facc");
	fN2 = getDoubleValue(config,"fN2");
	fpi = getDoubleValue(config,"fpi");
	fin = getDoubleValue(config,"fin");
	fp = getDoubleValue(config,"fp");
	ffsi = getDoubleValue(config,"ffsi");
	
	}
};
