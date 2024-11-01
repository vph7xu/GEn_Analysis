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
  double coin_time_mean = 0.0;
  double coin_time_L = 0.0;
  double coin_time_H = 0.0;
  double coin_time_width = 0.0;
  double coin_time_ac_L = 0.0;
  double coin_time_ac_H = 0.0;
  double W2_mean = 0.0;
  double W2_L = 0.0;
  double W2_H = 0.0;
  double W2_width = 0.0;
  double dy_L = 0.0;
  double dy_H = 0.0;
  double dx_L = 0.0;
  double dx_H = 0.0;
  double dy_p_L = 0.0;
  double dy_p_H = 0.0;
  double dx_p_L = 0.0;
  double dx_p_H = 0.0;
  double dy_ac_L = 0.0;
  double dy_ac_H = 0.0;
  double dx_C = 0.0;
  double dx_R = 0.0;
  double dy_C = 0.0;
  double dy_R = 0.0;
  double IHWP_flip = 0.0;
  double run_num_L = 0.0;
  double run_num_H = 0.0;
  double N_cut_cx = 0.0;
  double N_cut_cy = 0.0;
  double N_cut_rx = 0.0;
  double N_cut_ry = 0.0;

//function to get the assign cuts
  void parsecuts(const std::map<std::string, std::string> &config){

	coin_time_mean = getDoubleValue(config,"coin_time_mean");
  	coin_time_L = getDoubleValue(config,"coin_time_L");
	coin_time_H = getDoubleValue(config,"coin_time_H");
	coin_time_width = getDoubleValue(config,"coin_time_width");

	coin_time_ac_L = getDoubleValue(config,"coin_time_ac_L");
	coin_time_ac_H = getDoubleValue(config,"coin_time_ac_H");

	W2_mean = getDoubleValue(config,"W2_mean");
	W2_L = getDoubleValue(config,"W2_L");
	W2_H = getDoubleValue(config,"W2_H");
	W2_width = getDoubleValue(config,"W2_width");

        dy_L = getDoubleValue(config,"dy_L");
        dy_H = getDoubleValue(config,"dy_H");

        dx_L = getDoubleValue(config,"dx_L");
        dx_H = getDoubleValue(config,"dx_H");

        dy_p_L = getDoubleValue(config,"dy_p_L");
        dy_p_H = getDoubleValue(config,"dy_p_H");

        dx_p_L = getDoubleValue(config,"dx_p_L");
        dx_p_H = getDoubleValue(config,"dx_p_H");

        dy_ac_L = getDoubleValue(config,"dy_ac_L");
        dy_ac_H = getDoubleValue(config,"dy_ac_H");

        dx_C = getDoubleValue(config,"dx_C");
        dx_R = getDoubleValue(config,"dx_R");
        
        dy_C = getDoubleValue(config,"dy_C");
        dy_R = getDoubleValue(config,"dy_R");

        IHWP_flip = getDoubleValue(config,"IHWP_flip");

        run_num_L = getDoubleValue(config,"run_num_L");
        run_num_H = getDoubleValue(config,"run_num_H");

        N_cut_cx = getDoubleValue(config,"N_cut_cx");
        N_cut_cy = getDoubleValue(config,"N_cut_cy");
        N_cut_rx = getDoubleValue(config,"N_cut_rx");
        N_cut_ry = getDoubleValue(config,"N_cut_ry");

  }

};
