#include "models.h"
#include "parse.h"
#include "plotdxdy.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <TDatime.h>



void physics_asymmetry(const char* kin, bool flag_eHCAL_cut){
	std::map<std::string, std::string> corr = parseConfig(Form("DB/corrections/%s_corrections_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str())); //read corrections from file
	std::map<int, RunData> rundata = readRunDataFromFile(Form("txt/%s_raw_neutron_asymmetry_results_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str())); //read the Aexp from file
    std::map<std::string, std::string> avgpoldata = parseConfig(Form("txt/%s_average_polarization_results_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));
    
    //read corrections
	double Aacc = getDoubleValue(corr,"Aacc");
	double Api = getDoubleValue(corr,"Api");	
	double Ain = getDoubleValue(corr,"Ain");
	double Ap = getDoubleValue(corr,"Ap");
	double Afsi = getDoubleValue(corr,"Afsi");
	
	double facc = getDoubleValue(corr,"facc");
	double fN2 = getDoubleValue(corr,"fN2");
	double fpi = getDoubleValue(corr,"fpi");
	double fin = getDoubleValue(corr,"fin");
	double fp = getDoubleValue(corr,"fp");
	double ffsi = getDoubleValue(corr,"ffsi");

    double errAacc = getDoubleValue(corr,"errAacc");
    double errApi = getDoubleValue(corr,"errApi");  
    double errAin = getDoubleValue(corr,"errAin");
    double errAp = getDoubleValue(corr,"errAp");
    double errAfsi = getDoubleValue(corr,"errAfsi");

    double errfacc = getDoubleValue(corr,"errfacc");
    double errfN2 = getDoubleValue(corr,"errfN2");
    double errfpi = getDoubleValue(corr,"errfpi");
    double errfin = getDoubleValue(corr,"errfin");
    double errfp = getDoubleValue(corr,"errfp");
    double errffsi = getDoubleValue(corr,"errffsi");

    //read average polarization data
    double avg_He3pol = getDoubleValue(avgpoldata,"avg_He3pol");
    double avg_beampol = getDoubleValue(avgpoldata,"avg_beampol");
    double avg_Pn = getDoubleValue(avgpoldata,"avg_Pn");
    double err_avg_He3pol = getDoubleValue(avgpoldata,"err_avg_He3pol");
    double err_avg_beampol = getDoubleValue(avgpoldata,"err_avg_beampol");
    double err_avg_Pn = getDoubleValue(avgpoldata,"err_avg_Pn");

	double pol_n = 0.96;

	double Aphy_sum_numerator = 0.0;
	double Aphy_sum_denominator = 0.0;
	double P_sum_numerator = 0.0;
    double P_sum_denominator = 0.0;


	const double epsilon = 1e-10; //to handle division by zero

	//write physics asymmetry values to files
	std::ofstream outfile;
	outfile.open(Form("txt/%s_physics_neutron_asymmetry_results_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));

    //write sum physics asymmetry values to files
    std::ofstream outfile_sum;
    outfile_sum.open(Form("txt/%s_physics_neutron_asymmetry_summary_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));

    for (const auto& entry : rundata) {    
        int runNumber = entry.first;
        const RunData& runData = entry.second;

        // Check for NaN and very small values, skip if any are problematic
        if (std::isnan(runData.Aexp) || std::isnan(runData.errAexp) || 
            std::isnan(runData.pol_He3) || std::isnan(runData.pol_beam) ||
            runData.pol_He3 < epsilon || runData.pol_beam < epsilon /*|| runData.errAexp == 0*/) {
            std::cerr << "Skipping run " << runNumber << " due to NaN or small value(s)." << std::endl;
            continue;
        }

        /*std::cout << "Run_number: " << runNumber
                  << " | Aexp: " << runData.Aexp
                  << " | errAexp: " << runData.errAexp
                  << " | pol_He3: " << runData.pol_He3
                  << " | pol_beam: " << runData.pol_beam << std::endl;*/

        // Calculate denominator and check if it is safe to use
        double denominator = runData.pol_He3 * 0.01 * runData.pol_beam * 0.01 * pol_n * (1 - facc - fN2 - fpi - fin - fp - ffsi);
        if (std::fabs(denominator) < epsilon || runData.Aexp == 0) {
            std::cerr << "Skipping run " << runNumber << " due to small denominator, potential division by zero." << std::endl;
            continue;
        }

        // Calculate Aphy and errAphy_stat (statistical error)

        double Aphy = (runData.Aexp - facc * Aacc - fpi * Api - fin * Ain - fp * Ap - ffsi * Afsi) / denominator;
        double errAphy_stat = runData.errAexp / denominator;
        

        outfile << "Run_number: " << runNumber << " | Aphy: " << Aphy << " | errAphy: " << errAphy_stat << std::endl;

        if (errAphy_stat>epsilon){ //weight statistically
        	Aphy_sum_numerator += (Aphy / (errAphy_stat * errAphy_stat));
        	Aphy_sum_denominator += (1 / (errAphy_stat * errAphy_stat)); 
            P_sum_numerator += (runData.pol_He3 * 0.01 * runData.pol_beam * 0.01 * pol_n / (errAphy_stat * errAphy_stat));
            P_sum_denominator += (1 / (errAphy_stat * errAphy_stat));
    	}
    }


    // Compute weighted average if denominator is non-zero
    if (Aphy_sum_denominator > epsilon) {
        double Aphy_sum = Aphy_sum_numerator / Aphy_sum_denominator;
        double errAphy_sum = 1/sqrt(Aphy_sum_denominator);
        double P = P_sum_numerator/P_sum_denominator;
        double fn = (1 - facc - fN2 - fpi - fin - fp - ffsi);
        //get the systematic error
        double errAphy_sys = sqrt( (facc*facc*errAacc*errAacc + fpi*fpi*errApi*errApi + fin*fin*errAin*errAin + fp*fp*errAp*errAp + ffsi*ffsi*errAfsi*errAfsi)/(P*P*fn*fn)
            + pow(Aphy_sum*errfN2/fn,2) + pow((P*Aphy_sum-Aacc)*errfacc/(P*fn),2) + pow((P*Aphy_sum-Api)*errfpi/(P*fn),2) 
            + pow((P*Aphy_sum-Ain)*errfin/(P*fn),2) + pow((P*Aphy_sum-Ap)*errfp/(P*fn),2) + pow((P*Aphy_sum-Afsi)*errffsi/(P*fn),2)
            + Aphy_sum*Aphy_sum*((err_avg_He3pol*err_avg_He3pol) + ((err_avg_Pn/pol_n)*(err_avg_Pn/pol_n)) + ((err_avg_beampol/avg_beampol)*(err_avg_beampol/avg_beampol))) //sigmaHe3pol is taken as 0.05*P, 0,05 is hard coded  
            );


        std::cout <<"Aphy_sum_numerator : " << Aphy_sum_numerator << " Aphy_sum_denominator : "<< Aphy_sum_denominator <<" Aphy_sum: " << Aphy_sum 
        <<" errAphy_stat_sum : "<< errAphy_sum <<" errAphy_sys : "<< errAphy_sys<< std::endl;

        outfile_sum<<"Aphy_sum = "<<Aphy_sum<<endl;
        outfile_sum<<"err_Aphy_stat_sum = "<<errAphy_sum<<endl;
        outfile_sum<<"err_Aphy_sys = "<<errAphy_sys<<endl; 

    } else {
        std::cerr << "Error: Aphy_sum_denominator is zero or too small, cannot calculate Aphy_sum." << std::endl;
    }

    outfile.close();

}	