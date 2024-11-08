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


void physics_asymmetry(const char* kin){
	std::map<std::string, std::string> config = parseConfig(Form("DB/corrections/%s_corrections.txt",kin)); //read corrections from file

	std::map<int, RunData> rundata = readRunDataFromFile(Form("%s_raw_neutron_asymmetry_results.txt",kin)); //read the Aexp from file

	double Aacc = getDoubleValue(config,"Aacc");
	double Api = getDoubleValue(config,"Api");	
	double Ain = getDoubleValue(config,"Ain");
	double Ap = getDoubleValue(config,"Ap");
	double Afsi = getDoubleValue(config,"Afsi");
	
	double facc = getDoubleValue(config,"facc");
	double fN2 = getDoubleValue(config,"fN2");
	double fpi = getDoubleValue(config,"fpi");
	double fin = getDoubleValue(config,"fin");
	double fp = getDoubleValue(config,"fp");
	double ffsi = getDoubleValue(config,"ffsi");

	double pol_n = 0.96;

	double Aphy_sum_numerator = 0.0;
	double Aphy_sum_denominator = 0.0;
	
	const double epsilon = 1e-10; //to handle division by zero

	//write physics asymmetry values to files
	std::ofstream outfile;
	outfile.open(Form("%s_physics_neutron_asymmetry_results.txt",kin));


    for (const auto& entry : rundata) {    
        int runNumber = entry.first;
        const RunData& runData = entry.second;

        // Check for NaN and very small values, skip if any are problematic
        if (std::isnan(runData.Aexp) || std::isnan(runData.errAexp) || 
            std::isnan(runData.pol_He3) || std::isnan(runData.pol_beam) ||
            runData.pol_He3 < epsilon || runData.pol_beam < epsilon) {
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
        if (std::fabs(denominator) < epsilon) {
            std::cerr << "Skipping run " << runNumber << " due to small denominator, potential division by zero." << std::endl;
            continue;
        }

        // Calculate Aphy and errAphy_stat
        double Aphy = (runData.Aexp - facc * Aacc - fpi * Api - fin * Ain - fp * Ap - ffsi * Afsi) / denominator;
        double errAphy_stat = runData.errAexp / denominator;

        outfile << "Run_number: " << runNumber << " | Aphy: " << Aphy << " | errAphy: " << errAphy_stat << std::endl;

        if (errAphy_stat>epsilon){
        	Aphy_sum_numerator += (Aphy / (errAphy_stat * errAphy_stat));
        	Aphy_sum_denominator += (1 / (errAphy_stat * errAphy_stat)); 
    	}
    }

    // Compute weighted average if denominator is non-zero
    if (Aphy_sum_denominator > epsilon) {
        double Aphy_sum = Aphy_sum_numerator / Aphy_sum_denominator;
        double errAphy_sum = 1/sqrt(Aphy_sum_denominator);
        std::cout <<"Aphy_sum_numerator : " << Aphy_sum_numerator << " Aphy_sum_denominator : "<< Aphy_sum_denominator <<" Aphy_sum: " << Aphy_sum 
        <<" errAphy_sum : "<< errAphy_sum << std::endl;

    } else {
        std::cerr << "Error: Aphy_sum_denominator is zero or too small, cannot calculate Aphy_sum." << std::endl;
    }

    outfile.close();

}	