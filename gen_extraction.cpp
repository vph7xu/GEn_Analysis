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


struct FieldAngles {
    double horizontal;
    double vertical;
};

FieldAngles lookupFieldAngles(const std::string& kinematic, const std::string& filePath) {
    std::ifstream file(filePath);
    std::string line;
    FieldAngles angles = { -1, -1 };  // Default values if kinematic is not found

    // Skip the header line
    std::getline(file, line);

    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string kin, horizontalStr, verticalStr;

        // Read fields from the line
        std::getline(ss, kin, ',');
        std::getline(ss, horizontalStr, ',');
        std::getline(ss, verticalStr, ',');

        // Remove any leading/trailing whitespace
        kin.erase(0, kin.find_first_not_of(" \t"));
        kin.erase(kin.find_last_not_of(" \t") + 1);

        if (kin == kinematic) {
            angles.horizontal = std::stod(horizontalStr);
            angles.vertical = std::stod(verticalStr);
            break;
        }
    }

    if (angles.horizontal == -1 && angles.vertical == -1) {
        std::cerr << "Kinematic setting not found." << std::endl;
    }

    return angles;
}


void calAvgKin(const char* filename, const char* kin, bool flag_eHCAL_cut){

	TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

	//read the field angles
	FieldAngles angles = lookupFieldAngles(kin, "DB/Field_Meas.csv");

	//set target polarization direction
	TVector3 TargetPolDirection;
	TargetPolDirection.SetMagThetaPhi(1.0,angles.horizontal*TMath::Pi()/180,angles.vertical*TMath::Pi()/180);

	std::map<std::string, std::string> config1 = parseConfig(Form("cuts/cut_%s_He3.txt",kin)); //parse the cuts
	std::map<int, int> HelicityCheck = readCSVToMap("DB/Helicity_quality.csv");
	std::map<int, int> MollerQuality = readCSVToMap("DB/Moller_quality.csv");

	double coin_time_L = getDoubleValue(config1,"coin_time_L");
	double coin_time_H = getDoubleValue(config1,"coin_time_H");

	double W2_L = getDoubleValue(config1,"W2_L");
	double W2_H = getDoubleValue(config1,"W2_H");

	double dy_L = getDoubleValue(config1,"dy_L");
	double dy_H = getDoubleValue(config1,"dy_H");

	double dx_L = getDoubleValue(config1,"dx_L");
	double dx_H = getDoubleValue(config1,"dx_H");

	double eHCAL_L = getDoubleValue(config1,"eHCAL_L");

	double dy_p_L = getDoubleValue(config1,"dy_p_L");
	double dy_p_H = getDoubleValue(config1,"dy_p_H");

	double dx_p_L = getDoubleValue(config1,"dx_p_L");
	double dx_p_H = getDoubleValue(config1,"dx_p_H");

	double IHWP_flip = getDoubleValue(config1,"IHWP_flip");

	double run_num_L = getDoubleValue(config1,"run_num_L");
	double run_num_H = getDoubleValue(config1,"run_num_H");


	double Q2 = 0.0;
	double ebeam = 0.0;
	double trP = 0.0;
	double trPx = 0.0;
	double trPy = 0.0;
	double trPz = 0.0;
	double etheta =0.0;

	int runnum = 0.0;
	int helicity = 0.0;
	double eHCAL = 0.0;
	double coin_time = 0.0;
	double W2 = 0.0;
	double dx = 0.0;
	double dy = 0.0; 

	double m = 0.938;

	const int nexp = 6;
	double T_avg[nexp] = {0};
	double T_1_Q2_avg = 0;
	double T_1_tau_avg = 0;
	double T_1_epsilon_avg = 0;
	double T_1_Px_avg = 0;
	double T_1_Pz_avg = 0;

	double Q2_avg = 0;
	double tau_avg = 0;
	double epsilon_avg = 0;
	double Px_avg = 0;
	double Pz_avg = 0;

	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("ebeam",&ebeam);
	tree->SetBranchAddress("trP",&trP);
	tree->SetBranchAddress("trPx",&trPx);
	tree->SetBranchAddress("trPy",&trPy);
	tree->SetBranchAddress("trPz",&trPz);
	tree->SetBranchAddress("etheta",&etheta);

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("W2",&W2);

    if (flag_eHCAL_cut == false) eHCAL_L=0.0;

    std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;

	int nentries = tree->GetEntries();

	int count = 1;

	for (int i=0; i<nentries*1; i++){
		tree->GetEntry(i);

		bool goodMoller = (lookupValue(HelicityCheck,runnum)==1);
		bool goodHelicity = (lookupValue(MollerQuality,runnum)==1) and (abs(helicity)==1);
		bool goodrun = (run_num_L<runnum and runnum<run_num_H);
		bool goodQE = eHCAL>eHCAL_L and (coin_time_L<coin_time) and (coin_time<coin_time_H) and (W2_L<W2) and (W2<W2_H) and (dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H);

		if (goodMoller and goodHelicity and goodQE and goodrun){
			TLorentzVector Pe(0,0,ebeam,ebeam);
			TLorentzVector Peprime(trPx,trPy,trPz,trP);

			TLorentzVector q = Pe - Peprime;
			TVector3 Pe_vect = Pe.Vect();
			TVector3 Peprime_vect = Peprime.Vect();
			TVector3 q_vect = q.Vect();
			TVector3 normal = Pe_vect.Cross(Peprime_vect);
			normal = normal.Unit();
			q_vect = q_vect.Unit();

			double tau = Q2/(4*m*m);
			double tan_theta = tan(etheta/2);
			double epsilon = 1.0/(1+2*(1+tau)*pow(tan_theta,2));

			double Px = normal.Dot(q_vect.Cross(TargetPolDirection));
			double Pz = q_vect.Dot(TargetPolDirection);

			double B = -2 * sqrt(tau*(1+tau)) * tan_theta * Px;
			double C = -2 * tau * sqrt(1+tau+ pow((1+tau)*tan_theta,2) * tan_theta * Pz);
			double D = tau/epsilon;

			double T_0 = C/D;
			double T_1 = B/D;
		 	double T_2 = -1*C / (D*D);
	  		double T_3 = -1*B / (D*D);
	  		double T_4 = C / (D*D*D);
	  		double T_5 = B / (D*D*D);

	  	  	// Now we update the average for each value
			T_avg[0] += (T_0 - T_avg[0]) / count;
			T_avg[1] += (T_1 - T_avg[1]) / count;
			T_avg[2] += (T_2 - T_avg[2]) / count;
			T_avg[3] += (T_3 - T_avg[3]) / count;
			T_avg[4] += (T_4 - T_avg[4]) / count;
			T_avg[5] += (T_5 - T_avg[5]) / count;
			
			T_1_Q2_avg += (T_1*Q2 - T_1_Q2_avg) / count;
			T_1_tau_avg += (T_1*tau - T_1_tau_avg) / count;
			T_1_epsilon_avg += (T_1*epsilon - T_1_epsilon_avg) / count;
			T_1_Px_avg += (T_1*Px - T_1_Px_avg) / count;
			T_1_Pz_avg += (T_1*Pz - T_1_Pz_avg) / count;


			//Q2_avg = T_1_Q2_avg / T_avg[1];

			// We also get the average kinematic variables for later
			Q2_avg += (Q2 - Q2_avg) / count;
			tau_avg += (tau - tau_avg) / count;
			epsilon_avg += (epsilon - epsilon_avg) / count;
			Px_avg += (Px - Px_avg) / count;
			Pz_avg += (Pz - Pz_avg) / count;

			count+=1;
		}
        if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
        std::cout.flush();

	}

	std::ofstream outfile;
	outfile.open(Form("txt/%s_average_kinematic_values_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));

	outfile<<"Q2_avg = "<<T_1_Q2_avg<<endl;
	outfile<<"tau_avg = "<<T_1_tau_avg<<endl;
	outfile<<"epsilon_avg = "<<T_1_epsilon_avg<<endl;
	outfile<<"Px_avg = "<<T_1_Px_avg<<endl;
	outfile<<"Pz_avg = "<<T_1_Pz_avg<<endl;
	std::cout<<"count = "<<count<<endl;

}


void gen_extraction(const char* kin = "GEN3", bool flag_eHCAL_cut = false){

	std::map<std::string, std::string> Aphyresults = parseConfig(Form("txt/%s_He3_physics_neutron_asymmetry_summary_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));
	std::map<std::string, std::string> avg_kin_val = parseConfig(Form("txt/%s_average_kinematic_values_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));

	//physics asymmetries
	double Aphy_sum = getDoubleValue(Aphyresults,"Aphy_sum");
	double err_Aphy_stat_sum = getDoubleValue(Aphyresults,"err_Aphy_stat_sum");
	double err_Aphy_sys = getDoubleValue(Aphyresults,"err_Aphy_sys");

	//kinematic averages
	double Q2_avg = getDoubleValue(avg_kin_val,"Q2_avg");
	double tau_avg = getDoubleValue(avg_kin_val,"tau_avg");
	double epsilon_avg = getDoubleValue(avg_kin_val,"epsilon_avg");
	double Px_avg = getDoubleValue(avg_kin_val,"Px_avg");
	double Pz_avg = getDoubleValue(avg_kin_val,"Pz_avg");

	//define variables to solve for the ratio
	double A = (epsilon_avg/tau_avg)*Aphy_sum;
	double B = (sqrt((2*epsilon_avg*(1-epsilon_avg))/tau_avg))*Px_avg;
	double C = Aphy_sum + (sqrt(1-(epsilon_avg*epsilon_avg)))*Pz_avg;

	double lambda = (-B+sqrt(B*B-4*A*C))/(2*A);

	double sigma2_lambda_stat = (pow((C/(A*sqrt(B*B-4*A*C)) + lambda/A)*(epsilon_avg/tau_avg) + 1/sqrt((B*B -4*A*C))),2)*err_Aphy_stat_sum*err_Aphy_stat_sum;
	double sigma2_lambda_sys = (pow((C/(A*sqrt(B*B-4*A*C)) + lambda/A)*(epsilon_avg/tau_avg) + 1/sqrt((B*B -4*A*C))),2)*err_Aphy_sys*err_Aphy_sys;

	double sigma_lambda_stat = sqrt(sigma2_lambda_stat);
	double sigma_lambda_sys = sqrt(sigma2_lambda_sys);
	double sigma_lambda = sqrt(sigma2_lambda_stat + sigma2_lambda_sys);

	double GMn = -0.0189;
	double err_GMn = 0.0002;

	cout<<" lambda : "<< lambda << " sigma_lambda_stat : " <<sigma_lambda_stat << " sigma_lambda_sys : "<<sigma_lambda_sys <<" sigma_lambda : "<<sigma_lambda<<endl;

	std::ofstream outfile;
	outfile.open(Form("txt/%s_results_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));
	outfile<<"lambda = "<< lambda <<endl;
	outfile<<"sigma_lambda_stat = " <<sigma_lambda_stat <<endl;
	outfile<<"sigma_lambda_sys = "<<sigma_lambda_sys <<endl;
	outfile<<"sigma_lambda : "<<sigma_lambda<<endl;


}