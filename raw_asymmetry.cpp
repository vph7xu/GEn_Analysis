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

// Struct to hold each row of CSV data
struct beamPol {
    TDatime startTime;
    TDatime endTime;
    double beam_polarization;
    double error;
};

TDatime parseDateTime(const std::string &datetimeStr) {
    int year, month, day, hour, minute, second;
    sscanf(datetimeStr.c_str(), "%d-%d-%d %d:%d:%d", &year, &month, &day, &hour, &minute, &second);
    return TDatime(year, month, day, hour, minute, second);
}

// Function to load the CSV data into a vector of DataEntry
std::vector<beamPol> loadCSVData(const std::string &filename) {
    std::vector<beamPol> data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error opening CSV file!" << std::endl;
        return data;
    }

    std::string line;
    std::getline(file, line); // Skip header line

    // Read each line of the CSV file
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string startStr, endStr, polarizationStr, errorStr;

        // Parse each field in the row
        std::getline(ss, startStr, ',');
        std::getline(ss, endStr, ',');
        std::getline(ss, polarizationStr, ',');
        std::getline(ss, errorStr, ',');

        // Create a DataEntry with parsed data
        beamPol beampol;
        beampol.startTime = parseDateTime(startStr);
        beampol.endTime = parseDateTime(endStr);
        beampol.beam_polarization = std::stod(polarizationStr);
        beampol.error = std::stod(errorStr);

        data.push_back(beampol);
    }
    file.close();
    return data;
}


// Function to search for an entry in the data vector that matches the given lookupTime
std::pair<double, double> searchData(const std::vector<beamPol> &data, const TDatime &lookupTime) {
    for (const auto &entry : data) {
        if (lookupTime >= entry.startTime && lookupTime <= entry.endTime) {
            return {entry.beam_polarization, entry.error}; // Return the polarization value if a match is found
        }
    }
    return {-1.0,-1.0}; // Return -1.0 if no matching time range is found
}


void raw_asymmetry(const char* filename, const char* printfilename, const char* kin){

	std::map<std::string, std::string> config1 = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
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

	double dy_p_L = getDoubleValue(config1,"dy_p_L");
	double dy_p_H = getDoubleValue(config1,"dy_p_H");

	double dx_p_L = getDoubleValue(config1,"dx_p_L");
	double dx_p_H = getDoubleValue(config1,"dx_p_H");

	double IHWP_flip = getDoubleValue(config1,"IHWP_flip");

	double run_num_L = getDoubleValue(config1,"run_num_L");
	double run_num_H = getDoubleValue(config1,"run_num_H");

	//end of parsing cuts
	
	// Load CSV data
    std::string csvFile = "DB/Beam_pol.csv";
    std::vector<beamPol> data = loadCSVData(csvFile);
    if (data.empty()) {
        return 1; // Exit if loading CSV data failed
    }


	TFile* file = TFile::Open(filename);
	//TFile* sim_file = TFile::Open(sim_filename);

    TTree* tree = (TTree*)file->Get("Tout");
	//TTree* sim_tree = (TTree*)sim_file->Get("Tout");
	
	double dx = 0.0;
	double dy = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;
	double eHCAL = 0.0;

	//double coin_time = 0.0;

	int runnum = 0;
	int helicity = 0;
	int IHWP = 0;
	double He3Pol = 0.0;

	double Nplus = 0;
	double Nminus = 0;

	double Nplus_total = 0.0;
	double Nminus_total = 0.0;

	double Pplus = 0;
	double Pminus = 0;

	double Aexp = 0.0;
	double errAexp = 0.0;

	double Ap_exp = 0.0;
	double errAp_exp = 0.0;

	double avg_He3Pol = 0.0;
	double avg_beampol = 0.0;
	double total_QE = 0.0;

	double err_avg_beampol = 0.0;

	double runx=0;

	TDatime *datetime = nullptr;

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("IHWP",&IHWP);
	tree->SetBranchAddress("He3Pol",&He3Pol);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("datetime", &datetime);

	TH1I *h_IHWP = new TH1I("h_IHWP","IHWP",10,-2,2);
	TH1I *h_runnum = new TH1I("h_runnum","runnum",300,0,3000);

	TGraphErrors *gAsym = new TGraphErrors();
	TGraphErrors *gAp_sym = new TGraphErrors();

	//neutron spot cut
	//TCutG *cutg = CreateOvalCut("cutg",0,-0.025,0.35,0.2,100);

	TH2D* nspot_cut = new TH2D("nspot_cut","neutron spot",250,-2,2,250,-4,4);

	int nentries = tree->GetEntries();

	//write asymmetry values to files
	std::ofstream outfile;
	outfile.open(Form("%s_raw_neutron_asymmetry_results.txt",kin)); 

	std::ofstream outfile_p;
	outfile_p.open(Form("%s_raw_proton_asymmetry_results.txt",kin));

	std::ofstream outfile_pol_He3;
	outfile_pol_He3.open(Form("%s_He3_polarization_results.txt",kin));

	std::ofstream outfile_pol_beam;
	outfile_pol_beam.open(Form("%s_beam_polarization_results.txt",kin));

	std::ofstream outfile_n_asym;
	outfile_n_asym.open(Form("%s_raw_neutron_asymmetry_only_results.txt",kin));

	std::ofstream outfile_avg_polarizations;
	outfile_avg_polarizations.open(Form("%s_average_polarization_results.txt",kin));


	for (int i=0; i<nentries; i++){
		tree->GetEntry(i);
		
		std::string lookupDatenTime = std::to_string(datetime->GetYear())+"-"+std::to_string(datetime->GetMonth())+"-"+std::to_string(datetime->GetDay())+" "
		+std::to_string(datetime->GetHour())+":"+std::to_string(datetime->GetMinute())+":"+std::to_string(datetime->GetSecond());
		//std::cout<<"date : "<<datetime->AsString(); 

		TDatime lookupTime = *datetime;

		h_IHWP->Fill(IHWP);
		h_runnum->Fill(runnum);
		if (eHCAL>0.3 and (W2_L<W2 and W2<W2_H)and(coin_time_L<coin_time and coin_time<coin_time_H)){
			nspot_cut->Fill(dy,dx);
		}
		if (i==0){
			runx = runnum;
		}
		

		//if(abs(helicity)>0){
			if(runnum==runx and lookupValue(HelicityCheck,runnum)==1 and lookupValue(MollerQuality,runnum)==1){
				
				if (IHWP == 1) helicity = -IHWP*IHWP_flip*helicity;
	            else if (IHWP == -1) helicity = -IHWP*IHWP_flip*helicity;
	            else continue;

				if ( eHCAL>0.3 and (W2_L<W2 and W2<W2_H)and(coin_time_L<coin_time and coin_time<coin_time_H)){
					//if(cutg->IsInside(dy,dx)){
					if ((dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H)){

						if (helicity==1){
							Nplus=Nplus+1;
						}
						else if (helicity==-1){
							Nminus=Nminus+1;
						}
					}
					else if((dx_p_L<dx and dx<dx_p_H) and (dy_p_L<dy and dy<dy_p_H)){
						//std::cout<<"here"<<endl;
						if (helicity==1){
							Pplus=Pplus+1;
						}
						else if (helicity==-1){
							Pminus=Pminus+1;
						}					
					}
						//}
						//else{
						//	if (helicity==1){
	                                        //                Nminus=Nminus+1;
	                                        //        }
	                                        //        else if (helicity==-1){
	                                        //                Nplus=Nplus+1;
	                                        //        }
						//}
					
				}
			}

			else if(runnum!=runx and lookupValue(HelicityCheck,runx)==1 and lookupValue(MollerQuality,runnum)==1){
				Aexp = (Nplus-Nminus)*100/(Nplus+Nminus);
				errAexp = 2*100*sqrt((Nplus*Nminus)*(Nplus+Nminus))/((Nplus+Nminus)*(Nplus+Nminus));

				Ap_exp = (Pplus-Pminus)*100/(Pplus+Pminus);
				errAp_exp = 2*100*sqrt((Pplus*Pminus)*(Pplus+Pminus))/((Pplus+Pminus)*(Pplus+Pminus));

				auto [beam_polarization, error] = searchData(data, lookupTime);

				avg_He3Pol += (Nplus+Nminus)*He3Pol*0.01;
				avg_beampol += (Nplus+Nminus)* beam_polarization*0.01;
				total_QE +=(Nplus+Nminus);

				err_avg_beampol += (Nplus+Nminus)*error*0.01;

				Nplus_total+=Nplus;
				Nminus_total+=Nminus;

				gAsym->SetPoint(runx,runx,Aexp);
				gAsym->SetPointError(runx,0,errAexp);

				gAp_sym->SetPoint(runx,runx,Ap_exp);
				gAp_sym->SetPointError(runx,0,errAp_exp);

				//write asymmetries to files
				outfile << "Run_number: " << runx << " | Aexp: " << Aexp/100 << " | errAexp: " << errAexp/100 << " | pol_He3: " << He3Pol <<  " | pol_beam: " << beam_polarization << std::endl;
	            outfile_n_asym << "Run_number: " << runx << " Nplus : "<< Nplus<<" Nminus : "<<Nminus<<" | Aexp: " << Aexp/100 << " | errAexp: " << errAexp/100 << std::endl;
	            outfile_p << "Run_number: " << runx << " | Ap_exp: " << Ap_exp/100 << " | errAexp: " << errAp_exp/100 << std::endl;
	            outfile_pol_He3 << "Run_number: " << runx << " | pol_He3: " << He3Pol <<std::endl; //<< " | errAexp: " << errAp_exp << std::endl;
	            outfile_pol_beam<< "Run_number: " << runx << " | pol_beam: " << beam_polarization<<std::endl;
	            //write He3 polarizations to files


	            //write beam polarization to files



				std::cout<<"runnum : "<<runx<<" Nplus : "<<Nplus<<" Nminus : "<<Nminus<<" Aexp : "<<Aexp<<endl;
				std::cout<<"runnum : "<<runx<<" Pplus : "<<Pplus<<" Pminus : "<<Pminus<<" Ap_exp : "<<Ap_exp<<endl;
	            std::cout<< "Run_number: " << runx << " | pol_He3: " << He3Pol <<std::endl; //<< " | errAexp: " << errAp_exp << std::endl;
	            std::cout<< "Run_number: " << runx << " | pol_beam: " << beam_polarization<<std::endl;

				runx=runnum;
				Nplus=0.0;
				Nminus=0.0;
				Aexp=0.0;
				Pplus=0.0;
				Pminus=0.0;
				Ap_exp=0.0;
			}

			else if(runnum!=runx){
	                      	runx=runnum;
	                        Nplus=0.0;
	                        Nminus=0.0;
	                        Aexp=0.0;
	                        Pplus=0.0;
							Pminus=0.0;
							Ap_exp=0.0;
			}

		//}
	}


	avg_He3Pol = avg_He3Pol/total_QE;
	avg_beampol = avg_beampol/total_QE;
	err_avg_beampol = err_avg_beampol/total_QE;

	outfile_avg_polarizations<< "avg_He3pol = " << avg_He3Pol<<endl;
	outfile_avg_polarizations<< "avg_beampol = "<< avg_beampol<<endl;
	outfile_avg_polarizations<< "avg_Pn = "<< 0.96 <<endl;
	outfile_avg_polarizations<< "err_avg_He3pol = "<< 0.05<<endl;
	outfile_avg_polarizations<< "err_avg_beampol = "<< err_avg_beampol<<endl;
	outfile_avg_polarizations<< "err_Pn = "<< 0.005<<endl;

	std::cout<<"Nplus_total : "<<Nplus_total<<" Nminus_total : "<<Nminus_total<<" Aexp_total : "<<(Nplus_total - Nminus_total)/(Nplus_total + Nminus_total)<<
				" errAexp_total : " << 2*sqrt((Nplus_total*Nminus_total)*(Nplus_total+Nminus_total))/((Nplus_total+Nminus_total)*(Nplus_total+Nminus_total))<<endl;



	outfile.close();
	outfile_p.close();

	//TCutG * cutsq = CreateSquareCut(dy_L,dx_L,dy_H,dx_H);
	//TCutG * cutsq1 = CreateSquareCut(dy_p_L,dx_p_L,dy_p_H,dx_p_H);

	TCanvas *c = new TCanvas("c","c",3600,3000);
	TCanvas *c1 = new TCanvas("c1","c1",3600,900);
	
	c->Divide(2,2);
	c->cd(1);
	h_IHWP->Draw();
	c->cd(2);
	h_runnum->Draw();
	c->cd(3);
	nspot_cut->Draw("COLZ");
	CreateSquareCut(dy_L,dx_L,dy_H,dx_H);
	CreateSquareCut(dy_p_L,dx_p_L,dy_p_H,dx_p_H);
	//cutg->Draw("L");
	
	c1->Divide(1,1);
	c1->cd(1);
	gAsym->GetXaxis()->SetTitle("run number");
	gAsym->GetYaxis()->SetTitle("Raw Asymmetry");
	gAsym->GetXaxis()->SetLimits(run_num_L,run_num_H);
	gAsym->GetYaxis()->SetRangeUser(-10,10);
	gAsym->SetMarkerStyle(7);
	gAsym->SetMarkerSize(5.0);
	gAsym->SetMarkerColor(kBlue);
	gAsym->Draw("AP");

	gAp_sym->GetXaxis()->SetLimits(run_num_L,run_num_H);
	gAp_sym->GetYaxis()->SetRangeUser(-10,10);
	gAp_sym->SetMarkerStyle(7);
	gAp_sym->SetMarkerSize(5.0);
	gAp_sym->SetMarkerColor(kRed);
	gAp_sym->Draw("P SAME");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);  // Coordinates for the legend box (x1, y1, x2, y2)
    
    // Add entries to the legend
    legend->AddEntry(gAsym, "(e,e'n) events", "lp"); // "lp" means line and point in the legend
    legend->AddEntry(gAp_sym, "(e,e'p)events", "lp");

    // Set legend properties
    legend->SetTextSize(0.03);    // Set text size
    legend->SetBorderSize(0);     // No border around the legend box
    legend->SetFillStyle(0);      // Transparent background for the legend

    // Draw the legend
    legend->Draw();

	c->Print(Form("%s.pdf(",printfilename));
	c1->Print(Form("%s.pdf)",printfilename));	
	c1->SaveAs(Form("%s.png",printfilename));

}
