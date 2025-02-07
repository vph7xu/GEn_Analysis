
#include "models.h"
#include "parse.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <TDatime.h>

void proton_contamination(const char* filename, const char* printfilename, const char* kin, bool flag_eHCAL_cut){

	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	std::map<int, int> HelicityCheck = readCSVToMap("DB/Helicity_quality.csv");
	std::map<int, int> MollerQuality = readCSVToMap("DB/Moller_quality.csv");

	double coin_time_L = getDoubleValue(config,"coin_time_L");
	double coin_time_H = getDoubleValue(config,"coin_time_H");

    double coin_time_ac_L = getDoubleValue(config,"coin_time_ac_L");
    double coin_time_ac_H = getDoubleValue(config,"coin_time_ac_H");

	double W2_L = getDoubleValue(config,"W2_L");
	double W2_H = getDoubleValue(config,"W2_H");

	double dy_L = getDoubleValue(config,"dy_L");
	double dy_H = getDoubleValue(config,"dy_H");

	double dx_L = getDoubleValue(config,"dx_L");
	double dx_H = getDoubleValue(config,"dx_H");

	double eHCAL_L = getDoubleValue(config,"eHCAL_L");

	double dy_p_L = getDoubleValue(config,"dy_p_L");
	double dy_p_H = getDoubleValue(config,"dy_p_H");

	double dx_p_L = getDoubleValue(config,"dx_p_L");
	double dx_p_H = getDoubleValue(config,"dx_p_H");

	double dy_ac_L = getDoubleValue(config,"dy_ac_L");
	double dy_ac_H = getDoubleValue(config,"dy_ac_H");

	double P_sbs_L = getDoubleValue(config,"P_sbs_L");
	double P_sbs_H = getDoubleValue(config,"P_sbs_H");

	double run_num_L = getDoubleValue(config,"run_num_L");
	double run_num_H = getDoubleValue(config,"run_num_H");

	double IHWP_flip = getDoubleValue(config,"IHWP_flip");

    double coin_time_offset_L = coin_time_L+30.0;
    double coin_time_offset_H = coin_time_H+30.0;

	//end of parsing cuts
	
	TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");
	
	int runnum = 0;
	int helicity = 0;
	int IHWP = 0;
    double dx = 0.0;
    double dy = 0.0;
    double W2 = 0.0;
    double Q2 = 0.0;
	double coin_time = 0.0;
	double ntrack = 0;
	double ntrack_sbs = 0;
    double vz = 0.0;
    double eHCAL = 0.0;
    double ePS = 0.0;
    double trP_sbs = 0.0;

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("IHWP",&IHWP);
    tree->SetBranchAddress("dx",&dx);
    tree->SetBranchAddress("dy",&dy);
    tree->SetBranchAddress("W2",&W2);
    tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("ntrack", &ntrack);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
    tree->SetBranchAddress("vz", &vz);
    tree->SetBranchAddress("eHCAL", &eHCAL);
    tree->SetBranchAddress("ePS",&ePS);
    tree->SetBranchAddress("trP_sbs",&trP_sbs);


    int nentries = tree->GetEntries();
    double Nplus_total = 0.0;
    double Nminus_total = 0.0;
    double proton_events = 0.0;
    double QE_events = 0.0;

	if (flag_eHCAL_cut == false) eHCAL_L=0.0;

    std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;

    TH2D *h_dxdy_p_cut_W2_cointime = new TH2D("h_dxdy_p_cut_W2_cointime","dxdy (ntrack_sbs>0); dy (m); dx (m)",200,-4,4,200,-4,4);
    TH1D *h_dx_p_cut_W2_cointime = new TH1D("h_dx_p_cut_W2_cointime", "dx (ntrack_sbs>0); dx (m)", 200,-4,4);
    TH1D *h_dy_p_cut_W2_cointime = new TH1D("h_dy_p_cut_W2_cointime", "dy (ntrack_sbs>0); dy (m)", 200,-4,4);
    TH1D *h_runnum_sbs_tracks = new TH1D("h_runnum_sbs_tracks","sbs tracks across the kinematic", run_num_H-run_num_L+1,run_num_L,run_num_H);

    TH1D *h_trP_sbs = new TH1D("h_trP_sbs","tracking momentum of the hadron ; pN(GeV)",200,0,8);

    for (int i = 0; i<nentries; i++){
        tree->GetEntry(i);
        if(run_num_L<runnum and runnum<run_num_H){
	        if(lookupValue(HelicityCheck,runnum)==1 and lookupValue(MollerQuality,runnum)==1 and ( helicity == -1 or helicity == 1)){
	        	bool cut_eHCAL = (eHCAL)>eHCAL_L;
	        	bool cut_W2 = (W2_L<W2 and W2<W2_H);
	        	bool cut_coin = (coin_time_L<coin_time and coin_time<coin_time_H);
	        	bool cut_QE = (W2_L<W2 and W2<W2_H) and (coin_time_L<coin_time and coin_time<coin_time_H) and (dy_L<dy and dy<dy_H) and (dx_L<dx and dx<dx_H);
	        	bool cut_p = (W2_L<W2 and W2<W2_H) and (coin_time_L<coin_time and coin_time<coin_time_H) and (dy_p_L<dy and dy<dy_p_H) and (dx_p_L<dx and dx<dx_p_H);
	        	bool cut_sbs_track = (ntrack_sbs>0);
	        	bool cut_trP_sbs = P_sbs_L<trP_sbs and trP_sbs<P_sbs_H;

				helicity = -1*IHWP*IHWP_flip*helicity;

				if(cut_sbs_track){
					h_runnum_sbs_tracks->Fill(runnum);
				}

				if(cut_eHCAL and cut_W2 and cut_coin and cut_sbs_track ){
					h_trP_sbs->Fill(trP_sbs);

				}

				if(cut_eHCAL and cut_W2 and cut_coin and cut_sbs_track and cut_trP_sbs){	
					h_dxdy_p_cut_W2_cointime->Fill(dy,dx);
					h_dx_p_cut_W2_cointime->Fill(dx);
					h_dy_p_cut_W2_cointime->Fill(dy);
				}

				if (cut_eHCAL and cut_QE and cut_sbs_track and cut_trP_sbs){ //proton events under the neutron peak
					proton_events+=1;
				}

				if (cut_eHCAL and cut_QE){ //QE events in the neutron peak
					QE_events+=1;
				}

				if (cut_eHCAL and cut_p and cut_trP_sbs){ //events under the proton peak
					if (helicity == 1){
						Nplus_total+=1;
					}
					if (helicity == -1){
						Nminus_total+=1;
					}
				}
			}
		}
		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
        std::cout.flush();
	}

    TCanvas * c = new TCanvas("c","c",3600,3000);
    TCanvas * c1 = new TCanvas("c1","c1",3600,3000);
    c->Divide(2,2);
    c1->Divide(2,2);
		
    //Create a box
    TBox* box_dxdy_n = new TBox(dy_L,dx_L,dy_H,dx_H);
    box_dxdy_n->SetFillStyle(0);
    box_dxdy_n->SetLineColor(kRed);

    TBox* box_dxdy_p = new TBox(dy_p_L,dx_p_L,dy_p_H,dx_p_H);
    box_dxdy_p->SetFillStyle(0);
    box_dxdy_p->SetLineColor(kBlue);

    //Create lines
    TLine *line1 = new TLine(P_sbs_L,0.0,P_sbs_L,h_trP_sbs->GetMaximum());
	line1->SetLineColor(kBlue);
	line1->SetLineWidth(2);

    TLine *line2 = new TLine(P_sbs_H,0.0,P_sbs_H,h_trP_sbs->GetMaximum());
	line2->SetLineColor(kBlue);
	line2->SetLineWidth(2);

    c->cd(1);
    h_dxdy_p_cut_W2_cointime->Draw("COLZ");
    box_dxdy_n->Draw("SAME");
    box_dxdy_p->Draw("SAME");

    c->cd(2);
    h_runnum_sbs_tracks->Draw();
    h_runnum_sbs_tracks->SetXTitle("run number");

    c->cd(3);
    h_dx_p_cut_W2_cointime->Draw();

    c->cd(4);
    h_dy_p_cut_W2_cointime->Draw();

    c1->cd(1);
    h_trP_sbs->Draw();
    line1->Draw();
    line2->Draw();

    c->SaveAs(Form("plots/proton_plots_for_%s_eHCAL_cut_%s.pdf(",kin,std::to_string(flag_eHCAL_cut).c_str()));
    c1->SaveAs(Form("plots/proton_plots_for_%s_eHCAL_cut_%s.pdf)",kin,std::to_string(flag_eHCAL_cut).c_str()));
    //c->SaveAs(Form("plots/proton_plots_for_%s_eHCAL_cut_%s.pdf",kin,std::to_string(flag_eHCAL_cut).c_str()));
    //c->SaveAs(Form("plots/proton_plots_for_%s_eHCAL_cut_%s.pdf",kin,std::to_string(flag_eHCAL_cut).c_str()));

	double proton_frac = proton_events/QE_events;// this is not correct should remove other fractions before doing this
	double errproton_frac = (proton_events/QE_events)*sqrt((1/proton_events)+(1/QE_events)); // this is not correct should remove other fractions before doing this
	
	//get the proton asymmetry
	double Ap = (Nplus_total-Nminus_total)/(Nplus_total+Nminus_total);
	double errAp = 2*sqrt((Nplus_total*Nminus_total)*(Nplus_total+Nminus_total))/((Nplus_total+Nminus_total)*(Nplus_total+Nminus_total));;
		
	std::cout<<"Nplus : "<<Nplus_total<<endl;
	std::cout<<"Nminus : "<<Nminus_total<<endl;
	std::cout<<"A_p : "<<Ap<<endl;
	std::cout<<"err_A_p : "<<errAp<<endl;
	std::cout<<"proton_events : "<< proton_events<<endl;
	std::cout<<"QE_events : "<<QE_events<<endl;
	std::cout<<"proton fraction (double counting) : "<<proton_frac<<endl; //not correct , double counting
	std::cout<<"error proton fraction (double counting) : "<<errproton_frac<<endl; //not correct , double counting

	std::ofstream outfile;
    outfile.open(Form("txt/%s_proton_asymmetry_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));
	outfile<<"Nplus = "<<Nplus_total<<endl;
	outfile<<"Nminus = "<<Nminus_total<<endl;
	outfile<<"A_p = "<<Ap<<endl;
	outfile<<"err_A_p = "<<errAp<<endl;
	outfile<<"proton_events = "<< proton_events<<endl;
	outfile<<"QE_events = "<<QE_events<<endl;
	outfile<<"f_p = "<<proton_frac<<endl; //not correct , double counting
	outfile<<"err_f_p = "<<errproton_frac<<endl; //not correct , double counting

}