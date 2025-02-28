
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


TGraphErrors* CalculateProtonAsymmetry(std::vector<TH1D*>& Helicity_histograms,
                                 const char* printfilename,
                                 const char* kin,
                                 bool flag_eHCAL_cut)
{
    // Prepare arrays to store the bin centers, asymmetries, and errors
    int nBins = Helicity_histograms.size();
    std::vector<double> dx_bin_centers(nBins);
    std::vector<double> asymmetries(nBins);
    std::vector<double> errors(nBins);
    
    // Open a file to save the results
    std::ofstream outfile(Form("txt/%s_proton_asymmetry_graph_dx_binned_eHCAL_cut_%s.txt",
                               kin,
                               std::to_string(flag_eHCAL_cut).c_str()));
    outfile << "dx_bin_center, Proton asymmetry, Error\n";

    for (int i = 0; i < nBins; ++i) {
        TH1D* h = Helicity_histograms[i];
        
        // Get the bin content for helicity = +1 and helicity = -1
        int bin_plus  = h->FindBin( 1.0);  // Helicity +1
        int bin_minus = h->FindBin(-1.0);  // Helicity -1
        
        double N_plus  = h->GetBinContent(bin_plus);
        double N_minus = h->GetBinContent(bin_minus);
        double N_total = N_plus + N_minus;

        double asymmetry = 0.0;
        double error     = 0.0;

        if (N_total > 0.0) {
            asymmetry = 100.0 * (N_plus - N_minus) / N_total;
            // Binomial error (approx)
            error = 100.0 * std::sqrt((4.0 * N_plus * N_minus) / 
                                      (N_total * N_total * N_total));
        }

        // For demonstration, define the bin center
        // (You might want to base this on real cointime bin edges.)
        double dx_bin_low_edge = -4.0;
        double dx_bin_width    = 0.1;
        double dx_bin_center   = i * dx_bin_width 
                                     + dx_bin_low_edge 
                                     + 0.5 * dx_bin_width;

        // Store in vectors
        dx_bin_centers[i] = dx_bin_center;
        asymmetries[i]          = asymmetry;
        errors[i]               = error;

        // Write to file
        outfile << dx_bin_center << "," 
                << asymmetry << "," 
                << error << "\n";
    }
    
    outfile.close();

    // Create the TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nBins,
                                           &dx_bin_centers[0],
                                           &asymmetries[0],
                                           nullptr,  // no x errors
                                           &errors[0]);
    graph->SetTitle("Proton Asymmetry vs dx; dx (m); Asymmetry (%)");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    // Return the graph to be drawn elsewhere
    return graph;
}


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
    std::cout<<"run_num_L: "<< run_num_L <<" run_num_H: "<< run_num_H <<endl;
    std::cout<<"coin_time_L: "<< coin_time_L <<" coin_time_H: "<< coin_time_H <<endl;
    std::cout<<"W2_L: "<< W2_L <<" W2_H: "<< W2_H <<endl;
    std::cout<<"dx_L: "<< dx_L<<" dx_H: "<< dx_H<<endl;
    std::cout<<"dy_L: "<< dy_L<<" dy_H: "<< dy_H<<endl;


    TH2D *h_dxdy_p_cut_W2_cointime = new TH2D("h_dxdy_p_cut_W2_cointime","dxdy (ntrack_sbs>0); dy (m); dx (m)",200,-4,4,200,-4,4);
    TH1D *h_dx_p_cut_W2_cointime = new TH1D("h_dx_p_cut_W2_cointime", "dx (ntrack_sbs>0); dx (m)", 200,-4,4);
    TH1D *h_dx_p_cut_W2_cointime_dy = new TH1D("h_dx_p_cut_W2_cointime_dy", "dx (ntrack_sbs>0); dx (m)", 200,-4,4);
    TH1D *h_dy_p_cut_W2_cointime = new TH1D("h_dy_p_cut_W2_cointime", "dy (ntrack_sbs>0); dy (m)", 200,-4,4);
    TH1D *h_runnum_sbs_tracks = new TH1D("h_runnum_sbs_tracks","sbs tracks across the kinematic", run_num_H-run_num_L+1,run_num_L,run_num_H);

    TH1D *h_trP_sbs = new TH1D("h_trP_sbs","tracking momentum of the hadron ; pN(GeV)",200,0,8);


    // Setup bins for dx
    const double binMin   = -4.0;
    const double binMax   = 4.0;
    const double binWidth = 0.1;
    const int nBins = static_cast<int>((binMax - binMin) / binWidth) + 1;

    // Prepare histograms for helicity in each dx bin
    std::vector<TH1D*> Helicity_histograms;
    Helicity_histograms.reserve(nBins);
    for (int i = 0; i < nBins; ++i) {
        double binLowEdge  = binMin + i*binWidth;
        double binHighEdge = binLowEdge + binWidth;
        TH1D* hist = new TH1D(Form("hist_bin_%d", i),
                              Form("Helicity for bin [%.1f,%.1f)", 
                                   binLowEdge, binHighEdge),
                              5, -2.5, 2.5);
        Helicity_histograms.push_back(hist);
    }


    for (int i = 0; i<nentries; i++){
        tree->GetEntry(i);
        if(run_num_L<runnum and runnum<run_num_H){
	        if(lookupValue(HelicityCheck,runnum)==1 and lookupValue(MollerQuality,runnum)==1 and ( helicity == -1 or helicity == 1)){
	        	bool cut_eHCAL = (eHCAL)>eHCAL_L;
	        	bool cut_W2 = (W2_L<W2 and W2<W2_H);
	        	bool cut_coin = (coin_time_L<coin_time and coin_time<coin_time_H);
	        	bool cut_dy_p = ((dy_p_L<dy and dy<dy_p_H));
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

				if(cut_eHCAL and cut_W2 and cut_coin and cut_sbs_track and cut_trP_sbs and cut_dy_p){
					h_dx_p_cut_W2_cointime_dy->Fill(dx);
					// Fill the helicity histogram for dx bin
	                int binIndex = static_cast<int>((dx - binMin) / binWidth);
	                if (binIndex >= 0 && binIndex < nBins) {
	                    Helicity_histograms[binIndex]->Fill(helicity);
	                }
            	}
			}
		}
		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
        std::cout.flush();
	}


    TCanvas * c = new TCanvas("c","c",3600,3000);
    TCanvas * c1 = new TCanvas("c1","c1",3600,3000);
    TCanvas * c2 = new TCanvas("c2","c2",3600,3000);

    c->Divide(2,2);
    c1->Divide(2,2);
    c2->Divide(2,2);

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

    c2->cd(1);
    // top sub‐pad ( ~70% of vertical space )
    TPad* p2_top = new TPad("p2_top", "p2_top", 0.0, 0.3, 1.0, 1.0);
    p2_top->SetBottomMargin(0.02); // small bottom margin
    p2_top->Draw();
    p2_top->cd();

    // Optionally hide the x label here, since we'll put it on bottom
    h_dx_p_cut_W2_cointime_dy->GetXaxis()->SetTitle("");
    h_dx_p_cut_W2_cointime_dy->Draw("hist");
    //box_anti1->Draw("SAME");
    //box_anti2->Draw("SAME");

    p2_top->Update();

    // --------------------------------------------------
    // bottom sub‐pad ( ~30% ) for TGraphErrors
    // --------------------------------------------------
    c2->cd(1); // go back to main pad #2
    TPad* p2_bottom = new TPad("p2_bottom","p2_bottom", 0.0, 0.0, 1.0, 0.3);
    p2_bottom->SetTopMargin(0.02);
    p2_bottom->SetBottomMargin(0.25); // enough space for axis labels
    p2_bottom->Draw();
    p2_bottom->cd();

    // Build the TGraphErrors from Helicity_histograms
    TGraphErrors* gAsym = CalculateProtonAsymmetry(Helicity_histograms,
                                             printfilename,
                                             kin,
                                             flag_eHCAL_cut);

    // Draw the TGraph in the bottom pad
    gAsym->Draw("AP");
    gAsym->SetMarkerStyle(kFullCircle);
    gAsym->SetMarkerColor(kBlue);
    gAsym->SetLineColor(kBlue);

    // Force x range to match the histogram range (0 → 250)
    //gAsym->GetXaxis()->SetRangeUser(-4.0, 4.0);
    gAsym->GetXaxis()->SetLimits(-4.0, 4.0);

    // Create a horizontal line at y=0, spanning -4 to 4 on the x-axis
    TLine* lineZero = new TLine(-4.0, 0.0, 4.0, 0.0);
    lineZero->SetLineColor(kRed);
    lineZero->SetLineStyle(2);  // e.g. dashed line
    lineZero->SetLineWidth(1);
    lineZero->Draw("SAME");

    // Now set the bottom pad’s X axis label
    gAsym->GetXaxis()->SetTitle("dx (ns)");

    gAsym->SetMaximum(20);
    gAsym->SetMinimum(-20);

    p2_bottom->Update();


    c->SaveAs(Form("plots/proton_plots_for_%s_eHCAL_cut_%s.pdf(",kin,std::to_string(flag_eHCAL_cut).c_str()));
    c1->SaveAs(Form("plots/proton_plots_for_%s_eHCAL_cut_%s.pdf",kin,std::to_string(flag_eHCAL_cut).c_str()));
    c2->SaveAs(Form("plots/proton_plots_for_%s_eHCAL_cut_%s.pdf)",kin,std::to_string(flag_eHCAL_cut).c_str()));
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