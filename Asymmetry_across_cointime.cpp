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

void CalculateAsymmetry(std::vector<TH1D*>& Helicity_histograms, const char* printfilename, const char* kin, bool flag_eHCAL_cut) {
    // Prepare arrays to store the W2 bin centers, asymmetries, and errors
    int nBins = Helicity_histograms.size();
    std::vector<double> cointime_bin_centers(nBins);
    std::vector<double> asymmetries(nBins);
    std::vector<double> errors(nBins);
    
    // Open a file to save results
    std::ofstream outfile(Form("txt/%s_asymmetry_graph_cointime_binned_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));
    outfile << "cointime_bin_low_edge, Asymmetry, Error\n";

    for (int i = 0; i < nBins; ++i) {
        TH1D* h = Helicity_histograms[i];
        
        // Get the bin content for helicity = +1 and helicity = -1
        int bin_plus = h->FindBin(1.0);   // Helicity +1
        int bin_minus = h->FindBin(-1.0); // Helicity -1
        
        double N_plus = h->GetBinContent(bin_plus);
        double N_minus = h->GetBinContent(bin_minus);
        
        // Total counts in this bin
        double N_total = N_plus + N_minus;

        // Calculate asymmetry and its statistical error
        double asymmetry = 0.0;
        double error = 0.0;

        if (N_total > 0) {
            asymmetry = 100*(N_plus - N_minus) / N_total;
            
            // Error on asymmetry (binomial error)
            error = 100* sqrt((4 * N_plus * N_minus) / (N_total * N_total * N_total));
        }

        // Get the bin center (since we care about the middle of each W2 bin)
        double cointime_bin_low_edge = 0.0;
        double cointime_bin_width = 20;
        double cointime_bin_center = i*cointime_bin_width + cointime_bin_low_edge + 0.5 * cointime_bin_width;

        // Store results in arrays
        cointime_bin_centers[i] = cointime_bin_center;
        asymmetries[i] = asymmetry;
        errors[i] = error;

        // Write to file: W2 low edge, asymmetry, and error
        outfile << cointime_bin_low_edge << "," << asymmetry << "," << error << "\n";
    }
    
    outfile.close();

    // Create a TGraphErrors with the asymmetries and their errors
    TGraphErrors* graph = new TGraphErrors(nBins, &cointime_bin_centers[0], &asymmetries[0], 0, &errors[0]);
    graph->SetTitle("Helicity Asymmetry vs cointime; cointime (ns); Asymmetry (%)");
    
    // Customize the graph appearance if necessary
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    
    // Save the graph into a file or draw it on a canvas
    TCanvas* c = new TCanvas("c", "Asymmetry Graph", 800, 600);
    graph->Draw("AP");  // "A" draws axes, "P" draws points with error bars
    c->SaveAs(Form("plots/%s_asymmetry_graph_cointime_binned_eHCAL_cut_%s.pdf",kin,std::to_string(flag_eHCAL_cut).c_str())); // Save as PDF, or you can choose another format

    // Optionally save graph into a ROOT file
    //TFile* outFile = new TFile("asymmetry_graph.root", "RECREATE");
    //graph->Write();
    //outFile->Close();
}


void Asymmetry_across_cointime(const char* filename, const char* printfilename, const char* kin, bool flag_eHCAL_cut){


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

	double dy_ac_L = getDoubleValue(config,"dy_ac_L");
	double dy_ac_H = getDoubleValue(config,"dy_ac_H");

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
    double vz = 0.0;
    double eHCAL = 0.0;
    double ePS = 0.0;

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("IHWP",&IHWP);
    tree->SetBranchAddress("dx",&dx);
    tree->SetBranchAddress("dy",&dy);
    tree->SetBranchAddress("W2",&W2);
    tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("ntrack", &ntrack);
    tree->SetBranchAddress("vz", &vz);
    tree->SetBranchAddress("eHCAL", &eHCAL);
    tree->SetBranchAddress("ePS",&ePS);
    //TH1D *h_dx = new TH1D("h_dx","dx",200,-10,10);
    //TH1D *h_dy = new TH1D("h_dy","dy",200,-10,10);
    //TH1D *h_W2 = new TH1D("h_W2","W2",250,-4,8);
    //TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",100,dy_L,dy_H,100,-4,3);
   	TH1D *h_coin_time = new TH1D("h_coin_time","cointime",100,0,200);
    TH2D *h_dxdy_cut_W2_coin = new TH2D("h_dxdy_cut_W2_coin","dxdy with W2 and coin cuts",200,-4,4,200,-4,4);

    double N_plus_all = 0.0;
    double N_minus_all = 0.0;
    double accidental_events = 0.0;
    double QE_events = 0.0;

	//bins for cointime
	const double binMin = 0.0;
	const double binMax = 200.0;
	const double binWidth = 20;
	const int nBins = static_cast<int>((binMax-binMin)/binWidth)+1;

	std::vector<TH1D*>Helicity_histograms;
	for (int i = 0; i<nBins;++i){
		double binLowEdge = binMin +i*binWidth;
		double binHighEdge = binLowEdge + binWidth;
		Helicity_histograms.push_back(new TH1D(Form("hist_bin_%d", i),Form("Helicity for Bin[%f,%f)",binLowEdge, binHighEdge),5,-2.5,2.5));
	}


    int nentries = tree->GetEntries();

	std::cout<<"coin_L and H : "<<coin_time_L<<"  "<<coin_time_H<<endl;
	std::cout<<"dx_L and H : "<<dy_L <<" "<<dy_H <<endl;

    //eHCAL cut yes or no
    if (flag_eHCAL_cut == false) eHCAL_L=0.0;

    std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;

	//fill corrected helicity each bin from data
	for (int i = 0; i<nentries; i++){
        tree->GetEntry(i);
        if(lookupValue(HelicityCheck,runnum)==1 and lookupValue(MollerQuality,runnum)==1 and abs(vz)<0.27 and ePS>0.2 and( helicity == -1 or helicity == 1)){
        	if(/*(eHCAL>eHCAL_L)and*/( W2_L<W2 and W2<W2_H ) and (dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H)){
        		//Fill cointime for illustration
                if(eHCAL>eHCAL_L){
                    h_coin_time->Fill(coin_time);
                }
                //to get the asymmetry from accidentals 
                if (((coin_time_ac_L<coin_time) and (coin_time<coin_time_L)) or ((coin_time_H<coin_time) and (coin_time<coin_time_ac_H))){
                    //std::cout<<"here"<<endl;
                    if (-1*IHWP*IHWP_flip*helicity == 1){
                        N_plus_all+=1;
                    }
                    else if(-1*IHWP*IHWP_flip*helicity == -1){
                        N_minus_all+=1;
                    }
                }

                //to get the fraction
                //offset cut , no eHCAL cut for asymmetry cal
                if(eHCAL>eHCAL_L){
                    if (coin_time_offset_L<coin_time and coin_time<coin_time_offset_H){
                        accidental_events+=1;
                    }
                        //QE events in the window
                    if(coin_time_L<coin_time and coin_time<coin_time_H){
                        QE_events+=1;
                    }
                }
                
                //for the Asymmetry variation across cointime distribution
                int binIndex = static_cast<int>((coin_time-binMin)/binWidth);
				if (binIndex >= 0 && binIndex<nBins){
					Helicity_histograms[binIndex]->Fill(-1*IHWP*IHWP_flip*helicity); //check -1
				}
			}

            if ((eHCAL>eHCAL_L) and ( W2_L<W2 and W2<W2_H ) and (coin_time_L<coin_time and coin_time<coin_time_H)){
                h_dxdy_cut_W2_coin->Fill(dy,dx);

            }

        }
        if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
        std::cout.flush();
		
    }

    //fraction and asymmetry
    double f_acc = accidental_events/QE_events;
    double A_acc = (N_plus_all-N_minus_all)/(N_plus_all+N_minus_all);

    //stat error in fraction and asymmetry
    double err_f_acc = sqrt(accidental_events)/QE_events; //assuming poisson stat
    double err_A_acc = sqrt((4 * N_plus_all * N_minus_all) / pow(N_plus_all+N_minus_all,3));

    // Create a box and draw it inside the h_coin_time histogram
    double x1 = coin_time_L+30;  // Left edge of the box (in nS)
    double x2 = coin_time_H+30; // Right edge of the box (in nS)
    double y1 = -10;   // Bottom edge of the box (in counts)
    double y2 = h_coin_time->GetMaximum(); // Top edge of the box, 80% of the histogram's max

    // Create a box
    TBox* box_offset = new TBox(x1, y1, x2, y2);
    box_offset->SetFillColorAlpha(6, 0.3); // Red color with 30% transparency
    box_offset->SetLineColor(6);         // Outline color (optional)

    // Create a box
    TBox* box_anti1 = new TBox(60, y1, coin_time_L, y2);
    box_anti1->SetFillColorAlpha(kRed, 0.3); // Red color with 30% transparency
    box_anti1->SetLineColor(kRed);         // Outline color (optional)

    // Create a box
    TBox* box_anti2 = new TBox(coin_time_H, y1, 170, y2);
    box_anti2->SetFillColorAlpha(kRed, 0.3); // Red color with 30% transparency
    box_anti2->SetLineColor(kRed);         // Outline color (optional)

    //Create a box
    TBox* box_dxdy = new TBox(dy_L,dx_L,dy_H,dx_H);
    box_dxdy->SetFillStyle(0);
    box_dxdy->SetLineColor(kRed);

    TCanvas * ccoin = new TCanvas("ccoin","ccoin",3600,3000);
    ccoin->Divide(2,2);

    ccoin->cd(1);
    h_coin_time->SetXTitle("cointime (ns)");
    h_coin_time->Draw();
    box_offset->Draw("SAME");

    ccoin->cd(2);
    h_coin_time->SetXTitle("cointime (ns)");
    h_coin_time->Draw();
    box_anti1->Draw("SAME");
    box_anti2->Draw("SAME");

    ccoin->cd(3);
    h_dxdy_cut_W2_coin->SetYTitle("dx (m)");
    h_dxdy_cut_W2_coin->SetXTitle("dy (m)");
    h_dxdy_cut_W2_coin->Draw("COLZ");
    box_dxdy->Draw("SAME");

    ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%s.pdf",kin,std::to_string(flag_eHCAL_cut).c_str()));
    ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%s.png",kin,std::to_string(flag_eHCAL_cut).c_str()));
    ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%s.jpg",kin,std::to_string(flag_eHCAL_cut).c_str()));

    std::cout<<"N_plus_all : "<<N_plus_all<<endl;
    std::cout<<"N_minus_all : "<<N_minus_all<<endl;
    std::cout<<"A_acc : "<<A_acc<<"+-"<<err_A_acc<<endl;
    std::cout<<"accidental_events : "<<accidental_events<<endl;
    std::cout<<"QE_events : "<<QE_events<<endl;
    std::cout<<"f_acc : "<<f_acc<<"+-"<<err_f_acc<<endl;

    std::ofstream outfile;
    outfile.open(Form("txt/%s_accidental_asymmetry_eHCAL_cut_%s.txt",kin,std::to_string(flag_eHCAL_cut).c_str()));
    outfile<<"N_plus_all = "<<N_plus_all<<endl;
    outfile<<"N_minus_all = "<<N_minus_all<<endl;
    outfile<<"A_acc = "<<A_acc<<endl;
    outfile<<"err_A_acc = "<<err_A_acc<<endl;
    outfile<<"accidental_events = "<<accidental_events<<endl;
    outfile<<"QE_events = "<<QE_events<<endl;
    outfile<<"f_acc = "<<f_acc<<endl;
    outfile<<"err_f_acc = "<<err_f_acc<<endl;



    //Process Helicity histogram

    CalculateAsymmetry(Helicity_histograms,printfilename,kin,flag_eHCAL_cut);


}