#include "models.h"
#include "parse.h"
#include "plotdxdy.h"
#include <locale>
#include <TStyle.h>
#include <cmath>

void CalculateAsymmetry(std::vector<TH1D*>& Helicity_histograms, const char* printfilename) {
    // Prepare arrays to store the W2 bin centers, asymmetries, and errors
    int nBins = Helicity_histograms.size();
    std::vector<double> W2_bin_centers(nBins);
    std::vector<double> asymmetries(nBins);
    std::vector<double> errors(nBins);
    
    // Open a file to save results
    std::ofstream outfile(printfilename);
    outfile << "W2_bin_low_edge, Asymmetry, Error\n";

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
            error = 100*sqrt((4 * N_plus * N_minus) / (N_total * N_total * N_total));
        }

        // Get the bin center (since we care about the middle of each W2 bin)
        double W2_bin_low_edge = -4.0;
        double W2_bin_width = 0.25;
        double W2_bin_center = i*W2_bin_width + W2_bin_low_edge + 0.5 * W2_bin_width;

        // Store results in arrays
        W2_bin_centers[i] = W2_bin_center;
        asymmetries[i] = asymmetry;
        errors[i] = error;

        // Write to file: W2 low edge, asymmetry, and error
        outfile << W2_bin_low_edge << "," << asymmetry << "," << error << "\n";
    }
    
    outfile.close();

    // Create a TGraphErrors with the asymmetries and their errors
    TGraphErrors* graph = new TGraphErrors(nBins, &W2_bin_centers[0], &asymmetries[0], 0, &errors[0]);
    graph->SetTitle("Helicity Asymmetry vs W2; W2; Asymmetry(%)");
    
    // Customize the graph appearance if necessary
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);
    
    // Save the graph into a file or draw it on a canvas
    TCanvas* c = new TCanvas("c", "Asymmetry Graph", 800, 600);
    graph->Draw("AP");  // "A" draws axes, "P" draws points with error bars
    c->SaveAs("plots/asymmetry_graph_W2binned.pdf"); // Save as PDF, or you can choose another format

    // Optionally save graph into a ROOT file
    //TFile* outFile = new TFile("asymmetry_graph.root", "RECREATE");
    //graph->Write();
    //outFile->Close();
}


TGraphErrors* graphAsymmetry(std::vector<TH1D*> helicity_histograms, double binMin, double binWidth, int nBins) {
    std::vector<double> W2_bin_centers;
    std::vector<double> Asymmetry;
    std::vector<double> AsymmetryErrors;

    for (int i = 0; i < nBins; ++i) {
        double binLowEdge = binMin + i * binWidth;
        double binCenter = binLowEdge + binWidth / 2.0;
        W2_bin_centers.push_back(binCenter);

        int bin_plus = helicity_histograms[i]->FindBin(1.0);   // Helicity +1
        int bin_minus = helicity_histograms[i]->FindBin(-1.0); // Helicity -1

        double N_plus = helicity_histograms[i]->GetBinContent(bin_plus); // helicity +1
        double N_minus = helicity_histograms[i]->GetBinContent(bin_minus); // helicity -1

        if (N_plus + N_minus > 0) {
            // Calculate asymmetry
            double A = 100*(N_plus - N_minus) / (N_plus + N_minus);
            Asymmetry.push_back(A);

            // Error calculation using binomial error propagation
            double A_error = 2.0 *100* sqrt((N_plus * N_minus) / pow((N_plus + N_minus), 3));
            AsymmetryErrors.push_back(A_error);
        } else {
            Asymmetry.push_back(0);  // Handle bins with no entries
            AsymmetryErrors.push_back(0);
        }
    }

    // Create TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nBins, &W2_bin_centers[0], &Asymmetry[0], 0, &AsymmetryErrors[0]);
    return graph;
}

void Asymmetry_across_W2(const char* filename, const char* printfilename, const char* kin) {


	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	std::map<int, int> HelicityCheck = readCSVToMap("DB/Helicity_quality.csv");
	std::map<int, int> MollerQuality = readCSVToMap("DB/Moller_quality.csv");

	double coin_time_L = getDoubleValue(config,"coin_time_L");
	double coin_time_H = getDoubleValue(config,"coin_time_H");

	double W2_L = getDoubleValue(config,"W2_L");
	double W2_H = getDoubleValue(config,"W2_H");

	double dy_L = getDoubleValue(config,"dy_L");
	double dy_H = getDoubleValue(config,"dy_H");

	double dx_L = getDoubleValue(config,"dx_L");
	double dx_H = getDoubleValue(config,"dx_H");

	double dy_ac_L = getDoubleValue(config,"dy_ac_L");
	double dy_ac_H = getDoubleValue(config,"dy_ac_H");

	double run_num_L = getDoubleValue(config,"run_num_L");
	double run_num_H = getDoubleValue(config,"run_num_H");

	double IHWP_flip = getDoubleValue(config,"IHWP_flip");

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

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("IHWP",&IHWP);
    tree->SetBranchAddress("dx",&dx);
    tree->SetBranchAddress("dy",&dy);
    tree->SetBranchAddress("W2",&W2);
    tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("ntrack", &ntrack);

    //TH1D *h_dx = new TH1D("h_dx","dx",200,-10,10);
    //TH1D *h_dy = new TH1D("h_dy","dy",200,-10,10);
    TH1D *h_W2_dx_antidy = new TH1D("h_W2_dx_antidy","W^{2} distribution",60,-4,8);
    TH1D *h_W2_dx = new TH1D("h_W2_dx","W^{2} distribution",60,-4,8);
    TH1D *h_W2_dx_dy = new TH1D("h_W2_dx_dy","W^{2} distribution",60,-4,8);

    //TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",100,dy_L,dy_H,100,-4,3);
   	//TH1D *h_coin_time = new TH1D("h_coin_time","coin_time",1000,40,200);


	//bins for W2
	const double binMin = -4.0;
	const double binMax = 7.0;
	const double binWidth = 0.25;
	const int nBins = static_cast<int>((binMax-binMin)/binWidth)+1;

    //for asymmetry calculations
    /*double N_plus_all = 0.0;
    double N_minus_all = 0.0;
    double inelastic_events = 0.0;
    double QE_events = 0.0;
    */

	std::vector<TH1D*>Helicity_histograms;
    std::vector<TH1D*>Helicity_histograms_1;  
	for (int i = 0; i<nBins;++i){
		double binLowEdge = binMin +i*binWidth;
		double binHighEdge = binLowEdge + binWidth;
		Helicity_histograms.push_back(new TH1D(Form("hist_bin_%d", i),Form("Helicity for Bin[%f,%f)",binLowEdge, binHighEdge),5,-2.5,2.5));
        Helicity_histograms_1.push_back(new TH1D(Form("hist_bin_1_%d", i),Form("Helicity for Bin[%f,%f)",binLowEdge, binHighEdge),5,-2.5,2.5));
	}

    int nentries = tree->GetEntries();

	std::cout<<"coin_L and H : "<<coin_time_L<<"  "<<coin_time_H<<endl;
	std::cout<<"dx_L and H : "<<dy_L <<" "<<dy_H <<endl;

	//fill corrected helicity each bin from data
	for (int i = 0; i<nentries; i++){
        tree->GetEntry(i);
        if(lookupValue(HelicityCheck,runnum)==1 and lookupValue(MollerQuality,runnum)==1){


            /*if((coin_time_L<coin_time && coin_time<coin_time_H) && (dx_L<dx && dx<dx_H) && (dy_ac_L>dy || dy>dy_ac_H)){
                if (-1*IHWP*IHWP_flip*helicity == 1){
                    N_plus_all+=1;
                }
                else if(-1*IHWP*IHWP_flip*helicity == -1){
                    N_minus_all+=1;
                }
            }*/

            //for plots
        	if((coin_time_L<coin_time && coin_time<coin_time_H) && (dx_L<dx && dx<dx_H) && (dy_ac_L>dy || dy>dy_ac_H)){
        		h_W2_dx_antidy->Fill(W2);
                int binIndex = static_cast<int>((W2-binMin)/binWidth);
				if (binIndex >= 0 && binIndex<nBins){
					Helicity_histograms[binIndex]->Fill(-1*IHWP*IHWP_flip*helicity);//check -1
				}
			}
            if((coin_time_L<coin_time && coin_time<coin_time_H) && (dx_L<dx && dx<dx_H)){
                h_W2_dx->Fill(W2);
                int binIndex_1 = static_cast<int>((W2-binMin)/binWidth);
                if (binIndex_1 >= 0 && binIndex_1<nBins){
                    Helicity_histograms_1[binIndex_1]->Fill(-1*IHWP*IHWP_flip*helicity);//check -1
                }
            }
            if((coin_time_L<coin_time && coin_time<coin_time_H) && (dx_L<dx && dx<dx_H) && (dy_L<dy && dy<dy_H)){
                h_W2_dx_dy->Fill(W2);
            }

        }

        if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
        std::cout.flush();
		
    }

    h_W2_dx_antidy->SetLineColor(kBlue);
    h_W2_dx->SetLineColor(kBlack);
    h_W2_dx_dy->SetLineColor(kGreen);

    TCanvas *cW2 = new TCanvas("cW2","cW2",3600,3000);

    TLegend* legendW2 = new TLegend(0.1, 0.9, 0.5, 0.8);
    legendW2->AddEntry(h_W2_dx_antidy, "|dx|<0.4m and (dy<-1.5m or dy>1.0m)", "l");
    legendW2->AddEntry(h_W2_dx, "|dx|<0.4m", "l");
    legendW2->AddEntry(h_W2_dx_dy, "|dx|<0.4m and |dy|<0.4m ", "l");
 

    cW2->Divide(2,2);
    cW2->cd(1);
    h_W2_dx_antidy->Draw("HIST");
    h_W2_dx_antidy->SetXTitle("W^{2} (GeV^{2})");
    h_W2_dx->Draw("HIST SAME");
    h_W2_dx_dy->Draw("HIST SAME");
    legendW2->Draw();

    cW2->SaveAs("plots/W2_for_inelastics.pdf"); 

    //Process Helicity histograms

    //CalculateAsymmetry(Helicity_histograms,printfilename);
    TGraphErrors* g_asym_1 = graphAsymmetry(Helicity_histograms, binMin, binWidth, nBins);
    TGraphErrors* g_asym_2 = graphAsymmetry(Helicity_histograms_1, binMin, binWidth, nBins);

    // Set the style for both graphs
    g_asym_1->SetTitle("Asymmetry vs W^{2}");
    g_asym_1->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
    g_asym_1->GetYaxis()->SetTitle("Asymmetry (%)");

    g_asym_1->SetMinimum(-8);
    g_asym_1->SetMaximum(8);

    g_asym_1->SetMarkerStyle(21);
    g_asym_1->SetMarkerColor(kBlue);
    g_asym_1->SetLineColor(kBlue);
    g_asym_1->SetMarkerSize(5);

    g_asym_2->SetMarkerStyle(21);
    g_asym_2->SetMarkerColor(kBlack);
    g_asym_2->SetLineColor(kBlack);
    g_asym_2->SetMarkerSize(5);

    // Create a canvas and plot both graphs
    TCanvas* c = new TCanvas("c", "Asymmetry vs W^{2}", 3600, 3000);
    g_asym_1->Draw("AP");  // A for axis, P for points
    g_asym_2->Draw("P SAME");  // SAME to draw on the same canvas
    
    // Add legend to distinguish between the two graphs
    TLegend* legend = new TLegend(0.5, 0.8, 0.9, 0.9);
    legend->AddEntry(g_asym_1, "|dx|<0.4m and (dy<-1.5m or dy>1.0m)", "lp");
    legend->AddEntry(g_asym_2, "|dx|<0.4m", "lp");
    legend->Draw();

    // Save the canvas
    c->SaveAs(Form("plots/%s_asymmetry_across_W2.pdf",printfilename));

}