#pragma once

//#include "parse.h"
#include "cuts.h"
#include "plotdxdy.h"

void plotW2(const char* filename, const char* printfilename, const char *kin){
	
	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts

	cuts cutsobject;
        cutsobject.parsecuts(config);

        //double coin_time_L = 90;//getDoubleValue(config,"coin_time_L");
        //double coin_time_H = 115;//getDoubleValue(config,"coin_time_H");


	TFile* file = TFile::Open(filename);
        TTree* tree = (TTree*)file->Get("Tout");

        double dx = 0.0;
        double dy = 0.0;
        double W2 = 0.0;
        double Q2 = 0.0;
        double coin_time = 0.0;
	double xSH = 0.0;
	double ySH = 0.0;

        tree->SetBranchAddress("dx",&dx);
        tree->SetBranchAddress("dy",&dy);
        tree->SetBranchAddress("W2",&W2);
        tree->SetBranchAddress("Q2",&Q2);
        tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("xSH",&xSH);
	tree->SetBranchAddress("ySH",&ySH);

	//bins for the xSH
	const double binMin = -1.1;
	const double binMax = 1.1;
	const double binWidth = 0.1;
	const int nBins = static_cast<int>((binMax-binMin)/binWidth)+1;

	std::vector<TH1D*>W2_histograms;
	for (int i = 0; i<nBins; ++i){
		double binLowEdge = binMin +i*binWidth;
		double binHighEdge = binLowEdge + binWidth;
		W2_histograms.push_back(new TH1D(Form("hist_bin_%d", i),Form("W2 for Bin[%f,%f)",binLowEdge, binHighEdge),1000,-4,8));
	}

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; ++i){
		tree->GetEntry(i);
		if (coin_time<cutsobject.coin_time_L or coin_time>cutsobject.coin_time_H){
			int binIndex = static_cast<int>((xSH-binMin)/binWidth);
			if (binIndex >= 0 && binIndex<nBins){
				W2_histograms[binIndex]->Fill(W2);
			}
		}
                if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();
		
	}

	TFile outFile(Form("%s_histograms.root",printfilename),"RECREATE");

	for (auto hist : W2_histograms){
		hist->Write();
	}

	outFile.Close();

    	std::vector<double> binCenters, means, meanErrors, sigmas, sigmaErrors;

	const int canvasesNeeded = (nBins+15)/16;

	for (int canvasIndex = 0; canvasIndex < canvasesNeeded; ++canvasIndex){
		TCanvas *canvas = new TCanvas(Form("canvas_%d",canvasIndex),"Histograms",1200,1200);
		canvas->Divide(4,4);

		for (int histIndex = 0; histIndex < 16; ++histIndex){
			int overallIndex = canvasIndex * 16 + histIndex;
			if (overallIndex < nBins){
				canvas->cd(histIndex+1);
				W2_histograms[overallIndex]->Draw();

				TF1 *gausFit = new TF1("gausFit","gaus",-2,1.2);
				W2_histograms[overallIndex]->Fit(gausFit,"RQ");

				gausFit->Draw("SAME");

                		gPad->Update(); // Update the pad to make sure stats box is created
		               
			       	// Get fit parameters
                		double mean = gausFit->GetParameter(1);
                		double sigma = gausFit->GetParameter(2);
				double meanError = gausFit->GetParError(1);
				double sigmaError = gausFit->GetParError(2);

				//store fit parameter
				double binCenter = (binMin + overallIndex * binWidth) + binWidth / 2;
               			binCenters.push_back(binCenter);
                		means.push_back(mean);
                		meanErrors.push_back(meanError);
                		sigmas.push_back(sigma);
                		sigmaErrors.push_back(sigmaError);
                
               		 	// Draw the fit statistics on the canvas
                		TLatex latex;
                		latex.SetNDC();
                		latex.SetTextSize(0.05);
                		latex.SetTextColor(kRed);
                		latex.DrawLatex(0.6, 0.7, Form("Mean = %.2f", mean));
                		latex.DrawLatex(0.6, 0.6, Form("Sigma = %.2f", sigma));	
		
			}
		}
		if (canvasIndex==0) canvas->Print(Form("%s_W2_bigbite_acceptance.pdf(",printfilename));
		else if (canvasIndex==(canvasesNeeded-1)) canvas->Print(Form("%s_W2_bigbite_acceptance.pdf)",printfilename));
		else canvas->Print(Form("%s_W2_bigbite_acceptance.pdf",printfilename));

		TGraphErrors *meanGraph = new TGraphErrors(nBins, &binCenters[0], &means[0], nullptr, &meanErrors[0]);
    		TGraphErrors *sigmaGraph = new TGraphErrors(nBins, &binCenters[0], &sigmas[0], nullptr, &sigmaErrors[0]);

    		// Customize and draw the graphs
    		TCanvas *graphCanvas = new TCanvas("graphCanvas", "Fit Parameters", 1200, 600);
    		graphCanvas->Divide(2, 1);

    		graphCanvas->cd(1);
    		meanGraph->SetTitle("Mean of Gaussian Fits");
    		meanGraph->GetXaxis()->SetTitle("Bin Center");
    		meanGraph->GetYaxis()->SetTitle("Mean");
		meanGraph->GetYaxis()->SetRangeUser(0.8,1.2);
   	 	meanGraph->SetMarkerStyle(7);
    		meanGraph->Draw("AP");

    		graphCanvas->cd(2);
    		sigmaGraph->SetTitle("Sigma of Gaussian Fits");
    		sigmaGraph->GetXaxis()->SetTitle("Bin Center");
    		sigmaGraph->GetYaxis()->SetTitle("Sigma");
		sigmaGraph->GetYaxis()->SetRangeUser(0.1,0.4);
    		sigmaGraph->SetMarkerStyle(7);
    		sigmaGraph->Draw("AP");

    		graphCanvas->Print(Form("%s_W2_QE_peak.pdf",printfilename));

	}
	
	for (auto hist : W2_histograms){
		delete hist;
	}

}

void plotdxdy(const char* filename,const char* printfilename,const char* kin){
	
	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	cuts cutsobject;	
	cutsobject.parsecuts(config);

	std::cout<<"coin_time : "<<cutsobject.coin_time_L<<endl;
	
	TFile* file = TFile::Open(filename);
	TTree* tree = (TTree*)file->Get("Tout");

	double dx = 0.0; 
	double dy = 0.0; 
	double W2 = 0.0;
    double Q2 = 0.0;
	double coin_time = 0.0;

	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);

	TH1D *h_dx = new TH1D("h_dx","dx",1000,-10,10);
	TH1D *h_dy = new TH1D("h_dy","dy",1000,-10,10);
	TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
	TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",250,-2,2,250,-4,4);
	TH2D *h_dyW2 = new TH2D("h_dyW2","dy v W2",400,-1,3,800,-4,4);
	TH2D *h_dxW2 = new TH2D("h_dxW2","dx v W2",400,-1,3,400,-2,2);
	TH1D *h_coin_time = new TH1D("h_coin_time","coin_time",2000,0,200);


	TH1D *h_dx_W2_cut = new TH1D("h_dx_W2_cut","dx after W2 and Coin cut",1000,-10,10);
	TH1D *h_dy_W2_cut = new TH1D("h_dy_W2_cut","dy after W2 and Coin cut",1000,-10,10);
	//TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
	TH2D *h_dxdy_W2_cut = new TH2D("h_dxdy_W2_cut","dxdy after W2 and Coin cut",250,-2,2,250,-4,4);
	TH1D *h_coin_time_W2_cut = new TH1D("h_coin_time_W2_cut","coin_time after W2 cut",2000,0,200);
		
	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		
		//before adding a cut on W2
		h_dx->Fill(dx);
		h_dy->Fill(dy);
		h_dxdy->Fill(dy,dx);
		//h_dyW2->Fill(W2,dy);
		//h_dxW2->Fill(W2,dx);
		h_W2->Fill(W2);
		h_coin_time->Fill(coin_time);
               
		if (cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H){
			h_dyW2->Fill(W2,dy);
			h_dxW2->Fill(W2,dx);
		}

		//add a cut on W2 and coin time
		if ((cutsobject.W2_L<W2 and W2<cutsobject.W2_H) and (cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H)){
			h_dx_W2_cut->Fill(dx);
			h_dy_W2_cut->Fill(dy);
			h_dxdy_W2_cut->Fill(dy,dx);
			//h_dyW2->Fill(W2,dy);
			
		}	
		if((cutsobject.W2_L<W2 and W2<cutsobject.W2_H)){
			h_coin_time_W2_cut->Fill(coin_time);
		}	

	        if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();


	}
	
	//oval cut
	//TCutG *cutg = CreateOvalCut("cutg",N_cut_cx,N_cut_cy,N_cut_rx,N_cut_ry,100);

	//square cut
	TCutG *cutsq = CreateSquareCut(cutsobject.dy_L,cutsobject.dx_L,cutsobject.dy_H,cutsobject.dy_H);

	//add lines to W2 plot to show the cut
	
	double W2_Ymax = h_W2->GetMaximum();

	TLine *line1 = new TLine(cutsobject.W2_L,0.0,cutsobject.W2_L,W2_Ymax);
	line1->SetLineColor(kBlue);
	line1->SetLineWidth(2);
	
	TLine *line2 = new TLine(cutsobject.W2_H,0.0,cutsobject.W2_H,W2_Ymax);
	line2->SetLineColor(kBlue);
	line2->SetLineWidth(2);
	
	TCanvas* c = new TCanvas("c","c",1200,1600);
	TCanvas* c1 = new TCanvas("c1","c1",1200,1600);
	TCanvas* c2 = new TCanvas("c2","c2",1200,1600);
	TCanvas* c3 = new TCanvas("c3","c3",1800,2400);

	c->Divide(3,2);
	c->cd(1);
	h_dx->SetXTitle("HCAL_X(exp)-HCAL_X(act)");
	h_dx->Draw();
	c->cd(2);
	h_dy->Draw();
	h_dy->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	c->cd(3);
	gStyle->SetPalette(kRainBow);
	h_dxdy->Draw("COLZ");
	h_dxdy->SetYTitle("HCAL_X(exp)-HCAL_X(act)");
	h_dxdy->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	c->cd(4);
	h_W2->Draw();
	c->cd(5);
	h_coin_time->Draw();
	c->Print(Form("../plots/%s_dxdyplots.pdf",printfilename));


	c1->Divide(3,2);
	c1->cd(1);
	h_dx_W2_cut->Draw();
	h_dx_W2_cut->SetXTitle("HCAL_X(exp)-HCAL_X(act)");
	c1->cd(2);
	h_dy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	h_dy_W2_cut->Draw();
	c1->cd(3);
	gStyle->SetPalette(kRainBow);
	h_dxdy_W2_cut->Draw("COLZ");
	//cutg->Draw("L");
	//cutsq->Draw("L");
	h_dxdy_W2_cut->SetYTitle("HCAL_X(exp)-HCAL_X(act)");
	h_dxdy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	c->cd(4);
	c1->cd(4);
	h_W2->Draw();
	line1->Draw("same");
	line2->Draw("same");
	//c1->cd(5);
	//h_dyW2->Draw();
	c1->cd(5);
	h_coin_time_W2_cut->Draw();

	c2->Divide(2,1);
	c2->cd(1);
	h_dyW2->Draw();
	c2->cd(2);
	h_dxW2->Draw();

	c3->Divide(2,1);
	c3->cd(1);
	h_dxdy_W2_cut->Draw("COLZ");
	h_dxdy_W2_cut->SetYTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dxdy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	cutsq->Draw("L");

	c1->Print(Form("../plots/%s_dxdyplots_withW2Cut.pdf(",printfilename));
	c2->Print(Form("../plots/%s_dxdyplots_withW2Cut.pdf)",printfilename));
	c3->SaveAs(Form("../plots/%s_dxdyplots_withW2Cut.png",printfilename));
}
