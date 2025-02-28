#pragma once

//#include "parse.h"
#include "cuts.h"
#include "plotdxdy.h"
//#include "parse.h"

void plotW2(const char* filename, const char* printfilename, const char *kin, bool flag_eHCAL_cut){
	
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
	double eHCAL = 0.0;

    tree->SetBranchAddress("dx",&dx);
    tree->SetBranchAddress("dy",&dy);
    tree->SetBranchAddress("W2",&W2);
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("xSH",&xSH);
	tree->SetBranchAddress("ySH",&ySH);
	tree->SetBranchAddress("eHCAL",&eHCAL);

	//bins for the xSH
	const double binMin = -1.1;
	const double binMax = 1.1;
	const double binWidth = 0.1;
	const int nBins = static_cast<int>((binMax-binMin)/binWidth)+1;

	double eHCAL_L = cutsobject.eHCAL_L;
		
	if (flag_eHCAL_cut == false) eHCAL_L=0.0;

    std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;

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

	TFile outFile(Form("rootfiles/%s_histograms.root",printfilename),"RECREATE");

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

    		graphCanvas->Print(Form("plots/%s_W2_QE_peak_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	}
	
	for (auto hist : W2_histograms){
		delete hist;
	}

}

void plotdxdy(const char* filename,const char* printfilename,const char* kin, bool flag_eHCAL_cut){
	
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
	double eHCAL = 0.0;
	int runnum = 0;
	double ntrack_sbs =0.0;

	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);

	TH1D *h_dx = new TH1D("h_dx","dx",200,-10,10);
	TH1D *h_dy = new TH1D("h_dy","dy",200,-10,10);
	TH1D *h_W2 = new TH1D("h_W2","W2",100,-4,8);
	TH1D *h_W2_cut_cointime = new TH1D("h_W2_cut_cointime","W2 with cointime cut",100,-4,8);
	TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",250,-2,2,250,-4,4);
	TH2D *h_dyW2 = new TH2D("h_dyW2","dy v W2",200,-1,3,200,-4,4);
	TH2D *h_dxW2 = new TH2D("h_dxW2","dx v W2",200,-1,3,200,-2,2);
	TH1D *h_coin_time = new TH1D("h_coin_time","coin_time",200,0,200);
	TH1D *h_eHCAL_cut_QE = new TH1D("h_eHCAL_cut_QE","eHCAL with cointime, W2, dx and dy cuts;eHCAL(GeV)", 200, 0, 1.5);
	TH1D *h_eHCAL_cut_cointime_W2 = new TH1D("h_eHCAL_cut_cointime_W2","eHCAL with cointime and W2 cuts;eHCAL(GeV)",200,0,1.5);
	TH1D *h_eHCAL_cut_cointime = new TH1D("h_eHCAL_cut_cointime","eHCAL with cointime cut;eHCAL(GeV)",200,0,1.5);
	TH1D *h_eHCAL = new TH1D("h_eHCAL","eHCAL;eHCAL(GeV)",200,0,1.5);

	TH1D *h_dx_W2_cut = new TH1D("h_dx_W2_cut","dx after W2, cointime and dy cut",200,-10,10);
	TH1D *h_dy_W2_cut = new TH1D("h_dy_W2_cut","dy after W2 and cointime cut",200,-10,10);
	//TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
	TH2D *h_dxdy_W2_cut = new TH2D("h_dxdy_W2_cut","dxdy after W2 and cointime cut",250,-2,2,250,-4,4);
	TH1D *h_coin_time_W2_cut = new TH1D("h_coin_time_W2_cut","cointime after W2 cut",200,0,200);


	double eHCAL_L = cutsobject.eHCAL_L;
		
	if (flag_eHCAL_cut == false) eHCAL_L=0.0;

    std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;
    std::cout<<"run_num_L: "<<cutsobject.run_num_L<<" run_num_H: "<<cutsobject.run_num_H<<endl;
    std::cout<<"coin_time_L: "<<cutsobject.coin_time_L<<" coin_time_H: "<<cutsobject.coin_time_H<<endl;
    std::cout<<"W2_L: "<<cutsobject.W2_L<<" W2_H: "<<cutsobject.W2_H<<endl;
    std::cout<<"dx_L: "<<cutsobject.dx_L<<" dx_H: "<<cutsobject.dx_H<<endl;
    std::cout<<"dy_L: "<<cutsobject.dy_L<<" dy_H: "<<cutsobject.dy_H<<endl;

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries*1; i++){
		tree->GetEntry(i);
		if(cutsobject.run_num_L<runnum and runnum<cutsobject.run_num_H /*and ntrack_sbs>0*/){
			//before adding a cut on W2
			h_dx->Fill(dx);
			h_dy->Fill(dy);
			h_dxdy->Fill(dy,dx);
			h_eHCAL->Fill(eHCAL);
			//h_dyW2->Fill(W2,dy);
			//h_dxW2->Fill(W2,dx);
			if (eHCAL>eHCAL_L){
				h_W2->Fill(W2);
				h_coin_time->Fill(coin_time);
	        }

	        if (eHCAL>eHCAL_L and cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H /*and dy<cutsobject.dy_H and cutsobject.dy_L<dy*/){
	        	h_W2_cut_cointime->Fill(W2);
	        }


			if (eHCAL>eHCAL_L and cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H ){
				h_dyW2->Fill(W2,dy);
				h_dxW2->Fill(W2,dx);
				h_eHCAL_cut_cointime->Fill(eHCAL);
			}

			//add a cut on W2 and coin time
			if (eHCAL>eHCAL_L and (cutsobject.W2_L<W2 and W2<cutsobject.W2_H) and (cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H) ){
				if (dy<cutsobject.dy_H and cutsobject.dy_L<dy){
					h_dx_W2_cut->Fill(dx);
					//h_dxdy_W2_cut->Fill(dy,dx);
					if (dx<cutsobject.dx_H and cutsobject.dx_L<dx){
						h_eHCAL_cut_QE->Fill(eHCAL);
					}
				}

				h_dy_W2_cut->Fill(dy);
				h_dxdy_W2_cut->Fill(dy,dx);
				h_eHCAL_cut_cointime_W2->Fill(eHCAL);
				//h_dyW2->Fill(W2,dy);
				
			}	
			if(eHCAL>eHCAL_L and (cutsobject.W2_L<W2 and W2<cutsobject.W2_H) /*and dy<cutsobject.dy_H and cutsobject.dy_L<dy*/){
				h_coin_time_W2_cut->Fill(coin_time);
			}	
		}

	    if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();


	}
	
	//oval cut
	//TCutG *cutg = CreateOvalCut("cutg",N_cut_cx,N_cut_cy,N_cut_rx,N_cut_ry,100);

	//square cut
	//TCutG *cutsq = CreateSquareCut(cutsobject.dy_L,cutsobject.dx_L,cutsobject.dy_H,cutsobject.dy_H);

	//add lines to W2 plot to show the cut
	
	double W2_Ymax = h_W2_cut_cointime->GetMaximum();
	double cointime_Ymax = h_coin_time_W2_cut->GetMaximum();

	TLine *line1 = new TLine(cutsobject.W2_L,0.0,cutsobject.W2_L,W2_Ymax);
	line1->SetLineColor(kBlue);
	line1->SetLineWidth(2);
	
	TLine *line2 = new TLine(cutsobject.W2_H,0.0,cutsobject.W2_H,W2_Ymax);
	line2->SetLineColor(kBlue);
	line2->SetLineWidth(2);
	
	TLine *line3 = new TLine(cutsobject.coin_time_L,0.0,cutsobject.coin_time_L,cointime_Ymax);
	line3->SetLineColor(kBlue);
	line3->SetLineWidth(2);
	
	TLine *line4 = new TLine(cutsobject.coin_time_H,0.0,cutsobject.coin_time_H,cointime_Ymax);
	line4->SetLineColor(kBlue);
	line4->SetLineWidth(2);

	TLine *line5 = new TLine(cutsobject.eHCAL_L,0.0,cutsobject.eHCAL_L,h_eHCAL_cut_cointime->GetMaximum());
	line5->SetLineColor(kBlue);
	line5->SetLineColor(2);

	TBox* box_dxdy_n = new TBox(cutsobject.dy_L,cutsobject.dx_L,cutsobject.dy_H,cutsobject.dx_H);
    box_dxdy_n->SetFillStyle(0);
    box_dxdy_n->SetLineColor(kRed);
    box_dxdy_n->SetLineWidth(3);

	TCanvas* c = new TCanvas("c","c",3600,3000);
	TCanvas* c1 = new TCanvas("c1","c1",3600,3000);
	TCanvas* c2 = new TCanvas("c2","c2",3600,3000);
	TCanvas* c3 = new TCanvas("c3","c3",3600,3000);
	TCanvas* c4 = new TCanvas("c4","c4",3600,3000);

	c->Divide(3,2);
	c->cd(1);
	h_dx->SetXTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dx->Draw();
	c->cd(2);
	h_dy->Draw();
	h_dy->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	c->cd(3);
	gStyle->SetPalette(kRainBow);
	h_dxdy->Draw("COLZ");
	h_dxdy->SetYTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dxdy->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	c->cd(4);
	h_W2->Draw();
	c->cd(5);
	h_coin_time->Draw();
	h_coin_time->SetXTitle("cointime(ns)");
	c->cd(6);
	h_W2_cut_cointime->Draw();
	h_W2_cut_cointime->SetXTitle("W2(GeV2)");
	c->Print(Form("plots/%s_dxdyplots_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	c1->Divide(3,2);
	c1->cd(1);
	h_dx_W2_cut->Draw();
	h_dx_W2_cut->SetXTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	c1->cd(2);
	h_dy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	h_dy_W2_cut->Draw();
	c1->cd(3);
	gStyle->SetPalette(kRainBow);
	h_dxdy_W2_cut->Draw("COLZ");
	h_dxdy_W2_cut->SetYTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dxdy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	box_dxdy_n->Draw();
	c1->cd(4);
	h_W2->Draw();
	h_W2->SetXTitle("W2(GeV2)");
	line1->Draw("same");
	line2->Draw("same");
	c1->cd(5);
	h_coin_time_W2_cut->Draw();
	h_coin_time_W2_cut->SetXTitle("cointime(ns)");
	line3->Draw("same");
	line4->Draw("same");
	c1->cd(6);
	h_W2_cut_cointime->Draw();
	h_W2_cut_cointime->SetXTitle("W2(GeV2)");
	line1->Draw("same");
	line2->Draw("same");

	c2->Divide(2,2);
	c2->cd(1);
	h_dyW2->Draw();
	c2->cd(2);
	h_dxW2->Draw();

	c3->Divide(2,2);
	c3->cd(1);
	h_dxdy_W2_cut->Draw("COLZ");
	h_dxdy_W2_cut->SetYTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dxdy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	//cutsq->Draw("L");

	c4->Divide(2,2);
	c4->cd(1);
	h_eHCAL->Draw();
	line5->Draw();
	c4->cd(2);
	h_eHCAL_cut_cointime->Draw();
	line5->Draw();
	c4->cd(3);
	h_eHCAL_cut_cointime_W2->Draw();
	c4->cd(4);
	h_eHCAL_cut_QE->Draw();

	c1->Print(Form("plots/%s_dxdyplots_withW2Cut_eHCAL_cut_%s.pdf(",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c2->Print(Form("plots/%s_dxdyplots_withW2Cut_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c4->Print(Form("plots/%s_dxdyplots_withW2Cut_eHCAL_cut_%s.pdf)",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c3->SaveAs(Form("plots/%s_dxdyplots_withW2Cut_eHCAL_cut_%s.png",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c1->SaveAs(Form("plots/%s_dxdyplots_withW2Cut_eHCAL_cut_%s.png",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
}
