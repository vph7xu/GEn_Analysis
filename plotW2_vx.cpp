void plotW2_vx(const char* filename, const char* printfilename, const char *kin/*, const char *var*/){
	TFile* file = TFile::Open(filename);
        TTree* tree = (TTree*)file->Get("Tout");

        double dx = 0.0;
        double dy = 0.0;
        double W2 = 0.0;
        double Q2 = 0.0;
        double coin_time = 0.0;
	double xSH = 0.0;
	double ySH = 0.0;
	double BPMAx = 0.0;
	double BPMAy = 0.0;
	double vx = 0.0;
	double vy = 0.0;


        tree->SetBranchAddress("dx",&dx);
        tree->SetBranchAddress("dy",&dy);
        tree->SetBranchAddress("W2",&W2);
        tree->SetBranchAddress("Q2",&Q2);
        tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("xSH",&xSH);
	tree->SetBranchAddress("ySH",&ySH);
	tree->SetBranchAddress("BPMAx",&BPMAx);
	tree->SetBranchAddress("BPMAy",&BPMAy);
	tree->SetBranchAddress("vx",&vx);
	tree->SetBranchAddress("vy",&vy);

	//bins for the xSH
	const double binMin = -0.0035;
	const double binMax = 0.003;
	const double binWidth = 0.0001;
	const int nBins = static_cast<int>((binMax-binMin)/binWidth)+1;

	std::vector<TH1D*>W2_histograms;
	for (int i = 0; i<nBins; ++i){
		double binLowEdge = binMin +i*binWidth;
		double binHighEdge = binLowEdge + binWidth;
		W2_histograms.push_back(new TH1D(Form("hist_bin_%d", i),Form("W2 for Bin[%f,%f)",binLowEdge, binHighEdge),1000,-4,8));
	}
	
	//double Var = 0.0;

	//if (var == "BPMAx") tree->SetBranchAddress("Var",&BPMAx);
       	//else if(var == "BPMAy") tree->SetBranchAddress("Var",&BPMAy);	

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; ++i){
		tree->GetEntry(i);
		int binIndex = static_cast<int>((vy-binMin)/binWidth);
		if (binIndex >= 0 && binIndex<nBins){
			W2_histograms[binIndex]->Fill(W2);
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
		if (canvasIndex==0) canvas->Print(Form("%s_W2.pdf(",printfilename));
		else if (canvasIndex==(canvasesNeeded-1)) canvas->Print(Form("%s_W2.pdf)",printfilename));
		else canvas->Print(Form("%s_W2.pdf",printfilename));

		TGraphErrors *meanGraph = new TGraphErrors(nBins, &binCenters[0], &means[0], nullptr, &meanErrors[0]);
    		TGraphErrors *sigmaGraph = new TGraphErrors(nBins, &binCenters[0], &sigmas[0], nullptr, &sigmaErrors[0]);

    		// Customize and draw the graphs
    		TCanvas *graphCanvas = new TCanvas("graphCanvas", "Fit Parameters", 1200, 600);
    		graphCanvas->Divide(2, 1);

    		graphCanvas->cd(1);
    		meanGraph->SetTitle("Mean of Gaussian Fits");
    		meanGraph->GetXaxis()->SetTitle("Bin Center");
		meanGraph->GetXaxis()->SetLimits(-0.003,0.002);
    		meanGraph->GetYaxis()->SetTitle("Mean");
		meanGraph->GetYaxis()->SetRangeUser(0.86,0.91);
   	 	meanGraph->SetMarkerStyle(7);
    		meanGraph->Draw("AP");

    		graphCanvas->cd(2);
    		sigmaGraph->SetTitle("Sigma of Gaussian Fits");
    		sigmaGraph->GetXaxis()->SetTitle("Bin Center");
		sigmaGraph->GetXaxis()->SetLimits(-0.0031,0.002);
    		sigmaGraph->GetYaxis()->SetTitle("Sigma");
		sigmaGraph->GetYaxis()->SetRangeUser(0.13,0.18);
    		sigmaGraph->SetMarkerStyle(7);
    		sigmaGraph->Draw("AP");

    		graphCanvas->Print(Form("%s_W2_QE_peak.pdf",printfilename));

	}
	
	for (auto hist : W2_histograms){
		delete hist;
	}

}

