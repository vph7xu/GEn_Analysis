#pragma once

//#include "parse.h"
#include "cuts.h"
#include "plotdxdy.h"
//#include "parse.h"


// ----------------------------------------------------------------------------
//  1.  Helper: fit a TH1D in μ̂ ± n·σ̂ and return the TF1* (drawn on pad)
// ----------------------------------------------------------------------------
TF1* FitGaussNSigma(TH1D* h,
                    const char* tag,
                    double nSigma = 3,
                    bool    useMode  = false)          // set 4 for ±4 σ
{
    if (!h || h->GetEntries() < 50) return nullptr;      // low‑stat guard

    //const double mu    = h->GetMean();
    const double sigma = h->GetStdDev();

    double center;
    if (useMode) {
        int    imax   = h->GetMaximumBin();
        center        = h->GetXaxis()->GetBinCenter(imax);
    } else {
        center        = h->GetMean();
    }

    if (sigma <= 0) return nullptr;

    const double xmin  = std::max(center - nSigma * sigma, h->GetXaxis()->GetXmin());
    const double xmax  = std::min(center + nSigma * sigma, h->GetXaxis()->GetXmax());

    TF1 *g = new TF1(Form("g_%s",tag),"gaus", xmin, xmax);
    g->SetParameters(h->GetBinContent(h->GetMaximumBin()), center, sigma);
    g->SetLineColor(kRed);            // visible on the pad
    g->SetLineWidth(2);

    h->Fit(g, "RQ");                  // R = draw; Q = quiet console
    // Now g is drawn *and* lives in memory for later use
    return g;
}


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

void plotdxdy(const char* filename,const char* printfilename,const char* kin, bool flag_eHCAL_cut, bool veto, bool ETOF){
	
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
	double ntrack_sbs = 0.0;
	double ntrack = 0.0;
	double ePS = 0.0;
	double eSH = 0.0;
	double trP =0.0;
	double RFtime = 0.0;
	double RFtime_sbs = 0.0;
	double hcal_tdctime_tw = 0.0;
	double hcal_tdctime = 0.0;
	double trigbits = 0.0;
	double hodo_time_vec[100];
	double hodo_bar_id_vec[100];
	double hcal_clus_blk_id;
	double hcal_clus_blk_e;
	double hodo_etof_vec[100];
	int nhodo_clus = 0;

	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
	tree->SetBranchAddress("ntrack",&ntrack);
	tree->SetBranchAddress("ePS",&ePS);
	tree->SetBranchAddress("eSH",&eSH);
	tree->SetBranchAddress("trP",&trP);
	tree->SetBranchAddress("RFtime", &RFtime);
	tree->SetBranchAddress("RFtime_sbs", &RFtime_sbs);
	tree->SetBranchAddress("hcal_tdctime_tw", &hcal_tdctime_tw);
	tree->SetBranchAddress("hcal_tdctime",&hcal_tdctime);
	tree->SetBranchAddress("trigbits",&trigbits);
	tree->SetBranchAddress("hodo_time", &hodo_time_vec);
	tree->SetBranchAddress("hodo_bar_id", &hodo_bar_id_vec);
	tree->SetBranchAddress("hcal_clus_blk_id", &hcal_clus_blk_id);
	tree->SetBranchAddress("hcal_clus_blk_e", &hcal_clus_blk_e);
	tree->SetBranchAddress("hodo_etof", &hodo_etof_vec);
	tree->SetBranchAddress("nhodo_clus",&nhodo_clus);


	TH1D *h_dx = new TH1D("h_dx","dx",200,-6,3);
	TH1D *h_dy = new TH1D("h_dy","dy",200,-4,4);
	TH1D *h_W2 = new TH1D("h_W2","W2",100,0,4);
	TH1D *h_W2_cut_cointime = new TH1D("h_W2_cut_cointime","W2 with cointime cut",100,-4,8);
	TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",250,-2,2,250,-4,4);
	TH2D *h_dyW2 = new TH2D("h_dyW2","dy v W2",200,-1,3,200,-4,4);
	TH2D *h_dxW2 = new TH2D("h_dxW2","dx v W2",200,-1,3,200,-2,2);
	TH1D *h_coin_time = new TH1D("h_coin_time","coin_time",200,-50,50);
	TH1D *h_eHCAL_cut_QE = new TH1D("h_eHCAL_cut_QE","eHCAL with cointime, W2, dx and dy cuts;eHCAL(GeV)", 200, 0, 1.5);
	TH1D *h_eHCAL_cut_cointime_W2 = new TH1D("h_eHCAL_cut_cointime_W2","eHCAL with cointime and W2 cuts;eHCAL(GeV)",200,0,1.5);
	TH1D *h_eHCAL_cut_cointime = new TH1D("h_eHCAL_cut_cointime","eHCAL with cointime cut;eHCAL(GeV)",200,0,1.5);
	TH1D *h_eHCAL = new TH1D("h_eHCAL","eHCAL;eHCAL(GeV)",200,0,1.5);

	TH1D *h_hodo_meantime = new TH1D("h_hodo_meantime","h_hodo_meantime;bb.hodotdc.clus.bar.tdc.meantime[0]",200,-30,50);
	TH1D *h_hcal_tdctime = new TH1D("h_hcal_tdctime","h_hcal_tdctime;sbs.hcal.clus_blk.tdctime[0]",200,-30,50);
	TH1D *h_hcal_tdctime_tw = new TH1D("h_hcal_tdctime_tw","h_hcal_tdctime_tw;sbs.hcal.clus_blk.tdctime_tw[0]",200,-30,50);
	TH1D *h_cointime_new = new TH1D("h_cointime_new","h_cointime;hcal_tdc - hodo_tdc (ns)",200,-30,50);
	TH1D *h_cointime_new_tw = new TH1D("h_cointime_new_tw","h_cointime_tw;hcal_tdc_tw - hodo_tdc (ns)",200,-30,50);
	TH1D *h_cointime_new_tw_etof = new TH1D("h_cointime_new_tw_etof","h_cointime_tw_etof;hcal_tdc_tw_etof - hodo_tdc (ns)",200,-30,50);
	TH1D *h_cointime_new_RF = new TH1D("h_cointime_new_RF","h_cointime_RF;(hcal_tdc-sbs_rf) - (hodo_tdc-bb_rf) (ns)",200,-30,50);
	TH1D *h_cointime_new_RF_tw = new TH1D("h_cointime_new_RF_tw","h_cointime_tw_RF;(hcal_tdc_tw-sbs_rf) - (hodo_tdc-bb_rf) (ns)",200,-30,50);

	TH2D *h_hodo_meantime_vs_hodo_id = new TH2D("h_hodo_meantime_vs_hodo_id","h_hodo_meantime_vs_hodo_id",101,-0.5,100.5,200,-30,50);
	TH2D *h_cointime_new_vs_hodo_id = new TH2D("h_cointime_new_vs_hodo_id","h_cointime_new_vs_hodo_id",101,-0.5,100.5,200,-30,50);
	TH2D *h_cointime_new_tw_vs_hodo_id = new TH2D("h_cointime_new_tw_vs_hodo_id","h_cointime_new_tw_vs_hodo_id",101,-0.5,100.5,200,-30,50);

	TH2D *h_hcal_tdctime_vs_hcal_id = new TH2D("h_hcal_tdctime_vs_hcal_id","h_hcal_tdctime_vs_hcal_id",301,-0.5,300.5,200,-30,50);
	TH2D *h_hcal_tdctime_tw_vs_hcal_id = new TH2D("h_hcal_tdctime_tw_vs_hcal_id","h_hcal_tdctime_tw_vs_hcal_id",301,-0.5,300.5,200,-30,50);
	TH2D *h_cointime_new_vs_hcal_id = new TH2D("h_cointime_new_vs_hcal_id","h_cointime_new_vs_hcal_id",301,-0.5,300.5,200,-30,50);
	TH2D *h_cointime_new_tw_vs_hcal_id = new TH2D("h_cointime_new_tw_vs_hcal_id","h_cointime_new_tw_vs_hcal_id",301,-0.5,300.5,200,-30,50);

	TH2D *h_cointime_new_vs_dx = new TH2D("h_cointime_new_vs_dx","h_cointime_new_vs_dx",200,-4,3,200,-30,50);
	TH2D *h_cointime_new_tw_vs_dx = new TH2D("h_cointime_new_tw_vs_dx","h_cointime_new_tw_vs_dx",200,-4,3,200,-30,50);

	TH2D *h_cointime_new_vs_dy = new TH2D("h_cointime_new_vs_dy","h_cointime_new_vs_dy",200,-4,3,200,-30,50);
	TH2D *h_cointime_new_tw_vs_dy = new TH2D("h_cointime_new_tw_vs_dy","h_cointime_new_tw_vs_dy",200,-4,3,200,-30,50);

	TH1D *h_dx_W2_cut = new TH1D("h_dx_W2_cut","dx after W2, cointime and dy cut",200,-6,3);
	TH1D *h_dy_W2_cut = new TH1D("h_dy_W2_cut","dy after W2 and cointime cut",200,-4,4);
	//TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
	TH2D *h_dxdy_W2_cut = new TH2D("h_dxdy_W2_cut","dxdy after W2 and cointime cut",250,-2,2,250,-4,4);
	TH1D *h_coin_time_W2_cut = new TH1D("h_coin_time_W2_cut","cointime after W2 cut",200,-50,50);
	TH1D *h_W2_QE_cuts = new TH1D("h_W2_QE_cuts","W2",100,0,4);

	// Histogram array for time vs energy
	std::vector<TH2D*> hhcal_tdctime(289,nullptr);
	std::vector<TH2D*> hhcal_tdctime_tw(289,nullptr);

	// Histogram array for time
	std::vector<TH1D*> hhcal_tdctime_1D(289,nullptr);
	std::vector<TH1D*> hhcal_tdctime_tw_1D(289,nullptr);

	for(int i=0; i<289; ++i){
  		hhcal_tdctime[i] = new TH2D(Form("ht_blk%03d",i),Form("HCAL block %03d: E vs  hcalt (ns); E [GeV]; t [ns]", i),100, 0.0, 0.2, 200, -50.0, 50.0);
 	 	hhcal_tdctime_tw[i] = new TH2D(Form("htwCorr_blk%03d",i),Form("HCAL block %03d: E vs hcalt tw corr (ns); E [GeV]; t_{tw corr} [ns]", i),100, 0.0, 0.2,200, -50.0, 50.0);
		
		hhcal_tdctime_1D[i] = new TH1D(Form("ht_blk_1D%03d",i),Form("HCAL block %03d : t ; t [ns]  ", i),200,-50.0,50.0);
		hhcal_tdctime_tw_1D[i] = new TH1D(Form("htwCorr_blk_1D%03d",i),Form("HCAL block %03d : t tw corr ; t tw corr [ns]  ", i),200,-50.0,50.0);
	}


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

		if (veto and ntrack_sbs>0){
			continue;
		}

		if(nhodo_clus==0){
				continue;
		}

		if(cutsobject.run_num_L<runnum and runnum<cutsobject.run_num_H and ntrack>0 and ePS>0.2){
			//before adding a cut on W2
			//h_dx->Fill(dx);
			//h_dy->Fill(dy);
			h_dxdy->Fill(dy,dx);
			h_eHCAL->Fill(eHCAL);
			//h_dyW2->Fill(W2,dy);
			//h_dxW2->Fill(W2,dx);

			//double hodo_time = hodo_time_vec[0];
			//double hodo_bar_id = hodo_bar_id_vec[0];

			if (eHCAL>eHCAL_L){
				h_dx->Fill(dx);
				h_dy->Fill(dy);
				h_W2->Fill(W2);
				h_coin_time->Fill(coin_time);
	        }

	        if (eHCAL>eHCAL_L and cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H /*and dy<cutsobject.dy_H and cutsobject.dy_L<dy*/){
	        	h_W2_cut_cointime->Fill(W2);
	        }

	        if (trigbits==4 and eHCAL>eHCAL_L and cutsobject.W2_L <W2 and W2<cutsobject.W2_H and ePS>0.2 and (eSH+ePS)/trP>0.65){

	        	double hodo_time = hodo_time_vec[0];
				double hodo_bar_id = hodo_bar_id_vec[0];
				
				double hodo_etof = 0.0;

				if (ETOF == 1){
					hodo_etof = hodo_etof_vec[0];
				}
				else{
					hodo_etof = 0.0;
				}

	        	h_hodo_meantime->Fill(hodo_time);
	        	h_hcal_tdctime->Fill(hcal_tdctime);
	        	h_hcal_tdctime_tw->Fill(hcal_tdctime_tw);
	        	h_cointime_new->Fill(hcal_tdctime-(hodo_time-hodo_etof));
	        	h_cointime_new_tw->Fill(hcal_tdctime_tw-(hodo_time-hodo_etof));
	        	h_cointime_new_tw_etof->Fill(hcal_tdctime_tw-(hodo_time-hodo_etof));
	        	h_cointime_new_RF->Fill((hcal_tdctime - RFtime_sbs) - (hodo_time - RFtime - hodo_etof));
	        	h_cointime_new_RF_tw->Fill((hcal_tdctime_tw - RFtime_sbs) - (hodo_time - RFtime - hodo_etof));

	        	h_hodo_meantime_vs_hodo_id->Fill(hodo_bar_id,hodo_time);
	        	h_hcal_tdctime_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime);
	        	h_hcal_tdctime_tw_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime_tw);

	        	h_cointime_new_vs_hodo_id->Fill(hodo_bar_id,hcal_tdctime-(hodo_time-hodo_etof));
	        	h_cointime_new_tw_vs_hodo_id->Fill(hodo_bar_id,hcal_tdctime_tw-(hodo_time-hodo_etof));
	        	h_cointime_new_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime-(hodo_time-hodo_etof));
	        	h_cointime_new_tw_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime_tw-(hodo_time-hodo_etof));

	        	h_cointime_new_vs_dx->Fill(dx,hcal_tdctime-(hodo_time-hodo_etof));
	        	h_cointime_new_tw_vs_dx->Fill(dx,hcal_tdctime_tw-(hodo_time-hodo_etof));
	            h_cointime_new_vs_dy->Fill(dy,hcal_tdctime-(hodo_time-hodo_etof));
	        	h_cointime_new_tw_vs_dy->Fill(dy,hcal_tdctime_tw-(hodo_time-hodo_etof));

	       		if( hcal_clus_blk_id>= 0 &&  hcal_clus_blk_id< 289){
				    hhcal_tdctime[hcal_clus_blk_id]->Fill( hcal_clus_blk_e,  hcal_tdctime);
				    hhcal_tdctime_tw[hcal_clus_blk_id]->Fill( hcal_clus_blk_e, hcal_tdctime_tw);

				    hhcal_tdctime_1D[hcal_clus_blk_id]->Fill(hcal_tdctime);
				    hhcal_tdctime_tw_1D[hcal_clus_blk_id]->Fill(hcal_tdctime_tw);
				}

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
						//h_W2_QE_cuts->Fill(W2);
					}
				}
				//h_dx_W2_cut->Fill(dx);
				h_dy_W2_cut->Fill(dy);
				h_dxdy_W2_cut->Fill(dy,dx);
				h_eHCAL_cut_cointime_W2->Fill(eHCAL);
				//h_dyW2->Fill(W2,dy);
				
			}	

			if (eHCAL>eHCAL_L and (cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H) ){
				if ((dy<cutsobject.dy_H and cutsobject.dy_L<dy)){
					if ((dx<cutsobject.dx_H and cutsobject.dx_L<dx) /*or (dx<cutsobject.dx_p_H and cutsobject.dx_p_L<dx)*/   ){
						h_W2_QE_cuts->Fill(W2);
					}
				}
			}	


			if(eHCAL>eHCAL_L and (cutsobject.W2_L<W2 and W2<cutsobject.W2_H) /*and dy<cutsobject.dy_H and cutsobject.dy_L<dy*/){
				h_coin_time_W2_cut->Fill(coin_time);

			}	
		}

	    if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();


	}
	

	///////////////////trying to do a t0 correction to HCAL tdc time ////////////////////////////////
	int nBlk = 289;
	double t0[nBlk];
	double tw_t0[nBlk];

	std::fill_n(t0, nBlk, std::numeric_limits<double>::quiet_NaN());
	std::fill_n(tw_t0, nBlk, std::numeric_limits<double>::quiet_NaN());

	TH1D *h_t0 = new TH1D("h_t0",
	                      "HCAL per–block t_{0}; block ID ; t_{0}  [ns]",
	                      nBlk, -0.5, nBlk - 0.5);

	TH1D *h_tw_t0 = new TH1D("h_tw_t0",
	                      "HCAL per–block t_{0}; block ID ; t_{0}  [ns]",
	                      nBlk, -0.5, nBlk - 0.5);

	TF1 *g = new TF1("g", "gaus", -20., 20.);
	TF1 *g_tw = new TF1("g_tw", "gaus", -20., 20.);

	double nSigma = 0.5; 

	for (int i = 0; i<nBlk; ++i){
		TH1D *h = hhcal_tdctime_1D[i];

		if(!h || h->GetEntries()<50){
			continue;
		}

	    double mean  = h->GetMean();
	    double sigma = h->GetStdDev();     // same as RMS but unbiased

	    // Sanity guard: a pathological distribution (all bins in one place)
	    if (sigma <= 0) continue;

	    double lo = mean - nSigma * sigma;
	    double hi = mean + nSigma * sigma;

	    // Ensure the window lies inside the histogram range
	    lo = std::max(lo, h->GetXaxis()->GetXmin());
	    hi = std::min(hi, h->GetXaxis()->GetXmax());

	    TF1 *gfit = g;                    // re‑use the TF1 object
	    gfit->SetRange(lo, hi);

	    double amp   = h->GetBinContent(h->GetMaximumBin());
	    gfit->SetParameters(amp, mean, sigma);

	    int status = h->Fit(gfit, "RQ"); // Quiet, No‑draw, Range is from TF1
	    if (status == 0) {
	        t0[i] = gfit->GetParameter(1);           // fitted μ
	        h_t0->SetBinContent(i + 1, t0[i]);
	        h_t0->SetBinError  (i + 1, gfit->GetParError(1));
	    }

	}

	for (int i = 0; i<nBlk; ++i){
		TH1D *h_tw = hhcal_tdctime_tw_1D[i];

		if(!h_tw || h_tw->GetEntries()<50){
			continue;
		}

	    double mean_tw  = h_tw->GetMean();
	    double sigma_tw = h_tw->GetStdDev();     // same as RMS but unbiased

	    // Sanity guard: a pathological distribution (all bins in one place)
	    if (sigma_tw <= 0) continue;

	    double lo_tw = mean_tw - nSigma * sigma_tw;
	    double hi_tw = mean_tw + nSigma * sigma_tw;

	    // Ensure the window lies inside the histogram range
	    lo_tw = std::max(lo_tw, h_tw->GetXaxis()->GetXmin());
	    hi_tw = std::min(hi_tw, h_tw->GetXaxis()->GetXmax());

	    TF1 *gfit_tw = g_tw;                    // re‑use the TF1 object
	    gfit_tw->SetRange(lo_tw, hi_tw);

	    double amp_tw   = h_tw->GetBinContent(h_tw->GetMaximumBin());
	    gfit_tw->SetParameters(amp_tw, mean_tw, sigma_tw);

	    int status_tw = h_tw->Fit(gfit_tw, "RQ"); // Quiet, No‑draw, Range is from TF1
	    if (status_tw == 0) {
	        tw_t0[i] = gfit_tw->GetParameter(1);           // fitted μ
	        h_tw_t0->SetBinContent(i + 1, tw_t0[i]);
	        h_tw_t0->SetBinError  (i + 1, gfit_tw->GetParError(1));
	    }

	}

	/////////////////////////// t0 correction //////////////////////////////////////////

	TH2D *h_hcal_tdctime_t0_corr_vs_hcal_id = new TH2D("h_hcal_tdctime_t0_corr_vs_hcal_id","h_hcal_tdctime_t0_corr_vs_hcal_id",301,-0.5,300.5,200,-30,30);
	TH2D *h_hcal_tdctime_tw_t0_corr_vs_hcal_id = new TH2D("h_hcal_tdctime_tw_t0_corr_vs_hcal_id","h_hcal_tdctime_to_corr_tw_vs_hcal_id",301,-0.5,300.5,200,-30,30);
	TH2D *h_cointime_new_t0_corr_vs_hcal_id = new TH2D("h_cointime_new_t0_corr_vs_hcal_id","h_cointime_new_t0_corr_vs_hcal_id",301,-0.5,300.5,200,-30,30);
	TH2D *h_cointime_new_tw_t0_corr_vs_hcal_id = new TH2D("h_cointime_new_tw_t0_corr_vs_hcal_id","h_cointime_new_tw_t0_corr_vs_hcal_id",301,-0.5,300.5,200,-30,30);

	TH1D *h_hcal_tdctime_t0_corr = new TH1D("h_hcal_tdctime_t0_corr","h_hcal_tdctime_t0_corr;sbs.hcal.clus_blk.tdctime[0]",200,-30,30);
	TH1D *h_hcal_tdctime_tw_t0_corr = new TH1D("h_hcal_tdctime_tw_to_corr","h_hcal_tdctime_tw_t0_corr;sbs.hcal.clus_blk.tdctime_tw[0]",200,-30,30);
	TH1D *h_cointime_new_t0_corr = new TH1D("h_cointime_new_t0_corr","h_cointime_new_t0_corr;sbs.hcal.clus_blk.tdctime[0]-bb.hodotdc.clus.bar.tdc.meantime[0] (ns)",200,-30,30);
	TH1D *h_cointime_new_tw_t0_corr = new TH1D("h_cointime_new_tw_t0_corr","h_cointime_new_tw_t0_corr;sbs.hcal.clus_blk.tdctime_tw[0]-bb.hodotdc.clus.bar.tdc.meantime[0] (ns)",200,-30,30);
	TH1D *h_cointime_new_RF_t0_corr = new TH1D("h_cointime_new_RF_t0_corr","h_cointime_new_RF_t0_corr;sbs.hcal.clus_blk.tdctime[0]-bb.hodotdc.clus.bar.tdc.meantime[0] (ns)",200,-30,30);
	TH1D *h_cointime_new_RF_tw_t0_corr = new TH1D("h_cointime_new_RF_tw_t0_corr","h_cointime_new_RF_tw_t0_corr;sbs.hcal.clus_blk.tdctime_tw[0]-bb.hodotdc.clus.bar.tdc.meantime[0] (ns)",200,-30,30);

	// Histogram array for time
	std::vector<TH1D*> h_cointime_tdctime_t0_1D(289,nullptr);
	std::vector<TH1D*> h_cointime_tdctime_tw_t0_1D(289,nullptr);

	for(int i=0; i<289; ++i){
		
		h_cointime_tdctime_t0_1D[i] = new TH1D(Form("hcoin_t0_blk_1D%03d",i),Form("HCAL block %03d : cointime_t0 ; t [ns]  ", i),200,-30.0,30.0);
		h_cointime_tdctime_tw_t0_1D[i] = new TH1D(Form("hcointwCorr_t0_blk_1D%03d",i),Form("HCAL block %03d : cointime tw corr t0; t tw corr [ns]  ", i),200,-30.0,30.0);
	}


	for (int i = 0; i<nentries*1; i++){
		tree->GetEntry(i);

	    if (veto and ntrack_sbs>0){
			continue;
		}

		if(nhodo_clus==0){
				continue;
		}

	    if (trigbits==4 and eHCAL>eHCAL_L and cutsobject.W2_L <W2 and W2<cutsobject.W2_H and ePS>0.2 and (eSH+ePS)/trP>0.65){

	        double hodo_time_1 = hodo_time_vec[0];
			double hodo_bar_id_1 = hodo_bar_id_vec[0];
			double hodo_etof_1 = hodo_etof_vec[0]; 

			int hcal_blk = static_cast<int>(hcal_clus_blk_id);

			h_hcal_tdctime_t0_corr->Fill(hcal_tdctime-t0[hcal_blk]);
			h_hcal_tdctime_tw_t0_corr->Fill(hcal_tdctime_tw - tw_t0[hcal_blk]);
			h_cointime_new_t0_corr->Fill(hcal_tdctime-hodo_time_1-t0[hcal_blk]);
			h_cointime_new_tw_t0_corr->Fill(hcal_tdctime_tw - hodo_time_1 - tw_t0[hcal_blk] + hodo_etof_1);
			h_cointime_new_RF_t0_corr->Fill(hcal_tdctime-hodo_time_1-t0[hcal_blk]- RFtime_sbs + RFtime + hodo_etof_1);
			h_cointime_new_RF_tw_t0_corr->Fill(hcal_tdctime_tw - hodo_time_1 - tw_t0[hcal_blk]- RFtime_sbs + RFtime + hodo_etof_1);			


			h_hcal_tdctime_t0_corr_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime-t0[hcal_blk]);
			h_hcal_tdctime_tw_t0_corr_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime_tw - tw_t0[hcal_blk]);
			h_cointime_new_t0_corr_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime-hodo_time_1-t0[hcal_blk] + hodo_etof_1);
			h_cointime_new_tw_t0_corr_vs_hcal_id->Fill(hcal_clus_blk_id,hcal_tdctime_tw - hodo_time_1 - tw_t0[hcal_blk] + hodo_etof_1);

	       	if( hcal_blk>= 0 &&  hcal_blk< 289){
				h_cointime_tdctime_t0_1D[hcal_blk]->Fill(hcal_tdctime-hodo_time_1-t0[hcal_blk] + hodo_etof_1);
				h_cointime_tdctime_tw_t0_1D[hcal_blk]->Fill(hcal_tdctime_tw - hodo_time_1 - tw_t0[hcal_blk] + hodo_etof_1);
			}
		}

	    if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();
	}

	gStyle->SetOptStat(0);   // hide plain histogram stats
	gStyle->SetOptFit(111); 

	////////////////////////////////////////////////////////////////////////////////////
	TF1 *g_cointime_new =FitGaussNSigma(h_cointime_new,"h_cointime_new",0.4,1);
	TF1 *g_cointime_new_tw =FitGaussNSigma(h_cointime_new_tw,"h_cointime_new_tw",0.4,1);
	TF1 *g_cointime_new_tw_etof =FitGaussNSigma(h_cointime_new_tw_etof,"h_cointime_new_tw_etof",0.4,1);
	TF1 *g_cointime_new_RF =FitGaussNSigma(h_cointime_new_RF,"h_cointime_new_RF",0.4,1);
	TF1 *g_cointime_new_RF_tw =FitGaussNSigma(h_cointime_new_RF_tw,"h_cointime_new_RF_tw",0.4,1);


	TF1 *g_cointime_new_t0_corr =FitGaussNSigma(h_cointime_new_t0_corr,"h_cointime_new_t0_corr",0.4,1);
	TF1 *g_cointime_new_tw_t0_corr =FitGaussNSigma(h_cointime_new_tw_t0_corr,"h_cointime_new_tw_t0_corr",0.4,1);
	TF1 *g_cointime_new_RF_t0_corr =FitGaussNSigma(h_cointime_new_RF_t0_corr,"h_cointime_new_RF_t0_corr",0.4,1);
	TF1 *g_cointime_new_RF_tw_t0_corr =FitGaussNSigma(h_cointime_new_RF_tw_t0_corr,"h_cointime_new_RF_tw_t0_corr",0.4,1);

	////////////////////////////////////////////////////////////////////////////////////

	//TH1D *h_t0_corr = new TH1D("h_t0_corr",
	//                      "HCAL per–block t_{0}; block ID ; t_{0}  [ns]",
	//                      nBlk, -0.5, nBlk - 0.5);

	//TH1D *h_tw_t0_corr = new TH1D("h_tw_t0_corr",
	//                      "HCAL per–block t_{0}; block ID ; t_{0}  [ns]",
	//                      nBlk, -0.5, nBlk - 0.5);

	TF1 *g_t0_corr = new TF1("g_t0_corr", "gaus", -20., 20.);
	TF1 *g_tw_t0_corr = new TF1("g_tw_t0_corr", "gaus", -20., 20.);

	const double nSigma1 = 0.5; 

	for (int i = 0; i<nBlk; ++i){
		TH1D *h1 = h_cointime_tdctime_t0_1D[i];

		if(!h1 || h1->GetEntries()<50){
			continue;
		}

	    double mean1  = h1->GetMean();
	    double sigma1 = h1->GetStdDev();     // same as RMS but unbiased

	    // Sanity guard: a pathological distribution (all bins in one place)
	    if (sigma1 <= 0) continue;

	    double center1;
	    int    imax1   = h1->GetMaximumBin();
	    center1        = h1->GetXaxis()->GetBinCenter(imax1);


	    double lo1 = center1 - nSigma1 * sigma1;
	    double hi1 = center1 + nSigma1 * sigma1;

	    // Ensure the window lies inside the histogram range
	    lo1 = std::max(lo1, h1->GetXaxis()->GetXmin());
	    hi1 = std::min(hi1, h1->GetXaxis()->GetXmax());

	    TF1 *gfit_t0_corr = g_t0_corr;                    // re‑use the TF1 object
	    gfit_t0_corr->SetRange(lo1, hi1);

	    double amp1   = h1->GetBinContent(h1->GetMaximumBin());
	    gfit_t0_corr->SetParameters(amp1, mean1, sigma1);

	    int status1 = h1->Fit(gfit_t0_corr, "RQ"); // Quiet, No‑draw, Range is from TF1
	    //if (status1 == 0) {
	        //t0[i] = gfit->GetParameter(1);           // fitted μ
	        //h_t0->SetBinContent(i + 1, t0[i]);
	        //h_t0->SetBinError  (i + 1, gfit->GetParError(1));
	    //}

	}


	const double nSigma2 = 0.5; 

	for (int i = 0; i<nBlk; ++i){
		TH1D *h2 = h_cointime_tdctime_tw_t0_1D[i];

		if(!h2 || h2->GetEntries()<50){
			continue;
		}

	    double mean2  = h2->GetMean();
	    double sigma2 = h2->GetStdDev();     // same as RMS but unbiased

	    // Sanity guard: a pathological distribution (all bins in one place)
	    if (sigma2 <= 0) continue;

	    double center2;
	    int    imax2   = h2->GetMaximumBin();
	    center2        = h2->GetXaxis()->GetBinCenter(imax2);


	    double lo2 = center2 - nSigma2 * sigma2;
	    double hi2 = center2 + nSigma2 * sigma2;

	    // Ensure the window lies inside the histogram range
	    lo2 = std::max(lo2, h2->GetXaxis()->GetXmin());
	    hi2 = std::min(hi2, h2->GetXaxis()->GetXmax());

	    TF1 *gfit_tw_t0_corr = g_tw_t0_corr;                    // re‑use the TF1 object
	    gfit_tw_t0_corr->SetRange(lo2, hi2);

	    double amp2   = h2->GetBinContent(h2->GetMaximumBin());
	    gfit_tw_t0_corr->SetParameters(amp2, mean2, sigma2);

	    int status2 = h2->Fit(gfit_tw_t0_corr, "RQ"); // Quiet, No‑draw, Range is from TF1
	    //if (status1 == 0) {
	        //t0[i] = gfit->GetParameter(1);           // fitted μ
	        //h_t0->SetBinContent(i + 1, t0[i]);
	        //h_t0->SetBinError  (i + 1, gfit->GetParError(1));
	    //}

	}

    ////////////////////////////////////////////////////////////////////////////////////
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
	TCanvas* c5 = new TCanvas("c5","c5",3600,3000);
	TCanvas* c6 = new TCanvas("c6","c6",3600,3000);
	TCanvas* c7 = new TCanvas("c7","c7",3600,3000);
	TCanvas* c8 = new TCanvas("c8","c8",3600,3000);
	TCanvas* c9 = new TCanvas("c9","c9",3600,3000);
	TCanvas* c10 = new TCanvas("c10","c10",3600,3000);
	TCanvas* c11 = new TCanvas("c11","c11",3600,3000);
	TCanvas* c12 = new TCanvas("c12","c12",3600,3000);
	TCanvas* c13 = new TCanvas("c13","c13",3600,3000);
	TCanvas* c14 = new TCanvas("c14","c14",3600,3000);

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

	c5->Divide(2,2);
	c5->cd(1);
	h_dx->SetLineWidth(2);
	h_dx->Draw();
	h_dx_W2_cut->SetLineWidth(2);
	h_dx_W2_cut->SetLineColor(kGreen);
	h_dx_W2_cut->Draw("SAME");
	c5->cd(2);
	h_dy->SetLineWidth(2);
	h_dy->Draw();
	h_dy_W2_cut->SetLineWidth(2);
	h_dy_W2_cut->SetLineColor(kGreen);
	h_dy_W2_cut->Draw("SAME");
	c5->cd(3);
	h_W2->SetLineWidth(2);
	h_W2->Draw();
	h_W2_QE_cuts->SetLineWidth(2);
	h_W2_QE_cuts->SetLineColor(kGreen);
	h_W2_QE_cuts->Draw("SAME");

	c5->Print(Form("plots/%s_cuts_vs_nocuts_eHCAL_cut_%s.png",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	c6->Divide(2,2);
	c6->cd(1);
	h_hodo_meantime->Draw();
	c6->cd(2);
	h_hcal_tdctime->Draw();
	c6->cd(3);
	h_hcal_tdctime_tw->Draw();

	c7->Divide(2,2);
	c7->cd(1);
	h_cointime_new->Draw();
	c7->cd(2);
	h_cointime_new_tw->Draw();
	c7->cd(3);
	h_cointime_new_RF->Draw();
	c7->cd(4);
	h_cointime_new_RF_tw->Draw();

	c8->Divide(2,2);
	c8->cd(1);
	h_hodo_meantime_vs_hodo_id->Draw("COLZ");
	c8->cd(2);
	h_hcal_tdctime_vs_hcal_id->Draw("COLZ");
	c8->cd(3);
	h_hcal_tdctime_tw_vs_hcal_id->Draw("COLZ");

	c9->Divide(2,2);
	c9->cd(1);
	h_cointime_new_vs_hodo_id->Draw("COLZ");
	c9->cd(2);
	h_cointime_new_tw_vs_hodo_id->Draw("COLZ");
	c9->cd(3);
	h_cointime_new_vs_hcal_id->Draw("COLZ");
	c9->cd(4);
	h_cointime_new_tw_vs_hcal_id->Draw("COLZ");


	c10->Divide(2,2);
	c10->cd(1);
	h_cointime_new_vs_dx->Draw("COLZ");
	c10->cd(2);
	h_cointime_new_tw_vs_dx->Draw("COLZ");
	c10->cd(3);
	h_cointime_new_vs_dy->Draw("COLZ");
	c10->cd(4);
	h_cointime_new_tw_vs_dy->Draw("COLZ");

	c11->Divide(2,2);
	c11->cd(1);
	h_hcal_tdctime_t0_corr->Draw();
	c11->cd(2);
	h_cointime_new_t0_corr->Draw();
	c11->cd(3);
	h_hcal_tdctime_t0_corr_vs_hcal_id->Draw("COLZ");
	c11->cd(4);
	h_cointime_new_t0_corr_vs_hcal_id->Draw("COLZ");

	c12->Divide(2,2);
	c12->cd(1);
	h_hcal_tdctime_tw_t0_corr->Draw();
	c12->cd(2);
	h_cointime_new_tw_t0_corr->Draw();
	c12->cd(3);
	h_hcal_tdctime_tw_t0_corr_vs_hcal_id->Draw("COLZ");
	c12->cd(4);
	h_cointime_new_tw_t0_corr_vs_hcal_id->Draw("COLZ");

	c13->Divide(2,2);
	c13->cd(1);
	h_cointime_new_RF_t0_corr->Draw();
	c13->cd(2);
	h_cointime_new_RF_tw_t0_corr->Draw();

	c14->Divide(2,2);
	c14->cd(1);
	h_cointime_new_tw_etof->Draw();

	c6->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf(",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c7->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c8->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c9->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c10->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c11->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));	
	c12->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c13->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	c14->Print(Form("plots/%s_timing_eHCAL_cut_%s.pdf)",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
    // Create a canvas
	TCanvas *chcal = new TCanvas("chcal","HCAL blocks",1200,900);

	// Open the multi‐page PDF
	chcal->Print(Form("plots/HCAL_blk_tw_%s_eHCAL_cut_%s.pdf(",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	// Loop over all blocks: draw & print
	for(int i=0; i<289; ++i){
	   chcal->Clear();
	   chcal->Divide(2,1);

	   chcal->cd(1);
	   hhcal_tdctime[i]->Draw("COLZ");         // corrected-time

	   chcal->cd(2);
	   hhcal_tdctime_tw[i]->Draw("COLZ");          // raw-time

	   chcal->Print(Form("plots/HCAL_blk_tw_%s_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	}

	  // Close the multi‐page PDF
	chcal->Print(Form("plots/HCAL_blk_tw_%s_eHCAL_cut_%s.pdf)",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	/////////////////////////////////////////////////////////////////////////////////
    // Create a canvas
	TCanvas *chcal_1D = new TCanvas("chcal_1D","HCAL blocks 1D",1200,900);

	// Open the multi‐page PDF
	chcal_1D->Print(Form("plots/HCAL_blk_tw_1D_%s_eHCAL_cut_%s.pdf(",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	// Loop over all blocks: draw & print
	for(int i=0; i<289; ++i){
	   chcal_1D->Clear();
	   chcal_1D->Divide(2,1);

	   chcal_1D->cd(1);
	   hhcal_tdctime_1D[i]->Draw();         // corrected-time

	   chcal_1D->cd(2);
	   hhcal_tdctime_tw_1D[i]->Draw();          // raw-time

	   chcal_1D->Print(Form("plots/HCAL_blk_tw_1D_%s_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	}

	  // Close the multi‐page PDF
	chcal_1D->Print(Form("plots/HCAL_blk_tw_1D_%s_eHCAL_cut_%s.pdf)",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

    /////////////////////////////////////////////////////////////////////////////////
	TCanvas *ct0 = new TCanvas("ct0","ct0",1200,900);

	ct0->Divide(2,1);
	ct0->cd(1);
	h_t0->Draw();
	ct0->cd(2);
	h_tw_t0->Draw();
	ct0->Print(Form("plots/HCAL_t0_%s_eHCAL_cut_%s.pdf)",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	/////////////////////////////////////////////////////////////////////////////////

    // Create a canvas
	TCanvas *ccoin_1D = new TCanvas("ccoin_1D","HCAL blocks 1D",1200,900);

	// Open the multi‐page PDF
	ccoin_1D->Print(Form("plots/HCAL_blk_coin_tw_1D_%s_eHCAL_cut_%s.pdf(",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

	// Loop over all blocks: draw & print
	for(int i=0; i<289; ++i){
	   ccoin_1D->Clear();
	   ccoin_1D->Divide(2,1);

	   ccoin_1D->cd(1);
	   h_cointime_tdctime_t0_1D[i]->Draw();         // corrected-time

	   ccoin_1D->cd(2);
	   h_cointime_tdctime_tw_t0_1D[i]->Draw();          // raw-time

	   ccoin_1D->Print(Form("plots/HCAL_blk_coin_tw_1D_%s_eHCAL_cut_%s.pdf",printfilename,std::to_string(flag_eHCAL_cut).c_str()));
	}

	  // Close the multi‐page PDF
	ccoin_1D->Print(Form("plots/HCAL_blk_coin_tw_1D_%s_eHCAL_cut_%s.pdf)",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

}
