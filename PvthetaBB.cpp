#include <iostream>
#include <cmath>
#include "cuts.h"
#include "plotdxdy.h"

void PvthetaBB(const char* filename, const char* printfilename, const char *kin){

	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	cuts cutsobject;	
	cutsobject.parsecuts(config);

	//std::cout<<"coin_time : "<<cutsobject.coin_time_L<<endl;
	
	TFile* file = TFile::Open(filename);
	TTree* tree = (TTree*)file->Get("Tout");

	double dx = 0.0; 
	double dy = 0.0; 
	double W2 = 0.0;
    	double Q2 = 0.0;
	double coin_time = 0.0;
	double trP = 0.0;
	double etheta = 0.0;
	double eSH = 0.0;
	double ePS = 0.0;
	//double ptheta = 0.0;
	//double pphi = 0.0; 
	//double ptheta_cal = 0.0;
	//double pphi_cal = 0.0;
	//double pN_expect = 0.0;
	//double trP_sbs = 0.0;
	//double trx_sbs = 0.0;
	//double try_sbs = 0.0;
	//double trth_sbs = 0.0;
	//double trph_sbs = 0.0;

	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("trP",&trP);
	tree->SetBranchAddress("etheta",&etheta);
	tree->SetBranchAddress("eSH",&eSH);
	tree->SetBranchAddress("ePS",&ePS);

	//tree->SetBranchAddress("ptheta",&ptheta);
	//tree->SetBranchAddress("pphi",&pphi);
	//tree->SetBranchAddress("ptheta_cal",&ptheta_cal);
	//tree->SetBranchAddress("pphi_cal",&pphi_cal);
	//tree->SetBranchAddress("pN_expect",&pN_expect);
	//tree->SetBranchAddress("trP_sbs",&trP_sbs);
	//tree->SetBranchAddress("trx_sbs",&trx_sbs);
	//tree->SetBranchAddress("try_sbs",&try_sbs);
	//tree->SetBranchAddress("trth_sbs",&trth_sbs);
	//tree->SetBranchAddress("trph_sbs",&trph_sbs);

	TH1D *h_coin_time = new TH1D("h_coin_time","h_coin_time",500,40,140);
	TH1D *h_epratio = new TH1D("h_epratio","h_epratio",500,-1,4);

	TH2D *h_trPvetheta = new TH2D("h_trPvetheta","h_trPvetheta",200,0.3,0.75,200,1,5);
	TH2D *h_trPvetheta_cor = new TH2D("h_trPvetheta_cor","h_trPvetheta_cor",200,0.3,0.75,200,1,15);
	TH1D *h_trP = new TH1D("h_trP","h_trP",200,0,5);
	TH1D *h_trP_cor = new TH1D("h_trP_cor","h_trP_cor",200,0,10);	
	TH1D *h_etheta = new TH1D("h_etheta","h_etheta",200,0.4,0.75);

	bool coin_time_cut = false;

	int nentries = tree->GetEntries();
	
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		h_coin_time->Fill(coin_time);
		h_epratio->Fill((ePS+eSH)/trP);

		if (coin_time){
			h_trPvetheta->Fill(etheta,trP);
			h_etheta->Fill(etheta);
			h_trP->Fill(trP);

		}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();
	}

	double autotrPmin = h_trP->GetMean() - 6 * h_trP->GetRMS(); 
    double autotrPmax = h_trP->GetMean() + 3 * h_trP->GetRMS(); ;

	double autoethetamin = h_etheta->GetMean() - 3 * h_etheta->GetRMS();  
    double autoethetamax = h_etheta->GetMean() + 5 * h_etheta->GetRMS();

    std::cout<<"trPmin : "<<autotrPmin<<"trPmax : "<<autotrPmax<<endl;
    std::cout<<"ethetamin : "<<autoethetamin<<"autoethetamax : "<<autoethetamax<<endl;

	h_trPvetheta->GetXaxis()->SetRangeUser(autoethetamin, autoethetamax);
	h_trPvetheta->GetYaxis()->SetRangeUser(autotrPmin, autotrPmax);

	double correlation = h_trPvetheta->GetCorrelationFactor();

	std::cout<<"correlation factor : "<<correlation<<endl;

	TF1* line = new TF1("line", "pol1", h_trPvetheta->GetXaxis()->GetBinCenter(1), h_trPvetheta->GetXaxis()->GetBinCenter(h_trPvetheta->GetNbinsX()));
	h_trPvetheta->Fit(line);

	std::cout<<"slope : "<<line->GetParameter(0)<<"intercept : "<<line->GetParameter(1)<<endl;

	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);

		if (coin_time){
			h_trP_cor->Fill(trP+line->GetParameter(0)*etheta);
			h_trPvetheta_cor->Fill(etheta,trP+line->GetParameter(0)*etheta);
		}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();
	}


	h_trPvetheta_cor->GetXaxis()->SetRangeUser(autoethetamin, autoethetamax);
	h_trPvetheta_cor->GetYaxis()->SetRangeUser(autotrPmin+2, autotrPmax+2);

	h_trP_cor->GetXaxis()->SetRangeUser(autotrPmin+2, autotrPmax+3);

	TCanvas* c = new TCanvas("c","c",2400,1800);
	TCanvas* c1 = new TCanvas("c1","c1",2400,1800);
	TCanvas* c2 = new TCanvas("c2","c2",2400,1800);

	c->Divide(2,1);
	c1->Divide(2,2);
	c2->Divide(2,2);

	c->cd(1);
	TLine *line01 = new TLine(cutsobject.coin_time_mean-cutsobject.coin_time_width, 0, cutsobject.coin_time_mean-cutsobject.coin_time_width, 50000);
    TLine *line02 = new TLine(cutsobject.coin_time_mean+cutsobject.coin_time_width, 0, cutsobject.coin_time_mean+cutsobject.coin_time_width, 50000);
	h_coin_time->SetXTitle("cointime(ns)");
	h_coin_time->Draw();
    line01->SetLineColor(kRed);
    line02->SetLineColor(kRed);
    line01->SetLineWidth(2);
    line02->SetLineWidth(2);
	line01->Draw("same");
    line02->Draw("same");

	c->cd(2);
	h_epratio->SetXTitle("e/p");
	h_epratio->Draw();

	c1->cd(1);
	h_etheta->SetXTitle("etheta(rad)");
	h_etheta->Draw();

	c1->cd(2);
	h_trP->SetXTitle("trP(GeV)");
	h_trP->Draw();

	c1->cd(3);
	h_trPvetheta->SetXTitle("etheta(rad)");
	h_trPvetheta->SetYTitle("trP(GeV)");
	h_trPvetheta->Draw("COLZ");
	line->Draw("same");

	c2->cd(1);
	h_trP_cor->SetXTitle("trP corrected(GeV)");
	h_trP_cor->Draw();

	c2->cd(2);
	h_trPvetheta_cor->SetXTitle("etheta(rad)");
	h_trPvetheta_cor->SetYTitle("trP corrected(GeV)");
	h_trPvetheta_cor->Draw("COLZ");
	//line->Draw("same");	

	c->Print(Form("../plots/%s_bbtracking.pdf(",printfilename));
	c1->Print(Form("../plots/%s_bbtracking.pdf",printfilename));
	c2->Print(Form("../plots/%s_bbtracking.pdf)",printfilename));	

	c1->SaveAs(Form("../plots/%s_bbtracking.png",printfilename));
	//c1->SaveAs(Form("../plots/%s_pphi_sbstracking.png",printfilename));
}
