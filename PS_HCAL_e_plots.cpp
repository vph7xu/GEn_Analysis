#include <iostream>
#include <locale>
#include <TStyle.h>
void PS_HCAL_e_plot(const char* filename, const char* printfilename){

	TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

	double ePS = 0.0;
	double eHCAL = 0.0;
	double vz = 0.0;
	double ntrack = 0.0;
	double nhits = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;

	tree->SetBranchAddress("ePS",&ePS);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("ntrack",&ntrack);
	tree->SetBranchAddress("nhits",&nhits);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("coin_time",&coin_time);

	TH1D *h_ePS = new TH1D("h_ePS","Preshower Energy",1000,0,4);
	TH1D *h_eHCAL = new TH1D("h_eHCAL","HCAL Energy",1000,0,2);
	TH1D *h_trvz = new TH1D("h_trvz","vertex z",1000,-1,1);
	TH1D *h_ntrack = new TH1D("h_ntrack","ntrack",60,0,6);
	TH1D *h_nhits = new TH1D("h_nhits","nhits",60,0,6);
	TH1D *h_W2 = new TH1D("h_W2","W^{2}",1000,-2,4);
	TH1D *h_cointime = new TH1D("h_cointime","coincidence time",1000,40,180);

	TH1D *h_trvz_c = new TH1D("h_trvz_c","vertex z",1000,-1,1);

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		if(ePS>0 and eHCAL>0){
			h_ePS->Fill(ePS);
			h_eHCAL->Fill(eHCAL);
			h_trvz->Fill(vz);
			h_ntrack->Fill(ntrack);
			h_nhits->Fill(nhits);
			h_W2->Fill(W2);
			h_cointime->Fill(coin_time);
		}
		//if (ePS>0.2 and eHCAL>0){
		//	h_trvz_c->Fill(vz);
		//	h_eHCAL->Fill(eHCAL);
		//}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();

	}

	TCanvas* c = new TCanvas("c","c",1600,1200);
	c->Divide(2,2);
	c->cd(1);
	h_ePS->SetXTitle("Preshower Energy(GeV)");
	h_ePS->Draw();
	TLine *line00 = new TLine(0.2, 0, 0.2, h_ePS->GetMaximum());
    line00->SetLineColor(kRed);
    line00->SetLineWidth(2);
    line00->Draw("same");

	c->cd(2);
	h_cointime->SetXTitle("coincidence time(ns)");
	h_cointime->Draw();
	TLine *line01 = new TLine(90, h_cointime->GetMinimum(), 90, h_cointime->GetMaximum());
    TLine *line02 = new TLine(110, h_cointime->GetMinimum(), 110, h_cointime->GetMaximum());
    line01->SetLineColor(kRed);
    line02->SetLineColor(kRed);
    line01->SetLineWidth(2);
    line02->SetLineWidth(2);
    line01->Draw("same");
    line02->Draw("same");
	//h_eHCAL->SetXTitle("HCal Energy(GeV)");
	//h_eHCAL->Draw();

	c->cd(3);
	h_trvz->SetXTitle("Target Vertex(m)");
	h_trvz->Draw();
	TLine *line03 = new TLine(-0.27, 0, -0.27, h_trvz->GetMaximum());
    TLine *line04 = new TLine(0.27, 0, 0.27, h_trvz->GetMaximum());
    line03->SetLineColor(kRed);
    line04->SetLineColor(kRed);
    line03->SetLineWidth(2);
    line04->SetLineWidth(2);
    line03->Draw("same");
    line04->Draw("same");
	
	c->cd(4);
	h_W2->SetXTitle("W^{2}(GeV^{2})");
	h_W2->Draw();
	TLine *line05 = new TLine(0, 0, 0, h_W2->GetMaximum());
    TLine *line06 = new TLine(1.6, 0, 1.6, h_W2->GetMaximum());
    line05->SetLineColor(kRed);
    line06->SetLineColor(kRed);
    line05->SetLineWidth(2);
    line06->SetLineWidth(2);
    line05->Draw("same");
    line06->Draw("same");
	//h_ntrack->Draw();
	//c->cd(5);
	//h_nhits->Draw();
	//h_trvz_c->Draw();

	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->cd(1);
	c1->SetLogy();
	h_eHCAL->Draw();
	c1->Update();

	c->Print(Form("../plots/%s_ePS_eHCAL.pdf(",printfilename));	
	c->SaveAs(Form("../plots/%s_ePS_eHCAL.png",printfilename));

	c1->Print(Form("../plots/%s_ePS_eHCAL.pdf)",printfilename));

}
