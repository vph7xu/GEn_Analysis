#include <iostream>
#include <locale>
#include <TStyle.h>

#include "cuts.h"
#include "plotdxdy.h"

void PS_HCAL_e_plot(const char* filename, const char* printfilename, const char* kin){

	TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

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

	double dx_p_L = getDoubleValue(config,"dx_p_L");
	double dx_p_H = getDoubleValue(config,"dx_p_H");

	double eHCAL_L = getDoubleValue(config,"eHCAL_L");

	double dy_ac_L = getDoubleValue(config,"dy_ac_L");
	double dy_ac_H = getDoubleValue(config,"dy_ac_H");

	double run_num_L = getDoubleValue(config,"run_num_L");
	double run_num_H = getDoubleValue(config,"run_num_H");

	double IHWP_flip = getDoubleValue(config,"IHWP_flip");

	double ePS = 0.0;
	double eSH = 0.0;
	double trP = 0.0;
	double eHCAL = 0.0;
	double vz = 0.0;
	double ntrack = 0.0;
	double nhits = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;
	int runnum = 0.0;
	int helicity = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	//double eHCAL = 0.0;

	tree->SetBranchAddress("ePS",&ePS);
	tree->SetBranchAddress("eSH",&eSH);
	tree->SetBranchAddress("trP",&trP);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("ntrack",&ntrack);
	tree->SetBranchAddress("nhits",&nhits);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);

	//tree->SetBranchAddress("eHCAL",&eHCAL);

	TH1D *h_ePS = new TH1D("h_ePS","Preshower Energy",200,0,4);
	TH1D *h_eHCAL = new TH1D("h_eHCAL","HCAL Energy",100,0,1);
	TH1D *h_trvz = new TH1D("h_trvz","vertex z",200,-1,1);
	TH1D *h_ntrack = new TH1D("h_ntrack","ntrack",60,0,6);
	TH1D *h_nhits = new TH1D("h_nhits","nhits",60,0,6);
	TH1D *h_W2 = new TH1D("h_W2","W^{2}",200,-2,4);
	TH1D *h_cointime = new TH1D("h_cointime","Coincidence time",200,coin_time_L-50,coin_time_H+50);
	TH1D *h_eoverp = new TH1D("h_eoverp","E/p BBCAL ; E/p", 100,0.3,1.7);
	//TH1D *h_eHCAL = new TH1D("h_eHCAL","HCAL energy ; HCAL energy (GeV)",200,0,2);

	TH1D *h_trvz_c = new TH1D("h_trvz_c","vertex z",1000,-1,1);

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);

		double eoverp = (eSH+ePS)/trP;

	    bool goodHelicity = (lookupValue(HelicityCheck, runnum) == 1);
	    bool goodMoller   = (lookupValue(MollerQuality, runnum) == 1);
	    bool goodVz       = abs(vz) < 0.27;
	    bool goodPS       = (ePS > 0.2);
	    bool goodRunRange = (run_num_L < runnum && runnum < run_num_H);
	    bool goodEHCAL    = (eHCAL > eHCAL_L); 
	    bool validHel     = (helicity == -1 || helicity == 1);
	    //bool goodGrinch = (grinch_track == 0) && (grinch_clus_size>2);
	    //bool goodSbs_track = ntrack_sbs>0 && abs(vz_sbs)<0.27;
	    bool goodEoverp = abs(eoverp-1)<0.8;
	    bool goodW2 = W2_L<W2 && W2<W2_H;
	    bool goodCoin = coin_time_L<coin_time && coin_time<coin_time_H;
	    bool goodSpot = dy_L<dy && dy<dy_H && ((dx_L<dx && dx<dx_H) || (dx_p_L<dx && dx<dx_p_H));


		if(ePS>0 && eHCAL>0){
			if (goodEHCAL && goodPS){
				h_trvz->Fill(vz);
			}

			if(goodVz && goodEHCAL){
				h_ePS->Fill(ePS);
			}
			
			if(goodVz && goodEoverp && goodCoin){
				h_eHCAL->Fill(eHCAL);
			}
			
			h_ntrack->Fill(ntrack);
			h_nhits->Fill(nhits);

			if(goodVz && goodPS && goodEHCAL && goodEoverp && goodCoin){
				h_W2->Fill(W2);
			}
			if(goodVz && goodPS && goodEHCAL && goodEoverp && goodW2){
				h_cointime->Fill(coin_time);
			}

			if(goodVz && goodPS && goodEHCAL && goodEoverp && goodW2 && goodCoin && goodSpot){
				h_eoverp->Fill(eoverp);
			}
		
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
	TLine *line01 = new TLine(coin_time_L, h_cointime->GetMinimum(), coin_time_L, h_cointime->GetMaximum());
    TLine *line02 = new TLine(coin_time_H, h_cointime->GetMinimum(), coin_time_H, h_cointime->GetMaximum());
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
	TLine *line05 = new TLine(W2_L, 0, W2_L, h_W2->GetMaximum());
    TLine *line06 = new TLine(W2_H, 0, W2_H, h_W2->GetMaximum());
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
	c1->Divide(2,2);
	c1->cd(1);
	//c1->SetLogy();
	TLine *line09 = new TLine(eHCAL_L, 0, eHCAL_L, h_eHCAL->GetMaximum());
	line09->SetLineColor(kRed);
	h_eHCAL->SetLineWidth(2);
	h_eHCAL->Draw();
	line09->Draw("same");
	//c1->Update();
	c1->cd(2);
	h_eoverp->SetLineWidth(2);
	h_eoverp->Draw();
	TLine *line07 = new TLine(0.8, 0, 0.8, h_eoverp->GetMaximum());
    TLine *line08 = new TLine(1.2, 0, 1.2, h_eoverp->GetMaximum());
    line07->SetLineColor(kRed);
    line08->SetLineColor(kRed);
    line07->SetLineWidth(2);
    line08->SetLineWidth(2);
    line07->Draw("same");
    line08->Draw("same");


	//TCanvas* c2 = new TCanvas("c2","c2",1600,1200);
	//c2->cd(1);


	c->Print(Form("../plots/%s_ePS_eHCAL.pdf(",printfilename));	
	c->SaveAs(Form("../plots/%s_ePS_eHCAL.png",printfilename));
	c1->SaveAs(Form("../plots/%s_ePS_eHCAL_1.png",printfilename));
	c1->Print(Form("../plots/%s_ePS_eHCAL.pdf)",printfilename));


}
