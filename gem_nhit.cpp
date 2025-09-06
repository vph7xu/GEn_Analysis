#include <iostream>
#include <locale>
#include <TStyle.h>

#include "cuts.h"
#include "plotdxdy.h"

void gem_nhit(const char* filename, const char* printfilename, const char* kin){

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
	double track_nhits = 0.0;
	double nhits = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;
	int runnum = 0.0;
	int helicity = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	double track_chi2 = 0.0;

	//double eHCAL = 0.0;

	tree->SetBranchAddress("ePS",&ePS);
	tree->SetBranchAddress("eSH",&eSH);
	tree->SetBranchAddress("trP",&trP);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("ntrack_chi2",&track_chi2);
	tree->SetBranchAddress("ntrack_hits",&track_nhits);
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

	TH1D *h_dx_3hit = new TH1D("h_dx_3hit","delta-x distribution (3 hit tracks) (m); dx (m)", 100, -4, 3);
	TH1D *h_dx_4hit = new TH1D("h_dx_4hit","delta-x distribution (4 hit tracks) (m); dx (m)", 100, -4, 3);
	TH1D *h_dx_5hit = new TH1D("h_dx_5hit","delta-x distribution (5 hit tracks) (m); dx (m)", 100, -4, 3);
	TH1D *h_dx_all_hit = new TH1D("h_dx_all_hit","delta-x distribution (all hit tracks) (m); dx (m)", 100, -4, 3);
	TH1D *h_dx_4_5hit = new TH1D("h_dx_4_5hit","delta-x distribution (4 & 5 hit tracks) (m); dx (m)", 100, -4, 3);
	//TH1D *h_eHCAL = new TH1D("h_eHCAL","HCAL energy ; HCAL energy (GeV)",200,0,2);

	TH1D *h_W2_3hit = new TH1D("h_W2_3hit","W^{2} distribution (3 hit tracks) (GeV^{2}); W^{2} (GeV^{2})", 100, -3, 5);
	TH1D *h_W2_4hit = new TH1D("h_W2_4hit","W^{2} distribution (4 hit tracks) (GeV^{2}); W^{2} (GeV^{2})", 100, -3, 5);
	TH1D *h_W2_5hit = new TH1D("h_W2_5hit","W^{2} distribution (5 hit tracks) (GeV^{2}); W^{2} (GeV^{2})", 100, -3, 5);
	TH1D *h_W2_all_hit = new TH1D("h_W2_all_hit","W^{2} distribution (all hit tracks) (GeV^{2}); W^{2} (GeV^{2})", 100, -3, 5);
	TH1D *h_W2_4_5hit = new TH1D("h_W2_4_5hit","W^{2} distribution (4 & 5 hit tracks) (GeV^{2}); W^{2} (GeV^{2})", 100, -3, 5);

	TH1D *h_cointime_3hit = new TH1D("h_cointime_3hit","cointime distribution (3 hit tracks) (ns); cointime (ns)", 100, coin_time_L-15, coin_time_H+15);
	TH1D *h_cointime_4hit = new TH1D("h_cointime_4hit","cointime distribution (4 hit tracks) (ns); cointime (ns)", 100, coin_time_L-15, coin_time_H+15);
	TH1D *h_cointime_5hit = new TH1D("h_cointime_5hit","cointime distribution (5 hit tracks) (ns); cointime (ns)", 100, coin_time_L-15, coin_time_H+15);
	TH1D *h_cointime_all_hit = new TH1D("h_cointime_all_hit","cointime distribution (all hit tracks) (ns); cointime (ns)", 100, coin_time_L-15, coin_time_H+15);
	TH1D *h_cointime_4_5hit = new TH1D("h_cointime_4_5hit","cointime distribution (4 & 5 hit tracks) (ns); cointime (ns)", 100, coin_time_L-15, coin_time_H+15);

	TH1D *h_vz_3hit = new TH1D("h_vz_3hit","vz distribution (3 hit tracks) (m); vz (m)", 100, -0.3, 0.3);
	TH1D *h_vz_4hit = new TH1D("h_vz_4hit","vz distribution (4 hit tracks) (m); vz (m)", 100, -0.3, 0.3);
	TH1D *h_vz_5hit = new TH1D("h_vz_5hit","vz distribution (5 hit tracks) (m); vz (m)", 100, -0.3, 0.3);
	TH1D *h_vz_all_hit = new TH1D("h_vz_all_hit","vz distribution (all hit tracks) (m); vz (m)", 100, -0.3, 0.3);
	TH1D *h_vz_4_5hit = new TH1D("h_vz_4_5hit","vz distribution (4 & 5 hit tracks) (m); vz (m)", 100, -0.3, 0.3);

	TH1D *h_chi2_3hit = new TH1D("h_chi2_3hit","chi2 distribution (3 hit tracks) ; chi2 ", 100, 0, 60);
	TH1D *h_chi2_4hit = new TH1D("h_chi2_4hit","chi2 distribution (4 hit tracks) ; chi2 ", 100, 0, 60);
	TH1D *h_chi2_5hit = new TH1D("h_chi2_5hit","chi2 distribution (5 hit tracks) ; chi2 ", 100, 0, 60);
	TH1D *h_chi2_all_hit = new TH1D("h_chi2_all_hit","chi2 distribution (all hit tracks) ; chi2 ", 100, 0, 60);
	TH1D *h_chi2_4_5hit = new TH1D("h_chi2_4_5hit","chi2 distribution (4 & 5 hit tracks) ; chi2 ", 100, 0, 60);

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
	    bool goodEHCAL    = eHCAL > eHCAL_L; 
	    bool validHel     = (helicity == -1 || helicity == 1);
	    //bool goodGrinch = (grinch_track == 0) && (grinch_clus_size>2);
	    //bool goodSbs_track = ntrack_sbs>0 && abs(vz_sbs)<0.27;
	    bool goodEoverp = abs(eoverp-1)<0.8;
	    bool goodW2 = W2_L<W2 && W2<W2_H;
	    bool goodCoin = coin_time_L<coin_time && coin_time<coin_time_H;
	    bool goodDy = dy_L<dy && dy<dy_H;
	    bool goodDx = dx_L<dx && dx<dx_H;
	    bool goodSpot = dy_L<dy && dy<dy_H && ((dx_L<dx && dx<dx_H) || (dx_p_L<dx && dx<dx_p_H));
	    bool goodChi2 = track_chi2<15;


	    if(goodChi2){

		    if (goodHelicity && goodMoller && goodVz && goodPS && goodRunRange && 
		    	goodEHCAL && validHel && goodEoverp && goodCoin && goodDy && goodW2){

		    	h_dx_all_hit->Fill(dx);

		    	if(track_nhits == 3){
		    		h_dx_3hit->Fill(dx);
		    	}

		    	if(track_nhits == 4){
		    		h_dx_4hit->Fill(dx);
		    	}

		    	if(track_nhits == 5){
		    		h_dx_5hit->Fill(dx);
		    	}

		    	if(track_nhits >=4){
		    		h_dx_4_5hit->Fill(dx);
		    	}


		    }

		    if (goodHelicity && goodMoller && goodVz && goodPS && goodRunRange && 
		    	goodEHCAL && validHel && goodEoverp && goodCoin && goodSpot){

		    	h_W2_all_hit->Fill(W2);

		    	if(track_nhits == 3){
		    		h_W2_3hit->Fill(W2);
		    	}

		    	if(track_nhits == 4){
		    		h_W2_4hit->Fill(W2);
		    	}

		    	if(track_nhits == 5){
		    		h_W2_5hit->Fill(W2);
		    	}

		    	if(track_nhits >=4){
		    		h_W2_4_5hit->Fill(W2);
		    	}

		    }

		    if (goodHelicity && goodMoller && goodVz && goodPS && goodRunRange && 
		    	goodEHCAL && validHel && goodEoverp && goodSpot && goodW2 ){

		    	h_cointime_all_hit->Fill(coin_time);

		    	if(track_nhits == 3){
		    		h_cointime_3hit->Fill(coin_time);
		    	}

		    	if(track_nhits == 4){
		    		h_cointime_4hit->Fill(coin_time);
		    	}

		    	if(track_nhits == 5){
		    		h_cointime_5hit->Fill(coin_time);
		    	}

		    	if(track_nhits >=4){
		    		h_cointime_4_5hit->Fill(coin_time);
		    	}

		    }


		    if (goodHelicity && goodMoller && goodPS && goodRunRange && 
		    	goodEHCAL && validHel && goodEoverp && goodCoin && goodSpot && goodW2){

		    	h_vz_all_hit->Fill(vz);

		    	if(track_nhits == 3){
		    		h_vz_3hit->Fill(vz);
		    	}

		    	if(track_nhits == 4){
		    		h_vz_4hit->Fill(vz);
		    	}

		    	if(track_nhits == 5){
		    		h_vz_5hit->Fill(vz);
		    	}

		    	if(track_nhits >=4){
		    		h_vz_4_5hit->Fill(vz);
		    	}

		    }
		}

	    if (goodHelicity && goodMoller && goodVz && goodPS && goodRunRange && 
	    	goodEHCAL && validHel && goodEoverp && goodCoin && goodSpot && goodW2){

	    	h_chi2_all_hit->Fill(track_chi2);

	    	if(track_nhits == 3){
	    		h_chi2_3hit->Fill(track_chi2);
	    	}

	    	if(track_nhits == 4){
	    		h_chi2_4hit->Fill(track_chi2);
	    	}

	    	if(track_nhits == 5){
	    		h_chi2_5hit->Fill(track_chi2);
	    	}

	    	if(track_nhits >=4){
	    		h_chi2_4_5hit->Fill(track_chi2);
	    	}


	    }


	    
	    if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
	    std::cout.flush();

	}

	TCanvas *c = new TCanvas("c","c",2400,1500);
	TCanvas *c1 = new TCanvas("c1","c1",2400,1500);

	c->Divide(2,2);
	c->cd(1);
	h_dx_all_hit->SetLineColor(kAzure);
	h_dx_4_5hit->SetLineColor(kGreen);
	h_dx_3hit->SetLineColor(kOrange);
	h_dx_4hit->SetLineColor(kMagenta);
	h_dx_5hit->SetLineColor(kBlack);

	h_dx_all_hit->Draw();
	h_dx_4_5hit->Draw("same");
	h_dx_3hit->Draw("same");
	h_dx_4hit->Draw("same");
	h_dx_5hit->Draw("same");

	c->cd(2);
	// create legend in NDC coordinates (x1,y1,x2,y2)
	auto leg = new TLegend(0.2, 0.2, 0.9, 0.9);
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.04);

	// add entries: (hist, label, option)
	leg->AddEntry(h_dx_all_hit, "All hit tracks",        "l");
	leg->AddEntry(h_dx_4_5hit,    "4 & 5 hit tracks",     "l");
	leg->AddEntry(h_dx_3hit,      "3 hit tracks",         "l");
	leg->AddEntry(h_dx_4hit,      "4 hit tracks",         "l");
	leg->AddEntry(h_dx_5hit,      "5 hit tracks",         "l");

	// finally draw legend
	leg->Draw();

	c->cd(3);
	h_W2_all_hit->SetLineColor(kAzure);
	h_W2_4_5hit->SetLineColor(kGreen);
	h_W2_3hit->SetLineColor(kOrange);
	h_W2_4hit->SetLineColor(kMagenta);
	h_W2_5hit->SetLineColor(kBlack);

	h_W2_all_hit->Draw();
	h_W2_4_5hit->Draw("same");
	h_W2_3hit->Draw("same");
	h_W2_4hit->Draw("same");
	h_W2_5hit->Draw("same");

	c->cd(4);
	h_cointime_all_hit->SetLineColor(kAzure);
	h_cointime_4_5hit->SetLineColor(kGreen);
	h_cointime_3hit->SetLineColor(kOrange);
	h_cointime_4hit->SetLineColor(kMagenta);
	h_cointime_5hit->SetLineColor(kBlack);

	h_cointime_all_hit->Draw();
	h_cointime_4_5hit->Draw("same");
	h_cointime_3hit->Draw("same");
	h_cointime_4hit->Draw("same");
	h_cointime_5hit->Draw("same");


	c1->Divide(2,2);
	c1->cd(1);
	h_vz_all_hit->SetLineColor(kAzure);
	h_vz_4_5hit->SetLineColor(kGreen);
	h_vz_3hit->SetLineColor(kOrange);
	h_vz_4hit->SetLineColor(kMagenta);
	h_vz_5hit->SetLineColor(kBlack);

	h_vz_all_hit->Draw();
	h_vz_4_5hit->Draw("same");
	h_vz_3hit->Draw("same");
	h_vz_4hit->Draw("same");
	h_vz_5hit->Draw("same");

	c1->cd(2);
	h_chi2_all_hit->SetLineColor(kAzure);
	h_chi2_4_5hit->SetLineColor(kGreen);
	h_chi2_3hit->SetLineColor(kOrange);
	h_chi2_4hit->SetLineColor(kMagenta);
	h_chi2_5hit->SetLineColor(kBlack);

	h_chi2_all_hit->Draw();
	h_chi2_4_5hit->Draw("same");
	h_chi2_3hit->Draw("same");
	h_chi2_4hit->Draw("same");
	h_chi2_5hit->Draw("same");


	c->Print(Form("plots/distributions_for_nhits_%s.pdf(",printfilename));
	c1->Print(Form("plots/distributions_for_nhits_%s.pdf)",printfilename));

}
