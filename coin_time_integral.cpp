// this is quick script for background subtracted coin_time integrals (side band method)
#include "cuts.h"
#include "plotdxdy.h"

void coin_time_integral(const char* filename, const char* printfilename, const char* kin, double bL, double bH, double sbL, double sbH){
	
	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	cuts cutsobject;	
	cutsobject.parsecuts(config);

	std::cout<<"coin_time : "<<cutsobject.coin_time_L<<endl;
	
	TFile* file = TFile::Open(filename);
	TTree* tree = (TTree*)file->Get("Tout");

	double vz = 0.0;
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

	tree->SetBranchAddress("vz",&vz);
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

	TH1D* h_coin_time = new TH1D("h_coin_time" , "coincidence time (ns) ; coincidence time (ns)", 200, 150, 250);
	TH1D* h_dx = new TH1D("h_dx","delta-x distribution (m); dx (m)", 100, -4, 3);
	TH1D* h_W2 = new TH1D("h_W2","W^{2} distribution (GeV^{2}) ; W^{2} (GeV^{2})", 200, -3, 3);

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries*1; i++){
		tree->GetEntry(i);
		if ( abs(vz)>0.27 || ePS<0.2 || eHCAL<0.025 /*|| runnum<5570*/) {
			continue;
		}

		if (W2>3 || abs(dy)>1 ) continue;

		h_coin_time->Fill(coin_time);

		if (abs(dy)>0.5 || abs(coin_time-187)>7) continue;
		h_W2->Fill(W2);


		if(W2<-0.23 || W2>1.43 || abs(dy)>0.5 || abs(coin_time-184.5)>3.5 || eHCAL<0.325) continue;
		h_dx->Fill(dx);

	}

	double PeakIntegral = h_coin_time->Integral(h_coin_time->FindBin(bL),h_coin_time->FindBin(bH));
	double SideBandIntegral = h_coin_time->Integral(h_coin_time->FindBin(sbL),h_coin_time->FindBin(sbH));

	double QE_count = PeakIntegral - SideBandIntegral;


	TCanvas*c = new TCanvas("c","c",2400,1500);
	c->Divide(2,2);
	c->cd(1);
	TLine *linebL = new TLine(bL,0,bL,0.9*h_coin_time->GetMaximum());
	TLine *linebH = new TLine(bH,0,bH,0.9*h_coin_time->GetMaximum());

	TLine *linesbL = new TLine(sbL,0,sbL,0.9*h_coin_time->GetMaximum());
	TLine *linesbH = new TLine(sbH,0,sbH,0.9*h_coin_time->GetMaximum());

	linebL->SetLineColor(kBlue);
	linebH->SetLineColor(kBlue);

	linesbL->SetLineColor(kRed);
	linesbH->SetLineColor(kRed);

	h_coin_time->Draw();
	linebL->Draw("same");
	linebH->Draw("same");
	linesbL->Draw("same");
	linesbH->Draw("same");

	c->cd(2);
	// create a little white box in NDC coordinates:
	auto pt = new TPaveText(0.1, 0.6, 0.9, 0.9, "NDC");
	pt->SetFillColor(kWhite);
	pt->SetTextFont(42);
	pt->SetTextSize(0.05);
	pt->SetTextAlign(12);  // left/bottom alignment

	pt->AddText( Form("PeakIntegral     = %.1f", PeakIntegral) );
	pt->AddText( Form("SideBandIntegral = %.1f", SideBandIntegral) );
	pt->AddText( Form("QE count         = %.1f", QE_count) );

	pt->Draw();
	c->Update();

	c->cd(3);
	h_dx->Draw();

	c->cd(4);
	h_W2->Draw();

	c->Print(Form("%s.pdf",printfilename));

}