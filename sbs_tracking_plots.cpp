#include <iostream>
#include <cmath>
#include "cuts.h"
#include "plotdxdy.h"

void sbs_tracking_plots(const char* filename, const char* printfilename, const char *kin){

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
	double ptheta = 0.0;
	double pphi = 0.0; 
	double ptheta_cal = 0.0;
	double pphi_cal = 0.0;
	double pN_expect = 0.0;
	double trP_sbs = 0.0;
	double trx_sbs = 0.0;
	double try_sbs = 0.0;
	double trth_sbs = 0.0;
	double trph_sbs = 0.0;
	double vz = 0.0;
	double vz_sbs = 0.0;
	double ntrack_sbs = 0.0;

	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("ptheta",&ptheta);
	tree->SetBranchAddress("pphi",&pphi);
	tree->SetBranchAddress("ptheta_cal",&ptheta_cal);
	tree->SetBranchAddress("pphi_cal",&pphi_cal);
	tree->SetBranchAddress("pN_expect",&pN_expect);
	tree->SetBranchAddress("trP_sbs",&trP_sbs);
	tree->SetBranchAddress("trx_sbs",&trx_sbs);
	tree->SetBranchAddress("try_sbs",&try_sbs);
	tree->SetBranchAddress("trth_sbs",&trth_sbs);
	tree->SetBranchAddress("trph_sbs",&trph_sbs);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("vz_sbs",&vz_sbs);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);

	TH1D *h_delta_ptheta = new TH1D("h_delta_ptheta","h_delta_ptheta",200,-0.1,0.1);
	TH1D *h_delta_pphi = new TH1D("h_delta_pphi","h_delta_pphi",200,-0.1,0.1);
	TH1D *h_delta_pmag = new TH1D("h_delta_pmag","h_delta_pmag",200,-1,1);
	
	TH1D *h_pmag_ratio = new TH1D("h_pmag_ratio","h_pmag_ratio",200,0,2);
	
	TH2D *h_ptheta_corr = new TH2D("h_ptheta_corr","h_ptheta_corr",200,0.25,0.5,200,0.25,0.5);
	TH2D *h_pphi_corr = new TH2D("h_pphi_corr","h_pphi_corr",200,2.8,3.6,200,-0.3,0.45);
	TH2D *h_pmag_corr = new TH2D("h_pmag_corr","h_pmag_corr",200,4,5,200,3,5);

	TH2D *h_delta_ptheta_v_trthsbs = new TH2D("h_delta_ptheta_v_trthsbs","h_delta_ptheta_v_trthsbs",500,-0.23,0,500,-0.1,-0.1);
	TH2D *h_delta_ptheta_v_trphsbs = new TH2D("h_delta_ptheta_v_trphsbs","h_delta_ptheta_v_trphsbs",500,-0.04,0.06,500,-0.1,0.1);
	TH2D *h_delta_ptheta_v_trxsbs = new TH2D("h_delta_ptheta_v_trxsbs","h_delta_ptheta_v_trxsbs",500,-0.3,0.8,500,-0.1,0.1);
	TH2D *h_delta_ptheta_v_trysbs = new TH2D("h_delta_ptheta_v_trysbs","h_delta_ptheta_v_trysbs",500,-0.25,0.25,500,-0.1,0.1);

	TH2D *h_delta_pphi_v_trthsbs = new TH2D("h_delta_pphi_v_trthsbs","h_delta_pphi_v_trthsbs",500,-0.23,0.0,500,-0.1,0.1);
	TH2D *h_delta_pphi_v_trphsbs = new TH2D("h_delta_pphi_v_trphsbs","h_delta_pphi_v_trphsbs",500,-0.04,0.06,500,-0.1,0.1);
	TH2D *h_delta_pphi_v_trxsbs = new TH2D("h_delta_pphi_v_trxsbs","h_delta_pphi_v_trxsbs",500,-0.3,0.8,500,-0.1,0.1);
	TH2D *h_delta_pphi_v_trysbs = new TH2D("h_delta_pphi_v_trysbs","h_delta_pphi_v_trysbs",500,-0.25,0.25,500,-0.1,0.1);

	TH2D *h_delta_pmag_v_trthsbs = new TH2D("h_delta_pmag_v_trthsbs","h_delta_pmag_v_trthsbs",500,-0.23,0.0,500,-1,1);
	TH2D *h_delta_pmag_v_trphsbs = new TH2D("h_delta_pmag_v_trphsbs","h_delta_pmag_v_trphsbs",500,-0.04,0.06,500,-1,1);
	TH2D *h_delta_pmag_v_trxsbs = new TH2D("h_delta_pmag_v_trxsbs","h_delta_pmag_v_trxsbs",500,-0.3,0.8,500,-1,1);
	TH2D *h_delta_pmag_v_trysbs = new TH2D("h_delta_pmag_v_trysbs","h_delta_pmag_v_trysbs",500,-0.25,0.25,500,-1,1);

	TH2D *h_delta_vz_v_trthsbs = new TH2D("h_delta_vz_v_trthsbs","h_delta_vz_v_trthsbs",500,-0.23,0.0,500,-0.5,0.5);
	TH2D *h_delta_vz_v_trphsbs = new TH2D("h_delta_vz_v_trphsbs","h_delta_vz_v_trphsbs",500,-0.04,0.06,500,-0.5,0.5);
	TH2D *h_delta_vz_v_trxsbs = new TH2D("h_delta_vz_v_trxsbs","h_delta_vz_v_trxsbs",500,-0.3,0.8,500,-0.5,0.5);
	TH2D *h_delta_vz_v_trysbs = new TH2D("h_delta_vz_v_trysbs","h_delta_vz_v_trysbs",500,-0.25,0.25,500,-0.5,0.5);

	TH1D *h_coin_time = new TH1D("h_coin_time","h_coin_time",1000,75,120);
	TH1D *h_W2 = new TH1D("h_W2","h_W2",1000,-2,4);
	TH2D *h_dxdy = new TH2D("h_dxdy","h_dxdy",200,-2,2,200,-6,2);

	bool cutcoin_time = false;
	bool cutW2_coin_time = false;
	bool cutQE = false;

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		//no QE cuts for coin_time plot
		if(ntrack_sbs>0){
		h_coin_time->Fill(coin_time);

		//coin_time cut only
		cutcoin_time = abs(coin_time-95)<8;
		if (cutcoin_time) {
			h_W2->Fill(W2);
			h_delta_vz_v_trthsbs->Fill(trth_sbs,(vz_sbs-vz));
			h_delta_vz_v_trphsbs->Fill(trph_sbs,(vz_sbs-vz));
			h_delta_vz_v_trxsbs->Fill(trx_sbs,(vz_sbs-vz));
			h_delta_vz_v_trysbs->Fill(try_sbs,(vz_sbs-vz));
		};
		
		//coin_time and W2 cuts only
		cutW2_coin_time = abs(W2-0.88)<0.9&&abs(coin_time-95)<8;
		if (cutW2_coin_time) h_dxdy->Fill(dy,dx);

		cutQE = (abs(W2-0.88)<0.9&&abs(coin_time-95)<8&&(sqrt(pow((dx+1.55)/0.5,2)+pow((dy-0.1)/0.5,2)))<1);
		if (cutQE){
			h_delta_ptheta->Fill(ptheta-ptheta_cal);
			h_delta_pphi->Fill(pphi-pphi_cal+M_PI);
			h_delta_pmag->Fill(pN_expect-trP_sbs);
	
			h_pmag_ratio->Fill(pN_expect/trP_sbs);

			h_ptheta_corr->Fill(ptheta_cal,ptheta);
			h_pphi_corr->Fill(pphi_cal,pphi);
			h_pmag_corr->Fill(pN_expect,trP_sbs);

			h_delta_ptheta_v_trthsbs->Fill(trth_sbs,ptheta-ptheta_cal);
			h_delta_ptheta_v_trphsbs->Fill(trph_sbs,ptheta-ptheta_cal);
			h_delta_ptheta_v_trxsbs->Fill(trx_sbs,ptheta-ptheta_cal);
			h_delta_ptheta_v_trysbs->Fill(try_sbs,ptheta-ptheta_cal);

			h_delta_pphi_v_trthsbs->Fill(trth_sbs,pphi-pphi_cal+M_PI);
			h_delta_pphi_v_trphsbs->Fill(trph_sbs,pphi-pphi_cal+M_PI);
			h_delta_pphi_v_trxsbs->Fill(trx_sbs,pphi-pphi_cal+M_PI);
			h_delta_pphi_v_trysbs->Fill(try_sbs,pphi-pphi_cal+M_PI);

			h_delta_pmag_v_trthsbs->Fill(trth_sbs,pN_expect-trP_sbs);
			h_delta_pmag_v_trphsbs->Fill(trph_sbs,pN_expect-trP_sbs);
			h_delta_pmag_v_trxsbs->Fill(trx_sbs,pN_expect-trP_sbs);
			h_delta_pmag_v_trysbs->Fill(try_sbs,pN_expect-trP_sbs);

		}
		}
		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();

	}
	TCanvas* c = new TCanvas("c","c",2400,1800);
	TCanvas* c1 = new TCanvas("c1","c1",2400,1800);
	TCanvas* c2 = new TCanvas("c2","c2",2400,1800);
	TCanvas* c3 = new TCanvas("c3","c3",2400,1800);
	TCanvas* c4 = new TCanvas("c4","c4",2400,1800);
	TCanvas* c5 = new TCanvas("c5","c5",2400,1800);
	TCanvas* c6 = new TCanvas("c6","c6",2400,1800);
	TCanvas* c7 = new TCanvas("c7","c7",2400,1800);

	c->Divide(2,2);
	c1->Divide(2,2);
	c2->Divide(2,2);
	c3->Divide(2,2);
	c4->Divide(2,2);
	c5->Divide(2,2);
	c6->Divide(2,2);
	c7->Divide(2,2);

	c->cd(1);
	h_ptheta_corr->Draw();
	h_ptheta_corr->SetXTitle("ptheta_cal (rad)");
	h_ptheta_corr->SetYTitle("ptheta (rad)");

	c1->cd(1);
	h_pphi_corr->Draw();
	h_pphi_corr->SetXTitle("pphi_cal (rad)");
	h_pphi_corr->SetYTitle("pphi (rad)");

	c2->cd(1);
	h_pmag_corr->Draw();
	h_pmag_corr->SetXTitle("pN_expect (GeV)");
	h_pmag_corr->SetYTitle("trP_sbs (GeV)");

	c3->cd(1);
	h_coin_time->Draw();
	h_coin_time->SetXTitle("coin_time (ns)");
	TLine *line01 = new TLine(95-8, 0, 95-8, 275);
    TLine *line02 = new TLine(95+8, 0, 95+8, 275);
    line01->SetLineColor(kRed);
    line02->SetLineColor(kRed);
    line01->SetLineWidth(2);
    line02->SetLineWidth(2);
    line01->Draw("same");
    line02->Draw("same");

	c->cd(2);
	h_delta_ptheta->Draw();
	h_delta_ptheta->SetXTitle("ptheta-ptheta_cal (rad)");

	c1->cd(2);
	h_delta_pphi->Draw();
	h_delta_pphi->SetXTitle("pphi-pphi_cal+pi (rad)");

	c2->cd(2);
	h_delta_pmag->Draw();
	h_delta_pmag->SetXTitle("pN_expect-trP_sbs (GeV)");

	c3->cd(2);
	h_W2->Draw();
	h_W2->SetXTitle("W2(GeV2)");
    	TLine *line1 = new TLine(0.88-0.9, 0, 0.88-0.9, 300);
    	TLine *line2 = new TLine(0.88+0.9, 0, 0.88+0.9, 300);
    	line1->SetLineColor(kRed);
    	line2->SetLineColor(kRed);
    	line1->SetLineWidth(2);
    	line2->SetLineWidth(2);
    	line1->Draw("same");
    	line2->Draw("same");

	c2->cd(3);
	h_pmag_ratio->Draw();
	h_pmag_ratio->SetXTitle("pN_expect/trP_sbs");
	

	c3->cd(3);
	h_dxdy->Draw("COLZ");
	h_dxdy->SetXTitle("dy(m)");
	h_dxdy->SetYTitle("dx(m)");
	TEllipse *ellipse = new TEllipse(0.1, -1.55, 0.5, 0.5);
    ellipse->SetLineColor(kRed);
    ellipse->SetLineWidth(2);
    ellipse->SetFillStyle(0); // No fill, just the outline
    ellipse->Draw("same");

    //c4
    c4->cd(1);
    h_delta_ptheta_v_trthsbs->Draw("COLZ");
    h_delta_ptheta_v_trthsbs->SetXTitle("trth_sbs (rad)");
    h_delta_ptheta_v_trthsbs->SetYTitle("ptheta-ptheta_cal (rad)");

    c4->cd(2);
    h_delta_ptheta_v_trphsbs->Draw("COLZ");
    h_delta_ptheta_v_trphsbs->SetXTitle("trph_sbs (rad)");
    h_delta_ptheta_v_trphsbs->SetYTitle("ptheta-ptheta_cal (rad)");

    c4->cd(3);
    h_delta_ptheta_v_trxsbs->Draw("COLZ");
    h_delta_ptheta_v_trxsbs->SetXTitle("trx_sbs (m)");
    h_delta_ptheta_v_trxsbs->SetYTitle("ptheta-ptheta_cal (rad)");

    c4->cd(4);
    h_delta_ptheta_v_trysbs->Draw("COLZ");
    h_delta_ptheta_v_trysbs->SetXTitle("try_sbs (m)");
    h_delta_ptheta_v_trysbs->SetYTitle("ptheta-ptheta_cal (rad)");

    //c5
    c5->cd(1);
    h_delta_pphi_v_trthsbs->Draw("COLZ");
    h_delta_pphi_v_trthsbs->SetXTitle("trth_sbs (rad)");
    h_delta_pphi_v_trthsbs->SetYTitle("pphi-pphi_cal+pi (rad)");

    c5->cd(2);
    h_delta_pphi_v_trphsbs->Draw("COLZ");
    h_delta_pphi_v_trphsbs->SetXTitle("trph_sbs (rad)");
    h_delta_pphi_v_trphsbs->SetYTitle("pphi-pphi_cal+pi (rad)");
    
    c5->cd(3);
    h_delta_pphi_v_trxsbs->Draw("COLZ");
    h_delta_pphi_v_trxsbs->SetXTitle("trx_sbs (m)");
    h_delta_pphi_v_trxsbs->SetYTitle("pphi-pphi_cal+pi (rad)");

    c5->cd(4);
    h_delta_pphi_v_trysbs->Draw("COLZ");
    h_delta_pphi_v_trysbs->SetXTitle("try_sbs (m)");
    h_delta_pphi_v_trysbs->SetYTitle("pphi-pphi_cal+pi (rad)");


    //c6
    c6->cd(1);
    h_delta_pmag_v_trthsbs->Draw("COLZ");
    h_delta_pmag_v_trthsbs->SetXTitle("trth_sbs (rad)");
    h_delta_pmag_v_trthsbs->SetYTitle("pN_expect-trP_sbs (GeV)");

    c6->cd(2);
    h_delta_pmag_v_trphsbs->Draw("COLZ");
    h_delta_pmag_v_trphsbs->SetXTitle("trph_sbs (rad)");
    h_delta_pmag_v_trphsbs->SetYTitle("pN_expect-trP_sbs (GeV)");

    c6->cd(3);
    h_delta_pmag_v_trxsbs->Draw("COLZ");
    h_delta_pmag_v_trxsbs->SetXTitle("trx_sbs (m)");
    h_delta_pmag_v_trxsbs->SetYTitle("pN_expect-trP_sbs (GeV)");

    c6->cd(4);
    h_delta_pmag_v_trysbs->Draw("COLZ");
    h_delta_pmag_v_trysbs->SetXTitle("try_sbs (m)");
    h_delta_pmag_v_trysbs->SetYTitle("pN_expect-trP_sbs (GeV)");

    //c7
    c7->cd(1);
    h_delta_vz_v_trthsbs->Draw("COLZ");
    h_delta_vz_v_trthsbs->SetXTitle("trth_sbs (rad)");
    h_delta_vz_v_trthsbs->SetYTitle("vz_sbs - vz_bb (m)");
    
    c7->cd(2);
    h_delta_vz_v_trphsbs->Draw("COLZ");
    h_delta_vz_v_trphsbs->SetXTitle("trph_sbs (rad)");
    h_delta_vz_v_trphsbs->SetYTitle("vz_sbs - vz_bb (m)");

    c7->cd(3);
    h_delta_vz_v_trxsbs->Draw("COLZ");
    h_delta_vz_v_trxsbs->SetXTitle("trx_sbs (m)");
    h_delta_vz_v_trxsbs->SetYTitle("vz_sbs - vz_bb (m)");

    c7->cd(4);
    h_delta_vz_v_trysbs->Draw("COLZ");
    h_delta_vz_v_trysbs->SetXTitle("try_sbs (m)");
    h_delta_vz_v_trysbs->SetYTitle("vz_sbs - vz_bb (m)");

	c->Print(Form("../plots/%s_sbstracking.pdf(",printfilename));
	c1->Print(Form("../plots/%s_sbstracking.pdf",printfilename));
	c2->Print(Form("../plots/%s_sbstracking.pdf",printfilename));
	c3->Print(Form("../plots/%s_sbstracking.pdf",printfilename));
	c4->Print(Form("../plots/%s_sbstracking.pdf",printfilename));
	c5->Print(Form("../plots/%s_sbstracking.pdf",printfilename));
	c6->Print(Form("../plots/%s_sbstracking.pdf",printfilename));
	c7->Print(Form("../plots/%s_sbstracking.pdf)",printfilename));


	c->SaveAs(Form("../plots/%s_ptheta_sbstracking.png",printfilename));
	c1->SaveAs(Form("../plots/%s_pphi_sbstracking.png",printfilename));
	c2->SaveAs(Form("../plots/%s_pmag_sbstracking.png",printfilename));
	c3->SaveAs(Form("../plots/%s_QE_cuts_sbstracking.png",printfilename));
	c4->SaveAs(Form("../plots/%s_deltaptheta_sbstracking.png",printfilename));
	c5->SaveAs(Form("../plots/%s_deltapphi_sbstracking.png",printfilename));
	c6->SaveAs(Form("../plots/%s_deltapmag_sbstracking.png",printfilename));
	c7->SaveAs(Form("../plots/%s_deltavz_sbstracking.png",printfilename));
}
