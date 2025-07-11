#include <iostream>
#include <cmath>
#include "cuts.h"
#include "plotdxdy.h"

void missing_mass(const char* filename, const char* printfilename, const char *kin){

	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	cuts cutsobject;	
	cutsobject.parsecuts(config);

	//std::cout<<"coin_time : "<<cutsobject.coin_time_L<<endl;
	
	double coin_time_L = getDoubleValue(config,"coin_time_L");
	double coin_time_H = getDoubleValue(config,"coin_time_H");

    double coin_time_ac_L = getDoubleValue(config,"coin_time_ac_L");
    double coin_time_ac_H = getDoubleValue(config,"coin_time_ac_H");

	double W2_L = getDoubleValue(config,"W2_L");
	double W2_H = getDoubleValue(config,"W2_H");

	double dy_L = getDoubleValue(config,"dy_L");
	double dy_H = getDoubleValue(config,"dy_H");

	double dx_L = getDoubleValue(config,"dx_L");
	double dx_H = getDoubleValue(config,"dx_H");

	double dx_p_L = getDoubleValue(config,"dx_p_L");
	double dx_p_H = getDoubleValue(config,"dx_p_H");

	double dx_C = getDoubleValue(config,"dx_C");
	double dx_R = getDoubleValue(config,"dx_R");

	double dy_C = getDoubleValue(config,"dy_C");
	double dy_R = getDoubleValue(config,"dy_R");

	double dy_ac_L = getDoubleValue(config,"dy_ac_L");
	double dy_ac_H = getDoubleValue(config,"dy_ac_H");

	double run_num_L = getDoubleValue(config,"run_num_L");
	double run_num_H = getDoubleValue(config,"run_num_H");

	double IHWP_flip = getDoubleValue(config,"IHWP_flip");

	double eHCAL_L = getDoubleValue(config,"eHCAL_L");
	double P_sbs_L = getDoubleValue(config,"P_sbs_L");
	double P_sbs_H = getDoubleValue(config,"P_sbs_H");

    double coin_time_offset_L = coin_time_L+30.0;
    double coin_time_offset_H = coin_time_H+30.0;

    std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;
    std::cout<<"run_num_L: "<< run_num_L <<" run_num_H: "<< run_num_H <<endl;
    std::cout<<"coin_time_L: "<< coin_time_L <<" coin_time_H: "<< coin_time_H <<endl;
    std::cout<<"W2_L: "<< W2_L <<" W2_H: "<< W2_H <<endl;
    std::cout<<"dx_L: "<< dx_L<<" dx_H: "<< dx_H<<endl;
    std::cout<<"dy_L: "<< dy_L<<" dy_H: "<< dy_H<<endl;


	TFile* file = TFile::Open(filename);
	TTree* tree = (TTree*)file->Get("Tout");


	double ebeam = 0.0;
	double trP = 0.0;
	double trPx = 0.0;
	double trPy = 0.0;
	double trPz = 0.0;
	double etheta = 0.0;
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
	double trPx_sbs = 0.0;
	double trPy_sbs = 0.0;
	double trPz_sbs = 0.0;
	double trx_sbs = 0.0;
	double try_sbs = 0.0;
	double trth_sbs = 0.0;
	double trph_sbs = 0.0;
	double vz = 0.0;
	double vz_sbs = 0.0;
	double ntrack_sbs = 0.0;
	double ePS = 0.0;
	double eHCAL = 0.0;
	int runnum = 0.0;

	tree->SetBranchAddress("ebeam", &ebeam);
	tree->SetBranchAddress("trP", &trP);
	tree->SetBranchAddress("trPx", &trPx);
	tree->SetBranchAddress("trPy", &trPy);
	tree->SetBranchAddress("trPz", &trPz);
	tree->SetBranchAddress("etheta", &etheta);
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
	tree->SetBranchAddress("trPx_sbs", &trPx_sbs);
	tree->SetBranchAddress("trPy_sbs", &trPy_sbs);
	tree->SetBranchAddress("trPz_sbs", &trPz_sbs);
	tree->SetBranchAddress("trx_sbs",&trx_sbs);
	tree->SetBranchAddress("try_sbs",&try_sbs);
	tree->SetBranchAddress("trth_sbs",&trth_sbs);
	tree->SetBranchAddress("trph_sbs",&trph_sbs);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("vz_sbs",&vz_sbs);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
	tree->SetBranchAddress("ePS",&ePS);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("runnum",&runnum);

	TH1D* h_MM2 = new TH1D("h_MM2","Missing mass squared distribution ; M^{2} (GeV^{2})",200,-0.5,0.5);
	TH1D* h_MM2_cut_W2 = new TH1D("h_MM2_cut_W2","Missing mass squared distribution ; M^{2} (GeV^{2})",200,-0.5,0.5);
	TH1D* h_MM2_cut_W2_dx = new TH1D("h_MM2_cut_W2_dx","Missing mass squared distribution ; M^{2} (GeV^{2})",200,-0.5,0.5);

	TH1D* h_MM = new TH1D("h_MM","Missing mass distribution ; M (GeV)",200,-1.5,1.5);
	TH1D* h_MM_cut_W2 = new TH1D("h_MM_cut_W2","Missing mass distribution ; M (GeV)",200,-1.5,1.5);
	TH1D* h_MM_cut_W2_dx = new TH1D("h_MM_cut_W2_dx","Missing mass distribution ; M (GeV)",200,-1.5,1.5);

	TH1D* h_EM = new TH1D("h_EM","Missing Energy distribution ; Emiss (GeV)",200,-1.5,1.5);
	TH1D* h_EM_cut_W2 = new TH1D("h_EM_cut_W2","Missing Energy distribution ; Emiss (GeV)",200,-1.5,1.5);
	TH1D* h_EM_cut_W2_dx = new TH1D("h_EM_cut_W2_dx","Missing Energy distribution ; Emiss (GeV)",200,-1.5,1.5);

	TH1D* h_PM = new TH1D("h_PM","Missing Momentum distribution ; Pmiss (GeV)",200,-1.5,1.5);
	TH1D* h_PM_cut_W2 = new TH1D("h_PM_cut_W2","Missing Momentum distribution ; Pmiss (GeV)",200,-1.5,1.5);
	TH1D* h_PM_cut_W2_dx = new TH1D("h_PM_cut_W2_dx","Missing Momentum distribution ; Pmiss (GeV)",200,-1.5,1.5);

	TH1D* h_dx = new TH1D("h_dx","delta-x distribution ; dx (m)",200,-4,3);
	TH1D* h_W2_1 = new TH1D("h_W2_1","W^{2} distribution ; W^{2} (GeV^{2})",200,-2,5);

	TH1D *h_coin_time = new TH1D("h_coin_time","Coincidence time distribution ; coincidence time (ns)", 100,coin_time_L-30, coin_time_H+30);
	TH1D *h_W2 = new TH1D("h_W2","W^{2} distribution ; W^{2} (GeV^{2})",100,-2,4);
	TH2D *h_dxdy = new TH2D("h_dxdy","delta-x delta-y distribution ; dy (m) ; dx (m)",100,-4,3,100,-4,3);

	double Mp = 0.938;

	int nentries = tree->GetEntries();

	for (int i = 0; i<nentries; i++){

		tree->GetEntry(i);

		if(ntrack_sbs>0){

			if(ePS<0.2 || abs(vz)>0.27 || abs(vz_sbs)>0.27 ||  eHCAL<eHCAL_L) continue;
			h_W2->Fill(W2);

			if(W2<W2_H && dx_p_L<dx && dx<dx_p_H && dy_L<dy && dy<dy_H){
				h_coin_time->Fill(coin_time);
			}

			if(coin_time_L<coin_time && coin_time<coin_time_H && W2<W2_H){
				h_dxdy->Fill(dy,dx);
			}

			if(ePS<0.2 || abs(vz)>0.27 || abs(vz_sbs)>0.27 || 
				coin_time<coin_time_L || coin_time>coin_time_H || eHCAL<eHCAL_L 
				|| W2>6 || dy<dy_L || dy>dy_H ) continue;


	        TLorentzVector Pe(0,0,ebeam,ebeam), 
	        PeP(trPx,trPy,trPz,trP), 
	        Pn(0,0,0,Mp),
	        PnP(trPx_sbs,trPy_sbs,trPz_sbs,sqrt(pow(trP_sbs,2)+pow(Mp,2)));
	        
	        TLorentzVector q=Pe-PeP;
	        TLorentzVector Pmiss = Pe - PeP + Pn - PnP;

	        double MM2 = Pmiss.M2();
	        double MM = Pmiss.M();
	        double EM = Pmiss.E();
	        double PM = Pmiss.P();
			
			h_MM2->Fill(MM2);
			
			h_MM->Fill(MM);
			h_EM->Fill(EM);
			h_PM->Fill(PM);
			h_dx->Fill(dx);
			h_W2_1->Fill(W2);


			if(W2>W2_H) continue;
			h_MM_cut_W2->Fill(MM);
			h_MM2_cut_W2->Fill(MM2);

			h_EM_cut_W2->Fill(EM);
			h_PM_cut_W2->Fill(PM);

			if(dx>dx_p_H || dx<dx_p_L) continue;
			h_MM_cut_W2_dx->Fill(MM);
			h_MM2_cut_W2_dx->Fill(MM2);

			h_EM_cut_W2_dx->Fill(EM);
			h_PM_cut_W2_dx->Fill(PM);

		}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();

	}

	TCanvas* c = new TCanvas("c","c",2400,1500);
	TCanvas* c1 = new TCanvas("c1","c1",2400,1500);
	TCanvas* c2 = new TCanvas("c2","c2",2400,1500);

	c->Divide(2,2);
	c->cd(1);
	//h_MM2->SetMarkerStyle(8);
	//h_MM2_cut_W2->SetMarkerStyle(7);
	//h_MM2_cut_W2_dx->SetMarkerStyle(3);

	h_MM2->SetLineColor(kAzure);
	h_MM2_cut_W2->SetLineColor(kViolet);
	h_MM2_cut_W2_dx->SetLineColor(kBlack);

	h_MM2->Draw("hist");
	h_MM2_cut_W2->Draw("hist same");
	h_MM2_cut_W2_dx->Draw("hist same");

	auto leg = new TLegend(0.55, 0.65, 0.88, 0.85); // x1,y1,x2,y2 (NDC)
	leg->SetBorderSize(0);          // no frame
	leg->SetFillStyle(0);           // transparent
	leg->SetTextSize(0.04);         // adjust to taste

	leg->AddEntry(h_MM2,             "global + coin_time + dy cuts", "l");
	leg->AddEntry(h_MM2_cut_W2,      "+ W^{2} cut",            "l");
	leg->AddEntry(h_MM2_cut_W2_dx,   "+ W^{2} + dx cuts",      "l");

	leg->Draw();

	c->cd(2);
	h_MM->SetLineColor(kAzure);
	h_MM_cut_W2->SetLineColor(kViolet);
	h_MM_cut_W2_dx->SetLineColor(kBlack);

	h_MM->Draw("hist");
	h_MM_cut_W2->Draw("hist same");
	h_MM_cut_W2_dx->Draw("hist same");

	auto leg1 = new TLegend(0.55, 0.65, 0.88, 0.85); // x1,y1,x2,y2 (NDC)
	leg1->SetBorderSize(0);          // no frame
	leg1->SetFillStyle(0);           // transparent
	leg1->SetTextSize(0.04);         // adjust to taste

	leg1->AddEntry(h_MM,             "global + coin_time + dy cuts", "l");
	leg1->AddEntry(h_MM_cut_W2,      "+ W^{2} cut",            "l");
	leg1->AddEntry(h_MM_cut_W2_dx,   "+ W^{2} + dx cuts",      "l");

	leg1->Draw();

	c->cd(3);
	h_dx->Draw();
	
	c->cd(4);
	h_W2_1->Draw();

	c1->Divide(2,2);
	c1->cd(1);
	h_EM->SetLineColor(kAzure);
	h_EM_cut_W2->SetLineColor(kViolet);
	h_EM_cut_W2_dx->SetLineColor(kBlack);

	h_EM->Draw("hist");
	h_EM_cut_W2->Draw("hist same");
	h_EM_cut_W2_dx->Draw("hist same");

	auto leg2 = new TLegend(0.55, 0.65, 0.88, 0.85); // x1,y1,x2,y2 (NDC)
	leg2->SetBorderSize(0);          // no frame
	leg2->SetFillStyle(0);           // transparent
	leg2->SetTextSize(0.04);         // adjust to taste

	leg2->AddEntry(h_EM,             "global + coin_time + dy cuts", "l");
	leg2->AddEntry(h_EM_cut_W2,      "+ W^{2} cut",            "l");
	leg2->AddEntry(h_EM_cut_W2_dx,   "+ W^{2} + dx cuts",      "l");

	leg2->Draw();


	c1->cd(2);
	h_PM->SetLineColor(kAzure);
	h_PM_cut_W2->SetLineColor(kViolet);
	h_PM_cut_W2_dx->SetLineColor(kBlack);

	h_PM->Draw("hist");
	h_PM_cut_W2->Draw("hist same");
	h_PM_cut_W2_dx->Draw("hist same");

	auto leg3 = new TLegend(0.55, 0.65, 0.88, 0.85); // x1,y1,x2,y2 (NDC)
	leg3->SetBorderSize(0);          // no frame
	leg3->SetFillStyle(0);           // transparent
	leg3->SetTextSize(0.04);         // adjust to taste

	leg3->AddEntry(h_PM,             "global + coin_time + dy cuts", "l");
	leg3->AddEntry(h_PM_cut_W2,      "+ W^{2} cut",            "l");
	leg3->AddEntry(h_PM_cut_W2_dx,   "+ W^{2} + dx cuts",      "l");

	leg3->Draw();


	c2->Divide(2,2);
	c2->cd(1);
	TLine *linecoinL = new TLine(coin_time_L, 0, coin_time_L, h_coin_time->GetMaximum()*0.8);
	TLine *linecoinH = new TLine(coin_time_H, 0, coin_time_H, h_coin_time->GetMaximum()*0.8);
	h_coin_time->Draw();
	linecoinL->Draw("same");
	linecoinH->Draw("same");

	c2->cd(2);
	TLine *lineW2 = new TLine(W2_H, 0, W2_H, h_W2->GetMaximum()*0.8);
	h_W2->Draw();
	lineW2->Draw("same");


	c2->cd(3);
	TBox *boxdxdy = new TBox(dy_L,dx_p_L,dy_H,dx_p_H);
	boxdxdy->SetFillStyle(0);
	h_dxdy->Draw("COLZ");
	boxdxdy->Draw("same");

	c->SaveAs(Form("../plots/%s_missing_mass.pdf(",printfilename));
	c1->SaveAs(Form("../plots/%s_missing_mass.pdf",printfilename));
	c2->SaveAs(Form("../plots/%s_missing_mass.pdf)",printfilename));
}