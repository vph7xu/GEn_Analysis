#include "cuts.h"
#include "plotdxdy.h"
#include <cmath>
#include <fstream>


void Hcal_fsam(const char* filename, const char* printfilename, const char* kin){

	std::map<std::string, std::string> config = parseConfig(Form("cuts/calib/cut_%s.txt",kin)); //parse the cuts
    cuts cutsobject;
    cutsobject.parsecuts(config);

    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

	double eHCAL = 0.0;
	double xHCAL = 0.0;
	double yHCAL = 0.0;
	double Q2 = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;
	double trP_sbs = 0.0;
	double nblk_HCAL = 0.0;
	double ntrack_sbs = 0.0;
	double vz = 0.0;
	double vz_sbs = 0.0;
	double hcal_clus_id[1000];
	double hcal_clus_mem_id[1000];
	double hcal_clus_mem_e[1000];
	int runnum = 0;
	
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("xHCAL",&xHCAL);
	tree->SetBranchAddress("yHCAL",&yHCAL);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("trP_sbs",&trP_sbs);
	tree->SetBranchAddress("hcal_clus_id",&hcal_clus_id);
	tree->SetBranchAddress("hcal_clus_mem_id",&hcal_clus_mem_id);
	tree->SetBranchAddress("hcal_clus_mem_e",&hcal_clus_mem_e);
	tree->SetBranchAddress("nblk_HCAL",&nblk_HCAL);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("vz_sbs",&vz_sbs);
	tree->SetBranchAddress("runnum",&runnum);

	double mp = 0.938272;
	double fs = 0;
	bool elastic_cut = false;
	int nentries = tree->GetEntries();

	TH1D *hfs = new TH1D("hfs","Sampling fraction distribution",100,0,0.5);
	TH1D *her = new TH1D("her","Actual energy dep/expected energy dep",200,0.0,0.4);
	TH1D *hpN = new TH1D("hpN","Outgoing nucleon momentum",1000,-0.1,4);

	TH2D *hfsample = new TH2D("hfsample","Sampling fraction vs xHCAL distribution",200,-3,1.5,200,0,0.3);
	TH2D *hdxdy = new TH2D("hdxdy","dx vs dy distribution",400,-2,2,800,-4,4);
	TH2D *heratio = new TH2D("heratio", "Actual energy dep/expected energy dep vs xHCAL",200,-3,1.5,200,0.0,0.4);

	TH2D *hfsvblkid = new TH2D("hfsvblkid","sampling fraction vs blk id",300,0,300,200,0.025,0.3);

	TH1D *hcointime = new TH1D("hcointime", "cointime distribution",1000,40,200);
	TH1D *hW2 = new TH1D("hW2","W2 distribution",1000,-4,4);

	TH2D *hvertexcorr = new TH2D("hvertexcorr", "vertex correlation", 100,-0.5,0.5,100,-0.5,0.5);
	TH1D *hdeltavertex = new TH1D("hdeltavertex","delta vertex", 100, -0.5, 0.5);

	TGraphErrors *graph = new TGraphErrors(hfsample->GetNbinsX());
	
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		
		if (ntrack_sbs>0){
		
			elastic_cut = (abs(coin_time-cutsobject.coin_time_mean)<cutsobject.coin_time_width and abs(vz-vz_sbs)<0.1 /*and 0<W2 and W2<5*/);

			//elastic_cut = (((pow((dy-cutsobject.dy_C)/cutsobject.dy_R,2)+pow((dx-cutsobject.dx_C)/cutsobject.dx_R,2))<=2.5)
			//&&(abs(coin_time-cutsobject.coin_time_mean)<cutsobject.coin_time_width)&&(abs(W2-cutsobject.W2_mean)<cutsobject.W2_width));//? true:false;
			
			//std::cout<<"elastic cut : "<<elastic_cut<<endl;
			hcointime->Fill(coin_time);
			hvertexcorr->Fill(vz,vz_sbs);
			hdeltavertex->Fill(vz-vz_sbs);


			if (abs(coin_time-cutsobject.coin_time_mean)<cutsobject.coin_time_width and abs(vz-vz_sbs)<0.1 /*and 0<W2 and W2<5*/){

				hW2->Fill(W2);
			
			}


			if((abs(coin_time-cutsobject.coin_time_mean)<cutsobject.coin_time_width)and abs(vz-vz_sbs)<0.1 /*and 0<W2 and W2<5*/){//&&(abs(W2-cutsobject.W2_mean)<cutsobject.W2_width)){
				
				hdxdy->Fill(dy,dx);

			}

			if(elastic_cut){
				fs = eHCAL*2*mp/Q2;	
				//fs = eHCAL*2*mp/(trP_sbs*trP_sbs);
				hfs->Fill(fs);
				hfsample->Fill(xHCAL,fs);
				heratio->Fill(xHCAL,eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp));
				her->Fill(eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp));
				hpN->Fill(trP_sbs);

				for (int i = 0; i<nblk_HCAL; i++){
					//hfsvblkid->Fill(hcal_clus_id[i],eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp));
					//hfsvblkid->Fill(hcal_clus_mem_id[i],eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp));
					hfsvblkid->Fill(hcal_clus_mem_id[i],eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp));
				}

				hdxdy->Fill(dy,dx);
			}
		}
                if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();
	}




	//for (int i = 1; i<=(hfsample->GetNbinsX());i++){

		//TH1D *hfsample_yproj = hfsample->ProjectionY("_py",i,i);

		//TF1 * gaussFit = new TF1("gaussFit","gaus",0.02,0.15);
		//hfsample_yproj->Fit(gaussFit,"Q");
		

		//double meanY = /*gaussFit->GetParameter(1);*/hfsample_yproj->GetMean();
		//double stdY = /*gaussFit->GetParameter(2);*/hfsample_yproj->GetStdDev();
		
	//	double x = hfsample->GetXaxis()->GetBinCenter(i);

	//	graph->SetPoint(i-1,x,meanY);
	//	graph->SetPointError(i-1,0,stdY);

	//	delete hfsample_yproj;
	//}

	TProfile *profX = hfsample->ProfileX();
	TProfile *profX1 = heratio->ProfileX();
	TProfile *profX2 = hfsvblkid->ProfileX();

	// get the sampling fraction for each bin written to txt file with the index
	std::ofstream outFile(Form("txt/%s_sampling_fractions_each_blk.txt",kin));
	double sampling_fraction_means [300];

	for (int i = 1; i <=profX2->GetNbinsX();++i){

		sampling_fraction_means[i-1] = profX2->GetBinContent(i);
		outFile<<i-1<<" "<<sampling_fraction_means[i-1]<<"\n";
	}


	TCutG *cut_elipse = CreateOvalCut("cut_elipse",cutsobject.dy_C,cutsobject.dx_C,cutsobject.dy_R*sqrt(2.5),cutsobject.dx_R*sqrt(2.5),100); 

	TCanvas *c = new TCanvas("c","c",3200,2400);
	TCanvas *c1 = new TCanvas("c1","c1",3200,2400);
	TCanvas *c2 = new TCanvas("c2","c2",3200,2400);
	TCanvas *c3 = new TCanvas("c3","c3",3200,2400);
	TCanvas *c4 = new TCanvas("c4", "c4", 3200, 2400);
	
	c->Divide(1,2);
	c->cd(1);
	hfsample->SetXTitle("xHCAL (m)");
	hfsample->SetYTitle("eHCAL*2*mp/Q^2");
	hfsample->Draw("COLZ");
	profX->SetMarkerStyle(20);
	profX->SetLineColor(kBlack);
	profX->Draw("SAME");
	//graph->SetMarkerStyle(20); // Set marker style
    //graph->SetMarkerColor(kRed); // Set marker color
    //graph->SetLineColor(kRed); // Set line color
    //graph->Draw("P SAME");

	c->cd(2);
	hfs->SetXTitle("eHCAL*2*mp/Q^2");
	hfs->Draw();

	c1->Divide(2,2);
	
	c1->cd(1);
	hcointime->SetXTitle("cointime (ns)");
	hcointime->Draw();
	TLine *line01 = new TLine(cutsobject.coin_time_mean-cutsobject.coin_time_width, 0, cutsobject.coin_time_mean-cutsobject.coin_time_width, 2750);
    TLine *line02 = new TLine(cutsobject.coin_time_mean+cutsobject.coin_time_width, 0, cutsobject.coin_time_mean+cutsobject.coin_time_width, 2750);
    line01->SetLineColor(kRed);
    line02->SetLineColor(kRed);
    line01->SetLineWidth(2);
    line02->SetLineWidth(2);
    line01->Draw("same");
    line02->Draw("same");

	c1->cd(2);
	TLine *line1 = new TLine(cutsobject.W2_mean-cutsobject.W2_width, 0, cutsobject.W2_mean-cutsobject.W2_width, 500);
    TLine *line2 = new TLine(cutsobject.W2_mean+cutsobject.W2_width, 0, cutsobject.W2_mean+cutsobject.W2_width, 500);
	hW2->Draw();
	hW2->SetXTitle("W^2 (GeV^2)");
	line1->SetLineColor(kRed);
    line2->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line2->SetLineWidth(2);
    line1->Draw("same");
    line2->Draw("same");

	c1->cd(3);
	hdxdy->Draw("COLZ");
	hdxdy->SetXTitle("dy (m)");
	hdxdy->SetYTitle("dx (m)");
	cut_elipse->Draw("L");

	c1->cd(4);
	hpN->Draw();

	c2->Divide(1,2);
	c2->cd(1);
	heratio->SetYTitle("eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp)");
	heratio->SetXTitle("xHCAL (m)");
	heratio->Draw("COLZ");
	profX1->SetMarkerStyle(20);
	profX1->SetLineColor(kBlack);
	profX1->Draw("SAME");

	c2->cd(2);
	her->SetXTitle("eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp)");
	her->Draw();

	c3->Divide(2,2);
	c3->cd(1);
	hfsvblkid->Draw("COLZ");
	profX2->SetMarkerStyle(20);
	profX2->SetLineColor(kBlack);
	profX2->Draw("SAME");

	c4->Divide(2,2);
	c4->cd(1);
	hdeltavertex->SetXTitle("vz_bb-vz_sbs (m)");
	hdeltavertex->Draw();
	c4->cd(2);
	hvertexcorr->SetXTitle("vz_bb (m)");
	hvertexcorr->SetYTitle("vz_sbs (m)");
	hvertexcorr->Draw("COLZ");
	c4->cd(3);
	hcointime->Draw();
	c4->cd(4);
	hW2->Draw();



	c->SaveAs(Form("../plots/%s_test_fs.png",printfilename));
	c2->SaveAs(Form("../plots/%s_test_fs_sbstracking.png",printfilename));
	c4->SaveAs(Form("../plots/%s_test_fs_goodtrackscuts.png",printfilename));
	c->Print(Form("../plots/%s_test_fs.pdf(",printfilename));
	c1->Print(Form("../plots/%s_test_fs.pdf",printfilename));
	c2->Print(Form("../plots/%s_test_fs.pdf",printfilename));
	c3->Print(Form("../plots/%s_test_fs.pdf",printfilename));
	c4->Print(Form("../plots/%s_test_fs.pdf)",printfilename));
}