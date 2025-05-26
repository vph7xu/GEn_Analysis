#pragma once

//#include "parse.h"
#include "cuts.h"
#include "plotdxdy.h"
//#include "parse.h"

void sbs_veto(const char* filename,const char* printfilename,const char* kin, bool flag_eHCAL_cut, bool veto){
	
	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	cuts cutsobject;	
	cutsobject.parsecuts(config);
	
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
	double nhits_sbs = 0.0;
	double chi2ndf_sbs = 0.0;
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
	tree->SetBranchAddress("nhits_sbs", &nhits_sbs);
	tree->SetBranchAddress("chi2ndf_sbs", &chi2ndf_sbs);
	tree->SetBranchAddress("nhodo_clus",&nhodo_clus);


	TH1D *h_dx = new TH1D("h_dx","dx",200,-6,3);
	TH1D *h_dy = new TH1D("h_dy","dy",200,-4,4);
	TH1D *h_W2 = new TH1D("h_W2","W2",100,0,4);
	TH1D *h_W2_cut_cointime = new TH1D("h_W2_cut_cointime","W2 with cointime cut",100,-4,8);
	TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",250,-2,2,250,-4,4);
	TH2D *h_dyW2 = new TH2D("h_dyW2","dy v W2",200,-1,3,200,-4,4);
	TH2D *h_dxW2 = new TH2D("h_dxW2","dx v W2",200,-1,3,200,-2,2);
	TH1D *h_coin_time = new TH1D("h_coin_time","coin_time",200,150,250);
	TH1D *h_eHCAL_cut_QE = new TH1D("h_eHCAL_cut_QE","eHCAL with cointime, W2, dx and dy cuts;eHCAL(GeV)", 200, 0, 1.5);
	TH1D *h_eHCAL_cut_cointime_W2 = new TH1D("h_eHCAL_cut_cointime_W2","eHCAL with cointime and W2 cuts;eHCAL(GeV)",200,0,1.5);
	TH1D *h_eHCAL_cut_cointime = new TH1D("h_eHCAL_cut_cointime","eHCAL with cointime cut;eHCAL(GeV)",200,0,1.5);
	TH1D *h_eHCAL = new TH1D("h_eHCAL","eHCAL;eHCAL(GeV)",200,0,1.5);

	TH1D *h_dx_W2_cut = new TH1D("h_dx_W2_cut","dx after W2, cointime and dy cut",200,-6,3);
	TH1D *h_dy_W2_cut = new TH1D("h_dy_W2_cut","dy after W2 and cointime cut",200,-4,4);
	//TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
	TH2D *h_dxdy_W2_cut = new TH2D("h_dxdy_W2_cut","dxdy after W2 and cointime cut",250,-2,2,250,-4,4);
	TH1D *h_coin_time_W2_cut = new TH1D("h_coin_time_W2_cut","cointime after W2 cut",200,150,250);
	TH1D *h_W2_QE_cuts = new TH1D("h_W2_QE_cuts","W2",100,0,4);

	TH1D *h_dx_all = new TH1D("h_dx_all","dx distribution;dx(m)",200,-6,4);
	TH1D *h_dx_sbstracks = new TH1D("h_dx_sbstracks","dx distribution sbs tracks;dx(m)",200,-6,4);
	TH1D *h_dx_nosbstracks = new TH1D("h_dx_nosbstracks","dx distribution no sbs tracks;dx(m)",200,-6,4);

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
		
		if(nhodo_clus==0 or trigbits != 4){
			continue;
		}

		if(cutsobject.run_num_L<runnum and runnum<cutsobject.run_num_H and ePS>0.2){
			if (eHCAL>eHCAL_L and (cutsobject.W2_L<W2 and W2<cutsobject.W2_H) and (cutsobject.coin_time_L<coin_time and coin_time<cutsobject.coin_time_H) ){
				if (dy<cutsobject.dy_H and cutsobject.dy_L<dy){
					
					h_dx_all->Fill(dx); //fill all dx
					
					if(ntrack_sbs>0 and nhits_sbs>2 and chi2ndf_sbs<5){   //fill sbs tracks only
						h_dx_sbstracks->Fill(dx);
					}

					else{  //fill no sbs tracks
						h_dx_nosbstracks->Fill(dx);
					}

			}
		}
	}

		if (veto and ntrack_sbs>0){
			continue;
		}

		if(cutsobject.run_num_L<runnum and runnum<cutsobject.run_num_H and ePS>0.2){
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
    
	c6->Divide(1,1);
	c6->cd(1);
	h_dx_all->SetLineColor(kRed);
	h_dx_sbstracks->SetLineColor(kBlue);
	h_dx_nosbstracks->SetLineColor(kGreen);
	h_dx_all->SetLineWidth(3);
	h_dx_sbstracks->SetLineWidth(3);
	h_dx_nosbstracks->SetLineWidth(3);
	h_dx_all->Draw();
	h_dx_sbstracks->Draw("same");
	h_dx_nosbstracks->Draw("same");

	auto leg = new TLegend(0.65, 0.50, 0.90, 0.70);
	leg->AddEntry(h_dx_all,         "All ",           "l");
	leg->AddEntry(h_dx_sbstracks,   "Only SBS tracks",      "l");
	leg->AddEntry(h_dx_nosbstracks, "Without SBS tracks",   "l");

	leg->SetBorderSize(0);
	leg->SetFillStyle(0);
	leg->SetTextSize(0.03);

	leg->Draw();

	c6->Print(Form("plots/%s_sbs_veto_comparison%s.png",printfilename,std::to_string(flag_eHCAL_cut).c_str()));

   }
