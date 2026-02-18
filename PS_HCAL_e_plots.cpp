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

	double dy_p_L = getDoubleValue(config,"dy_p_L");
	double dy_p_H = getDoubleValue(config,"dy_p_H");

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
	double vy = 0.0;
	double vx = 0.0;
	double ntrack = 0.0;
	double track_nhits = 0.0;
	double track_chi2 = 0.0;
	double nhits = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;
	int runnum = 0.0;
	int helicity = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	double thtgt = 0.0;
	double thetabend = 0.0;
	double ebeam = 0.0;
	double etheta = 0.0;
	double BPMAy = 0.0;
	double BPMAx = 0.0;
	double rastery = 0.0;
	double rasterx = 0.0;
	double raster2y = 0.0;
	double raster2x = 0.0;
	//double eHCAL = 0.0;

	tree->SetBranchAddress("ePS",&ePS);
	tree->SetBranchAddress("eSH",&eSH);
	tree->SetBranchAddress("trP",&trP);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("vy",&vy);
	tree->SetBranchAddress("vx",&vx);
	tree->SetBranchAddress("ntrack",&ntrack);
	tree->SetBranchAddress("ntrack_chi2",&track_chi2);
	tree->SetBranchAddress("ntrack_hits",&track_nhits);
	tree->SetBranchAddress("nhits",&nhits);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("thtgt",&thtgt);
	tree->SetBranchAddress("thetabend",&thetabend);
	tree->SetBranchAddress("ebeam",&ebeam);
	tree->SetBranchAddress("etheta",&etheta);
	tree->SetBranchAddress("BPMAy",&BPMAy);
	tree->SetBranchAddress("BPMAx",&BPMAx);
	tree->SetBranchAddress("Rastery",&rastery);
	tree->SetBranchAddress("Rasterx",&rasterx);
	tree->SetBranchAddress("Raster2y",&raster2y);
	tree->SetBranchAddress("Raster2x",&raster2x);
	//tree->SetBranchAddress("eHCAL",&eHCAL);

	TH1D *h_ePS = new TH1D("h_ePS","Pre-shower energy",200,0,2.5);
	TH1D *h_eHCAL = new TH1D("h_eHCAL","HCal energy",100,0,1);
	TH1D *h_trvz = new TH1D("h_trvz","Vertex z",200,-0.8,0.8);
	TH1D *h_ntrack = new TH1D("h_ntrack","ntrack",60,0,6);
	TH1D *h_nhits = new TH1D("h_nhits","nhits",60,0,6);
	TH1D *h_W2 = new TH1D("h_W2","W^{2}",200,-2,4);
	TH1D *h_W2_0 = new TH1D("h_W2_0","W^{2}; W^{2} (GeV^{2})",50,-1,3);
	TH1D *h_W2_1 = new TH1D("h_W2_1","W^{2}",200,-2,4);
	TH1D *h_W2_2 = new TH1D("h_W2_2","W^{2}",200,-2,4);
	TH1D *h_W2_3 = new TH1D("h_W2_3","W^{2}",200,-2,4);
	TH1D *h_dx = new TH1D("h_dx","dx; dx(m)", 100,-4,3);
	TH1D *h_dy = new TH1D("h_dy","dy; dy(m)", 100,-4,3);

	TH1D *h_dx_0 = new TH1D("h_dx_0","dx",100,-4,3);
	TH1D *h_dx_1 = new TH1D("h_dx_1","dx",100,-4,3);
	TH1D *h_dx_2 = new TH1D("h_dx_2","dx",100,-4,3);
	TH1D *h_dx_3 = new TH1D("h_dx_3","dx",100,-4,3);

	TH1D *h_W2_0_0 = new TH1D("h_W2_0_0","W^{2}",200,-2,4);
	TH1D *h_W2_1_1 = new TH1D("h_W2_1_0","W^{2}",200,-2,4);
	TH1D *h_W2_2_2 = new TH1D("h_W2_2_0","W^{2}",200,-2,4);
	TH1D *h_W2_3_3 = new TH1D("h_W2_3_0","W^{2}",200,-2,4);

	TH1D *h_dx_0_0 = new TH1D("h_dx_0_0","dx",100,-4,3);
	TH1D *h_dx_1_1 = new TH1D("h_dx_1_0","dx",100,-4,3);
	TH1D *h_dx_2_2 = new TH1D("h_dx_2_0","dx",100,-4,3);
	TH1D *h_dx_3_3 = new TH1D("h_dx_3_0","dx",100,-4,3);

	TH1D *h_cointime = new TH1D("h_cointime","Coincidence time",200,coin_time_L-50,coin_time_H+50);
	TH1D *h_eoverp = new TH1D("h_eoverp","E/p BBCal ; E/p", 100,0.3,1.7);
	TH1D *h_chi2_all = new TH1D("h_chi2_all","#chi^{2} ; #chi^{2}",100,0,60);
	TH1D *h_chi2_3hit = new TH1D("h_chi2_3hit","#chi^{2} ; #chi^{2}",100,0,60);
	TH1D *h_chi2_4hit = new TH1D("h_chi2_4hit","#chi^{2} ; #chi^{2}",100,0,60);
	TH1D *h_chi2_5hit = new TH1D("h_chi2_5hit","#chi^{2} ; #chi^{2}",100,0,60);
	//TH1D *h_eHCAL = new TH1D("h_eHCAL","HCAL energy ; HCAL energy (GeV)",200,0,2);

	TH1D *h_trvz_c = new TH1D("h_trvz_c","vertex z",1000,-1,1);

	TH2D *h_dxdy = new TH2D("h_dxdy","dx vs dy ; dy(m) ; dx(m)",100,-4,3,100,-4,3);

	TH1D *h_W2_new = new TH1D("h_W2_new","W^{2} new ; W^{2} (GeV^{2})",50,-1,3);
	TH1D *h_W2_new1 = new TH1D("h_W2_new1","W^{2} new1 ; W^{2} (GeV^{2})",50,-1,3);

	TH1D *h_W2_new_1 = new TH1D("h_W2_new_1","W^{2} new ; W^{2} (GeV^{2})",50,-1,3);
	TH1D *h_W2_new_1_1 = new TH1D("h_W2_new_1_1","W^{2} new1 ; W^{2} (GeV^{2})",50,-1,3);

	TH1D *h_trP = new TH1D("h_trP","trP",100,0,7);

	TH1D *h_trP_new = new TH1D("h_trP_new","trP new",100,0,7);
	TH1D *h_trP_new1 = new TH1D("h_trP_new1","trP new1",100,0,7);

	TH2D *h_W2_v_thetatgt = new TH2D("h_W2_v_thetatgt","W2 vs thetatgt",100,-0.3,0.3,100,-2,5);
	TH2D *h_W2timesthetabend_v_thetatgt = new TH2D("h_W2timesthetabend_v_thetatgt","W2*thetabend vs thetatgt",100,-0.3,0.3,100,-0.1,0.6);

	TH2D *h_p_eltimesthetabend_v_thetatgt = new TH2D("h_p_eltimesthetabend_v_thetatgt","p_el*thetabend vs thetatgt;thetatgt(rad);p_el(GeV)*thetabend(rad)", 100,-0.2,0.1,100,0.1,0.5);

	TH2D *h_trPtimesthetabend_v_thetatgt = new TH2D("h_trPtimesthetabend_v_thetatgt","trP*thetabend vs thetatgt;thetatgt(rad);trP(GeV)*thetabend(rad)", 100,-0.3,0.3,100,-0.1,0.6);

	TH2D *h_p_el_v_BPMAy = new TH2D("h_p_el_v_BPMAy","#Delta p vs BPMAy;BPMAy(m);p/p_el - 1",50,-0.003,0.003,50,-0.1,0.3);
	TH2D *h_p_el_v_BPMAx = new TH2D("h_p_el_v_BPMAx","#Delta p vs BPMAx;BPMAx(m);p/p_el - 1",50,-0.003,0.003,50,-0.1,0.3);

	TH2D *h_p_el_v_vy = new TH2D("h_p_el_v_vy","#Delta p vs vy;vy(m);p/p_el - 1",50,-0.003,0.003,50,-0.1,0.1);
	TH2D *h_p_el_v_vx = new TH2D("h_p_el_v_vx","#Delta p vs vx;vx(m);p/p_el - 1",50,-0.003,0.003,50,-0.1,0.1);

	TH2D *h_p_el_v_rastery = new TH2D("h_p_el_v_rastery","#Delta p vs rastery;rastery;p/p_el - 1",50,34000,54000,50,-0.1,0.3);
	TH2D *h_p_el_v_rasterx = new TH2D("h_p_el_v_rasterx","#Delta p vs rasterx;rasterx;p/p_el - 1",50,34000,54000,50,-0.1,0.3);

	TH2D *h_p_el_v_raster2y = new TH2D("h_p_el_v_raster2y","#Delta p vs raster2y;raster2y;p/p_el - 1",50,34000,54000,50,-0.1,0.3);
	TH2D *h_p_el_v_raster2x = new TH2D("h_p_el_v_raster2x","#Delta p vs raster2x;raster2x;p/p_el - 1",50,34000,54000,50,-0.1,0.3);

	TH1D *h_dp = new TH1D("h_dp","#Delta p (Sean's);p/p_el - 1",100,-0.25,0.25);
	TH1D *h_dp_new = new TH1D("h_dp_new","#Delta p (Sean's vy corr);p/p_el - 1",100,-0.25,0.25);
	TH1D *h_dp_new1 = new TH1D("h_dp_new1","#Delta p (New vy corr);p/p_el - 1",100,-0.25,0.25);

	TH1D *h_ebeam = new TH1D("h_ebeam","beam energy; beam energy(GeV)",100,8.4,8.5);

	// Sensitivity study histograms
	TH2D *h_dW2_vs_dpf = new TH2D("h_dW2_vs_dpf", "#Delta W^{2} vs #Delta p;#Delta p;#Delta W^{2} (GeV^{2})",
	                               60, -0.06, 0.06, 120, -1.0, 1.0);
	TProfile *p_dW2_vs_dpf = new TProfile("p_dW2_vs_dpf", "#Delta W^{2} vs #Delta p;#Delta p;#LT #Delta W^{2} #GT",
	                                      60, -0.06, 0.06);
	TH1D *h_S_analytic = new TH1D("h_S_analytic", "Analytic sensitivity S;S (GeV^{2});Counts", 120, -20.0, -10.0);
	TH1D *h_S_empirical = new TH1D("h_S_empirical", "Empirical S = #Delta W^{2}/#Delta p;S (GeV^{2});Counts", 120, -20.0, -10.0);


	double A1 = 0.26056;//0.28640;//0.28661;//0.287186374;
	double B1 = 0.96542;//0.87122905;//0.839712501;//0.798543608;

	double A = 0.26056;//0.263077985;//0.287186374;
	double B = 0.96542;//0.999237626;//0.798543608;

	double m = 0.938;

	double D = -5.96546*2.65;//-4.43313*2.55;//-15.893095*2.4;//3.19086*2.4;//-0.00184*3.1;
	double C = -0.01193*2.65;//-0.06490*2.55;//-0.0147632089*2.4;//0.09312*2.4;//-5.17872*3.1;


	double x0 = 0.5*(dx_p_L + dx_p_H);  // center in x (dx)
	double y0 = 0.5*(dy_p_L + dy_p_H);                           // center in y (dy)
	double rx = 0.22;                           // semi-axis in x
	double ry = 0.18;                           // semi-axis in y

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);

		double eoverp = (eSH+ePS)/trP;

		double ee = A*(1+B*thtgt)/thetabend - (vy*D + C);
		double ee1 = A1*(1+B1*thtgt)/thetabend - (vy*D + C);

		double ee_vycorr = A*(1+B*thtgt)/thetabend - (vy*D + C);
		double ee1_vycorr = A1*(1+B1*thtgt)/thetabend - (vy*D + C);

		double Q2_new = 2*ebeam*ee_vycorr*(1-cos(etheta));
		double Q2_new1 = 2*ebeam*ee1_vycorr*(1-cos(etheta));

		double nu_new = (ebeam - ee_vycorr);
		double nu_new1 = (ebeam - ee1_vycorr);

		double W2_new = m*m+2*m*nu_new-Q2_new;
		double W2_new1 = m*m+2*m*nu_new1-Q2_new1;

		double p_el = ebeam/(1+ebeam*(1-cos(etheta))/m);

		double delta_p = trP/p_el -1 ;

		// Build W2 from a given momentum (keeping ebeam, etheta fixed)
		auto W2_from_p = [&](double p)->double{
		  // electron final energy ~ p (ultrarelativistic); M = nucleon mass
		  double Q2 = 2.0*ebeam*p*(1.0 - cos(etheta));
		  double nu = ebeam - p;
		  return m*m + 2.0*m*nu - Q2;
		};

		// Reference (at true elastic momentum) and reconstructed
		double W2_el   = W2_from_p(p_el);     // using ideal elastic p_el
		double W2_reco = W2_from_p(ee1);      // using reconstructed momentum

		double dW2 = W2_reco - W2_el;         // shift due to momentum error
		double dpfrac = ee1/p_el - 1;//delta_p;              // your Δp definition

		// Fill 2D and profile for linear fit
		if (std::abs(dpfrac) < 0.05) {        // keep to small fractional errors
		  h_dW2_vs_dpf->Fill(dpfrac, dW2);
		  p_dW2_vs_dpf->Fill(dpfrac, dW2);
		  if (std::abs(dpfrac) > 1e-6) h_S_empirical->Fill(dW2/dpfrac);
		}

		// Analytic sensitivity S = (dW2/dp)*p_el  with dW2/dp = -2M - 4 Ei sin^2(theta/2)
		double sin2 = sin(0.5*etheta); sin2 *= sin2;
		double dW2_dp = -2.0*m - 4.0*ebeam*sin2;
		double S_analytic = dW2_dp * p_el;
		h_S_analytic->Fill(S_analytic);

	    bool goodHelicity = (lookupValue(HelicityCheck, runnum) == 1);
	    bool goodMoller   = (lookupValue(MollerQuality, runnum) == 1);
	    bool goodVz       = abs(vz) < 0.27;
	    bool goodPS       = (ePS > 0.2);
	    bool goodRunRange = (run_num_L < runnum && runnum < run_num_H);
	    bool goodEHCAL    = (eHCAL > 0.025);//eHCAL_L); 
	    bool validHel     = (helicity == -1 || helicity == 1);
	    //bool goodGrinch = (grinch_track == 0) && (grinch_clus_size>2);
	    //bool goodSbs_track = ntrack_sbs>0 && abs(vz_sbs)<0.27;
	    bool goodEoverp = abs(eoverp-1)<0.2;
	    bool goodW2 = W2_L<W2 && W2<W2_H;
	    bool goodCoin = coin_time_L<coin_time && coin_time<coin_time_H;
	    bool goodDy = -0.2<dy && dy<0.2;
	    bool goodDx = ((pow((dy-y0)/ry,2)+pow((dx-x0)/rx,2))<=1);
	    bool goodpSpot = ((pow((dy-y0)/ry,2)+pow((dx-x0)/rx,2))<=1);
	    /*(-0.2<dx && dx<0.2) ||*/ /*(-0.2+(dx_p_L+dx_p_H)/2<dx && dx<0.2+(dx_p_L+dx_p_H)/2);*/
	    bool goodSpot = dy_L<dy && dy<dy_H && ((dx_L<dx && dx<dx_H) || (dx_p_L<dx && dx<dx_p_H));


		if(ePS>0 && eHCAL>0){
			//if (goodEHCAL && goodPS && goodEoverp && goodCoin && goodW2 && goodSpot ){
				h_trvz->Fill(vz);
			//}

			if(goodVz && goodEHCAL && goodCoin){
				h_ePS->Fill(ePS);
			}
			
			if(goodVz && goodCoin && goodPS){
				h_eHCAL->Fill(eHCAL);
			}
			
			h_ntrack->Fill(ntrack);
			h_nhits->Fill(nhits);

			if(goodVz && goodPS  && goodEoverp && goodEHCAL && goodCoin /*&& goodSpot*/){
				h_W2->Fill(W2);
			}

			if(goodVz && goodPS  && goodEoverp && goodEHCAL && goodW2 &&goodCoin){
				h_dy->Fill(dy);
				h_dxdy->Fill(dy,dx);
				if(goodDy) h_dx->Fill(dx);
			}

			if(goodVz && goodPS  && goodEoverp && goodEHCAL /*&& goodDy*/ &&goodCoin){

				if(goodpSpot){
					h_W2_0->Fill(W2);
					h_W2_new->Fill(W2_new);
					h_W2_new1->Fill(W2_new1);

					h_trP->Fill(trP);
					h_trP_new->Fill(ee);
					h_trP_new1->Fill(ee1);

					h_W2_v_thetatgt->Fill(thtgt,W2);
					h_W2timesthetabend_v_thetatgt->Fill(thtgt,W2*thetabend);

					h_p_eltimesthetabend_v_thetatgt->Fill(thtgt,p_el*thetabend);
					h_trPtimesthetabend_v_thetatgt->Fill(thtgt,trP*thetabend);

					h_p_el_v_BPMAy->Fill(BPMAy,ee/p_el - 1);
					h_p_el_v_BPMAx->Fill(BPMAx,ee/p_el - 1);

					h_p_el_v_vy->Fill(vy,ee/p_el - 1);
					h_p_el_v_vx->Fill(vx,ee/p_el - 1);

					h_dp->Fill(trP/p_el - 1);
					h_dp_new->Fill(ee/p_el - 1);
					h_dp_new1->Fill(ee1/p_el - 1);

					h_p_el_v_rastery->Fill(rastery,ee/p_el - 1);
					h_p_el_v_rasterx->Fill(rasterx,ee/p_el - 1);

					h_p_el_v_raster2y->Fill(raster2y,ee/p_el - 1);
					h_p_el_v_raster2x->Fill(raster2x,ee/p_el - 1);

					h_ebeam->Fill(ebeam);

				}
				//cout<<"W2_old : "<<W2<<"   W2 new : "<<W2_new<<endl;

				if(goodW2) h_dx_0->Fill(dx);


				if(coin_time_L<coin_time && coin_time<coin_time_H-1){
					h_W2_1->Fill(W2);

					if(goodW2) h_dx_1->Fill(dx);	
				}

				if(coin_time_L<coin_time && coin_time<coin_time_H-2){
					h_W2_2->Fill(W2);
					if(goodW2) h_dx_2->Fill(dx);	
				}

				if(coin_time_L<coin_time && coin_time<coin_time_H-3){
					h_W2_3->Fill(W2);
					if(goodW2) h_dx_3->Fill(dx);	
				}

			}

			if(goodVz && goodPS  && goodEoverp && goodDy &&goodCoin && eHCAL>eHCAL_L-0.3){

				h_W2_0_0->Fill(W2);

				if(goodW2) h_dx_0_0->Fill(dx);

				if(eHCAL>eHCAL_L-0.2){
					h_W2_1_1->Fill(W2);
					if(goodW2) h_dx_1_1->Fill(dx);	
				}

				if(eHCAL>eHCAL_L-0.1){
					h_W2_2_2->Fill(W2);
					if(goodW2) h_dx_2_2->Fill(dx);	
				}

				if(eHCAL>eHCAL_L-0.0){
					h_W2_3_3->Fill(W2);
					if(goodW2) h_dx_3_3->Fill(dx);	
				}

			}

			if(goodVz && goodPS && goodEHCAL && goodEoverp && goodW2 && goodSpot){
				h_cointime->Fill(coin_time);
			}

			if(goodVz && goodPS && goodEHCAL && goodW2 && goodCoin /*&& goodSpot*/){
				h_eoverp->Fill(eoverp);
			}

			if(goodVz && goodPS && goodEoverp && goodEHCAL && goodW2 && goodCoin && goodSpot){
				
				h_chi2_all->Fill(track_chi2);

				if(track_nhits == 3){
					h_chi2_3hit->Fill(track_chi2);
				}
				else if(track_nhits == 4){
					h_chi2_4hit->Fill(track_chi2);
				}
				else if(track_nhits == 5){
					h_chi2_5hit->Fill(track_chi2);
				}
			}
		}
		//if (ePS>0.2 and eHCAL>0){
		//	h_trvz_c->Fill(vz);
		//	h_eHCAL->Fill(eHCAL);
		//}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();

	}

	//gStyle->SetOptStat(0);

	TCanvas* c = new TCanvas("c","c",1600,1200);

	c->Divide(2,2);
	c->cd(1);
	h_ePS->SetXTitle("Preshower energy(GeV/c)");
	h_ePS->SetFillStyle(1001);
	h_ePS->SetFillColorAlpha(kViolet+2,0.3);
	h_ePS->SetLineColor(kBlack);
	h_ePS->SetLineWidth(2);
	h_ePS->Draw();
	TLine *line00 = new TLine(0.2, 0, 0.2, h_ePS->GetMaximum());
    line00->SetLineColor(kRed);
    line00->SetLineWidth(2);
    line00->Draw("same");

	c->cd(2);
	h_cointime->SetXTitle("Coincidence time(ns)");
	h_cointime->SetFillStyle(1001);
	h_cointime->SetFillColorAlpha(kViolet+2,0.3);
	h_cointime->SetLineColor(kBlack);
	h_cointime->SetLineWidth(2);
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
	h_trvz->SetXTitle("Vertex z(m)");
	h_trvz->SetFillStyle(1001);
	h_trvz->SetFillColorAlpha(kViolet+2,0.3);
	h_trvz->SetLineColor(kBlack);
	h_trvz->SetLineWidth(2);
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
	h_W2->SetXTitle("W^{2}((GeV/c)^{2})");
	h_W2->SetFillStyle(1001);
	h_W2_1->SetFillStyle(1001);
	h_W2_2->SetFillStyle(1001);
	h_W2_3->SetFillStyle(1001);

	h_W2->SetFillColorAlpha(TColor::GetColor("#f0f0f0"),0.3);
	h_W2->SetFillColorAlpha(TColor::GetColor("#bdbdbd"),0.3);
	h_W2->SetFillColorAlpha(TColor::GetColor("#969696"),0.3);
	h_W2->SetFillColorAlpha(TColor::GetColor("#636363"),0.3);

	h_W2->SetLineColor(kBlack);
	h_W2_1->SetLineColor(kBlack);
	h_W2_2->SetLineColor(kBlack);
	h_W2_3->SetLineColor(kBlack);

	h_W2->SetLineWidth(2);
	h_W2_1->SetLineWidth(2);
	h_W2_2->SetLineWidth(2);
	h_W2_3->SetLineWidth(2);

	h_W2->Draw();
	h_W2_1->Draw("same");
	h_W2_2->Draw("same");
	h_W2_3->Draw("same");

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
	//h_eHCAL->SetLineWidth(2);
	h_eHCAL->SetXTitle("HCal energy (GeV/c)");
	h_eHCAL->SetFillStyle(1001);
	h_eHCAL->SetFillColorAlpha(kViolet+2,0.3);
	h_eHCAL->SetLineColor(kBlack);
	h_eHCAL->SetLineWidth(2);
	h_eHCAL->Draw();
	line09->Draw("same");
	//c1->Update();
	c1->cd(2);
	//h_eoverp->SetLineWidth(2);
	h_eoverp->SetFillStyle(1001);
	h_eoverp->SetFillColorAlpha(kViolet+2,0.3);
	h_eoverp->SetLineColor(kBlack);
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

    c1->cd(3);
    h_chi2_all->SetFillStyle(1001);
    h_chi2_3hit->SetFillStyle(1001);
    h_chi2_4hit->SetFillStyle(1001);
    h_chi2_5hit->SetFillStyle(1001);

    h_chi2_all->SetLineColor(kBlack);
    h_chi2_3hit->SetLineColor(kBlack);
	h_chi2_4hit->SetLineColor(kBlack);
	h_chi2_5hit->SetLineColor(kBlack);

	h_chi2_all->SetLineWidth(2);
	h_chi2_3hit->SetLineWidth(2);
	h_chi2_4hit->SetLineWidth(2);
	h_chi2_5hit->SetLineWidth(2);

	h_chi2_all->SetFillColorAlpha(kGray,0.3);
	h_chi2_3hit->SetFillColorAlpha(kGray+1,0.3);
	h_chi2_4hit->SetFillColorAlpha(kGray+2,0.3);
	h_chi2_5hit->SetFillColorAlpha(kGray+3,0.3);

    h_chi2_all->Draw();
    h_chi2_3hit->Draw("same");
    h_chi2_4hit->Draw("same");
    h_chi2_5hit->Draw("same");

	auto leg = new TLegend(0.62, 0.60, 0.88, 0.88); // x1,y1,x2,y2 — tweak as needed
	leg->SetBorderSize(0);
	leg->SetFillStyle(0);           // transparent legend box
	leg->SetTextFont(42);
	leg->SetTextSize(0.035);
	//leg->SetNColumns(2);            // optional: 2 columns for 4 entries

	leg->AddEntry(h_chi2_all,  "3,4,5-hit tracks",    "f");
	leg->AddEntry(h_chi2_3hit, "3-hit tracks","f");
	leg->AddEntry(h_chi2_4hit, "4-hit tracks","f");
	leg->AddEntry(h_chi2_5hit, "5-hit tracks","f");

	leg->Draw();
	gPad->Update();


	TCanvas* c2 = new TCanvas("c2","c2",1600,1200);

	h_W2_0->SetLineColor(kBlack);
	h_W2_1->SetLineColor(kBlack);
	h_W2_2->SetLineColor(kBlack);
	h_W2_3->SetLineColor(kBlack);

	h_dx_0->SetLineColor(kBlack);
	h_dx_1->SetLineColor(kBlack);
	h_dx_2->SetLineColor(kBlack);
	h_dx_3->SetLineColor(kBlack);

	h_W2_0->SetLineWidth(2);
	h_W2_1->SetLineWidth(2);
	h_W2_2->SetLineWidth(2);
	h_W2_3->SetLineWidth(2);

	h_dx_0->SetLineWidth(2);
	h_dx_1->SetLineWidth(2);
	h_dx_2->SetLineWidth(2);
	h_dx_3->SetLineWidth(2);

	//h_W2_0->SetFillStyle(1001);
	h_W2_1->SetFillStyle(1001);
	h_W2_2->SetFillStyle(1001);
	h_W2_3->SetFillStyle(1001);

	h_dx_0->SetFillStyle(1001);
	h_dx_1->SetFillStyle(1001);
	h_dx_2->SetFillStyle(1001);
	h_dx_3->SetFillStyle(1001);

	//h_W2_0->SetFillColorAlpha(kGray,0.3);
	h_W2_1->SetFillColorAlpha(kGray+1,0.3);
	h_W2_2->SetFillColorAlpha(kGray+2,0.3);
	h_W2_3->SetFillColorAlpha(kGray+3,0.3);

	h_dx_0->SetFillColorAlpha(kGray,0.3);
	h_dx_1->SetFillColorAlpha(kGray+1,0.3);
	h_dx_2->SetFillColorAlpha(kGray+2,0.3);
	h_dx_3->SetFillColorAlpha(kGray+3,0.3);


	h_W2_0_0->SetLineColor(kBlack);
	h_W2_1_1->SetLineColor(kBlack);
	h_W2_2_2->SetLineColor(kBlack);
	h_W2_3_3->SetLineColor(kBlack);

	h_dx_0_0->SetLineColor(kBlack);
	h_dx_1_1->SetLineColor(kBlack);
	h_dx_2_2->SetLineColor(kBlack);
	h_dx_3_3->SetLineColor(kBlack);

	h_W2_0_0->SetLineWidth(2);
	h_W2_1_1->SetLineWidth(2);
	h_W2_2_2->SetLineWidth(2);
	h_W2_3_3->SetLineWidth(2);

	h_dx_0_0->SetLineWidth(2);
	h_dx_1_1->SetLineWidth(2);
	h_dx_2_2->SetLineWidth(2);
	h_dx_3_3->SetLineWidth(2);

	h_W2_0_0->SetFillStyle(1001);
	h_W2_1_1->SetFillStyle(1001);
	h_W2_2_2->SetFillStyle(1001);
	h_W2_3_3->SetFillStyle(1001);

	h_dx_0_0->SetFillStyle(1001);
	h_dx_1_1->SetFillStyle(1001);
	h_dx_2_2->SetFillStyle(1001);
	h_dx_3_3->SetFillStyle(1001);

	h_W2_0_0->SetFillColorAlpha(kGray,0.3);
	h_W2_1_1->SetFillColorAlpha(kGray+1,0.3);
	h_W2_2_2->SetFillColorAlpha(kGray+2,0.3);
	h_W2_3_3->SetFillColorAlpha(kGray+3,0.3);

	h_dx_0_0->SetFillColorAlpha(kGray,0.3);
	h_dx_1_1->SetFillColorAlpha(kGray+1,0.3);
	h_dx_2_2->SetFillColorAlpha(kGray+2,0.3);
	h_dx_3_3->SetFillColorAlpha(kGray+3,0.3);

	c2->Divide(2,2);

	c2->cd(1);
	h_dx_0->Draw();
	h_dx_1->Draw("same");
	h_dx_2->Draw("same");
	h_dx_3->Draw("same");

	c2->cd(2);
	h_W2_0->Draw();
	h_W2_1->Draw("same");
	h_W2_2->Draw("same");
	h_W2_3->Draw("same");

	c2->cd(3);
	h_dx_0_0->Draw();
	h_dx_1_1->Draw("same");
	h_dx_2_2->Draw("same");
	h_dx_3_3->Draw("same");

	c2->cd(4);
	h_W2_0_0->Draw();
	h_W2_1_1->Draw("same");
	h_W2_2_2->Draw("same");
	h_W2_3_3->Draw("same");

	TLine *linedx1 = new TLine(dx_L, h_dx->GetMinimum(), dx_L, h_dx->GetMaximum());
    TLine *linedx2 = new TLine(dx_H, h_dx->GetMinimum(), dx_H, h_dx->GetMaximum());

	TLine *linedy1 = new TLine(dy_L, h_dy->GetMinimum(), dy_L, h_dy->GetMaximum());
    TLine *linedy2 = new TLine(dy_H, h_dy->GetMinimum(), dy_H, h_dy->GetMaximum());

	TCanvas* c3 = new TCanvas("c3","c3",1600,1200);
	c3->Divide(2,2);
	c3->cd(1);
	h_dx->SetFillStyle(1001);
	h_dx->SetFillColorAlpha(kViolet+2,0.3);
	h_dx->SetLineColor(kBlack);
	h_dx->SetLineWidth(3);
	linedx1->SetLineWidth(3);
	linedx2->SetLineWidth(3);
	linedx1->SetLineColor(kRed);
	linedx2->SetLineColor(kRed);
	h_dx->Draw();
	linedx1->Draw("same");
	linedx2->Draw("same");

	c3->cd(2);
	h_dy->SetFillStyle(1001);
	h_dy->SetFillColorAlpha(kViolet+2,0.3);
	h_dy->SetLineColor(kBlack);
	h_dy->SetLineWidth(3);
	linedy1->SetLineWidth(3);
	linedy2->SetLineWidth(3);
	linedy1->SetLineColor(kRed);
	linedy2->SetLineColor(kRed);
	h_dy->Draw();
	linedy1->Draw("same");
	linedy2->Draw("same");

	c3->cd(3);

	gStyle->SetPalette(kRainBow);
	h_dxdy->Draw("COLZ");

	auto el = new TEllipse(y0, x0, ry, rx);    // last args: start/end/theta use defaults
	el->SetFillStyle(0);                       // hollow
	el->SetLineColor(kRed);
	el->SetLineWidth(3);
	el->Draw("same");                          // overlay on current pad
	gPad->Update();

	TCanvas *c4 =  new TCanvas("c4","c4",1600,1200);

	h_W2_new->SetLineWidth(3);
	h_W2_new->SetLineColor(kBlue);
	h_W2_new1->SetLineWidth(3);
	h_W2_new1->SetLineColor(kRed);

	//h_W2_new->SetFillStyle(1001);
	//h_W2_new->SetFillColorAlpha(kGray,0.3);

	c4->Divide(2,2);
	c4->cd(1);
	//h_W2_new1->Draw();
	h_W2_0->Draw();
	h_W2_0->SetTitle("W^{2} (Sean's)");
	//h_W2_new->Draw("same");

	auto legW2 = new TLegend(0.62, 0.60, 0.88, 0.88); // x1,y1,x2,y2 — tweak as needed
	legW2->SetBorderSize(0);
	legW2->SetFillStyle(0);           // transparent legend box
	legW2->SetTextFont(42);
	legW2->SetTextSize(0.035);
	//leg->SetNColumns(2);            // optional: 2 columns for 4 entries

	legW2->AddEntry(h_W2_0,  "Sean's",    "l");
	legW2->AddEntry(h_W2_new, "Sean's vy corr","l");
	legW2->AddEntry(h_W2_new1, "New vy corr","l");


	//legW2->Draw();
	gPad->Update();

	c4->cd(2);
	h_W2_new1->SetTitle("W^{2} (new)");
	h_W2_new1->Draw();
	

	c4->cd(3);
	h_trP->SetLineWidth(3);
	h_trP_new->SetLineWidth(3);
	h_trP_new1->SetLineWidth(3);
	h_trP->SetLineColor(kBlack);
	h_trP_new->SetLineColor(kBlue);
	h_trP_new1->SetLineColor(kRed);
	h_trP->Draw();
	h_trP_new->Draw();
	h_trP_new1->Draw();

	c4->cd(4);
	h_trP_new->Draw();

	TCanvas *c5 =  new TCanvas("c5","c5",1600,1200);
	c5->Divide(2,2);
	c5->cd(1);
	h_W2_v_thetatgt->Draw("COLZ");
	c5->cd(2);
	h_W2timesthetabend_v_thetatgt->Draw("COLZ");
	c5->cd(3);
	h_p_eltimesthetabend_v_thetatgt->Draw("COLZ");
	// parameters you can tweak
	double q = 0.25;      // 20th percentile traces the lower edge
	int    nmin = 5;     // minimum entries per x-bin to keep a point
	int    rebinx = 1;    // optional x rebin to stabilize statistics
	if (rebinx>1) h_p_eltimesthetabend_v_thetatgt->RebinX(rebinx);

	// build lower-edge points (one per x-bin)
	std::vector<double> vectorx, vectory;
	for (int ix = 1; ix <= h_p_eltimesthetabend_v_thetatgt->GetNbinsX(); ++ix) {
	  std::unique_ptr<TH1D> py(h_p_eltimesthetabend_v_thetatgt->ProjectionY(Form("py_%d",ix), ix, ix, "e"));
	  if (!py || py->GetEntries() < nmin) continue;

	  // compute the q-quantile of y in this x-bin
	  double prob[1] = { q };
	  double qval[1] = { 0.0 };
	  py->GetQuantiles(1, qval, prob);

	  vectorx.push_back(h_p_eltimesthetabend_v_thetatgt->GetXaxis()->GetBinCenter(ix));
	  vectory.push_back(qval[0]);
	}

	// make a graph of the lower edge and fit y = c + m x
	TGraph g(vectorx.size(), vectorx.data(), vectory.data());
	g.SetMarkerStyle(20);
	g.SetMarkerSize(0.9);
	g.SetMarkerColor(kBlack);
	g.Draw("P same");

	TF1 fline("fline","pol1",-0.14,-0.005);
	g.Fit(&fline, "Q","",-0.14,-0.005);         // quiet fit

	fline.SetLineColor(kRed+1);
	fline.SetLineWidth(3);
	fline.Draw("same");

	// ... after doing g.Fit(&fline,"Q"); and drawing the line
	double cc   = fline.GetParameter(0);
	double dc  = fline.GetParError(0);
	double mm   = fline.GetParameter(1);
	double dm  = fline.GetParError(1);
	double chi2 = fline.GetChisquare();
	int    ndf  = fline.GetNDF();

	// Add a nice box on the plot (NDC coords)
	TPaveText *box = new TPaveText(0.14, 0.80, 0.44, 0.93, "NDC");
	box->SetFillStyle(0);            // transparent background (or use SetFillColorAlpha(kWhite,0.8))
	box->SetLineColor(kBlack);
	box->SetTextFont(42);
	box->SetTextSize(0.035);
	box->SetTextAlign(13);           // left/top
	box->AddText(Form("m = %.5f #pm %.5f", mm, dm));
	box->AddText(Form("c = %.5f #pm %.5f", cc, dc));
	//box->AddText(Form("#chi^{2}/ndf = %.1f / %d", chi2, ndf));
	box->Draw();

	gPad->Modified(); gPad->Update();


	c5->cd(4);
	h_trPtimesthetabend_v_thetatgt->Draw("COLZ");

	TCanvas *c6 =  new TCanvas("c6","c6",1600,1200);
	c6->Divide(2,2);
	c6->cd(1);
	h_p_el_v_BPMAy->Draw("COLZ");
	c6->cd(2);
	h_p_el_v_BPMAx->Draw("COLZ");
	c6->cd(3);
	h_p_el_v_vy->Draw("COLZ");

	// parameters you can tweak
	int    nmin_vy   = 5;     // minimum entries per x-bin to keep a point
	int    rebinx_vy = 1;     // optional x rebin to stabilize statistics
	if (rebinx_vy > 1) h_p_el_v_vy->RebinX(rebinx_vy);

	// build "middle-50%" points (per x-bin)
	std::vector<double> vectorx_vy, vectory_vy, erry_vy;
	vectorx_vy.reserve(h_p_el_v_vy->GetNbinsX());
	vectory_vy.reserve(h_p_el_v_vy->GetNbinsX());
	erry_vy.reserve(h_p_el_v_vy->GetNbinsX());

	for (int ix = 1; ix <= h_p_el_v_vy->GetNbinsX(); ++ix) {
	  std::unique_ptr<TH1D> py_vy(h_p_el_v_vy->ProjectionY(Form("py_vy_%d", ix), ix, ix, "e"));
	  if (!py_vy || py_vy->GetEntries() < nmin_vy) continue;

	  // Compute Q1 (25%), Median (50%), Q3 (75%)
	  double probs[3] = {0.25, 0.50, 0.75};
	  double qvals[3] = {0.0,  0.0,  0.0};
	  py_vy->GetQuantiles(3, qvals, probs);

	  double x  = h_p_el_v_vy->GetXaxis()->GetBinCenter(ix);
	  double y  = qvals[1];                        // use the median as the central y
	  double IQR = qvals[2] - qvals[0];            // Q3 - Q1
	  double ey = (IQR > 0) ? (0.5 * IQR) : 0.0;   // robust scale for weighting (semi-IQR)

	  vectorx_vy.push_back(x);
	  vectory_vy.push_back(y);
	  erry_vy.push_back(ey);
	}

	// graph with robust errors -> weighted fit prefers tighter (less spread) bins
	TGraphErrors g_vy(vectorx_vy.size(), vectorx_vy.data(), vectory_vy.data(),
	                  /*ex*/ nullptr, erry_vy.data());
	g_vy.SetMarkerStyle(20);
	g_vy.SetMarkerSize(0.9);
	g_vy.SetMarkerColor(kBlack);
	g_vy.Draw("P same");

	// choose a clean fit range from available x points (or set your own)
	double xmin_fit = *std::min_element(vectorx_vy.begin(), vectorx_vy.end());
	double xmax_fit = *std::max_element(vectorx_vy.begin(), vectorx_vy.end());

	TF1 fline_vy("fline_vy", "pol1", xmin_fit, xmax_fit);

	// "W" uses the provided Y errors for chi2 weighting; "R" respects range; "Q" quiet
	g_vy.Fit(&fline_vy, "Q R W");

	fline_vy.SetLineColor(kRed+1);
	fline_vy.SetLineWidth(3);
	fline_vy.Draw("same");

	// annotate results
	double cc_vy   = fline_vy.GetParameter(0);
	double dc_vy   = fline_vy.GetParError(0);
	double mm_vy   = fline_vy.GetParameter(1);
	double dm_vy   = fline_vy.GetParError(1);
	double chi2_vy = fline_vy.GetChisquare();
	int    ndf_vy  = fline_vy.GetNDF();

	TPaveText *box_vy = new TPaveText(0.14, 0.80, 0.44, 0.93, "NDC");
	box_vy->SetFillStyle(0);
	box_vy->SetLineColor(kBlack);
	box_vy->SetTextFont(42);
	box_vy->SetTextSize(0.035);
	box_vy->SetTextAlign(13);
	box_vy->AddText(Form("m = %.5f #pm %.5f", mm_vy, dm_vy));
	box_vy->AddText(Form("c = %.5f #pm %.5f", cc_vy, dc_vy));
	// box_vy->AddText(Form("#chi^{2}/ndf = %.1f / %d", chi2_vy, ndf_vy));
	box_vy->Draw();

	gPad->Modified(); gPad->Update();



	c6->cd(4);
	h_p_el_v_vx->Draw("COLZ");


	// === Gaussian fits for Δp histograms (range: -0.3 to 0.3) ===================
	double fitMin = -0.03, fitMax = 0.03;

	// Create separate functions so ROOT keeps all three
	TF1 *g_dp      = new TF1("g_dp",      "gaus", fitMin, fitMax);
	TF1 *g_dp_new  = new TF1("g_dp_new",  "gaus", fitMin, fitMax);
	TF1 *g_dp_new1 = new TF1("g_dp_new1", "gaus", fitMin, fitMax);

	// Seed with sensible initial parameters
	g_dp     ->SetParameters(h_dp     ->GetMaximum(), h_dp     ->GetMean(), h_dp     ->GetRMS());
	g_dp_new ->SetParameters(h_dp_new ->GetMaximum(), h_dp_new ->GetMean(), h_dp_new ->GetRMS());
	g_dp_new1->SetParameters(h_dp_new1->GetMaximum(), h_dp_new1->GetMean(), h_dp_new1->GetRMS());

	// Quiet fits, respect range
	h_dp     ->Fit(g_dp,      "RQ");
	h_dp_new ->Fit(g_dp_new,  "RQ");
	h_dp_new1->Fit(g_dp_new1, "RQ");

	// Match fit line colors to the hist styling used later
	g_dp     ->SetLineColor(kBlack);
	g_dp_new ->SetLineColor(kBlue);
	g_dp_new1->SetLineColor(kRed);
	g_dp     ->SetLineWidth(3);
	g_dp_new ->SetLineWidth(3);
	g_dp_new1->SetLineWidth(3);

	// Helper to add a mean/sigma box to the current pad
	auto add_mu_sigma_box = [](TF1* f, double x1=0.58, double y1=0.70, double x2=0.88, double y2=0.88){
	  double mu   = f->GetParameter(1);
	  double dmu  = f->GetParError(1);
	  double sig  = f->GetParameter(2);
	  double dsig = f->GetParError(2);
	  TPaveText *box = new TPaveText(x1,y1,x2,y2,"NDC");
	  box->SetFillStyle(0);
	  box->SetLineColor(kBlack);
	  box->SetTextFont(42);
	  box->SetTextSize(0.035);
	  box->SetTextAlign(13);
	  box->AddText(Form("#mu = %.5f #pm %.5f",   mu,  dmu));
	  box->AddText(Form("#sigma = %.5f #pm %.5f", sig, dsig));
	  box->Draw();
	  gPad->Modified(); gPad->Update();
	};


	TCanvas *c7 =  new TCanvas("c7","c7",1600,1200);
	c7->Divide(2,2);

	// Panel 1: h_dp (Sean's)
	c7->cd(1);
	h_dp->Draw();              // already fit above; drawing will overlay the fit automatically
	g_dp->Draw("same");        // ensure fit line is visible
	add_mu_sigma_box(g_dp);

	// Panel 2: h_dp_new (Sean's vy corr)
	c7->cd(2);
	h_dp_new->Draw();
	g_dp_new->Draw("same");
	add_mu_sigma_box(g_dp_new);

	// Panel 3: h_dp_new1 (New vy corr)
	c7->cd(3);
	h_dp_new1->Draw();
	g_dp_new1->Draw("same");
	add_mu_sigma_box(g_dp_new1);

	// Panel 4: overlay all three + their fit curves
	c7->cd(4);
	h_dp->SetLineWidth(3);
	h_dp_new->SetLineWidth(3);
	h_dp_new1->SetLineWidth(3);
	h_dp->SetLineColor(kBlack);
	h_dp_new->SetLineColor(kBlue);
	h_dp_new1->SetLineColor(kRed);

	// draw in a consistent order so curves are visible
	h_dp_new1->Draw();         // red
	h_dp->Draw("same");        // black
	h_dp_new->Draw("same");    // blue

	// overlay the three fit curves on the combined plot
	g_dp     ->Draw("same");
	g_dp_new ->Draw("same");
	g_dp_new1->Draw("same");



	auto legdeltaP = new TLegend(0.62, 0.60, 0.88, 0.88); // x1,y1,x2,y2 — tweak as needed
	legdeltaP->SetBorderSize(0);
	legdeltaP->SetFillStyle(0);           // transparent legend box
	legdeltaP->SetTextFont(42);
	legdeltaP->SetTextSize(0.035);
	//leg->SetNColumns(2);            // optional: 2 columns for 4 entries

	legdeltaP->AddEntry(h_dp,  "Sean's",    "l");
	legdeltaP->AddEntry(h_dp_new, "Sean's vy corr","l");
	legdeltaP->AddEntry(h_dp_new1, "New vy corr","l");

	legdeltaP->Draw();
	gPad->Update();

	TCanvas *c8 =  new TCanvas("c8","c8",1600,1200);
	c8->Divide(2,2);
	c8->cd(1);
	h_p_el_v_rastery->Draw("COLZ");
	c8->cd(2);
	h_p_el_v_rasterx->Draw("COLZ");
	c8->cd(3);
	h_p_el_v_raster2y->Draw("COLZ");
	c8->cd(4);
	h_p_el_v_raster2x->Draw("COLZ");

	TCanvas *c9 =  new TCanvas("c9","c9",1600,1200);
	c9->Divide(2,2);
	c9->cd(1);
	h_ebeam->SetLineColor(kBlack);
	h_ebeam->SetLineWidth(3);
	h_ebeam->Draw();


	TCanvas *c10 = new TCanvas("c10","c10: W2 sensitivity to momentum",1600,1200);
	c10->Divide(2,2);

	// (a) 2D scatter ΔW2 vs Δp
	c10->cd(1);
	h_dW2_vs_dpf->Draw("COLZ");

	// (b) Profile + linear fit (slope ~ average sensitivity)
	c10->cd(2);
	p_dW2_vs_dpf->SetMarkerStyle(20);
	p_dW2_vs_dpf->SetMarkerSize(0.9);
	p_dW2_vs_dpf->Draw("E1");
	TF1 *fS = new TF1("fS","pol1",-0.03,0.03); // restrict to small |Δp|
	p_dW2_vs_dpf->Fit(fS,"RQ");
	fS->SetLineWidth(3);
	fS->SetLineColor(kRed);
	fS->Draw("same");

	// annotate slope
	{
	  double slope = fS->GetParameter(1);
	  double eslope = fS->GetParError(1);
	  TPaveText *pt = new TPaveText(0.15,0.78,0.55,0.90,"NDC");
	  pt->SetFillStyle(0); pt->SetTextFont(42); pt->SetTextSize(0.035);
	  pt->AddText(Form("Slope = dW^{2}/d(\\#Delta p) = %.4f #pm %.4f GeV^{2}", slope, eslope));
	  pt->Draw();
	}

	// (c) Compare analytic S distribution
	c10->cd(3);
	h_S_analytic->SetLineColor(kBlue+1);
	h_S_analytic->SetLineWidth(2);
	h_S_analytic->Draw();

	// (d) Empirical S = ΔW2/Δp distribution (sanity check)
	c10->cd(4);
	h_S_empirical->SetLineColor(kGreen+2);
	h_S_empirical->SetLineWidth(2);
	h_S_empirical->Draw();

	c->Print(Form("../plots/%s_ePS_eHCAL.pdf(",printfilename));	
	c1->Print(Form("../plots/%s_ePS_eHCAL.pdf)",printfilename));
	c->SaveAs(Form("../plots/%s_ePS_eHCAL.png",printfilename));
	c1->SaveAs(Form("../plots/%s_ePS_eHCAL_1.png",printfilename));
	c2->SaveAs(Form("../plots/%s_ePS_eHCAL_2.png",printfilename));
	c3->SaveAs(Form("../plots/%s_ePS_eHCAL_3.png",printfilename));
	c4->Print(Form("../plots/%s_W2_test.pdf(",printfilename));
	c5->Print(Form("../plots/%s_W2_test.pdf",printfilename));
	c6->Print(Form("../plots/%s_W2_test.pdf",printfilename));
	c7->Print(Form("../plots/%s_W2_test.pdf",printfilename));
	c8->Print(Form("../plots/%s_W2_test.pdf",printfilename));
	c9->Print(Form("../plots/%s_W2_test.pdf",printfilename));
	c10->Print(Form("../plots/%s_W2_test.pdf)",printfilename));


}
