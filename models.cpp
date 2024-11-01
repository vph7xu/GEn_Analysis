//this is to return models from the simulation and for the bkg

#include "models.h"
#include "parse.h"

TH1D* gExistingHist = nullptr;

Double_t customFunction(Double_t *x, Double_t *par) {
    if (gExistingHist == nullptr) return 0.0;
    Int_t bin = gExistingHist->FindBin(x[0]);
    return gExistingHist->GetBinContent(bin);
}

TH1D* bkg_model(TH2D *hist2D, double xmin, double xmax, const char* kin){

	//std::map<std::string, std::string> config = parseConfig(Form("cut_%s.txt",kin));

	TH2D *h2_anticut = (TH2D*)hist2D->Clone("h2_anticut");
   	h2_anticut->Reset();
	
	int nx = hist2D->GetNbinsX();
	int ny = hist2D->GetNbinsY();
    	
	for (int ix = 1; ix <= nx; ++ix) {
        	for (int iy = 1; iy <= ny; ++iy) {
            		double x = hist2D->GetXaxis()->GetBinCenter(ix);
            		double y = hist2D->GetYaxis()->GetBinCenter(iy);

            		if (!(x > xmin && x < xmax)) {
                		double content = hist2D->GetBinContent(ix, iy);
              	  		h2_anticut->SetBinContent(ix, iy, content);
         	   	}
       		}
	}

	TH1D *projY = h2_anticut->ProjectionX();
        TH1D *projX = h2_anticut->ProjectionY();

	TCanvas *c = new TCanvas("c","c",1600,1200);
	c->Divide(2,2);
	c->cd(1);
	projY->SetTitle("dy projection");
	projY->Draw("HIST");
	c->cd(2);
	projX->SetTitle("dx projection");
	projX->Draw("HIST");
	c->cd(3);
	h2_anticut->Draw("COLZ");
	c->cd(4);
	hist2D->GetXaxis()->SetRangeUser(-4,4);
        hist2D->GetYaxis()->SetRangeUser(-4,4);
	hist2D->Draw("COLZ");

	c->Print(Form("../plots/bkg_%s.pdf",kin));
	c->SaveAs(Form("bkg_%s.png",kin));

	gExistingHist = projY;

	TF1 *model = new TF1("model", customFunction,projY->GetXaxis()->GetXmin(),projY->GetXaxis()->GetXmax(),0);

	return projX;

}

std::pair<TH1D*, TH1D*> sim_hist(const char* sim_filename, double W2_L_sim, double W2_H_sim, double dy_L_sim, double dy_H_sim){
	TFile* sim_file = TFile::Open(sim_filename);

	TTree* sim_tree = (TTree*)sim_file->Get("Tout");

	double dx = 0.0;
	double dy = 0.0;
	double W2 = 0.0;
	double Q2 = 0.0;
	double weight = 0.0;
	double fnucl = 0.0;

	sim_tree->SetBranchAddress("dx",&dx);
	sim_tree->SetBranchAddress("dy",&dy);
	sim_tree->SetBranchAddress("W2",&W2);
	sim_tree->SetBranchAddress("Q2",&Q2);
	sim_tree->SetBranchAddress("weight",&weight);
	sim_tree->SetBranchAddress("fnucl",&fnucl);

        TH1D *h_dx_n = new TH1D("h_dx_n","dx neutrons",100,-4,3);
        TH1D *h_dx_p = new TH1D("h_dx_p","dx_protons",100,-4,3);

	int nentries = sim_tree->GetEntries();

	for (int i = 0; i<nentries; i++){
                sim_tree->GetEntry(i);

		if (W2_L_sim<W2 and W2<W2_H_sim and dy_L_sim<dy and dy<dy_H_sim){
                        if (fnucl == 0.0){
                                h_dx_n->Fill(dx,weight);
                        }
                        else if (fnucl == 1.0){
                                h_dx_p->Fill(dx,weight);
                        }
                }

                if (i %1000 == 0 ) std::cout<< (i * 100.0/ nentries)<<"sim file" << "% \r";
                std::cout.flush();
	}

	//TCanvas* ctest = new TCanvas("ctest","ctest",800,1200);
	//h_dx_n->Draw();

	return std::make_pair(h_dx_n,h_dx_p);
}


void models(const char* filename,const char* sim_filename,const char* printfilename, const char* kin){
        
	//read the sampling fraction for each blk
	std::ifstream inFile("sampling_fractions_each_blk.txt");	

	std::vector<double> sampling_fractions_each_blk;

    	std::string line;
    	while (std::getline(inFile, line)) {
        	std::stringstream ss(line);
        	double iblk, sampling_fraction;

        	// Extract the two columns (index and Y mean)
        	ss >> iblk >> sampling_fraction;

        	// Add the Y mean to the vector
        	sampling_fractions_each_blk.push_back(sampling_fraction);
    	}	

    	inFile.close();

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

	double dy_ac_L = getDoubleValue(config,"dy_ac_L");
	double dy_ac_H = getDoubleValue(config,"dy_ac_H");

	double run_num_L = getDoubleValue(config,"run_num_L");
	double run_num_H = getDoubleValue(config,"run_num_H");

	double IHWP_flip = getDoubleValue(config,"IHWP_flip");

	//end of parsing cuts
	
	TFile* file = TFile::Open(filename);
	//TFile* sim_file = TFile::Open(sim_filename);

        TTree* tree = (TTree*)file->Get("Tout");
	//TTree* sim_tree = (TTree*)sim_file->Get("Tout");

	//p and n dist
	auto sim_histograms = sim_hist(sim_filename,W2_L,W2_H,dy_L,dy_H);
	hist_p = sim_histograms.second;
	hist_n = sim_histograms.first;

	int runnum = 0;
	int helicity = 0;
	int IHWP = 0;
        double dx = 0.0;
        double dy = 0.0;
        double W2 = 0.0;
        double Q2 = 0.0;
	double coin_time = 0.0;
	double ntrack = 0;
	double theta_pq = 0.0;
	double pN_expect = 0.0;
	double eHCAL = 0.0;
	double nblk_HCAL = 0.0;
	double hcal_clus_e[1000];
	double hcal_clus_id[1000];
	double hcal_clus_mem_e[1000];
	double hcal_clus_mem_id[1000];
	double trP_sbs = 0.0;
	double ntrack_sbs = 0.0;
	double vz = 0.0;	

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("IHWP",&IHWP);
        tree->SetBranchAddress("dx",&dx);
        tree->SetBranchAddress("dy",&dy);
        tree->SetBranchAddress("W2",&W2);
        tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("ntrack", &ntrack);
	tree->SetBranchAddress("theta_pq",&theta_pq);
	tree->SetBranchAddress("pN_expect",&pN_expect);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("hcal_clus_e",&hcal_clus_e);
	tree->SetBranchAddress("hcal_clus_mem_e",&hcal_clus_mem_e);
	tree->SetBranchAddress("hcal_clus_id",&hcal_clus_id);
	tree->SetBranchAddress("hcal_clus_mem_id",&hcal_clus_mem_id);
	tree->SetBranchAddress("nblk_HCAL",&nblk_HCAL);
	tree->SetBranchAddress("trP_sbs",&trP_sbs);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
	tree->SetBranchAddress("vz",&vz);


        TH1D *h_dx = new TH1D("h_dx","dx",200,-10,10);
        TH1D *h_dy = new TH1D("h_dy","dy",200,-10,10);
        TH1D *h_W2 = new TH1D("h_W2","W2",250,-4,8);
        TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",100,dy_L,dy_H,100,-4,3);
        TH1D *h_coin_time = new TH1D("h_coin_time","coin_time",1000,40,200);

        TH1D *h_dx_W2_cut = new TH1D("h_dx_W2_cut","dx dist : data/sim comparison",100,-4,3);
        TH1D *h_dx_W2_cut_plotting = new TH1D("h_dx_W2_cut_plotting","dx dist : data/sim comparison",100,-4,3);
        TH1D *h_dy_W2_cut = new TH1D("h_dy_W2_cut","dy",100,-4,3);
        //TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
        TH2D *h_dxdy_W2_cut = new TH2D("h_dxdy_W2_cut","dx v dy",100,dy_L,dy_H,100,-4,3);
        TH2D *h_dxdy_bkg = new TH2D("h_dxdy_bkg","dxdy (anticut shaded)",100,-4,4,100,-4,3);

	TH1D *h_dx_sim_n_bkg = new TH1D("h_dx_sim_n_bkg","dx distribution : data/sim comparison",100,-4,3);

        int nentries = tree->GetEntries();
        double Nplus = 0.0;
        double Nminus = 0.0;

	std::cout<<"coin_L and H : "<<coin_time_L<<"  "<<coin_time_H<<endl;
	std::cout<<"dx_L and H : "<<dy_L <<" "<<dy_H <<endl;

	//get the bkg from data
	for (int i = 0; i<nentries; i++){
                tree->GetEntry(i);
                if(lookupValue(HelicityCheck,runnum)==1 and lookupValue(MollerQuality,runnum)==1){

 			double KinE = 0.0;

			for (int j = 0; j<nblk_HCAL; j++){
			//KinE += hcal_clus_e[j]/sampling_fractions_each_blk[hcal_clus_id[j]];
				KinE += hcal_clus_mem_e[j]/sampling_fractions_each_blk[hcal_clus_mem_id[j]];

				//if (i%1000 == 0){ 
				//std::cout<<"event : "<<i<<"blk id: "<<hcal_clus_mem_id[j]<<" blk energy : "<<hcal_clus_mem_e[j]<<" blk sf : " <<sampling_fractions_each_blk[hcal_clus_mem_id[j]]<<endl;
				//}
			}
			//if (i%1000 == 0){ 
			//std::cout<<"event : "<<i<<" kinE : "<<KinE<<endl;
			//}
			//double Pperp = theta_pq * pN_expect;
			double realPperp = (theta_pq * sqrt(-pow(0.938,2)+pow((KinE+0.938),2))); // assuming mN and mP is equal
			double realPpar = (cos(theta_pq) * sqrt(-pow(0.938,2)+pow((KinE+0.938),2)))-pN_expect;
			double Pmiss = sqrt(pow(realPperp,2)+pow(realPpar,2));          

			//before adding a cut on W2

                	helicity = IHWP*IHWP_flip*helicity;

                	h_dx->Fill(dx);
                	h_dy->Fill(dy);
                	h_dxdy->Fill(dy,dx);
                	h_W2->Fill(W2);

                	bool cut_Ppar = abs(realPpar)<1.5;
                	bool cut_eHCAL = (eHCAL)>0.3;

                	if (W2_L<W2 and W2<W2_H and cut_eHCAL){
                		h_coin_time->Fill(coin_time);
                	}
	        
			//add a cut on W2
                	if (cut_eHCAL and (W2_L<W2 and W2<W2_H) and (coin_time_L<coin_time and coin_time<coin_time_H) and dy_L<dy and dy<dy_H){
                		h_dx_W2_cut_plotting->Fill(dx);
                        	h_dx_W2_cut->Fill(dx);
                        	h_dy_W2_cut->Fill(dy);
                        	h_dxdy_W2_cut->Fill(dy,dx);
                	}

                	if(cut_eHCAL and (W2_L<W2 and W2<W2_H) and (coin_time_L<coin_time and coin_time<coin_time_H)){

                		h_dxdy_bkg->Fill(dy,dx);
                	}

                	//inelastic asymmetry using anticut, no W2, Ppar cut here 
                	if((coin_time_L<coin_time and coin_time<coin_time_H) and (dx_L<dx and dx<dx_H) and (dy_L>dy or dy>dy_H)){
                		if (helicity == 1){
                			Nplus += 1;
                		}
                		else if (helicity == -1){
                			Nminus += 1;
                		}
                	}
                	
        	}
                if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();

        }
        //add lines to W2 plot to show the cut
	//
	
	hist_bkg=bkg_model(h_dxdy_bkg, dy_ac_L, dy_ac_H, kin);//background dist, hard coded for GEN3
	
	//scale everything
	double scale_data = h_dx_W2_cut->Integral();
	h_dx_W2_cut->Scale(1.0/h_dx_W2_cut->Integral());
	hist_p->Scale(1.0/hist_p->Integral());
	hist_n->Scale(1.0/hist_n->Integral());
	hist_bkg->Scale(1.0/hist_bkg->Integral());

	//bkg_model(h_dxdy_W2_cut,-1.5,1.5, kin);
	TBox *box1 = new TBox(1.0,-4,2,3.05);
	box1->SetLineColor(kRed);
	box1->SetFillStyle(3002);
	box1->SetLineWidth(3);

	TBox *box2 = new TBox(-1.5,-4,-2,3.05);
	box2->SetLineColor(kRed);
	box2->SetFillStyle(3002);	
	box2->SetLineWidth(3);

	TCanvas *c1 = new TCanvas("c1","c1",3600,3000);
	c1->Divide(2,2);
	c1->cd(1);
	h_dxdy_W2_cut->Draw("COLZ");
	c1->cd(2);
	h_dxdy_bkg->SetStats(0);
	h_dxdy_bkg->GetXaxis()->SetRangeUser(-2,2);
	h_dxdy_bkg->SetYTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dxdy_bkg->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	h_dxdy_bkg->Draw("COLZ");
	box1->Draw("same");
	box2->Draw("same");
	c1->SaveAs(Form("%s_dxdy.png",kin));

	TF1 *fit_data = new TF1("fit_data", fit_sim_n_bkg, h_dx_W2_cut->GetXaxis()->GetXmin(), h_dx_W2_cut->GetXaxis()->GetXmax(),3);

	// Set better initial guesses
	double init_param0 = 1;
	double init_param1 = 1;
	double init_param2 = 1;

	//fit_data->SetParameters(init_param0, init_param1, init_param2);

	// Apply constraints
	fit_data->SetParLimits(0, 0.1, 100);
	fit_data->SetParLimits(1, 0.1, 100);
	fit_data->SetParLimits(2, 0.0, 100);

	// Fit the data with log-likelihood
	h_dx_W2_cut->Fit(fit_data, "RL");

	hist_p->Scale(fit_data->GetParameter(0));
	hist_n->Scale(fit_data->GetParameter(1)*fit_data->GetParameter(0));
	hist_bkg->Scale(fit_data->GetParameter(2)*fit_data->GetParameter(0));

	hist_p->SetLineColor(6);
	hist_n->SetLineColor(9);
	hist_bkg->SetLineColor(7);
	h_dx_sim_n_bkg->SetLineColor(3);

	hist_p->SetLineWidth(3);
	hist_n->SetLineWidth(3);
	hist_bkg->SetLineWidth(3);
	h_dx_sim_n_bkg->SetLineWidth(3);

	hist_p->Scale(scale_data);
	hist_n->Scale(scale_data);
	hist_bkg->Scale(scale_data);
	h_dx_W2_cut->Scale(scale_data);
	
	h_dx_W2_cut->GetXaxis()->SetRangeUser(-6,6);
	h_dx_W2_cut->SetMarkerStyle(20);
	h_dx_W2_cut->SetMarkerColor(kBlack);

	h_dx_W2_cut_plotting->GetXaxis()->SetRangeUser(-6,6);
	h_dx_W2_cut_plotting->SetMarkerStyle(20);
	h_dx_W2_cut_plotting->SetMarkerColor(kBlack);

	hist_p->SetFillColorAlpha(6,0.5);
	hist_p->SetFillStyle(3004);
	hist_n->SetFillColorAlpha(9,0.5);
	hist_n->SetFillStyle(3005);
	hist_bkg->SetFillColorAlpha(7,0.5);
	hist_bkg->SetFillStyle(3003);
	h_dx_sim_n_bkg->SetFillColorAlpha(19,0.1);
	h_dx_sim_n_bkg->SetFillStyle(3009);

	//marker style and color
	//h_dx_W2_cut->SetMarkerStyle(20);
	//h_dx_W2_cut->SetMarkerColor(kBlack);

	//h_dx_sim_n_bkg->SetMarkerStyle(30);
	//h_dx_sim_n_bkg->SetMarkerColor(kRed);


	//fill up the total model histogram

	int nBins = h_dx_sim_n_bkg->GetNbinsX();

	for (int i = 1; i<=nBins; ++i){
		h_dx_sim_n_bkg->SetBinContent(i,(hist_p->GetBinContent(i)+hist_n->GetBinContent(i)+hist_bkg->GetBinContent(i)));	
	}

	//get the inelastic fraction
	double inelastic_events = 0.0; 
	double QE_events = 0.0;
	double inelastic_frac = 0.0;
	double errinelastic_frac = 0.0;

	inelastic_events = hist_bkg->Integral(hist_bkg->FindBin(dx_L),hist_bkg->FindBin(dx_H)); //integrate background hist to get inelastic events under the neutron peak
	QE_events = h_dx_W2_cut->Integral(h_dx_W2_cut->FindBin(dx_L),h_dx_W2_cut->FindBin(dx_H)); //integrate data to get the QE neutrons within the dx cut
	inelastic_frac = inelastic_events/QE_events;// this is not correct should remove other fractions before doing this
	errinelastic_frac = (inelastic_events/QE_events)*sqrt((1/inelastic_events)+(1/QE_events)); // this is not correct should remove other fractions before doing this
	
	//get the inelastic asymmetry
	double Ain = (Nplus-Nminus)/(Nplus+Nminus);
	double errAin = 2*sqrt((Nplus*Nminus)*(Nplus+Nminus))/((Nplus+Nminus)*(Nplus+Nminus));;

        TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);
        legend->AddEntry(h_dx_W2_cut_plotting,"Data","p");
        legend->AddEntry(hist_p,"sim proton","lf");
        legend->AddEntry(hist_n,"sim neutron","lf");
        legend->AddEntry(hist_bkg,"bkg data","lf");
        legend->AddEntry(h_dx_sim_n_bkg,"full model","lf");
        //legend->Draw();

	TLine *line1 = new TLine(dx_L,0.0,dx_L,2000);
	line1->SetLineColor(kRed);
	line1->SetLineWidth(2);
	
	TLine *line2 = new TLine(dx_H,0.0,dx_H,2000);
	line2->SetLineColor(kRed);
	line2->SetLineWidth(2);

	TCanvas* cfit = new TCanvas("cfit","Fit with models",3600,3000);
	cfit->Divide(2,2);
	cfit->cd(1);
	h_dx_sim_n_bkg->SetXTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dx_sim_n_bkg->Draw("HIST");
	h_dx_W2_cut_plotting->Draw("P SAME");
	hist_p->Draw("HIST SAME");
	hist_n->Draw("HIST SAME");
	hist_bkg->Draw("HIST SAME");
	line1->Draw("SAME");
	line2->Draw("SAME");

	legend->Draw();
	//h_dx_sim_n_bkg->Add(hist_p);
	//h_dx_sim_n_bkg->Add(hist_n);
	//h_dx_sim_n_bkg->Add(hist_bkg);

	cfit->cd(2);
	h_dx_sim_n_bkg->Draw();

	cfit->cd(3);
	h_coin_time->Draw();



	std::cout<<" Neutrons : "<<hist_n->Integral()<<endl;


	std::cout<<"Nplus : "<<Nplus<<endl;
	std::cout<<"Nminus : "<<Nminus<<endl;
	std::cout<<"Ain : "<<Ain<<endl;
	std::cout<<"errAin : "<<errAin<<endl;
	std::cout<<"inelastic_events : "<< inelastic_events<<endl;
	std::cout<<"QE_events : "<<QE_events<<endl;
	std::cout<<"inelastic fraction (double counting) : "<<inelastic_frac<<endl; //not correct , double counting
	std::cout<<"error inelastic fraction (double counting) : "<<errinelastic_frac<<endl; //not correct , double counting
	//fit_data->Draw();

	cfit->Print(Form("models_%s.pdf",kin));
	cfit->SaveAs(Form("models_%s.png",kin));
	cfit->SaveAs(Form("models_%s.jpg",kin));
} 

