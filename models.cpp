//this is to return models from the simulation and for the bkg

#include "models.h"
#include "parse.h"

TH1D* gExistingHist = nullptr;

Double_t customFunction(Double_t *x, Double_t *par) {
    if (gExistingHist == nullptr) return 0.0;
    Int_t bin = gExistingHist->FindBin(x[0]);
    return gExistingHist->GetBinContent(bin);
}

TH1D* bkg_model(TH2D *hist2D, double xmin, double xmax, const char* kin, bool flag_eHCAL_cut, double eHCAL_L){

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

	c->Print(Form("plots/bkg_%s_eHCAL_cut_%s_%f.pdf",kin,std::to_string(flag_eHCAL_cut).c_str(),eHCAL_L));
	c->SaveAs(Form("plots/bkg_%s_eHCAL_cut_%s_%f.png",kin,std::to_string(flag_eHCAL_cut).c_str(),eHCAL_L));

	gExistingHist = projY;

	TF1 *model = new TF1("model", customFunction,projY->GetXaxis()->GetXmin(),projY->GetXaxis()->GetXmax(),0);

	return projX;

}

std::pair<TH1D*, TH1D*> sim_hist(const char* sim_filename, double W2_L_sim, double W2_H_sim, double dy_L_sim, double dy_H_sim, double eHCAL_L_sim){
	TFile* sim_file = TFile::Open(sim_filename);

	TTree* sim_tree = (TTree*)sim_file->Get("Tout");

	double dx = 0.0;
	double dy = 0.0;
	double W2 = 0.0;
	double Q2 = 0.0;
	double weight = 0.0;
	double fnucl = 0.0;
	double eHCAL = 0.0;

	sim_tree->SetBranchAddress("dx",&dx);
	sim_tree->SetBranchAddress("dy",&dy);
	sim_tree->SetBranchAddress("W2",&W2);
	sim_tree->SetBranchAddress("Q2",&Q2);
	sim_tree->SetBranchAddress("weight",&weight);
	sim_tree->SetBranchAddress("fnucl",&fnucl);
	sim_tree->SetBranchAddress("eHCAL",&eHCAL);


        TH1D *h_dx_n = new TH1D("h_dx_n","dx neutrons",100,-4,3);
        TH1D *h_dx_p = new TH1D("h_dx_p","dx_protons",100,-4,3);

	int nentries = sim_tree->GetEntries();

	for (int i = 0; i<nentries; i++){
                sim_tree->GetEntry(i);

		if (W2_L_sim<W2 and W2<W2_H_sim and dy_L_sim<dy and dy<dy_H_sim and eHCAL_L_sim<eHCAL){
                        if (fnucl == 0.0){
                                h_dx_n->Fill(dx,weight);
                        }
                        else if (fnucl == 1.0){
                                h_dx_p->Fill(dx+0.05,weight); // change this later
                        }
                }

                if (i %1000 == 0 ) std::cout<< (i * 100.0/ nentries)<<"sim file" << "% \r";
                std::cout.flush();
	}

	//TCanvas* ctest = new TCanvas("ctest","ctest",800,1200);
	//h_dx_n->Draw();

	return std::make_pair(h_dx_n,h_dx_p);
}


void models(const char* filename,const char* sim_filename,const char* printfilename, const char* kin, bool flag_eHCAL_cut){
        
	//read the sampling fraction for each blk
	std::ifstream inFile(Form("txt/%s_sampling_fractions_each_blk.txt",kin));	

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

	double eHCAL_L = getDoubleValue(config,"eHCAL_L");

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
	auto sim_histograms = sim_hist(sim_filename,W2_L,W2_H,dy_L,dy_H,eHCAL_L);
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
    	double grinch_clus_size = 0.0;
    	double grinch_track = 0.0;
    	double ePS       = 0.0;

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
    	tree->SetBranchAddress("grinch_track",  &grinch_track);
    	tree->SetBranchAddress("grinch_clus_size",   &grinch_clus_size);
	tree->SetBranchAddress("ePS",       &ePS);

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
	std::cout<<"dy_L and H : "<<dy_L <<" "<<dy_H <<endl;
	std::cout<<"dx_L and H : "<<dx_L <<" "<<dx_H <<endl;
	std::cout<<"W2_L and H : "<<W2_L <<" "<<W2_H <<endl;

	if (flag_eHCAL_cut == false) eHCAL_L=0.0;

    	std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;

	//get the bkg from data
	for (int i = 0; i<nentries; i++){
                tree->GetEntry(i);

	        // some conditions from your DB checks:
	        bool goodHelicity = (lookupValue(HelicityCheck, runnum) == 1);
	        bool goodMoller   = (lookupValue(MollerQuality, runnum) == 1);
	        bool goodVz       = abs(vz) < 0.27;
	        bool goodPS       = (ePS > 0.2);
	        bool goodRunRange = (run_num_L < runnum && runnum < run_num_H);
	        bool goodEHCAL    = (eHCAL > eHCAL_L); 
	        bool validHel     = (helicity == -1 || helicity == 1);
	        bool goodGrinch = (grinch_track == 0) && (grinch_clus_size>2);

                if(goodHelicity && goodMoller && goodPS && validHel && goodRunRange && goodEHCAL && goodVz && goodGrinch){

 			//double KinE = 0.0;

			//for (int j = 0; j<nblk_HCAL; j++){
			//KinE += hcal_clus_e[j]/sampling_fractions_each_blk[hcal_clus_id[j]];
			//	KinE += hcal_clus_mem_e[j]/sampling_fractions_each_blk[hcal_clus_mem_id[j]];

				//if (i%1000 == 0){ 
				//std::cout<<"event : "<<i<<"blk id: "<<hcal_clus_mem_id[j]<<" blk energy : "<<hcal_clus_mem_e[j]<<" blk sf : " <<sampling_fractions_each_blk[hcal_clus_mem_id[j]]<<endl;
				//}
			//}
			//if (i%1000 == 0){ 
			//std::cout<<"event : "<<i<<" kinE : "<<KinE<<endl;
			//}
			//double Pperp = theta_pq * pN_expect;
			//double realPperp = (theta_pq * sqrt(-pow(0.938,2)+pow((KinE+0.938),2))); // assuming mN and mP is equal
			//double realPpar = (cos(theta_pq) * sqrt(-pow(0.938,2)+pow((KinE+0.938),2)))-pN_expect;
			//double Pmiss = sqrt(pow(realPperp,2)+pow(realPpar,2));          

			//before adding a cut on W2

                	helicity = -1*IHWP*IHWP_flip*helicity;

                	h_dx->Fill(dx);
                	h_dy->Fill(dy);
                	h_dxdy->Fill(dy,dx);
                	h_W2->Fill(W2);

                	//bool cut_Ppar = abs(realPpar)<1.5;

                	if (W2_L<W2 && W2<W2_H){
                		h_coin_time->Fill(coin_time);
                	}
	        
			//add a cut on W2
                	if ((W2_L<W2 && W2<W2_H) && (coin_time_L<coin_time && coin_time<coin_time_H) && dy_L<dy && dy<dy_H){
                		h_dx_W2_cut_plotting->Fill(dx);
                        	h_dx_W2_cut->Fill(dx);
                        	h_dy_W2_cut->Fill(dy);
                        	h_dxdy_W2_cut->Fill(dy,dx);
                	}

                	if((W2_L<W2 && W2<W2_H) && (coin_time_L<coin_time && coin_time<coin_time_H)){

                		h_dxdy_bkg->Fill(dy,dx);
                	}

                	//inelastic asymmetry using anticut, no W2, Ppar cut here 
                	if((coin_time_L<coin_time && coin_time<coin_time_H) && (dx_L<dx && dx<dx_H) && (dy_L>dy && dy>dy_H)){
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
	
	hist_bkg=bkg_model(h_dxdy_bkg, dy_ac_L, dy_ac_H, kin, flag_eHCAL_cut, eHCAL_L);//background dist, hard coded for GEN3
	
	//scale everything
	double scale_data = h_dx_W2_cut->Integral();
	h_dx_W2_cut->Scale(1.0/h_dx_W2_cut->Integral());
	hist_p->Scale(1.0/hist_p->Integral());
	hist_n->Scale(1.0/hist_n->Integral());
	hist_bkg->Scale(1.0/hist_bkg->Integral());

	//bkg_model(h_dxdy_W2_cut,-1.5,1.5, kin);
	TBox *box1 = new TBox(dy_ac_H,-4,4,3.05);
	box1->SetLineColor(kRed);
	box1->SetFillStyle(3002);
	box1->SetLineWidth(3);

	TBox *box2 = new TBox(dy_ac_L,-4,-4,3.05);
	box2->SetLineColor(kRed);
	box2->SetFillStyle(3002);	
	box2->SetLineWidth(3);

	TCanvas *c1 = new TCanvas("c1","c1",3600,3000);
	c1->Divide(2,2);
	c1->cd(1);
	h_dxdy_W2_cut->Draw("COLZ");
	c1->cd(2);
	h_dxdy_bkg->SetStats(0);
	h_dxdy_bkg->GetXaxis()->SetRangeUser(-4,4);
	h_dxdy_bkg->SetYTitle("HCAL_X(exp)-HCAL_X(act) (m)");
	h_dxdy_bkg->SetXTitle("HCAL_Y(exp)-HCAL_Y(act) (m)");
	h_dxdy_bkg->Draw("COLZ");
	box1->Draw("same");
	box2->Draw("same");
	c1->SaveAs(Form("plots/%s_dxdy.png",kin));

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
	double N_events = 0.0;
	double inelastic_frac = 0.0;
	double errinelastic_frac = 0.0;

	inelastic_events = hist_bkg->Integral(hist_bkg->FindBin(dx_L),hist_bkg->FindBin(dx_H)); //integrate background hist to get inelastic events under the neutron peak
	QE_events = h_dx_W2_cut->Integral(h_dx_W2_cut->FindBin(dx_L),h_dx_W2_cut->FindBin(dx_H)); //integrate data to get the QE neutrons within the dx cut
	N_events = hist_n->Integral(hist_n->FindBin(dx_L),hist_n->FindBin(dx_H));
	inelastic_frac = inelastic_events/QE_events;// this is not correct should remove other fractions before doing this
	errinelastic_frac = (inelastic_events/QE_events)*sqrt((1/inelastic_events)+(1/QE_events)); // this is not correct should remove other fractions before doing this
	
	//get the inelastic asymmetry
	double Ain = (Nplus-Nminus)/(Nplus+Nminus);
	double errAin = 2*sqrt((Nplus*Nminus)*(Nplus+Nminus))/((Nplus+Nminus)*(Nplus+Nminus));;

        TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);
        legend->AddEntry(h_dx_W2_cut_plotting,"Data","p");
        legend->AddEntry(hist_p,"simulated protons","lf");
        legend->AddEntry(hist_n,"simulated neutrons","lf");
        legend->AddEntry(hist_bkg,"background from data","lf");
        legend->AddEntry(h_dx_sim_n_bkg,"N(sim_p+R*sim_n+Nbg*bkg)","lf");
        //legend->Draw();

	TLine *line1 = new TLine(dx_L,0.0,dx_L,hist_n->GetMaximum());
	line1->SetLineColor(kRed);
	line1->SetLineWidth(2);
	
	TLine *line2 = new TLine(dx_H,0.0,dx_H,hist_n->GetMaximum());
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

	cfit->cd(4);
	h_dx_W2_cut_plotting->Draw();


	std::cout<<" Neutrons : "<<hist_n->Integral()<<endl;


	std::cout<<"Nplus : "<<Nplus<<endl;
	std::cout<<"Nminus : "<<Nminus<<endl;
	std::cout<<"Ain : "<<Ain<<endl;
	std::cout<<"errAin : "<<errAin<<endl;
	std::cout<<"inelastic_events : "<< inelastic_events<<endl;
	std::cout<<"QE_events : "<<QE_events<<endl;
	std::cout<<"N_events : "<<N_events;
	std::cout<<"inelastic fraction (double counting) : "<<inelastic_frac<<endl; //not correct , double counting
	std::cout<<"error inelastic fraction (double counting) : "<<errinelastic_frac<<endl; //not correct , double counting
	//fit_data->Draw();

    	std::ofstream outfile;
    	outfile.open(Form("txt/%s_inelastic_asymmetry_eHCAL_cut_%s_%f.txt",kin,std::to_string(flag_eHCAL_cut).c_str(),eHCAL_L));
    	outfile<<"N_plus = "<<Nplus<<endl;
    	outfile<<"N_minus = "<<Nminus<<endl;
    	outfile<<"A_in = "<<Ain<<endl;
    	outfile<<"err_A_in= "<<errAin<<endl;
    	outfile<<"inelastic_events = "<<inelastic_events<<endl;
    	outfile<<"QE_events = "<<QE_events<<endl;
    	outfile<<"N_events = "<<N_events<<endl;
    	outfile<<"f_bkg = "<<inelastic_frac<<endl;
    	outfile<<"err_f_bkg = "<<errinelastic_frac<<endl;

	cfit->Print(Form("plots/models_%s_eHCAL_cut_%s_%f.pdf",kin,std::to_string(flag_eHCAL_cut).c_str(),eHCAL_L));
	cfit->SaveAs(Form("plots/models_%s_eHCAL_cut_%s_%f.png",kin,std::to_string(flag_eHCAL_cut).c_str(),eHCAL_L));
	cfit->SaveAs(Form("plots/models_%s_eHCAL_cut_%s_%f.jpg",kin,std::to_string(flag_eHCAL_cut).c_str(),eHCAL_L));

} 

