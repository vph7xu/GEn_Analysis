void plotdxdy_sim(const char* filename,const char* printfilename){
	
	TFile* file = TFile::Open(filename);
	TTree* tree = (TTree*)file->Get("Tout");

	double dx = 0.0; 
	double dy = 0.0; 
	double W2 = 0.0;
        double Q2 = 0.0;
	double weight = 0.0;
	double fnucl = 0.0;

	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("weight",&weight);
	tree->SetBranchAddress("fnucl",&fnucl);

	TH1D *h_dx = new TH1D("h_dx","dx",500,-10,10);
	TH1D *h_dy = new TH1D("h_dy","dy",500,-10,10);
	TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
	TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",250,-2,2,250,-4,4);

	TH1D *h_dx_W2_cut = new TH1D("h_dx_W2_cut","dx after W2 cut",1000,-10,10);
	TH1D *h_dy_W2_cut = new TH1D("h_dy_W2_cut","dy after W2 cut",1000,-10,10);
	//TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
	TH2D *h_dxdy_W2_cut = new TH2D("h_dxdy_W2_cut","dxdy after W2 cut",250,-2,2,250,-4,4);
	
	TH1D *h_dx_n = new TH1D("h_dx_n","dx neutrons",1000,-10,10);
	TH1D *h_dx_p = new TH1D("h_dx_p","dx_protons",1000,-10,10);

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		
		//before adding a cut on W2
		h_dx->Fill(dx,weight);
		h_dy->Fill(dy,weight);
		h_dxdy->Fill(dy,dx,weight);
		h_W2->Fill(W2);
                
		//add a cut on W2
		if (0.0<W2 and W2<1.6){
			h_dx_W2_cut->Fill(dx,weight);
			h_dy_W2_cut->Fill(dy,weight);
			h_dxdy_W2_cut->Fill(dy,dx,weight);
			
			if (fnucl == 0){
				h_dx_n->Fill(dx,weight);
			}
			else if (fnucl == 1){
				h_dx_p->Fill(dx,weight);
			}
		}		

	        if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
		std::cout.flush();


	}
	//add lines to W2 plot to show the cut
	
	double W2_Ymax = h_W2->GetMaximum();

	TLine *line1 = new TLine(0.0,0.0,0.0,W2_Ymax);
	line1->SetLineColor(kBlue);
	line1->SetLineWidth(2);
	
	TLine *line2 = new TLine(1.6,0.0,1.6,W2_Ymax);
	line2->SetLineColor(kBlue);
	line2->SetLineWidth(2);

	//gaus fit

	TCanvas* c = new TCanvas("c","c",1600,1200);
	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	TCanvas* c2 = new TCanvas("c2","c2",1600,1200);

	c->Divide(2,2);
	c->cd(1);
	h_dx->SetXTitle("HCAL_X(exp)-HCAL_X(act)");
	h_dx->Draw("HIST");
	c->cd(2);
	h_dy->Draw("HIST");
	h_dy->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	c->cd(3);
	gStyle->SetPalette(kRainBow);
	h_dxdy->Draw("COLZ");
	h_dxdy->SetYTitle("HCAL_X(exp)-HCAL_X(act)");
	h_dxdy->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	c->cd(4);
	h_W2->Draw();
	c->Print(Form("../plots/%s_dxdyplots.pdf",printfilename));


	c1->Divide(2,2);
	c1->cd(1);
	h_dx_W2_cut->Draw("HIST");
	h_dx_W2_cut->SetXTitle("HCAL_X(exp)-HCAL_X(act)");
	c1->cd(2);
	h_dy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	h_dy_W2_cut->Draw("HIST");
	c1->cd(3);
	gStyle->SetPalette(kRainBow);
	h_dxdy_W2_cut->Draw("COLZ");
	h_dxdy_W2_cut->SetYTitle("HCAL_X(exp)-HCAL_X(act)");
	h_dxdy_W2_cut->SetXTitle("HCAL_Y(exp)-HCAL_Y(act)");
	c->cd(4);
	c1->cd(4);
	h_W2->Draw();
	line1->Draw("same");
	line2->Draw("same");
	c1->Print(Form("../plots/%s_dxdyplots_withW2Cut.pdf",printfilename));

	c2->Divide(1,1);
	c2->cd(1);
	h_dx_n->SetLineColor(kBlue);
	h_dx_p->SetLineColor(kRed);
	
	h_dx_p->Draw("HIST");
	h_dx_n->Draw("HIST SAME");

	c2->Print(Form("../plots/%s_p_n.pdf",printfilename));

}
