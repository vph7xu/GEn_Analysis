void simulation_plots(const char* simfilename,const char* printfilename, const char *kin){
  	TFile* sim_file = TFile::Open(simfilename);

	TTree* sim_tree = (TTree*)sim_file->Get("Tout");

	double dx = 0.0;
	double dy = 0.0;
	double W2 = 0.0;
	double Q2 = 0.0;
	double weight = 0.0;
	double fnucl = 0.0;
	double eHCAL = 0.0;
	double ePS = 0.0;
	double eSH = 0.0;
	double trP = 0.0;
	double vz = 0.0;

	sim_tree->SetBranchAddress("dx",&dx);
	sim_tree->SetBranchAddress("dy",&dy);
	sim_tree->SetBranchAddress("W2",&W2);
	sim_tree->SetBranchAddress("Q2",&Q2);
	sim_tree->SetBranchAddress("weight",&weight);
	sim_tree->SetBranchAddress("fnucl",&fnucl);
	sim_tree->SetBranchAddress("eHCAL",&eHCAL);
	sim_tree->SetBranchAddress("ePS",&ePS);
	sim_tree->SetBranchAddress("eSH",&eSH);
	sim_tree->SetBranchAddress("trP",&trP);
	sim_tree->SetBranchAddress("vz",&vz);

	TH1D *h_dx_bkg = new TH1D("h_dx_bkg","dx bkg;dx (m)",100,-4,4);
    TH1D *h_dx_bkg_n = new TH1D("h_dx_bkg_n","dx bkg (neutrons);dx (m) ",100,-4,4);
    TH1D *h_dx_bkg_p = new TH1D("h_dx_bkg_p","dx bkg (protons);dx (m)",100,-4,4);

	TH1D *h_W2_bkg = new TH1D("h_W2_bkg","W^{2} bkg;W^{2} (GeV^{2})",100,-1,5);
    TH1D *h_W2_bkg_n = new TH1D("h_W2_bkg_n","W^{2} bkg (neutrons);W^{2} (GeV^{2}) ",100,-1,5);
    TH1D *h_W2_bkg_p = new TH1D("h_W2_bkg_p","W^{2} bkg (protons);W^{2} (GeV^{2})",100,-1,5);


	int nentries = sim_tree->GetEntries();

	for (int i = 0; i<nentries; i++){
                sim_tree->GetEntry(i);

                double eoverp = (eSH+ePS)/trP;

	        // some conditions from your DB checks:
	        bool goodVz       = abs(vz) < 0.27;
	        bool goodPS       = (ePS > 0.2);
	        bool goodEHCAL    = (eHCAL > 0.025); 
	        bool goodEoverp = abs(eoverp-1)<0.2;
			bool gooddy = abs(dy)<0.5;
			bool gooddx = abs(dx)<0.5 || abs(dx+2.7)<0.5;
			bool goodW2 = -1<W2 and W2<2;


		if (/*W2_L_sim<W2 and W2<W2_H_sim and  gooddy and*/ goodEHCAL and goodVz and goodPS and goodEoverp){
                        
						if (fnucl == 0.0 || fnucl == 1.0){
							h_dx_bkg->Fill(dx,weight);
							//h_W2_bkg->Fill(W2,weight);
						}
                        if (fnucl == 0.0){
                            h_dx_bkg_n->Fill(dx,weight);
							//h_W2_bkg_n->Fill(W2,weight);
                        }
                        if (fnucl == 1.0){
                            h_dx_bkg_p->Fill(dx,weight);
							//h_W2_bkg_p->Fill(W2,weight);
                        }
        }

		if (gooddx and gooddy and goodEHCAL and goodVz and goodPS and goodEoverp){
                        
						if (fnucl == 0.0 || fnucl == 1.0){
							//h_dx_bkg->Fill(dx,weight);
							h_W2_bkg->Fill(W2,weight);
						}
                        if (fnucl == 0.0){
                            //h_dx_bkg_n->Fill(dx,weight);
							h_W2_bkg_n->Fill(W2,weight);
                        }
                        if (fnucl == 1.0){
                            //h_dx_bkg_p->Fill(dx,weight);
							h_W2_bkg_p->Fill(W2,weight);
                        }
        }

                if (i %1000 == 0 ) std::cout<< (i * 100.0/ nentries)<<"sim file" << "% \r";
                std::cout.flush();
	}	

	double dx_integral = h_dx_bkg->Integral();
	double W2_integral = h_W2_bkg->Integral();

	TCanvas *C = new TCanvas("C","C",2000,800);

	// Legend for dx
	auto leg_dx = new TLegend(0.60, 0.70, 0.88, 0.88);
	leg_dx->SetBorderSize(0);
	leg_dx->SetFillStyle(0);
	leg_dx->SetTextSize(0.03);
	leg_dx->AddEntry(h_dx_bkg,   "Total background", "l");
	leg_dx->AddEntry(h_dx_bkg_n, "Neutron background", "l");
	leg_dx->AddEntry(h_dx_bkg_p, "Proton background", "l");
	

	// Legend for W2
	auto leg_W2 = new TLegend(0.15, 0.70, 0.43, 0.88);
	leg_W2->SetBorderSize(0);
	leg_W2->SetFillStyle(0);
	leg_W2->SetTextSize(0.03);
	leg_W2->AddEntry(h_W2_bkg,   "Total background", "l");
	leg_W2->AddEntry(h_W2_bkg_n, "Neutron background", "l");
	leg_W2->AddEntry(h_W2_bkg_p, "Proton background", "l");

	C->Divide(2,1);
	gStyle->SetOptStat(0);
	
	C->cd(1);
	h_dx_bkg->SetLineWidth(3);
	h_dx_bkg_n->SetLineWidth(3);
	h_dx_bkg_p->SetLineWidth(3);

	h_dx_bkg->SetLineColor(kBlack);
	h_dx_bkg_n->SetLineColor(kRed);
	h_dx_bkg_p->SetLineColor(kBlue);

	h_dx_bkg->Scale(1/dx_integral);
	h_dx_bkg_n->Scale(1/dx_integral);
	h_dx_bkg_p->Scale(1/dx_integral);

	h_dx_bkg->Draw();
	h_dx_bkg_n->Draw("same");
	h_dx_bkg_p->Draw("same");

	leg_dx->Draw("same");

	C->cd(2);
	h_W2_bkg->SetLineWidth(3);
	h_W2_bkg_n->SetLineWidth(3);
	h_W2_bkg_p->SetLineWidth(3);

	h_W2_bkg->SetLineColor(kBlack);
	h_W2_bkg_n->SetLineColor(kRed);
	h_W2_bkg_p->SetLineColor(kBlue);

	h_W2_bkg->Scale(1/W2_integral);
	h_W2_bkg_n->Scale(1/W2_integral);
	h_W2_bkg_p->Scale(1/W2_integral);

	h_W2_bkg->Draw();
	h_W2_bkg_n->Draw("same");
	h_W2_bkg_p->Draw("same");

	leg_W2->Draw("same");

	C->SaveAs(Form("plots/inelastic_dx_%s.png",kin));

} 
