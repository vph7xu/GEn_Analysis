void PS_HCAL_e_plot(const char* filename, const char* printfilename){

	TFile* file = TFile::Open(filename);
        TTree* tree = (TTree*)file->Get("Tout");

	double ePS = 0.0;
	double eHCAL = 0.0;
	double vz = 0.0;
	double ntrack = 0.0;
	double nhits = 0.0;

	tree->SetBranchAddress("ePS",&ePS);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("vz",&vz);
	tree->SetBranchAddress("ntrack",&ntrack);
	tree->SetBranchAddress("nhits",&nhits);

	TH1D *h_ePS = new TH1D("h_ePS","Preshower Energy",4000,0,4);
	TH1D *h_eHCAL = new TH1D("h_eHCAL","HCAL Energy",2000,0,2);
	TH1D *h_trvz = new TH1D("h_trvz","vertex z",2000,-1,1);
	TH1D *h_ntrack = new TH1D("h_ntrack","ntrack",60,0,6);
	TH1D *h_nhits = new TH1D("h_nhits","nhits",60,0,6);

	TH1D *h_trvz_c = new TH1D("h_trvz_c","vertex z",2000,-1,1);

	int nentries = tree->GetEntries();
	for (int i = 0; i<nentries; i++){
		tree->GetEntry(i);
		if(ePS>0 and eHCAL>0){
			h_ePS->Fill(ePS);
			h_eHCAL->Fill(eHCAL);
			h_trvz->Fill(vz);
			h_ntrack->Fill(ntrack);
			h_nhits->Fill(nhits);
		}
		//if (ePS>0.2 and eHCAL>0){
		//	h_trvz_c->Fill(vz);
		//	h_eHCAL->Fill(eHCAL);
		//}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();

	}

	TCanvas* c = new TCanvas("c","c",1600,1200);
	c->Divide(3,2);
	c->cd(1);
	h_ePS->Draw();
	c->cd(2);
	h_eHCAL->Draw();
	c->cd(3);
	h_trvz->Draw();
	c->cd(4);
	h_ntrack->Draw();
	c->cd(5);
	h_nhits->Draw();
	//h_trvz_c->Draw();

	TCanvas* c1 = new TCanvas("c1","c1",1600,1200);
	c1->cd(1);
	c1->SetLogy();
	h_eHCAL->Draw();
	c1->Update();

	c->Print(Form("../plots/%s_ePS_eHCAL.pdf(",printfilename));	
	c1->Print(Form("../plots/%s_ePS_eHCAL.pdf)",printfilename));

}
