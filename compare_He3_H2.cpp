void compare_He3_H2(const char* filename1, const char* filename2, const char* printfilename){
	TFile* file1 = TFile::Open(filename1);
	TFile* file2 = TFile::Open(filename2);

	TTree* tree1 = (TTree*)file1->Get("Tout");
	TTree* tree2 = (TTree*)file2->Get("Tout");

	double dx_He3 = 0.0;
	double dy_He3 = 0.0;
	double W2_He3 = 0.0;
	//double Q2_He3 = 0.0;

	double dx_H2 = 0.0;
	double dy_H2 = 0.0;
	double W2_H2 = 0.0;
	//double Q2_H2 = 0.0;

	tree1->SetBranchAddress("dx",&dx_He3);
	tree1->SetBranchAddress("dy",&dy_He3);
	tree1->SetBranchAddress("W2",&W2_He3);

	tree2->SetBranchAddress("dx",&dx_H2);
	tree2->SetBranchAddress("dy",&dy_H2);
	tree2->SetBranchAddress("W2",&W2_H2);
	
	TH1D *h_dx_W2_cut_He3 = new TH1D("h_dx_W2_cut_He3","dx after W2 cut",1000,-10,10);
        TH1D *h_dy_W2_cut_He3 = new TH1D("h_dy_W2_cut_He3","dy after W2 cut",1000,-10,10); 

        TH1D *h_dx_W2_cut_H2 = new TH1D("h_dx_W2_cut_H2","dx after W2 cut",1000,-10,10);
        TH1D *h_dy_W2_cut_H2 = new TH1D("h_dy_W2_cut_H2","dy after W2 cut",1000,-10,10);


	int nentries1 = tree1->GetEntries();
	for (int i = 0; i<nentries1; i++){
		tree1->GetEntry(i);
		if (0.0<W2_He3 and W2_He3<1.6){
			h_dx_W2_cut_He3->Fill(dx_He3);
			h_dy_W2_cut_He3->Fill(dy_He3);
		}
		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries1) << "% \r";
                std::cout.flush();

	}


	int nentries2 = tree2->GetEntries();
        for (int i = 0; i<nentries2; i++){
                tree2->GetEntry(i);
                if (0.0<W2_H2 and W2_H2<1.6){
                        h_dx_W2_cut_H2->Fill(dx_H2);
                        h_dy_W2_cut_H2->Fill(dy_H2);
                }
		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries2) << "% \r";
                std::cout.flush();

        }

	TCanvas* c = new TCanvas ("c","c",1600,1200);

	h_dx_W2_cut_He3->Scale(1.0/h_dx_W2_cut_He3->Integral());
	h_dx_W2_cut_H2->Scale(1.0/h_dx_W2_cut_H2->Integral());

	h_dx_W2_cut_He3->SetLineColor(kBlue);
	h_dx_W2_cut_H2->SetLineColor(kRed);

	h_dx_W2_cut_H2->Draw("HIST");
	h_dx_W2_cut_He3->Draw("HIST SAME");

	c->Print(Form("../plots/%s_compare_He3_H2.pdf",printfilename));
}
