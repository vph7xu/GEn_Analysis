void simulation_plots(const char* sim_filename){
	TFile* sim_file = TFile::Open(sim_filename);

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

	//TH1D *h_dx_bkg = new TH1D("h_dx_bkg","dx bkg",100,-4,3);
        //TH1D *h_dx_bkg_n = new TH1D("h_dx_bkg_n","dx bkg (neutrons)",100,-4,3);
        //TH1D *h_dx_bkg_p = new TH1D("h_dx_bkg_p","dx bkg (protons)",100,-4,3);

	//TH1D *h_weight = new TH1D();

	int nentries = sim_tree->GetEntries();

	double sum_weights = 0.0;
	double sum_weights_1 =0.0;

	for (int i = 0; i<nentries; i++){
                sim_tree->GetEntry(i);

                double eoverp = (eSH+ePS)/trP;

	        // some conditions from your DB checks:
	        bool goodVz       = abs(vz) < 0.27;
	        bool goodPS       = (ePS > 0.2);
	        bool goodEHCAL    = (eHCAL > 0.125); 
	        bool goodEoverp = abs(eoverp-1)<0.2;
		bool goodW2 = W2<1.5;
		bool gooddx = abs(dx)<0.5 || abs(dx+1.73)<0.5;
		bool gooddy = abs(dy)<0.5;

		sum_weights_1 += weight;

		if (goodVz && goodPS && goodEHCAL && goodEoverp && goodW2 && gooddx && gooddy){
			sum_weights += weight;
		}	

                if (i %1000 == 0 ) std::cout<< (i * 100.0/ nentries)<<"sim file" << "% \r";
                std::cout.flush();


	}
	cout<<"number of entries : "<<nentries<<endl;
	cout<<"sum of weights 1 : "<<sum_weights_1<<endl;
	cout<<"sum of weights : "<<sum_weights<<endl;

}
