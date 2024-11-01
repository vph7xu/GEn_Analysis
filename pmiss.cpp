#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "parse.h"
#include "plotdxdy.h"
#include <cmath>

void pmiss(const char* filename, const char* printfilename, const char *kin){
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

    double coin_time_L = getDoubleValue(config,"coin_time_L");
    double coin_time_H = getDoubleValue(config,"coin_time_H");

    double W2_L = getDoubleValue(config,"W2_L");
    double W2_H = getDoubleValue(config,"W2_H");

    double dy_L = getDoubleValue(config,"dy_L");
    double dy_H = getDoubleValue(config,"dy_H");

    double dx_L = getDoubleValue(config,"dx_L");
    double dx_H = getDoubleValue(config,"dx_H");

    double IHWP_flip = getDoubleValue(config,"IHWP_flip");

    double run_num_L = getDoubleValue(config,"run_num_L");
    double run_num_H = getDoubleValue(config,"run_num_H");

    double N_cut_cx = getDoubleValue(config,"N_cut_cx");
    double N_cut_cy = getDoubleValue(config,"N_cut_cy");
    double N_cut_rx = getDoubleValue(config,"N_cut_rx");
    double N_cut_ry = getDoubleValue(config,"N_cut_ry");


    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");
	
	double theta_pq = 0.0;
	double dx = 0.0;
	double dy = 0.0;
	double Q2 = 0.0;
	double pN_expect = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;
	double eHCAL = 0.0;
	double nblk_HCAL = 0.0;
	double hcal_clus_e[1000];
	double hcal_clus_id[1000];
	double hcal_clus_mem_e[1000];
	double hcal_clus_mem_id[1000];
	double trP_sbs = 0.0;
	double ntrack_sbs = 0.0;
	double vz = 0.0;

	tree->SetBranchAddress("theta_pq",&theta_pq);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("pN_expect",&pN_expect);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("hcal_clus_e",&hcal_clus_e);
	tree->SetBranchAddress("hcal_clus_mem_e",&hcal_clus_mem_e);
	tree->SetBranchAddress("hcal_clus_id",&hcal_clus_id);
	tree->SetBranchAddress("hcal_clus_mem_id",&hcal_clus_mem_id);
	tree->SetBranchAddress("nblk_HCAL",&nblk_HCAL);
	tree->SetBranchAddress("trP_sbs",&trP_sbs);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
	tree->SetBranchAddress("vz",&vz);

	//cut histograms
	TH1D *h_coin_time = new TH1D("h_coin_time","cointime (no cuts)", 1000,40,140);
	TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,4);
	TH2D *h_dxdy_nocuts = new TH2D("h_dxdy_nocuts","dxdy (no cuts)",200,-4,4,200,-4,4); 
	TH2D *h_dxdy = new TH2D("h_dxdy","dxdy (cointime cut included)",200,-4,4,200,-4,4);
	TH2D *h_dxdy_cut_cointime_W2 = new TH2D("h_dxdy_cut_cointime_W2","dxdy with cointime and W2 cuts",200,-4,4,200,-4,4);

	//cutting on Ppar
	TH1D *h_W2_cut_Ppar = new TH1D("h_W2_cut_Ppar","W2 distribution only with Ppar cut", 100,-2,4);
	TH1D *h_W2_cut_Cointime_Ppar = new TH1D("h_W2_cut_Cointime_Ppar","W2 distribution with Cointime and Ppar cut", 100,-2,4);
	TH1D *h_W2_cut_QE = new TH1D("h_W2_cut_QE","W2 distribution with coin, dx, and dy", 100,-2,4);
	TH1D *h_W2_cut_QEnPpar = new TH1D("h_W2_cut_QEnPpar","W2 distribution with coin, dx, dy and Ppar cuts", 100,-2,4);

	TH1D *h_coin_time_cut_Ppar = new TH1D("h_coin_time_cut_Ppar","cointime distribution only with Ppar cut",1000, 40, 140);
	TH2D *h_dxdy_cut_Ppar = new TH2D("h_dxdy_cut_Ppar","dxdy only with Ppar cut",200,-4,4,200,-4,4);
	TH2D *h_dxdy_cut_Ppar_cointime = new TH2D("h_dxdy_cut_Ppar_cointime","dxdy Ppar and cointime cut ",200,-4,4,200,-4,4);
	TH2D *h_dxdy_cut_Ppar_cointime_W2 = new TH2D("h_dxdy_cut_Ppar_cointime_W2","dxdy Ppar, cointime and W2 cut ",200,-4,4,200,-4,4);

	TH1D *h_dx_cut_coin_time_W2_dy = new TH1D("h_dx_cut_coin_time_W2","dx with cointime, W2, and dy cuts", 200, -4, 4);
	TH1D *h_dx_cut_Ppar_coin_time_W2_dy = new TH1D("h_dx_cut_Ppar_coin_time_W2","dx with Ppar, cointime, W2, and dy cuts", 200, -4, 4);

	//no cuts histogram
	TH2D *h_Ppar_coin_time = new TH2D("h_Ppar_coin_time","Ppar v cointime distribution",1000,40,140,200,-5,5);

	//cutting on eHCAL
	TH1D *h_W2_cut_eHCAL = new TH1D("h_W2_cut_eHCAL","W2 distribution only with eHCAL cut", 100,-2,4);
	TH1D *h_coin_time_cut_eHCAL = new TH1D("h_coin_time_cut_eHCAL","cointime distribution only with eHCAL cut",1000, 40, 140);
	TH2D *h_dxdy_cut_eHCAL = new TH2D("h_dxdy_cut_eHCAL","dxdy only with eHCAL cut",200,-4,4,200,-4,4);

	//KE
	TH1D *h_W2_cut_KE = new TH1D("h_W2_cut_KE","W2 distribution only with KE cut", 100,-2,4);
	TH1D *h_coin_time_cut_KE = new TH1D("h_coin_time_cut_KE","cointime distribution only with KE cut",1000, 40, 140);
	TH2D *h_dxdy_cut_KE = new TH2D("h_dxdy_cut_KE","dxdy only with KE cut",200,-4,4,200,-4,4);

	//correlation plot
	TH2D *h_Ppar_pNexpect_cut_cointime = new TH2D("h_Ppar_pNexpect_cut_cointime","Ppar v |q| only with coin cut ", 200, 0, 4, 200,-5, 5);
	TH2D *h_Ppar_pNexpect_cut_cointime_dxdy = new TH2D("h_Ppar_pNexpect_cut_cointime_dxdy","Ppar v |q| with coin, dx and dy cuts", 200, 0, 4, 200,-5, 5);
	TH2D *h_Ppar_vz_cut_cointime = new TH2D("h_Ppar_vz_cut_cointime","Ppar v vz only with coin cut", 200,-0.3,0.3,200,-5,5);
	TH2D *h_Ppar_vz_cut_cointime_dxdy = new TH2D("h_Ppar_vz_cut_cointime_dxdy","Ppar v vz with coin, dx and dy cuts", 200,-0.3,0.3,200,-5,5);
	TH2D *h_Ppar_eHCAL_cut_cointime = new TH2D("h_Ppar_eHCAL_cut_cointimee","Ppar v eHCAL (KinE) only with coin cut",200,0,10,200,-5,5);
	TH2D *h_Ppar_eHCAL_cut_cointime_dxdy = new TH2D("h_Ppar_eHCAL_cut_cointime_dxdy","Ppar v eHCAL (KinE) with coin, dx and dy cuts",200,0,10,200,-5,5);
	TH2D *h_Ppar_thetapq_cut_cointime = new TH2D("h_Ppar_thetapq_cut_cointime","Ppar v theta_pq only with coin cut",200,-0.1,0.5,200,-5,5);
	TH2D *h_Ppar_thetapq_cut_cointime_dxdy = new TH2D("h_Ppar_thetapq_cut_cointime_dxdy","Ppar v theta_pq with coin, dx and dy cut",200,-0.1,0.1,200,-5,5);	

	//Pmiss correlation plot
	TH2D *h_Pmiss_pNexpect_cut_cointime = new TH2D("h_Pmiss_pNexpect_cut_cointime","Pmiss v |q| only with coin cut ", 200, 0, 4, 200,0, 6);
	TH2D *h_Pmiss_pNexpect_cut_cointime_dxdy = new TH2D("h_Pmiss_pNexpect_cut_cointime_dxdy","Pmiss v |q| with coin, dx and dy cuts", 200, 0, 4, 200, 0, 6);
	TH2D *h_Pmiss_vz_cut_cointime = new TH2D("h_Pmiss_vz_cut_cointime","Pmiss v vz only with coin cut", 200,-0.3,0.3,200,0,6);
	TH2D *h_Pmiss_vz_cut_cointime_dxdy = new TH2D("h_Pmiss_vz_cut_cointime_dxdy","Pmiss v vz with coin, dx and dy cuts", 200,-0.3,0.3,200,0,6);
	TH2D *h_Pmiss_eHCAL_cut_cointime = new TH2D("h_Pmiss_eHCAL_cut_cointimee","Pmiss v eHCAL (KinE) only with coin cut",200,0,10,200,0,6);
	TH2D *h_Pmiss_eHCAL_cut_cointime_dxdy = new TH2D("h_Pmiss_eHCAL_cut_cointime_dxdy","Pmiss v eHCAL (KinE) with coin, dx and dy cuts",200,0,10,200,0,6);
	TH2D *h_Pmiss_thetapq_cut_cointime = new TH2D("h_Pmiss_thetapq_cut_cointime","Pmiss v theta_pq only with coin cut",200,-0.1,0.5,200,0,6);
	TH2D *h_Pmiss_thetapq_cut_cointime_dxdy = new TH2D("h_Pmiss_thetapq_cut_cointime_dxdy","Pmiss v theta_pq with coin, dx and dy cut",200,-0.1,0.1,200,0,6);	

	//pmiss histograms
	TH1D *h_pmiss = new TH1D("h_pmiss","Pmiss distribution",200, 0, 6);
	TH1D *h_pmiss_cut_coin = new TH1D("h_pmiss_cut_coin","Pmiss distribution only with cointime cut", 200, 0, 6);
	TH1D *h_pmiss_cut_coin_dxdy = new TH1D("h_pmiss_cut_coin_dxdy","Pmiss distribution with cointime, dx and dy cuts", 200, 0, 6);

	TH1D *h_Ppar_cut_hiW2 = new TH1D("h_Ppar_cut_hiW2","Ppar distribution for W2<0", 200 , -5, 5);
	TH1D *h_Ppar_cut_loW2 = new TH1D("h_Ppar_cut_loW2","Ppar distribution for W2>=0", 200 , -5, 5);

	TH1D *h_theta_pq = new TH1D("h_theta_pq","theta_pq only with cointime cut",200,0,2 );
	TH1D *h_Pperp = new TH1D("h_Pperp","Pperp only with cointime cut",200,0,2);
	TH1D *h_realPperp = new TH1D("h_realPperp","realPperp only with cointime cut",200,0,2);
	TH1D *h_realPpar = new TH1D("h_realPpar","realPpar only with cointime cut",200,-5,5);

	TH2D *h_dx_Pperp = new TH2D("h_dx_Pperp","dx v Pperp dist only with cointime cut",100,0,1,250,-3,3);
	TH2D *h_dx_realPperp = new TH2D("h_dx_realPperp","dx v realPperp dist only with cointime cut",100,0,1,250,-3,3);
	TH2D *h_dx_realPpar = new TH2D("h_dx_realPpar","dx v realPpar dist only with cointime cut",100,-5,5,250,-3,3);


	TH2D *h_Pperp_W2 = new TH2D("h_Pperp_W2","Pperp v W2 dist only with cointime cut",1200,0,3,100,0,1);
	TH2D *h_realPperp_W2 = new TH2D("h_realPperp_W2","realPperp v W2 dist only with cointime cut",1200,0,3,100,0,1);
	TH2D *h_realPpar_W2 = new TH2D("h_realPpar_W2","realPpar v W2 dist only with cointime cut",1200,0,3,100,-5,5);

	TH1D *h_theta_pq_cut_coin_time = new TH1D("h_theta_pq_cut_coin_time","theta_pq only with coin cut",100,0,1);
    TH1D *h_Pperp_cut_coin_time = new TH1D("h_Pperp_cut_coin_time","Pperp only with coin cut",100,0,1);
    TH1D *h_realPperp_cut_coin_time = new TH1D("h_realPperp_cut_coin_time","realPperp only with coin cut",200,-0.5,2);
    TH1D *h_realPpar_cut_coin_time = new TH1D("h_realPpar_cut_coin_time","realPpar only with coin cut",200,-5,5);
 
    TH2D *h_dx_Pperp_cut_coin_time = new TH2D("h_dx_Pperp_cut_coin_time","dx v Pperp dist only with coin cut",100,0,1,250,-4,3);
    TH2D *h_dx_realPperp_cut_coin_time = new TH2D("h_dx_realPperp_cut_coin_time","dx v realPperp dist only with coin cut",100,0,1,250,-4,3);
    TH2D *h_dx_realPpar_cut_coin_time = new TH2D("h_dx_realPpar_cut_coin_time","dx v realPpar dist only with coin cut",100,0,1,250,-4,3);	
	TH2D *h_Pperp_W2_cut_coin_time = new TH2D("h_Pperp_W2_cut_coin_time","Pperp v W2 dist only with coin cut",1200,0,3,100,0,1);
	TH2D *h_realPperp_W2_cut_coin_time = new TH2D("h_realPperp_W2_cut_coin_time","realPperp v W2 dist only with coin cut",1200,-2,4,200,-0.5,2);
	TH2D *h_realPpar_W2_cut_coin_time = new TH2D("h_realPpar_W2_cut_coin_time","realPpar v W2 dist only with coin cut",1200,-2,4,200,-5,5);

	TH1D *h_theta_pq_cut_coin_time_and_dx = new TH1D("h_theta_pq_cut_coin_time_and_dx","theta_pq coin, dx and dy cuts",100,0,1);
    TH1D *h_Pperp_cut_coin_time_and_dx = new TH1D("h_Pperp_cut_coin_time_and_dx","Pperp coin, dx and dy cuts",100,0,1);
    TH1D *h_realPperp_cut_coin_time_and_dx = new TH1D("h_realPperp_cut_coin_time_and_dx","realPperp coin, dx and dy cuts",200,-0.5,0.5);
    TH1D *h_realPpar_cut_coin_time_and_dx = new TH1D("h_realPpar_cut_coin_time_and_dx","realPpar coin, dx and dy cuts",200,-5,5);
 
    TH2D *h_dx_Pperp_cut_coin_time_and_dx = new TH2D("h_dx_Pperp_cut_coin_time_and_dx","dx v Pperp dist coin, dx and dy cuts",100,0,1,250,-4,3);
    TH2D *h_dx_realPperp_cut_coin_time_and_dx = new TH2D("h_dx_realPperp_cut_coin_time_and_dx","dx v realPperp dist coin, dx and dy cuts",100,0,1,250,-4,3);
    TH2D *h_dx_realPpar_cut_coin_time_and_dx = new TH2D("h_dx_realPpar_cut_coin_time_and_dx","dx v realPpar dist coin, dx and dy cuts",100,0,1,250,-4,3);	
	TH2D *h_Pperp_W2_cut_coin_time_and_dx = new TH2D("h_Pperp_W2_cut_coin_time_and_dx","Pperp v W2 dist coin, dx and dy cuts",1200,0,3,100,0,1);
	TH2D *h_realPperp_W2_cut_coin_time_and_dx = new TH2D("h_realPperp_W2_cut_coin_time_and_dx","realPperp v W2 dist coin, dx and dy cuts",1200,-2,4,200,-0.5,0.5);
	TH2D *h_realPpar_W2_cut_coin_time_and_dx = new TH2D("h_realPpar_W2_cut_coin_time_and_dx","realPpar v W2 dist coin, dx and dy cuts ",1200,-2,4,200,-5,5);

	TH1D *h_theta_pq_cut_QE = new TH1D("h_theta_pq_cut_QE","theta_pq with QE cuts",30,0,0.3);
	TH1D *h_Pperp_cut_QE = new TH1D("h_Pperp_cut_QE","Pperp with QE cuts",40,0,0.4);
	TH1D *h_realPperp_cut_QE = new TH1D("h_realPperp_cut_QE","realPperp with QE cuts",40,0,0.4);
	TH1D *h_realPpar_cut_QE = new TH1D("h_realPpar_cut_QE","realPpar with QE cuts",40,-5,5);
    TH2D *h_dx_Pperp_cut_QE = new TH2D("h_dx_Pperp_cut_QE","dx v Pperp dist with QE cuts",100,0,1,250,-4,3);
    TH2D *h_dx_realPperp_cut_QE = new TH2D("h_dx_realPperp_cut_QE","dx v realPperp dist with QE cuts",100,0,1,250,-4,3);
    TH2D *h_dx_realPpar_cut_QE = new TH2D("h_dx_realPpar_cut_QE","dx v realPpar dist with QE cuts",100,0,1,250,-4,3);
    TH2D *h_Pperp_W2_cut_QE = new TH2D("h_Pperp_W2_cut_QE","Pperp v W2 dist with QE cuts",1200,0,W2_H,100,0,1);
    TH2D *h_realPperp_W2_cut_QE = new TH2D("h_realPperp_W2_cut_QE","real Pperp v W2 dist with QE cuts",1200,0,W2_H,100,0,1);
    TH2D *h_realPpar_W2_cut_QE = new TH2D("h_realPpar_W2_cut_QE","real Ppar v W2 dist with QE cuts",1200,0,W2_H,100,-5,5);
	TH2D *h_dx_W2_cut_coin_time = new TH2D("h_dx_W2_cut_coin_time","dx v W2 dist with coin time cut", 100,W2_L,3,250,-4,3);

	TH1D *h_eHCAL = new TH1D("h_eHCAL","eHCAL distr",200,0,1);
	TH2D *h_Ppar_eHCAL = new TH2D("h_Ppar_eHCAL","Ppar vs eHCAL distribution (no cuts)",200,0,1,200,-5,4);

	TH1D *h_kinE = new TH1D("h_kinE","kinetic energy", 200,-1,6);
	TH1D *h_kinE_cut_coin_time = new TH1D("h_kinE_cut_coin_time","kinetic energy with cointime cut",200,-1,10);
	TH1D *h_kinE_cut_W2 = new TH1D("h_kinE_cut_W2","kintic energy with cointime and W2 cuts",200,-1,10);
	TH1D *h_kinE_cut_QE = new TH1D("h_kinE_cut_QE","kinetic energy with QE cuts",200,-1,10);

	int nentries = tree->GetEntries();

	for (int i = 0; i<nentries; i++){	
		tree->GetEntry(i);
		//if (ntrack_sbs>0){

		double KinE = 0.0;

		for (int j = 0; j<nblk_HCAL; j++){
			//KinE += hcal_clus_e[j]/sampling_fractions_each_blk[hcal_clus_id[j]];
			KinE += hcal_clus_mem_e[j]/sampling_fractions_each_blk[hcal_clus_mem_id[j]];

			//if (i%1000 == 0){ 
			//	std::cout<<"event : "<<i<<"blk id: "<<hcal_clus_mem_id[j]<<" blk energy : "<<hcal_clus_mem_e[j]<<" blk sf : " <<sampling_fractions_each_blk[hcal_clus_mem_id[j]]<<endl;
			//}
		}
		//if (i%1000 == 0){ 
		//	std::cout<<"event : "<<i<<" kinE : "<<KinE<<endl;
		//}
		double Pperp = theta_pq * pN_expect;
		double realPperp = (theta_pq * sqrt(-pow(0.938,2)+pow((KinE+0.938),2))); // assuming mN and mP is equal
		double realPpar = (cos(theta_pq) * sqrt(-pow(0.938,2)+pow((KinE+0.938),2)))-pN_expect;
		double Pmiss = sqrt(pow(realPperp,2)+pow(realPpar,2)); 
		
		//Filling cut histograms
		h_coin_time->Fill(coin_time);
		h_W2->Fill(W2);
		h_dxdy_nocuts->Fill(dy,dx);

		//Filling Ppar vs cointime
		h_Ppar_coin_time->Fill(coin_time,realPpar);


		h_pmiss->Fill(Pmiss);
		h_kinE->Fill(KinE);
		h_eHCAL->Fill(eHCAL);
		h_Ppar_eHCAL->Fill(eHCAL,realPpar);


		//kinE cut
		if (KinE>2.5){
			h_coin_time_cut_KE->Fill(coin_time);
			h_dxdy_cut_KE->Fill(dy,dx);
		}

		//eHCAL cut
		if (eHCAL>0.3){
			h_coin_time_cut_eHCAL->Fill(coin_time);
			h_dxdy_cut_eHCAL->Fill(dy,dx);
		}

		//applying a W2<0 cut only to see Ppar dist
		if (W2<0){
			h_Ppar_cut_hiW2->Fill(realPpar);
		}
		else{
			h_Ppar_cut_loW2->Fill(realPpar);

		}

		 //applying a Ppar cut
        if (abs(realPpar)<1.5){
			h_W2_cut_Ppar->Fill(W2);
			h_coin_time_cut_Ppar->Fill(coin_time);
			h_dxdy_cut_Ppar->Fill(dy,dx);
		}

		//applying cointime cut
		if(coin_time_L<coin_time and coin_time<coin_time_H){
			h_pmiss_cut_coin->Fill(Pmiss);

			//Filling dxdy histogram to illustrate the cut
			h_dxdy->Fill(dy,dx);

			//Filling Ppar correlation plots
			h_Ppar_pNexpect_cut_cointime->Fill(pN_expect,realPpar);
			h_Ppar_vz_cut_cointime->Fill(vz,realPpar);
			h_Ppar_eHCAL_cut_cointime->Fill(KinE,realPpar);
			h_Ppar_thetapq_cut_cointime->Fill(theta_pq,realPpar);

			//Filling Pmiss correlation plots
			h_Pmiss_pNexpect_cut_cointime->Fill(pN_expect,Pmiss);
			h_Pmiss_vz_cut_cointime->Fill(vz,Pmiss);
			h_Pmiss_eHCAL_cut_cointime->Fill(KinE,Pmiss);
			h_Pmiss_thetapq_cut_cointime->Fill(theta_pq,Pmiss);

			//Filling Pmiss histograms
			h_theta_pq->Fill(theta_pq);
			h_Pperp ->Fill(Pperp);
			h_realPperp->Fill(realPperp);
			h_realPpar->Fill(realPpar);

			h_dx_Pperp->Fill(Pperp,dx);
			h_Pperp_W2->Fill(W2,Pperp);

			h_dx_realPperp->Fill(realPperp,dx);
			h_realPperp_W2->Fill(W2,realPperp);
		
			h_dx_realPpar->Fill(realPpar,dx);
			h_realPpar_W2->Fill(W2,realPpar);	

			//
			h_theta_pq_cut_coin_time->Fill(theta_pq);
            h_Pperp_cut_coin_time->Fill(Pperp);
            h_dx_Pperp_cut_coin_time->Fill(Pperp,dx);
            h_Pperp_W2_cut_coin_time->Fill(W2,Pperp);

            h_realPperp_cut_coin_time->Fill(realPperp);
            h_dx_realPperp_cut_coin_time->Fill(realPperp,dx);
            h_realPperp_W2_cut_coin_time->Fill(W2,realPperp);

            h_realPpar_cut_coin_time->Fill(realPpar);
            h_dx_realPpar_cut_coin_time->Fill(realPpar,dx);
            h_realPpar_W2_cut_coin_time->Fill(W2,realPpar);

			h_dx_W2_cut_coin_time->Fill(W2,dx);
			h_kinE_cut_coin_time->Fill(KinE);

			if(W2<2){
				h_dxdy_cut_cointime_W2->Fill(dy,dx);
				if(dy_L<dy and dy<dy_H){
					h_dx_cut_coin_time_W2_dy->Fill(dx);
				}
			}


			if (abs(realPpar)<1.5){
				h_W2_cut_Cointime_Ppar->Fill(W2);
				h_dxdy_cut_Ppar_cointime->Fill(dy,dx);
				if(W2<2){
					h_dxdy_cut_Ppar_cointime_W2->Fill(dy,dx);
					if(dy_L<dy and dy<dy_H){
						h_dx_cut_Ppar_coin_time_W2_dy->Fill(dx);
					}
				}
			}
		}

		//applying cointime, dx and dy cuts
		if(coin_time_L<coin_time and coin_time<coin_time_H and dx_L<dx and dx<dx_H and dy_L<dy and dy<dy_H){
			
			h_pmiss_cut_coin_dxdy->Fill(Pmiss);

			//Filling Ppar correlation plots
			h_Ppar_pNexpect_cut_cointime_dxdy->Fill(pN_expect,realPpar);
			h_Ppar_vz_cut_cointime_dxdy->Fill(vz,realPpar);
			h_Ppar_eHCAL_cut_cointime_dxdy->Fill(KinE,realPpar);
			h_Ppar_thetapq_cut_cointime_dxdy->Fill(theta_pq,realPpar);

			//Filling Pmiss correlation plots
			h_Pmiss_pNexpect_cut_cointime_dxdy->Fill(pN_expect,Pmiss);
			h_Pmiss_vz_cut_cointime_dxdy->Fill(vz,Pmiss);
			h_Pmiss_eHCAL_cut_cointime_dxdy->Fill(KinE,Pmiss);
			h_Pmiss_thetapq_cut_cointime_dxdy->Fill(theta_pq,Pmiss);

			//

			h_theta_pq_cut_coin_time_and_dx->Fill(theta_pq);
            h_Pperp_cut_coin_time_and_dx->Fill(Pperp);
            h_dx_Pperp_cut_coin_time_and_dx->Fill(Pperp,dx);
            h_Pperp_W2_cut_coin_time_and_dx->Fill(W2,Pperp);

            h_realPperp_cut_coin_time_and_dx->Fill(realPperp);
            h_dx_realPperp_cut_coin_time_and_dx->Fill(realPperp,dx);
            h_realPperp_W2_cut_coin_time_and_dx->Fill(W2,realPperp);

            h_realPpar_cut_coin_time_and_dx->Fill(realPpar);
            h_dx_realPpar_cut_coin_time_and_dx->Fill(realPpar,dx);
            h_realPpar_W2_cut_coin_time_and_dx->Fill(W2,realPpar);

            h_W2_cut_QE->Fill(W2);
            //applying a Ppar cut
            if (abs(realPpar)<1.5){
            	h_W2_cut_QEnPpar->Fill(W2);
        	}
            //applying W2 cut
			if(W2_L<W2 and W2<W2_H){

				h_theta_pq_cut_QE->Fill(theta_pq);
				h_Pperp_cut_QE->Fill(Pperp);
				h_dx_Pperp_cut_QE->Fill(Pperp,dx);
				h_Pperp_W2_cut_QE->Fill(W2,Pperp);

				h_realPperp_cut_QE->Fill(realPperp);
				h_dx_realPperp_cut_QE->Fill(realPperp,dx);
				h_realPperp_W2_cut_QE->Fill(W2,realPperp);

				h_realPpar_cut_QE->Fill(realPpar);
				h_dx_realPpar_cut_QE->Fill(realPpar,dx);
				h_realPpar_W2_cut_QE->Fill(W2,realPpar);

				h_kinE_cut_QE->Fill(KinE);								
				//h_dx_W2_cut_coin_time->Fill(W2,dx);
				}
		}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();	
	}

	//square cut
	//TCutG *cutsq = CreateSquareCut(dy_L,dx_L,dy_H,dx_H);
	//TCutG *cutsq1 = CreateSquareCut(-10,dx_L,10,dx_H);

	TCanvas* ccuts = new TCanvas("ccuts","ccuts",3200,3600);
	TCanvas* c = new TCanvas("c","c",3200,3600);
	TCanvas* c1 = new TCanvas("c1","c1",3200,3600);
	TCanvas* c2 = new TCanvas("c2","c2",3200,3600);
	TCanvas* c3 = new TCanvas("c3","c3",3200,3600);
	TCanvas* c4 = new TCanvas("c4","c4",3200,3600);
	TCanvas* c5 = new TCanvas("c5","c5",3200,3600);
	TCanvas* c51 = new TCanvas("c51","c51",3200,3600);
	TCanvas* c6 = new TCanvas("c6","c6",3200,3600);
	TCanvas* c7 = new TCanvas("c7","c7",3200,3600);
	TCanvas* c8 = new TCanvas("c8","c8",3200,3600);
	TCanvas* c81 = new TCanvas("c81","c81",3200,3600);
	TCanvas* c9 = new TCanvas("c9","c9",3200,3600);
	TCanvas* c10 = new TCanvas("c10","c10",3200,3600);
	TCanvas* c11 = new TCanvas("c11","c11",3200,3600);
	TCanvas* c12 = new TCanvas("c12","c12",3200,3600);
	TCanvas* c13 = new TCanvas("c13","c13",3200,3600);
	TCanvas* c14 = new TCanvas("c14","c14",3200,3600);
	TCanvas* c15 = new TCanvas("c15","c15",3200,3600);
	TCanvas* c16 = new TCanvas("c16","c16",3200,3600);
	TCanvas* c17 = new TCanvas("c17","c17",3200,3600);
	TCanvas* c18 = new TCanvas("c18","c18",3200,3600);
	TCanvas* c19 = new TCanvas("c19","c19",3200,3600);
	TCanvas* c20 = new TCanvas("c20","c20",3200,3600);
	TCanvas* c21 = new TCanvas("c21","c21",3200,3600);
	TCanvas* c22 = new TCanvas("c22","c22",3200,3600);

	//c
    c->Divide(2,2);
	c->cd(1);
	h_theta_pq->Draw();
	c->cd(2);
	h_Pperp->Draw();
	c->cd(3);
	h_dx_Pperp->Draw("COLZ");
	c->cd(4);
	h_Pperp_W2->Draw("COLZ");


	//ccuts
	ccuts->Divide(2,2);
	ccuts->cd(1);
	h_coin_time->Draw();
	TLine *line01 = new TLine(coin_time_L, 0, coin_time_L, 2750);
    TLine *line02 = new TLine(coin_time_H, 0, coin_time_H, 2750);
    line01->SetLineColor(kRed);
    line02->SetLineColor(kRed);
    line01->SetLineWidth(2);
    line02->SetLineWidth(2);
    line01->Draw("same");
    line02->Draw("same");
	ccuts->cd(2);
	h_W2->Draw();
	ccuts->cd(3);
	h_dxdy->Draw("COLZ");
	//cutsq->Draw("same");

	//c1
	c1->Divide(2,2);
	c1->cd(1);
	h_theta_pq_cut_QE->SetXTitle("theta_pq (rad)");
	h_theta_pq_cut_QE->Draw();
    c1->cd(2);
	h_Pperp_cut_QE->SetXTitle("Pseudo Pmiss perp (GeV)");
    h_Pperp_cut_QE->Draw();
    c1->cd(3);
	h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
	h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c1->cd(4);
	h_Pperp_W2_cut_QE->SetXTitle("W2 (GeV)");
	h_Pperp_W2_cut_QE->SetYTitle("Pseudo Pmiss perp (GeV)");
    h_Pperp_W2_cut_QE->Draw("COLZ");

    //c2
    c2->Divide(2,2);
    c2->cd(1);
    h_theta_pq_cut_coin_time_and_dx->SetXTitle("theta_pq (rad)");
    h_theta_pq_cut_coin_time_and_dx->Draw();
    c2->cd(2);
    h_Pperp_cut_coin_time_and_dx->SetXTitle("Pseudo Pmiss perp (GeV)");
    h_Pperp_cut_coin_time_and_dx->Draw();
    c2->cd(3);
    h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
    h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c2->cd(4);
    h_Pperp_W2_cut_coin_time_and_dx->SetXTitle("W2 (GeV)");
    h_Pperp_W2_cut_coin_time_and_dx->SetYTitle("Pseudo Pmiss perp (GeV)");
    h_Pperp_W2_cut_coin_time_and_dx->Draw("COLZ");

    //c3
	c3->Divide(2,2);
	c3->cd(1);
	h_realPperp->SetXTitle("Pmiss perp (GeV)");
	h_realPperp->Draw();
	c3->cd(2);
	h_dx_realPperp->Draw("COLZ");
	c3->cd(3);
	h_realPperp_W2->Draw("COLZ");

	//c4
	c4->Divide(2,2);
	c4->cd(1);
	h_theta_pq_cut_QE->SetXTitle("theta_pq (rad)");
	h_theta_pq_cut_QE->Draw();
    c4->cd(2);
	h_realPperp_cut_QE->SetXTitle("Pmiss perp (GeV)");
    h_realPperp_cut_QE->Draw();
    c4->cd(3);
	h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
	h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c4->cd(4);
	h_realPperp_W2_cut_QE->SetXTitle("W2 (GeV)");
	h_realPperp_W2_cut_QE->SetYTitle("Pmiss perp (GeV)");
    h_realPperp_W2_cut_QE->Draw("COLZ");

    //c5
    c5->Divide(2,2);
    c5->cd(1);
    h_theta_pq_cut_coin_time_and_dx->SetXTitle("theta_pq (rad)");
    h_theta_pq_cut_coin_time_and_dx->Draw();
    c5->cd(3);
    h_realPperp_cut_coin_time_and_dx->SetXTitle("Pmiss perp (GeV)");
    h_realPperp_cut_coin_time_and_dx->Draw();
    c5->cd(2);
    h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
    h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c5->cd(4);
    h_realPperp_W2_cut_coin_time_and_dx->SetXTitle("W2 (GeV)");
    h_realPperp_W2_cut_coin_time_and_dx->SetYTitle("Pmiss perp (GeV)");
    h_realPperp_W2_cut_coin_time_and_dx->Draw("COLZ");

    c51->Divide(2,2);
    c51->cd(1);
    h_theta_pq_cut_coin_time->SetXTitle("theta_pq (rad)");
    h_theta_pq_cut_coin_time->Draw();
    c51->cd(3);
    h_realPperp_cut_coin_time->SetXTitle("Pmiss perp (GeV)");
    h_realPperp_cut_coin_time->Draw();
    c51->cd(2);
    h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
    h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c51->cd(4);
    h_realPperp_W2_cut_coin_time->SetXTitle("W2 (GeV)");
    h_realPperp_W2_cut_coin_time->SetYTitle("Pmiss perp (GeV)");
    h_realPperp_W2_cut_coin_time->Draw("COLZ");


    //c6
	c6->Divide(2,2);
	c6->cd(1);
	h_realPpar->SetXTitle("Pmiss par (GeV)");
	h_realPpar->Draw();
	c6->cd(2);
	h_dx_realPpar->Draw("COLZ");
	c6->cd(3);
	h_realPpar_W2->Draw("COLZ");

	//c7
	c7->Divide(2,2);
	c7->cd(1);
	h_theta_pq_cut_QE->SetXTitle("theta_pq (rad)");
	h_theta_pq_cut_QE->Draw();
    c7->cd(2);
	h_realPpar_cut_QE->SetXTitle("Pmiss par (GeV)");
    h_realPpar_cut_QE->Draw();
    c7->cd(3);
	h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
	h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c7->cd(4);
	h_realPpar_W2_cut_QE->SetXTitle("W2 (GeV)");
	h_realPpar_W2_cut_QE->SetYTitle("Pmiss par (GeV)");
    h_realPpar_W2_cut_QE->Draw("COLZ");

    //c8
    c8->Divide(2,2);
    c8->cd(1);
    h_theta_pq_cut_coin_time_and_dx->SetXTitle("theta_pq (rad)");
    h_theta_pq_cut_coin_time_and_dx->Draw();
    c8->cd(3);
    h_realPpar_cut_coin_time_and_dx->SetXTitle("Pmiss par (GeV)");
    h_realPpar_cut_coin_time_and_dx->Draw();
    c8->cd(2);
    h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
    h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c8->cd(4);
    h_realPpar_W2_cut_coin_time_and_dx->SetXTitle("W2 (GeV)");
    h_realPpar_W2_cut_coin_time_and_dx->SetYTitle("Pmiss par (GeV)");
    h_realPpar_W2_cut_coin_time_and_dx->Draw("COLZ");

    c81->Divide(2,2);
    c81->cd(1);
    h_theta_pq_cut_coin_time->SetXTitle("theta_pq (rad)");
    h_theta_pq_cut_coin_time->Draw();
    c81->cd(3);
    h_realPpar_cut_coin_time->SetXTitle("Pmiss par (GeV)");
    h_realPpar_cut_coin_time->Draw();
    c81->cd(2);
    h_dx_W2_cut_coin_time->SetXTitle("W2 (GeV)");
    h_dx_W2_cut_coin_time->SetYTitle("dx (m)");
    h_dx_W2_cut_coin_time->Draw("COLZ");
    c81->cd(4);
    h_realPpar_W2_cut_coin_time->SetXTitle("W2 (GeV)");
    h_realPpar_W2_cut_coin_time->SetYTitle("Pmiss par (GeV)");
    h_realPpar_W2_cut_coin_time->Draw("COLZ");

    //c9
    c9->Divide(2,2);
    c9->cd(1);
    h_kinE->SetXTitle("kinetic energy (GeV)");
    h_kinE->Draw();
    c9->cd(2);
    h_kinE_cut_coin_time->SetXTitle("kinetic energy (GeV)");
    h_kinE_cut_coin_time->Draw();
    c9->cd(3);
    h_kinE_cut_QE->SetXTitle("kinetic energy (GeV)");
    h_kinE_cut_QE->Draw();

    //c10
    c10->Divide(2,2);
    c10->cd(1);
    h_W2_cut_QE->SetXTitle("W2 (GeV2)");
    h_W2_cut_QE->Draw();
    c10->cd(2);
    h_W2_cut_QEnPpar->SetXTitle("W2 (GeV2)");
    h_W2_cut_QEnPpar->Draw();
    c10->cd(3);
    h_W2_cut_Ppar->SetXTitle("W2 (GeV2)");
    h_W2_cut_Ppar->Draw();
    c10->cd(4);
    h_W2_cut_Cointime_Ppar->SetXTitle("W2 (GeV2)");
    h_W2_cut_Cointime_Ppar->Draw();

    //c11
    c11->Divide(2,2);
    c11->cd(1);
    h_coin_time_cut_Ppar->SetXTitle("cointime (ns)");
    h_coin_time_cut_Ppar->Draw();
    c11->cd(2);
    h_dxdy_cut_Ppar->SetXTitle("dy (m)");
    h_dxdy_cut_Ppar->SetYTitle("dx (m)");
    h_dxdy_cut_Ppar->Draw("COLZ");
    c11->cd(3);
    h_Ppar_coin_time->SetXTitle("cointime (ns)");
    h_Ppar_coin_time->SetYTitle("Ppar (GeV)");
    h_Ppar_coin_time->Draw("COLZ");
    //c11->cd(4);
    //h_Ppar_cut_W2->SetXTitle("Ppar (GeV)");
    //h_Ppar_cut_W2->Draw();

    //c12
    c12->Divide(2,2);
    c12->cd(1);
    h_Ppar_pNexpect_cut_cointime->SetXTitle("|q| (GeV)");
    h_Ppar_pNexpect_cut_cointime->SetYTitle("Ppar (GeV)");
    h_Ppar_pNexpect_cut_cointime->Draw("COLZ");
    c12->cd(2);
    h_Ppar_pNexpect_cut_cointime_dxdy->SetXTitle("|q| (GeV)");
    h_Ppar_pNexpect_cut_cointime_dxdy->SetYTitle("Ppar (GeV)");
    h_Ppar_pNexpect_cut_cointime_dxdy->Draw("COLZ");
    c12->cd(3);
    h_Ppar_vz_cut_cointime->SetXTitle("vz (m)");
    h_Ppar_vz_cut_cointime->SetYTitle("Ppar (GeV)");
    h_Ppar_vz_cut_cointime->Draw("COLZ");
    c12->cd(4);
    h_Ppar_vz_cut_cointime_dxdy->SetXTitle("vz (m)");
    h_Ppar_vz_cut_cointime_dxdy->SetYTitle("Ppar (GeV)");   
    h_Ppar_vz_cut_cointime_dxdy->Draw("COLZ");

    //c13
    c13->Divide(2,2);
    c13->cd(1);
    h_Ppar_eHCAL_cut_cointime->SetXTitle("KinE (GeV)");
    h_Ppar_eHCAL_cut_cointime->SetYTitle("Ppar (GeV)");
    h_Ppar_eHCAL_cut_cointime->Draw("COLZ");
    c13->cd(2);
    h_Ppar_eHCAL_cut_cointime_dxdy->SetXTitle("KinE (GeV)");
    h_Ppar_eHCAL_cut_cointime_dxdy->SetYTitle("Ppar (GeV)");
    h_Ppar_eHCAL_cut_cointime_dxdy->Draw("COLZ");
    c13->cd(3);
    h_Ppar_thetapq_cut_cointime->SetXTitle("theta_pq (rad)");
    h_Ppar_thetapq_cut_cointime->SetYTitle("Ppar (GeV)");
    h_Ppar_thetapq_cut_cointime->Draw("COLZ");
    c13->cd(4);
    h_Ppar_thetapq_cut_cointime_dxdy->SetXTitle("theta_pq (rad)");
    h_Ppar_thetapq_cut_cointime_dxdy->SetYTitle("Ppar (GeV)");
    h_Ppar_thetapq_cut_cointime_dxdy->Draw("COLZ");

    //c14
    c14->Divide(2,2);
    c14->cd(1);
    h_dxdy_cut_Ppar->SetXTitle("dy (m)");
    h_dxdy_cut_Ppar->SetYTitle("dx (m)");
    h_dxdy_cut_Ppar->Draw("COLZ");
    c14->cd(2);
    h_dxdy_cut_Ppar_cointime->SetXTitle("dy (m)");
    h_dxdy_cut_Ppar_cointime->SetYTitle("dx (m)");
    h_dxdy_cut_Ppar_cointime->Draw("COLZ");
    c14->cd(3);
    h_dxdy_cut_Ppar_cointime_W2->SetXTitle("dy (m)");
    h_dxdy_cut_Ppar_cointime_W2->SetYTitle("dx (m)");
    h_dxdy_cut_Ppar_cointime_W2->Draw("COLZ");  
    c14->cd(4);
    h_dxdy_cut_cointime_W2->SetXTitle("dy (m)");
    h_dxdy_cut_cointime_W2->SetYTitle("dx (m)");
    h_dxdy_cut_cointime_W2->Draw("COLZ");

    //c15
    c15->Divide(2,2);
    c15->cd(1);
    h_pmiss->SetXTitle("Pmiss (GeV)");
    h_pmiss->Draw();
    c15->cd(2);
    h_pmiss_cut_coin->SetXTitle("Pmiss (GeV)");
    h_pmiss_cut_coin->Draw();    
    c15->cd(3);
    h_pmiss_cut_coin_dxdy->SetXTitle("Pmiss (GeV)");
    h_pmiss_cut_coin_dxdy->Draw();

    //c16
    c16->Divide(2,2);
    c16->cd(1);
    h_Ppar_cut_hiW2->SetXTitle("Ppar (GeV)");
    h_Ppar_cut_hiW2->Draw();
    c16->cd(2);
    h_Ppar_cut_loW2->SetXTitle("Ppar (GeV)");
    h_Ppar_cut_loW2->Draw();

    //c17
    c17->Divide(2,2);
    c17->cd(1);
    h_Pmiss_pNexpect_cut_cointime->SetXTitle("|q| (GeV)");
    h_Pmiss_pNexpect_cut_cointime->SetYTitle("Pmiss (GeV)");
    h_Pmiss_pNexpect_cut_cointime->Draw("COLZ");
    c17->cd(2);
    h_Pmiss_pNexpect_cut_cointime_dxdy->SetXTitle("|q| (GeV)");
    h_Pmiss_pNexpect_cut_cointime_dxdy->SetYTitle("Pmiss (GeV)");
    h_Pmiss_pNexpect_cut_cointime_dxdy->Draw("COLZ");
    c17->cd(3);
    h_Pmiss_vz_cut_cointime->SetXTitle("vz (m)");
    h_Pmiss_vz_cut_cointime->SetYTitle("Pmiss (GeV)");
    h_Pmiss_vz_cut_cointime->Draw("COLZ");
    c17->cd(4);
    h_Pmiss_vz_cut_cointime_dxdy->SetXTitle("vz (m)");
    h_Pmiss_vz_cut_cointime_dxdy->SetYTitle("Pmiss (GeV)");   
    h_Pmiss_vz_cut_cointime_dxdy->Draw("COLZ");

    //c18
    c18->Divide(2,2);
    c18->cd(1);
    h_Pmiss_eHCAL_cut_cointime->SetXTitle("KinE (GeV)");
    h_Pmiss_eHCAL_cut_cointime->SetYTitle("Pmiss (GeV)");
    h_Pmiss_eHCAL_cut_cointime->Draw("COLZ");
    c18->cd(2);
    h_Pmiss_eHCAL_cut_cointime_dxdy->SetXTitle("KinE (GeV)");
    h_Pmiss_eHCAL_cut_cointime_dxdy->SetYTitle("Pmiss (GeV)");
    h_Pmiss_eHCAL_cut_cointime_dxdy->Draw("COLZ");
    c18->cd(3);
    h_Pmiss_thetapq_cut_cointime->SetXTitle("theta_pq (rad)");
    h_Pmiss_thetapq_cut_cointime->SetYTitle("Pmiss (GeV)");
    h_Pmiss_thetapq_cut_cointime->Draw("COLZ");
    c18->cd(4);
    h_Pmiss_thetapq_cut_cointime_dxdy->SetXTitle("theta_pq (rad)");
    h_Pmiss_thetapq_cut_cointime_dxdy->SetYTitle("Pmiss (GeV)");
    h_Pmiss_thetapq_cut_cointime_dxdy->Draw("COLZ");

    //c19
    c19->Divide(2,2);
    c19->cd(1);
    h_coin_time->SetXTitle("cointime (ns)");
    h_coin_time->Draw();
    c19->cd(2);
    h_dxdy_nocuts->SetXTitle("dy (m)");
    h_dxdy_nocuts->SetYTitle("dx (m)");
    h_dxdy_nocuts->Draw("COLZ");
    c19->cd(3);
    h_coin_time_cut_Ppar->SetXTitle("cointime (ns)");
    h_coin_time_cut_Ppar->Draw();
    c19->cd(4);
    h_dxdy_cut_Ppar->SetXTitle("dy (m)");
    h_dxdy_cut_Ppar->SetYTitle("dx (m)");
    h_dxdy_cut_Ppar->Draw("COLZ");

    //c20
    c20->Divide(2,2);
    c20->cd(1);
    h_coin_time_cut_KE->SetXTitle("cointime (ns)");
    h_coin_time_cut_KE->Draw();
    c20->cd(2);
    h_dxdy_cut_KE->SetXTitle("dy (m)");
    h_dxdy_cut_KE->SetYTitle("dx (m)");
    h_dxdy_cut_KE->Draw("COLZ");
    c20->cd(3);
    h_coin_time_cut_eHCAL->SetXTitle("cointime (ns)");
    h_coin_time_cut_eHCAL->Draw();
    c20->cd(4);
    h_dxdy_cut_eHCAL->SetXTitle("dy (m)");
    h_dxdy_cut_eHCAL->SetYTitle("dx (m)");
    h_dxdy_cut_eHCAL->Draw("COLZ");  

    c21->Divide(2,2);
    c21->cd(1);
    h_dx_cut_coin_time_W2_dy->SetXTitle("dx");
    h_dx_cut_coin_time_W2_dy->Draw();
    c21->cd(2);
    h_dx_cut_Ppar_coin_time_W2_dy->SetXTitle("dx");
    h_dx_cut_Ppar_coin_time_W2_dy->Draw();

    c22->Divide(2,2);
    c22->cd(1);
    h_eHCAL->SetXTitle("eHCAL (GeV)");
    h_eHCAL->Draw();
    c22->cd(2);
    h_kinE->SetXTitle("kinetic energy (GeV)");
    h_kinE->Draw();
    c22->cd(3);
    h_Ppar_eHCAL->SetXTitle("eHCAL(GeV)");
    h_Ppar_eHCAL->SetYTitle("Ppar(GeV)");
    h_Ppar_eHCAL->Draw("COLZ");



    ccuts->Print(Form("%s_Pperp.pdf(",kin));
	c->Print(Form("%s_Pperp.pdf",kin));
	c1->Print(Form("%s_Pperp.pdf",kin));
	c2->Print(Form("%s_Pperp.pdf",kin));
	c3->Print(Form("%s_Pperp.pdf",kin));
	c4->Print(Form("%s_Pperp.pdf",kin));
	c5->Print(Form("%s_Pperp.pdf",kin));
	c51->Print(Form("%s_Pperp.pdf",kin));
	c6->Print(Form("%s_Pperp.pdf",kin));
	c7->Print(Form("%s_Pperp.pdf",kin));
	c8->Print(Form("%s_Pperp.pdf",kin));
	c81->Print(Form("%s_Pperp.pdf",kin));
	c9->Print(Form("%s_Pperp.pdf",kin));
	c10->Print(Form("%s_Pperp.pdf",kin));
	c11->Print(Form("%s_Pperp.pdf",kin));
	c12->Print(Form("%s_Pperp.pdf",kin));
	c13->Print(Form("%s_Pperp.pdf",kin));
	c14->Print(Form("%s_Pperp.pdf",kin));
	c15->Print(Form("%s_Pperp.pdf",kin));
	c16->Print(Form("%s_Pperp.pdf",kin));
	c17->Print(Form("%s_Pperp.pdf",kin));
	c18->Print(Form("%s_Pperp.pdf",kin));
	c19->Print(Form("%s_Pperp.pdf",kin));
	c20->Print(Form("%s_Pperp.pdf",kin));
	c21->Print(Form("%s_Pperp.pdf",kin));
	c22->Print(Form("%s_Pperp.pdf)",kin));
}


