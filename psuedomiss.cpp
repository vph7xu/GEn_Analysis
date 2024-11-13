#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "parse.h"
#include "plotdxdy.h"
#include <cmath>

void psuedomiss(const char* filename, const char* printfilename, const char *kin){
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
	double trP_sbs = 0.0;
	double ntrack_sbs = 0.0;

	tree->SetBranchAddress("theta_pq",&theta_pq);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("pN_expect",&pN_expect);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("eHCAL",&eHCAL);
	tree->SetBranchAddress("hcal_clus_e",&hcal_clus_e);
	tree->SetBranchAddress("hcal_clus_id",&hcal_clus_id);
	tree->SetBranchAddress("nblk_HCAL",&nblk_HCAL);
	tree->SetBranchAddress("trP_sbs",&trP_sbs);
	tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);

	TH1D *h_theta_pq = new TH1D("h_theta_pq","theta_pq",200,0,2 );
	TH1D *h_Pperp = new TH1D("h_Pperp","Pperp",200,0,2);
	TH1D *h_realPperp = new TH1D("h_realPperp","realPperp",200,0,2);
	TH1D *h_realPpar = new TH1D("h_realPpar","realPpar",200,-5,5);

	TH2D *h_dx_Pperp = new TH2D("h_dx_Pperp","dx v Pperp dist",100,0,1,250,-3,3);
	TH2D *h_dx_realPperp = new TH2D("h_dx_realPperp","dx v realPperp dist",100,0,1,250,-3,3);
	TH2D *h_dx_realPpar = new TH2D("h_dx_realPpar","dx v realPpar dist",100,-5,5,250,-3,3);


	TH2D *h_Pperp_W2 = new TH2D("h_Pperp_W2","Pperp v W2 dist",1200,0,3,100,0,1);
	TH2D *h_realPperp_W2 = new TH2D("h_realPperp_W2","realPperp v W2 dist",1200,0,3,100,0,1);
	TH2D *h_realPpar_W2 = new TH2D("h_realPpar_W2","realPpar v W2 dist",1200,0,3,100,-5,5);

	TH1D *h_theta_pq_cut_coin_time_and_dx = new TH1D("h_theta_pq_cut_coin_time_and_dx","theta_pq only with coin cut",100,0,1);
    TH1D *h_Pperp_cut_coin_time_and_dx = new TH1D("h_Pperp_cut_coin_time_and_dx","Pperp only with coin cut",100,0,1);
    TH1D *h_realPperp_cut_coin_time_and_dx = new TH1D("h_realPperp_cut_coin_time_and_dx","realPperp only with coin cut",200,-0.5,2);
    TH1D *h_realPpar_cut_coin_time_and_dx = new TH1D("h_realPpar_cut_coin_time_and_dx","realPpar only with coin cut",200,-5,5);
 
    TH2D *h_dx_Pperp_cut_coin_time_and_dx = new TH2D("h_dx_Pperp_cut_coin_time_dy_and_dx","dx v Pperp dist only with coin cut",100,0,1,250,-4,3);
    TH2D *h_dx_realPperp_cut_coin_time_and_dx = new TH2D("h_dx_realPperp_cut_coin_time_dy_and_dx","dx v realPperp dist only with coin cut",100,0,1,250,-4,3);
    TH2D *h_dx_realPpar_cut_coin_time_and_dx = new TH2D("h_dx_realPpar_cut_coin_time_dy_and_dx","dx v realPpar dist only with coin cut",100,0,1,250,-4,3);	
	TH2D *h_Pperp_W2_cut_coin_time_and_dx = new TH2D("h_Pperp_W2_cut_coin_time_dy and_dx","Pperp v W2 dist only with coin",1200,0,3,100,0,1);
	TH2D *h_realPperp_W2_cut_coin_time_and_dx = new TH2D("h_realPperp_W2_cut_coin_time_dy and_dx","realPperp v W2 dist only with coin cut",1200,-2,4,200,-0.5,2);
	TH2D *h_realPpar_W2_cut_coin_time_and_dx = new TH2D("h_realPpar_W2_cut_coin_time_dy_and_dx","realPpar v W2 dist only with coin cut ",1200,-2,4,200,-5,5);

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

	TH1D *h_kinE = new TH1D("h_kinE","kinetic energy", 200,-1,10);
	TH1D *h_kinE_cut_coin_time = new TH1D("h_kinE_cut_coin_time","kinetic energy with cointime cut",200,-1,10);
	TH1D *h_kinE_cut_W2 = new TH1D("h_kinE_cut_W2","kintic energy with cointime and W2 cuts",200,-1,10);
	TH1D *h_kinE_cut_QE = new TH1D("h_kinE_cut_QE","kinetic energy with QE cuts",200,-1,10);

	int nentries = tree->GetEntries();

	for (int i = 0; i<nentries; i++){	
		tree->GetEntry(i);
		//if (ntrack_sbs>0){

		//if (W2_L<W2 and W2<W2_H){
			double KinE = 0.0;

			for (int j = 0; j<nblk_HCAL; j++){
				KinE += hcal_clus_e[j]/sampling_fractions_each_blk[hcal_clus_id[j]];

				if (i%1000 == 0){ 
					std::cout<<"event : "<<i<<"blk id: "<<hcal_clus_id[j]<<" blk energy : "<<hcal_clus_e[j]<<" blk sf : " <<sampling_fractions_each_blk[hcal_clus_id[j]]<<endl;
				}
			}
			if (i%1000 == 0){ 
				std::cout<<"event : "<<i<<" kinE : "<<KinE<<endl;
			}
			double Pperp = theta_pq * pN_expect;
			double realPperp = (theta_pq * sqrt(-pow(0.938,2)+pow((KinE+0.938),2))); // this probably has to change depending on whether is a P or N 
			double realPpar = (cos(theta_pq) * sqrt(-pow(0.938,2)+pow((KinE+0.938),2)))-pN_expect; 
			
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

			h_kinE->Fill(KinE);

		//}
		if(coin_time_L<coin_time and coin_time<coin_time_H){ //and dx_L<dx and dx<dx_H and W2_L<W2 and W2<W2_H){
			
			h_dx_W2_cut_coin_time->Fill(W2,dx);
			h_kinE_cut_coin_time->Fill(KinE);

			//if (dx_L<dx and dx<dx_H and dy_L<dy and dy<dy_H){
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
			//}
		}
		//}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();	

	}

	//square cut
	//TCutG *cutsq = CreateSquareCut(W2_L,dx_L,W2_H,dx_H);
	//TCutG *cutsq1 = CreateSquareCut(-10,dx_L,10,dx_H);

	
	TCanvas* c = new TCanvas("c","c",1200,1600);
	TCanvas* c1 = new TCanvas("c1","c1",1200,1600);
	TCanvas* c2 = new TCanvas("c2","c2",1200,1600);
	TCanvas* c3 = new TCanvas("c3","c3",1200,1600);
	TCanvas* c4 = new TCanvas("c4","c4",1200,1600);
	TCanvas* c5 = new TCanvas("c5","c5",1200,1600);
	TCanvas* c6 = new TCanvas("c6","c6",1200,1600);
	TCanvas* c7 = new TCanvas("c7","c7",1200,1600);
	TCanvas* c8 = new TCanvas("c8","c8",1200,1600);
	TCanvas* c9 = new TCanvas("c9","c9",1200,1600);
	
    c->Divide(2,2);
	c->cd(1);
	h_theta_pq->Draw();
	c->cd(2);
	h_Pperp->Draw();
	c->cd(3);
	h_dx_Pperp->Draw();
	c->cd(4);
	h_Pperp_W2->Draw();

	TCutG *cutsq = CreateSquareCut(-2,dx_L,2,dx_H);

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
    h_dx_W2_cut_coin_time->Draw();
	cutsq->Draw("same");
    c1->cd(4);
	h_Pperp_W2_cut_QE->SetXTitle("W2 (GeV)");
	h_Pperp_W2_cut_QE->SetYTitle("Pseudo Pmiss perp (GeV)");
    h_Pperp_W2_cut_QE->Draw();

	TCutG *cutsq1 = CreateSquareCut(-10,dx_L,10,dx_H);

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
    h_dx_W2_cut_coin_time->Draw();
    cutsq1->Draw("L");
    c2->cd(4);
    h_Pperp_W2_cut_coin_time_and_dx->SetXTitle("W2 (GeV)");
    h_Pperp_W2_cut_coin_time_and_dx->SetYTitle("Pseudo Pmiss perp (GeV)");
    h_Pperp_W2_cut_coin_time_and_dx->Draw();

	c3->Divide(2,2);
	c3->cd(1);
	h_realPperp->SetXTitle("Pmiss perp (GeV)");
	h_realPperp->Draw();
	c3->cd(2);
	h_dx_realPperp->Draw();
	c3->cd(3);
	h_realPperp_W2->Draw();

	TCutG *cutsq2 = CreateSquareCut(-2,dx_L,2,dx_H);

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
    h_dx_W2_cut_coin_time->Draw();
	cutsq2->Draw("same");
    c4->cd(4);
	h_realPperp_W2_cut_QE->SetXTitle("W2 (GeV)");
	h_realPperp_W2_cut_QE->SetYTitle("Pmiss perp (GeV)");
    h_realPperp_W2_cut_QE->Draw();

	//TCutG *cutsq3 = CreateSquareCut(-10,dx_L,10,dx_H);

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
    h_dx_W2_cut_coin_time->Draw();
    //cutsq3->Draw("L");
    c5->cd(4);
    h_realPperp_W2_cut_coin_time_and_dx->SetXTitle("W2 (GeV)");
    h_realPperp_W2_cut_coin_time_and_dx->SetYTitle("Pmiss perp (GeV)");
    h_realPperp_W2_cut_coin_time_and_dx->Draw();


	c6->Divide(2,2);
	c6->cd(1);
	h_realPpar->SetXTitle("Pmiss par (GeV)");
	h_realPpar->Draw();
	c6->cd(2);
	h_dx_realPpar->Draw();
	c6->cd(3);
	h_realPpar_W2->Draw();

	//TCutG *cutsq2 = CreateSquareCut(-2,dx_L,2,dx_H);

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
    h_dx_W2_cut_coin_time->Draw();
	//cutsq2->Draw("same");
    c7->cd(4);
	h_realPpar_W2_cut_QE->SetXTitle("W2 (GeV)");
	h_realPpar_W2_cut_QE->SetYTitle("Pmiss par (GeV)");
    h_realPpar_W2_cut_QE->Draw();

	TCutG *cutsq3 = CreateSquareCut(-10,dx_L,-10,dx_H);

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
    h_dx_W2_cut_coin_time->Draw();
    cutsq3->Draw("L");
    c8->cd(4);
    h_realPpar_W2_cut_coin_time_and_dx->SetXTitle("W2 (GeV)");
    h_realPpar_W2_cut_coin_time_and_dx->SetYTitle("Pmiss par (GeV)");
    h_realPpar_W2_cut_coin_time_and_dx->Draw();

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


	c->Print(Form("plots/%s_Pperp.pdf(",kin));
	c1->Print(Form("plots/%s_Pperp.pdf",kin));
	c2->Print(Form("plots/%s_Pperp.pdf",kin));
	c3->Print(Form("plots/%s_Pperp.pdf",kin));
	c4->Print(Form("plots/%s_Pperp.pdf",kin));
	c5->Print(Form("plots/%s_Pperp.pdf",kin));
	c6->Print(Form("plots/%s_Pperp.pdf",kin));
	c7->Print(Form("plots/%s_Pperp.pdf",kin));
	c8->Print(Form("plots/%s_Pperp.pdf",kin));
	c9->Print(Form("plots/%s_Pperp.pdf)",kin));
}
