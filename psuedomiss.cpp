#include "parse.h"
#include "plotdxdy.h"

void psuedomiss(const char* filename, const char* printfilename, const char *kin){
	
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
	double W2 = 0.0;
	double coin_time = 0.0;

	tree->SetBranchAddress("theta_pq",&theta_pq);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("W2",&W2);
	tree->SetBranchAddress("coin_time",&coin_time);

	TH1D *h_theta_pq = new TH1D("h_theta_pq","theta_pq",200,0,2 );
	TH1D *h_Pperp = new TH1D("h_Pperp","Pperp",200,0,2);
	TH2D *h_dx_Pperp = new TH2D("h_dx_Pperp","dx v Pperp dist",100,0,1,250,-3,3);
	TH2D *h_Pperp_W2 = new TH2D("h_Pperp_W2","Pperp v W2 dist",1200,0,3,100,0,1);

	TH1D *h_theta_pq_cut_coin_time_and_dx = new TH1D("h_theta_pq_cut_coin_time_and_dx","theta_pq with coin and dx cuts",100,0,1);
        TH1D *h_Pperp_cut_coin_time_and_dx = new TH1D("h_Pperp_cut_coin_time_and_dx","Pperp with coin and dx cuts",100,0,1);
        TH2D *h_dx_Pperp_cut_coin_time_and_dx = new TH2D("h_dx_Pperp_cut_coin_time_and_dx","dx v Pperp dist with coin and dx cuts",100,0,1,250,-4,3);
        TH2D *h_Pperp_W2_cut_coin_time_and_dx = new TH2D("h_Pperp_W2_cut_coin_time_and_dx","Pperp v W2 dist with coin and dx cuts",1200,0,3,100,0,1);


	TH1D *h_theta_pq_cut_QE = new TH1D("h_theta_pq_cut_QE","theta_pq with QE cuts",30,0,0.3);
	TH1D *h_Pperp_cut_QE = new TH1D("h_Pperp_cut_QE","Pperp with QE cuts",40,0,0.4);
        TH2D *h_dx_Pperp_cut_QE = new TH2D("h_dx_Pperp_cut_QE","dx v Pperp dist with QE cuts",100,0,1,250,-4,3);
        TH2D *h_Pperp_W2_cut_QE = new TH2D("h_Pperp_W2_cut_QE","Pperp v W2 dist with QE cuts",1200,0,W2_H,100,0,0.3);
	TH2D *h_dx_W2_cut_coin_time = new TH2D("h_dx_W2_cut_coin_time","dx v W2 dist with coin time cut", 100,W2_L,3,250,-4,3);

	int nentries = tree->GetEntries();

	for (int i = 0; i<nentries; i++){	

		tree->GetEntry(i);
		//if (W2_L<W2 and W2<W2_H){
			double Pperp = theta_pq * sqrt(Q2);

			h_theta_pq->Fill(theta_pq);
			h_Pperp ->Fill(Pperp);
			h_dx_Pperp->Fill(Pperp,dx);
			h_Pperp_W2->Fill(W2,Pperp);
		//}
		if(coin_time_L<coin_time and coin_time<coin_time_H){ //and dx_L<dx and dx<dx_H and W2_L<W2 and W2<W2_H){
			
			h_dx_W2_cut_coin_time->Fill(W2,dx);

			if (dx_L<dx and dx<dx_H){
				h_theta_pq_cut_coin_time_and_dx->Fill(theta_pq);
                                h_Pperp_cut_coin_time_and_dx->Fill(Pperp);
                                h_dx_Pperp_cut_coin_time_and_dx->Fill(Pperp,dx);
                                h_Pperp_W2_cut_coin_time_and_dx->Fill(W2,Pperp);

				if(W2_L<W2 and W2<W2_H){
					h_theta_pq_cut_QE->Fill(theta_pq);
					h_Pperp_cut_QE->Fill(Pperp);
					h_dx_Pperp_cut_QE->Fill(Pperp,dx);
					h_Pperp_W2_cut_QE->Fill(W2,Pperp);
					//h_dx_W2_cut_coin_time->Fill(W2,dx);
				}
			}
		}

		if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();	

	}

	//square cut
	TCutG *cutsq = CreateSquareCut(W2_L,dx_L,W2_H,dx_H);
	//TCutG *cutsq1 = CreateSquareCut(-10,dx_L,10,dx_H);

	
	TCanvas* c = new TCanvas("c","c",1200,1600);
	TCanvas* c1 = new TCanvas("c1","c1",1200,1600);
	TCanvas* c2 = new TCanvas("c2","c2",1200,1600);
	
        c->Divide(2,2);
	c->cd(1);
	h_theta_pq->Draw();
	c->cd(2);
	h_Pperp->Draw();
	c->cd(3);
	h_dx_Pperp->Draw();
	c->cd(4);
	h_Pperp_W2->Draw();

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
	cutsq->Draw("L");
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

	c->Print(Form("%s_Pperp.pdf(",kin));
	c1->Print(Form("%s_Pperp.pdf",kin));
	c2->Print(Form("%s_Pperp.pdf)",kin));
	
}
