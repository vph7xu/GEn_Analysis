#include "models.h"
#include "parse.h"
#include "plotdxdy.h"
#include <cmath>

void raw_asymmetry(const char* filename, const char* printfilename, const char* kin){

	std::map<std::string, std::string> config1 = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	
	double coin_time_L = getDoubleValue(config1,"coin_time_L");
	double coin_time_H = getDoubleValue(config1,"coin_time_H");

	double W2_L = getDoubleValue(config1,"W2_L");
	double W2_H = getDoubleValue(config1,"W2_H");

	double dy_L = getDoubleValue(config1,"dy_L");
	double dy_H = getDoubleValue(config1,"dy_H");

	double dx_L = getDoubleValue(config1,"dx_L");
	double dx_H = getDoubleValue(config1,"dx_H");

	double IHWP_flip = getDoubleValue(config1,"IHWP_flip");

	double run_num_L = getDoubleValue(config1,"run_num_L");
	double run_num_H = getDoubleValue(config1,"run_num_H");

	//end of parsing cuts
	
	TFile* file = TFile::Open(filename);
	//TFile* sim_file = TFile::Open(sim_filename);

        TTree* tree = (TTree*)file->Get("Tout");
	//TTree* sim_tree = (TTree*)sim_file->Get("Tout");
	
	double dx = 0.0;
	double dy = 0.0;
	double W2 = 0.0;
	double coin_time = 0.0;
	//double coin_time = 0.0;

	int runnum = 0;
	int helicity = 0;
	int IHWP = 0;
	double He3Pol = 0.0;

	double Nplus = 0;
	double Nminus = 0;

	double Aexp = 0.0;
	double errAexp = 0.0;

	double runx=0;

	tree->SetBranchAddress("runnum",&runnum);
	tree->SetBranchAddress("helicity",&helicity);
	tree->SetBranchAddress("IHWP",&IHWP);
	tree->SetBranchAddress("He3Pol",&He3Pol);
	tree->SetBranchAddress("dx",&dx);
	tree->SetBranchAddress("dy",&dy);
	tree->SetBranchAddress("coin_time",&coin_time);
	tree->SetBranchAddress("W2",&W2);

	TH1I *h_IHWP = new TH1I("h_IHWP","IHWP",10,-2,2);
	TH1I *h_runnum = new TH1I("h_runnum","runnum",300,0,3000);

	TGraphErrors *gAsym = new TGraphErrors();

	//neutron spot cut
	TCutG *cutg = CreateOvalCut("cutg",-0.05,0,0.5,0.35,100);
	TH2D* nspot_cut = new TH2D("nspot_cut","neutron spot",250,-2,2,250,-4,4);

	int nentries = tree->GetEntries();

	for (int i=0; i<nentries; i++){
		tree->GetEntry(i);
		
		h_IHWP->Fill(IHWP);
		h_runnum->Fill(runnum);

		if (i==0){
			runx = runnum;
		}
		
		if(runnum==runx){
			
			if (IHWP == 1) helicity *= -1*IHWP_flip;
                       	else if (IHWP == -1) helicity *= -1*IHWP_flip;
                        else continue;

			if ((W2_L<W2 and W2<W2_H)and(coin_time_L<coin_time and coin_time<coin_time_H)){
				if(cutg->IsInside(dy,dx)){
				//if ((dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H)){
				nspot_cut->Fill(dy,dx);

					if (helicity==1){
						Nplus=Nplus+1;
					}
					else if (helicity==-1){
						Nminus=Nminus+1;
					}
					//}
					//else{
					//	if (helicity==1){
                                        //                Nminus=Nminus+1;
                                        //        }
                                        //        else if (helicity==-1){
                                        //                Nplus=Nplus+1;
                                        //        }
					//}
				}
			}
		}

		else if(runnum!=runx){
			Aexp = (Nplus-Nminus)*100/(Nplus+Nminus);
			errAexp = 2*100*sqrt((Nplus*Nminus)*(Nplus+Nminus))/((Nplus+Nminus)*(Nplus+Nminus));
			gAsym->SetPoint(runnum,runnum,Aexp);
			gAsym->SetPointError(runnum,0,errAexp);

			std::cout<<"runnum : "<<runnum<<" Nplus : "<<Nplus<<" Nminus : "<<Nminus<<" Aexp : "<<Aexp<<endl;

			runx=runnum;
			Nplus=0;
			Nminus=0;
		}



	}
	
	TCanvas *c = new TCanvas("c","c",800,1200);
	TCanvas *c1 = new TCanvas("c1","c1",800,1200);
	
	c->Divide(2,2);
	c->cd(1);
	h_IHWP->Draw();
	c->cd(2);
	h_runnum->Draw();
	c->cd(3);
	nspot_cut->Draw("COLZ");
	cutg->Draw("L");
	
	c1->Divide(1,1);
	c1->cd(1);
	gAsym->GetXaxis()->SetLimits(run_num_L,run_num_H);
	gAsym->GetYaxis()->SetRangeUser(-15,15);
	gAsym->SetMarkerStyle(7);
	gAsym->SetMarkerSize(3.0);
	gAsym->SetMarkerColor(kBlue);
	gAsym->Draw("AP");

	c1->Print(Form("%s.pdf",printfilename));	

}
