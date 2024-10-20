#include "models.h"
#include "parse.h"
#include "plotdxdy.h"
#include <cmath>

void accidental_asymmetry(const char* filename, const char* printfilename, const char* kin){

	std::map<std::string, std::string> config1 = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	std::map<int, int> HelicityCheck = readCSVToMap("DB/Helicity_quality.csv");
	std::map<int, int> MollerQuality = readCSVToMap("DB/Moller_quality.csv");
	
	double coin_time_L = getDoubleValue(config1,"coin_time_L");
	double coin_time_H = getDoubleValue(config1,"coin_time_H");

	double W2_L = getDoubleValue(config1,"W2_L");
	double W2_H = getDoubleValue(config1,"W2_H");

	double dy_L = getDoubleValue(config1,"dy_L");
	double dy_H = getDoubleValue(config1,"dy_H");

	double dx_L = getDoubleValue(config1,"dx_L");
	double dx_H = getDoubleValue(config1,"dx_H");

	double dy_p_L = getDoubleValue(config1,"dy_p_L");
	double dy_p_H = getDoubleValue(config1,"dy_p_H");

	double dx_p_L = getDoubleValue(config1,"dx_p_L");
	double dx_p_H = getDoubleValue(config1,"dx_p_H");

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

	double Nplus_total = 0.0;
	double Nminus_total = 0.0;

	double Aexp = 0.0;
	double errAexp = 0.0;

	double QE_events = 0.0;
	double acc_events = 0.0;

	double QE_events_total = 0.0;
	double acc_events_total = 0.0;

	double acc_frac = 0.0;

	double runx=0;

    // Varying coin time cuts
    //double step_size = 10;  // Adjust this step size as needed
    //int num_steps = (coin_time_H - coin_time_L) / step_size;

	//int helicity_new = 0;

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
	TH1D *h_cointime = new TH1D("h_cointime","cointime (ns)", 500,0,200);
	TH1D *h_cointime_anticut = new TH1D("h_cointime_anticut","cointime (ns) - anticut", 500,0,200);
	TH1D *h_cointime_offset = new TH1D("h_cointime_offset","cointime (ns) - offset", 500,0,200); 

	TGraphErrors *gAsym = new TGraphErrors();
	TGraphErrors *gAsymCut = new TGraphErrors();
	//TGraphErrors *gAp_sym = new TGraphErrors();

	//neutron spot cut
	//TCutG *cutg = CreateOvalCut("cutg",0,-0.025,0.35,0.2,100);

	TH2D* nspot_cut = new TH2D("nspot_cut","neutron spot",250,-2,2,250,-4,4);

	int nentries = tree->GetEntries();

	//write asymmetry values to files
	std::ofstream outfile;
	outfile.open("accidental_asymmetry_results.txt"); 

	std::ofstream outfile_p;
	outfile_p.open("accidental_fraction_results.txt");


	/*for (int step = 0; step <= num_steps; step++) {
        double current_cut_L = coin_time_L + step * step_size;
        double current_cut_H = current_cut_L + step_size;

        double Nplus_new = 0.0;
        double Nminus_new = 0.0;

        for (int i=0; i<nentries; i++) {
            tree->GetEntry(i);

            helicity_new = IHWP*IHWP_flip*helicity;

            if ((W2_L < W2 && W2 < W2_H) && (dx_L < dx && dx < dx_H) && (dy_L < dy && dy < dy_H) && (current_cut_L < coin_time && coin_time < current_cut_H)) {
                if (helicity_new == 1) {
                   	Nplus_new++;
                } else if (helicity_new == -1) {
                    Nminus_new++;
                }
            }

            if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r"; 
			std::cout.flush();
        }
        	// Calculate asymmetry and error
        double Aexp_new = (Nplus_new - Nminus_new) * 100 / (Nplus_new + Nminus_new);
        double errAexp_new = 2 * 100 * sqrt((Nplus_new * Nminus_new) * (Nplus_new + Nminus_new)) / ((Nplus_new + Nminus_new) * (Nplus_new + Nminus_new));

        	// Store the results in the graph
        gAsymCut->SetPoint(step, current_cut_L, Aexp_new);
        gAsymCut->SetPointError(step, 0, errAexp_new);

        std::cout << "Cut: [" << current_cut_L << ", " << current_cut_H << "] Asymmetry: " << Aexp_new << " ± " << errAexp_new << std::endl;
    }*/


	for (int i=0; i<nentries; i++){
		tree->GetEntry(i);
		
		h_IHWP->Fill(IHWP);
		h_runnum->Fill(runnum);

		//plots
        if ((W2_L<W2 and W2<W2_H)and(dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H)){
            h_cointime->Fill(coin_time);
        }

        if ((W2_L<W2 and W2<W2_H)and(dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H) and (coin_time_L+30<coin_time and coin_time<coin_time_H+30)){
            h_cointime_offset->Fill(coin_time);
        }

		if ((W2_L<W2 and W2<W2_H)and(dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H) and (coin_time_L>coin_time or coin_time>coin_time_H)){
			if ((40<coin_time and coin_time<180)){ //removing edging effect, for now hardcoded for GEN2
				h_cointime_anticut->Fill(coin_time);
			}
		}

		if ((W2_L<W2 and W2<W2_H)and(coin_time_L<coin_time and coin_time<coin_time_H)){
			nspot_cut->Fill(dy,dx);
		}

		if (i==0){
			runx = runnum;
		}
		
		if(runnum==runx and lookupValue(HelicityCheck,runnum)==1 and lookupValue(MollerQuality,runnum)==1){
			
			if (IHWP == 1) helicity = IHWP*IHWP_flip*helicity;
            else if (IHWP == -1) helicity = IHWP*IHWP_flip*helicity; //if this is only to show the effect of the half wave plats, right now it has no effect
            else continue;

            //count QE events
            if ((W2_L<W2 and W2<W2_H)and(dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H) and (coin_time_L<coin_time and coin_time<coin_time_H)){
            	QE_events+=1;
            }

            //count accidental events
            if ((W2_L<W2 and W2<W2_H)and(dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H) and (coin_time_L+30<coin_time and coin_time<coin_time_H+30)){
            	acc_events+=1;
            }

			if ((W2_L<W2 and W2<W2_H)and(dx_L<dx and dx<dx_H) and (dy_L<dy and dy<dy_H) and (coin_time_L>coin_time or coin_time>coin_time_H)){

				if ((40<coin_time and coin_time<180)){ //removing edging effect, for now hardcoded for GEN2

					if (helicity==1){
						Nplus=Nplus+1;
					}
					else if (helicity==-1){
						Nminus=Nminus+1;
					}
				}
			}
		}

		else if(runnum!=runx and lookupValue(HelicityCheck,runx)==1 and lookupValue(MollerQuality,runnum)==1){
			Aexp = (Nplus-Nminus)*100/(Nplus+Nminus);
			errAexp = 2*100*sqrt((Nplus*Nminus)*(Nplus+Nminus))/((Nplus+Nminus)*(Nplus+Nminus));
			acc_frac = acc_events/QE_events;

			gAsym->SetPoint(runx,runx,Aexp);
			gAsym->SetPointError(runx,0,errAexp);

			outfile << "Run: " << runx << " Aexp: " << Aexp << " ± " << errAexp << std::endl;
            
            outfile_p << "Run: " << runx <<" acc_frac : "<<acc_frac<<endl;

			std::cout<<"runnum : "<<runx<<" Nplus : "<<Nplus<<" Nminus : "<<Nminus<<" Aexp : "<<Aexp<<endl;
			std::cout<<"runnum : "<<runx<<" QE_events : "<<QE_events<<" acc_events : "<<acc_events<<" acc_frac : "<<acc_frac<<endl;

			Nplus_total+=Nplus;
			Nminus_total+=Nminus;
			QE_events_total+=QE_events;
            acc_events_total+=acc_events;


			runx=runnum;
			Nplus=0.0;
            Nminus=0.0;
            Aexp=0.0;
            QE_events=0.0;
            acc_events=0.0;

		}

		else if(runnum!=runx){
			runx=runnum;
            Nplus=0.0;
            Nminus=0.0;
            Aexp=0.0;
            QE_events=0.0;
            acc_events=0.0;
		}

	}
	outfile.close();
	outfile_p.close();

	//TCutG * cutsq = CreateSquareCut(dy_L,dx_L,dy_H,dx_H);
	//TCutG * cutsq1 = CreateSquareCut(dy_p_L,dx_p_L,dy_p_H,dx_p_H);

	std::cout<<"Nplus_total : "<<Nplus_total<<endl;
	std::cout<<"Nminus_total : "<<Nminus_total<<endl;
	std::cout<<"Aacc : " <<((Nplus_total-Nminus_total)/(Nplus_total+Nminus_total))<<endl;
	std::cout<<"errAacc : "<<2*sqrt((Nplus_total*Nminus_total)*(Nplus_total+Nminus_total))/((Nplus_total+Nminus_total)*(Nplus_total+Nminus_total))<<endl;
	std::cout<<"QE_events : "<<QE_events_total<<endl;
	std::cout<<"acc_events : "<<acc_events_total<<endl;
	std::cout<<"acc_frac : "<<(acc_events_total/QE_events_total)<<endl;
	std::cout<<"erracc_frac : "<<(acc_events_total/QE_events_total)*sqrt((1/acc_events_total)+(1/QE_events_total)) <<endl;


	TCanvas *c = new TCanvas("c","c",3600,3000);
	TCanvas *ccoin = new TCanvas("ccoin","ccoin",3600,3000);
	TCanvas *c1 = new TCanvas("c1","c1",3600,900);
	
	c->Divide(2,2);
	c->cd(1);
	h_IHWP->Draw();
	c->cd(2);
	h_runnum->Draw();
	c->cd(3);
	nspot_cut->Draw("COLZ");
	CreateSquareCut(dy_L,dx_L,dy_H,dx_H);
	CreateSquareCut(dy_p_L,dx_p_L,dy_p_H,dx_p_H);
	//cutg->Draw("L");
	
	ccoin->Divide(2,2);
	ccoin->cd(1);
	h_cointime->Draw();
	ccoin->cd(2);
	h_cointime_offset->Draw();
	ccoin->cd(3);
	h_cointime_anticut->Draw();

	c1->Divide(1,1);
	c1->cd(1);
	gAsym->GetXaxis()->SetTitle("run number");
	gAsym->GetYaxis()->SetTitle("Accidental Asymmetry");
	gAsym->GetXaxis()->SetLimits(run_num_L,run_num_H);
	gAsym->GetYaxis()->SetRangeUser(-10,10);
	gAsym->SetMarkerStyle(7);
	gAsym->SetMarkerSize(5.0);
	gAsym->SetMarkerColor(kBlue);
	gAsym->Draw("AP");

    TLegend *legend = new TLegend(0.7, 0.8, 0.9, 0.9);  // Coordinates for the legend box (x1, y1, x2, y2)
    
    // Add entries to the legend
    legend->AddEntry(gAsym, "(e,e'n) events", "lp"); // "lp" means line and point in the legend

    // Set legend properties
    legend->SetTextSize(0.03);    // Set text size
    legend->SetBorderSize(0);     // No border around the legend box
    legend->SetFillStyle(0);      // Transparent background for the legend

    // Draw the legend
    legend->Draw();

    /*TCanvas *cCut = new TCanvas("cCut", "Asymmetry vs Coin Time Cut", 1200, 900);
    gAsymCut->SetTitle("Asymmetry vs Coin Time Cut");
    gAsymCut->GetXaxis()->SetTitle("Coin Time Lower Bound (ns)");
    gAsymCut->GetYaxis()->SetTitle("Asymmetry (%)");
    gAsymCut->SetMarkerStyle(7);
    gAsymCut->SetMarkerSize(1.0);
    gAsymCut->SetMarkerColor(kBlue);
    gAsymCut->Draw("AP");

    cCut->SaveAs(Form("%s_varying_cut.png", printfilename));
    cCut->Print(Form("%s_varying_cut.pdf", printfilename));

	c->Print(Form("%s.pdf(",printfilename));
	ccoin->Print(Form("%s.pdf",printfilename));
	c1->Print(Form("%s.pdf)",printfilename));	
	c1->SaveAs(Form("%s.png",printfilename));*/
}
