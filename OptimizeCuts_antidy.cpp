// OptimizeCuts.C — complete ROOT macro for cut optimization
// ------------------------------------------------------------
// * Loads data & simulation trees via TTree::GetEntry loops
// * Auto-scales simulated neutrons & protons to data
// * For each trial rectangle, builds background, data, sim histograms
// * Composes model = simN + simP + background
// * Computes χ²(data, model) as figure of merit
// * Uses ROOT::Minuit2 to find optimal cuts on dx, dy, W2, eHCAL, coin_time
// ------------------------------------------------------------
#include "models.h"

#include <TF1.h> 
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TStyle.h>  // for extern TStyle *gStyle
#include <TROOT.h>   // if you also use gROOT anywhere
#include <TCanvas.h>


struct Store {
    std::vector<double> vz, ePS, dx, dy, W2, ct, eHCAL, w;
    std::vector<int> helicity, IHWP; 
    std::vector<char>   origin;  // 0=nMC,1=pMC,2=data
};

void OptimizeCuts(const char* dataFile,
                  const char* simFile,
                  const char* kin,
                  int    nbinsDx = 100,
                  double wN_user = -1.0,
                  double wP_user = -1.0)
{
    Store ev;

    std::ofstream outfile(Form("test_FOM_%s.txt",kin));
    std::ofstream outfile_all(Form("test_FOM_all_%s.txt",kin));

    gStyle->SetOptFit(0);

    // 1. Load DATA
    {
        TFile fData(dataFile, "READ");
        TTree* tData = (TTree*)fData.Get("Tout");
        double good_total_events = 0.0;

        double vz, ePS, eSH, trP, dx, dy, W2, ct, eHCAL;
        int helicity, IHWP;
        
        tData->SetBranchAddress("vz",&vz);
        tData->SetBranchAddress("ePS",&ePS);
        tData->SetBranchAddress("eSH",&eSH);
        tData->SetBranchAddress("trP",&trP);    
        tData->SetBranchAddress("dx", &dx);
        tData->SetBranchAddress("dy", &dy);
        tData->SetBranchAddress("W2", &W2);
        tData->SetBranchAddress("coin_time", &ct);
        tData->SetBranchAddress("eHCAL", &eHCAL);
        tData->SetBranchAddress("helicity",&helicity);
        tData->SetBranchAddress("IHWP",&IHWP);
        Long64_t nData = tData->GetEntries();
        ev.vz.reserve(nData+1000000);
        ev.ePS.reserve(nData+1000000);
        ev.dx.reserve(nData+1000000);
        ev.dy.reserve(nData+1000000);
        ev.W2.reserve(nData+1000000);
        ev.ct.reserve(nData+1000000);
        ev.eHCAL.reserve(nData+1000000);
        ev.helicity.reserve(nData+1000000);
        ev.IHWP.reserve(nData+1000000);
        ev.w.reserve(nData+1000000);
        ev.origin.reserve(nData+1000000);
        for(Long64_t i=0; i<1*nData; ++i) {

            tData->GetEntry(i);
            if(abs(vz)>0.27 && ePS<0.2 && abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL<0.025 && abs(helicity)!=1) continue; // you can add basic cuts here if needed
            ++good_total_events;
            ev.vz.push_back(vz);
            ev.ePS.push_back(ePS);
            ev.dx.push_back(dx);
            ev.dy.push_back(dy);
            ev.W2.push_back(W2);
            ev.ct.push_back(ct);
            ev.eHCAL.push_back(eHCAL);
            ev.helicity.push_back(helicity);
            ev.IHWP.push_back(IHWP);
            ev.w.push_back(1.0);
            ev.origin.push_back(2);
            if (i %1000 == 0 ) std::cout<< (i * 100.0/ nData)<<" data file" << "% \r";
            std::cout.flush();
        }
        std::cout<<"Loaded " << ev.dx.size() << " total events (data only)" << std::endl;
        std::cout<<"Good total events in data: " <<good_total_events<< std::endl;
    }



    // 2. Load SIM
    bool simHasTime = false;
    {
        TFile fSim(simFile, "READ");
        TTree* t = (TTree*)fSim.Get("Tout");
        if(t->GetBranch("coin_time")) simHasTime = true;
    }
    {
        TFile fSim(simFile, "READ");
        TTree* tSim = (TTree*)fSim.Get("Tout");
        double vz, ePS, eSH, trP, dx, dy, W2, ctSim, eHCAL, wt, fnucl;
        tSim->SetBranchAddress("vz", &vz);
        tSim->SetBranchAddress("ePS", &ePS);
        tSim->SetBranchAddress("eSH", &eSH);
        tSim->SetBranchAddress("trP", &trP);
        tSim->SetBranchAddress("dx", &dx);
        tSim->SetBranchAddress("dy", &dy);
        tSim->SetBranchAddress("W2", &W2);
        TBranch* bCT = tSim->GetBranch("coin_time");
        if(bCT) tSim->SetBranchAddress("coin_time", &ctSim);
        tSim->SetBranchAddress("eHCAL", &eHCAL);
        tSim->SetBranchAddress("weight", &wt);
        tSim->SetBranchAddress("fnucl", &fnucl);
        Long64_t nSim = tSim->GetEntries();
        ev.dx.reserve(ev.dx.size()+nSim);
        for(Long64_t i=0; i<nSim; ++i) {
            tSim->GetEntry(i);
            // only accept fnucl 0 or 1
            if(fnucl != 0.0 && fnucl != 1.0) continue;
            if(abs(vz)>0.27 && ePS<0.2 && abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL<0.025) continue; // you can add basic cuts here if needed
            ev.vz.push_back(vz);
            ev.ePS.push_back(ePS);
            ev.dx.push_back(dx);
            ev.dy.push_back(dy);
            ev.W2.push_back(W2);
            ev.ct.push_back(bCT? ctSim : 0.0);
            ev.eHCAL.push_back(eHCAL);
            ev.helicity.push_back(0);
            ev.IHWP.push_back(0);
            ev.w.push_back(wt);
            ev.origin.push_back(fnucl==0.0? 0 : 1);
            if (i %1000 == 0 ) std::cout<< (i * 100.0/ nSim)<<" sim file" << "% \r";
            std::cout.flush();

        }
    }

    std::cout << "Loaded " << ev.dx.size() << " total events (" 
              << "data+sim)" << std::endl;

   // 3. Scale factors
    double sumD=0, sumN=0, sumP=0; long nN=0, nP=0;
    for(size_t i=0;i<ev.dx.size();++i){
        int o=ev.origin[i]; double wt=ev.w[i];
        if(o==2) sumD+=wt;
        else if(o==0){ sumN+=wt; ++nN; }
        else          { sumP+=wt; ++nP; }
    }
    sumN=std::max(sumN,1e-9); sumP=std::max(sumP,1e-9);
    double wN=(wN_user>0? wN_user: sumD/sumN);
    double wP=(wP_user>0? wP_user: sumD/sumP);
    std::cout<<"MC neutrons="<<nN<<" protons="<<nP<<"\n";
    std::cout<<"wN="<<wN<<" wP="<<wP<<"\n";

    // 4. Define grid ranges and steps
    const int nSteps=5; // adjust grid resolution 
    std::vector<double> dyL_vals, dyH_vals,
                        W2L_vals, W2H_vals, eL_vals,
                        tL_vals, tH_vals;
    // example: dxL from -1 to 0
    //for(int i=0;i<nSteps;++i) dxL_vals.push_back(-1.0 + i*(1.0/nSteps));
    //for(int i=0;i<nSteps;++i) dxH_vals.push_back( 0.0 + i*(1.0/nSteps));
    for(int i=0;i<nSteps-1;++i) dyL_vals.push_back(-0.5 + i*(0.5/nSteps));
    for(int i=0;i<nSteps-1;++i) dyH_vals.push_back( 0.5 - i*(0.5/nSteps));
    for(int i=0;i<nSteps;++i) W2L_vals.push_back(-2.0 + i*(1.0/nSteps));
    for(int i=0;i<nSteps;++i) W2H_vals.push_back( 1.4 + i*(0.5/nSteps));
    for(int i=0;i<nSteps;++i) eL_vals.push_back(0.025 + i*(0.5/nSteps));
    //for(int i=0;i<nSteps;++i) eH_vals.push_back(0.5 + i*(2.0/nSteps));
    for(int i=1;i<4;++i) tL_vals.push_back(120 - i*(6/3));
    for(int i=1;i<4;++i) tH_vals.push_back(120 + i*(6/3));

    // 5. Loop over grid, compute FOM
    double bestFOM = 0.0001;
    std::array<double,10> bestCuts{};
    size_t total = dyL_vals.size() * dyH_vals.size() * W2L_vals.size()* W2H_vals.size() * eL_vals.size() * tL_vals.size() * tH_vals.size();
    size_t count = 0;


    //for(double dxL:dxL_vals) for(double dxH:dxH_vals)
    for(double dyL:dyL_vals) for(double dyH:dyH_vals)
    for(double W2L:W2L_vals) for(double W2H:W2H_vals)
    for(double eL:eL_vals)   //for(double eH:eH_vals)
    for(double tL:tL_vals)   for(double tH:tH_vals) {
        // build histograms
        
        if(dyL==dyH || W2L==W2H || tL==tH || abs(tL-tH)<3){
            continue;
        }

        TH1D *h_bkg = new TH1D("b",";dx;",nbinsDx,-4,3);
        TH1D *h_data = new TH1D("d",";dx;",nbinsDx,-4,3);
        TH1D *h_n = new TH1D("n",";dx;",nbinsDx,-4,3);
        TH1D *h_p = new TH1D("p",";dx;",nbinsDx,-4,3);

        double Nplus = 0.0;
        double Nminus = 0.0;
        // fill
        for(size_t i=0;i<ev.dx.size();++i){
            bool base = //ev.dx[i]>dxL && ev.dx[i]<dxH &&
                        abs(ev.vz[i])<0.27&&
                        ev.ePS[i]>0.2&&
                        ev.W2[i]>W2L && ev.W2[i]<W2H &&
                        ev.eHCAL[i]>eL;
            bool inct = ev.ct[i]>tL    && ev.ct[i]<tH;
            bool inDy = ev.dy[i]>dyL&&ev.dy[i]<dyH;
            bool antiDy = ev.dy[i]>1.0||ev.dy[i]<-1.1;
            bool inDx = ev.dx[i]>-0.4 && ev.dx[i]<0.4;
            bool poshel = ev.IHWP[i]*ev.helicity[i] == 1;
            bool neghel = ev.IHWP[i]*ev.helicity[i] == -1;

            int o=ev.origin[i]; double wt=ev.w[i];
            if(o==2){ if(base && inct && antiDy) h_bkg->Fill(ev.dx[i]);
                      if(base && inct && inDy)  h_data->Fill(ev.dx[i]);
                      if(base && inct && inDy && inDx && poshel) ++Nplus;
                      if(base && inct && inDy && inDx && neghel) ++Nminus;
                    }
            if(o==0 && base&&inDy) h_n->Fill(ev.dx[i], wt);
            if(o==1 && base&&inDy) h_p->Fill((ev.dx[i]/*+0.05*/), wt);
        }
        // normalize
        double scale_data = h_data->Integral();

        //std::cout<<"scale_data : "<<scale_data<<endl;
        //std::cout<<"n : "<<h_n->Integral()<<endl;
        //std::cout<<"p : "<<h_p->Integral()<<endl;
        //std::cout<<"bkg : "<<h_bkg->Integral()<<endl;


        if(h_bkg->Integral()>0) h_bkg->Scale(1.0/h_bkg->Integral());
        if(h_n->Integral()>0)   h_n->Scale(1.0/h_n->Integral());
        if(h_p->Integral()>0)   h_p->Scale(1.0/h_p->Integral());
        if(h_data->Integral()>0) h_data->Scale(1.0/h_data->Integral());

        hist_bkg = h_bkg;
        hist_p = h_p;
        hist_n = h_n;
        // compose
        //TH1D h_model = h_n; h_model.Add(&h_p); h_model.Add(&h_bkg);
        
        //std::cout<<"data : "<<h_data->Integral()<<endl;
        //std::cout<<"n : "<<h_n->Integral()<<endl;
        //std::cout<<"p : "<<h_p->Integral()<<endl;
        //std::cout<<"bkg : "<<h_bkg->Integral()<<endl;

        TF1 *fit_data = new TF1("fit_data", fit_sim_n_bkg, h_data->GetXaxis()->GetXmin(), h_data->GetXaxis()->GetXmax(),3);

        // Set better initial guesses
        double init_param0 = 1;
        double init_param1 = 1;
        double init_param2 = 1;

        //fit_data->SetParameters(init_param0, init_param1, init_param2);

        // Apply constraints
        fit_data->SetParLimits(0, 0.1, 100);
        fit_data->SetParLimits(1, 0.1, 100);
        fit_data->SetParLimits(2, 0.0, 100);

        // Fit the data with log-likelihood
        h_data->Fit(fit_data, "RLQ");

        hist_p->Scale(fit_data->GetParameter(0));
        hist_n->Scale(fit_data->GetParameter(1)*fit_data->GetParameter(0));
        hist_bkg->Scale(fit_data->GetParameter(2)*fit_data->GetParameter(0));

        hist_p->Scale(scale_data);
        hist_n->Scale(scale_data);
        hist_bkg->Scale(scale_data);
        h_data->Scale(scale_data);

        double nN = hist_n->Integral(hist_n->FindBin(-0.4),hist_n->FindBin(0.4));
        double nData = h_data->Integral(h_data->FindBin(-0.4),h_data->FindBin(0.4));
        double nBkg = hist_bkg->Integral(hist_bkg->FindBin(-0.4),hist_bkg->FindBin(0.4));

        // compute FOM
        double FOM = nN/sqrt(nN+nBkg);

        double FOM_1 = nN/(nData*sqrt(nData));

        double FOM_2 = nN/sqrt(Nplus*Nminus*(Nplus+Nminus)); 

        //std::cout<<"FOM : "<<FOM<<endl;
        //std::cout<<"dy_H :"<<dyH<<endl;
        outfile_all<<"FOM : "<<FOM<<" FOM_1 : "<<FOM_1<<" FOM_2 : "<<FOM_2<<" {dyL,dyH,W2L,W2H,eL,tL,tH}" <<"{"<<dyL<<","<<dyH<<","
        <<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<endl;


        //if (count % 10 == 0) {
        std::cout<< "Scan fit " << count << " / " << total<< " : " << (100.0*count/total) <<endl;
            //<< std::flush;
        //}

        count++;   

        if(FOM > bestFOM && 
            (h_data->Integral(h_data->FindBin(-0.4),h_data->FindBin(0.4)) > 
            hist_n->Integral(hist_n->FindBin(-0.4),hist_n->FindBin(0.4))) &&
            (h_data->Integral(h_data->FindBin(-1.2),h_data->FindBin(-0.8)) > 
            hist_p->Integral(hist_p->FindBin(-1.2),hist_n->FindBin(-0.8))) //some fail safes
            ) {
            bestFOM = FOM;
            bestCuts = {dyL,dyH,W2L,W2H,eL,tL,tH};
            //std::cout<<"NEW BEST FOM : "<<" {dyL,dyH,W2L,W2H,eL,tL,tH}" <<"{"<<dyL<<","<<dyH<<","
            //<<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<endl;

            outfile<<"NEW BEST FOM : "<<bestFOM<<" FOM_1 : "<<FOM_1<<" FOM_2 : "<<FOM_2<<" {dyL,dyH,W2L,W2H,eL,tL,tH}" <<"{"<<dyL<<","<<dyH<<","
            <<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<endl;

            hist_p->SetFillColorAlpha(6,0.5);
            hist_p->SetFillStyle(3004);
            hist_n->SetFillColorAlpha(9,0.5);
            hist_n->SetFillStyle(3005);
            hist_bkg->SetFillColorAlpha(7,0.5);
            hist_bkg->SetFillStyle(3003);
            h_data->SetMarkerStyle(20);
            h_data->SetMarkerColor(kBlack);

            TCanvas* cfit = new TCanvas("cfit","Fit with models",3600,3000);
            cfit->Divide(2,2);
            cfit->cd(1);
            //h_dx_sim_n_bkg->SetXTitle("HCAL_X(exp)-HCAL_X(act) (m)");
            //h_dx_sim_n_bkg->Draw("HIST");
            h_data->SetXTitle("HCAL_X(exp)-HCAL_X(act) (m)");
            h_data->Draw("P");
            hist_p->Draw("HIST SAME");
            hist_n->Draw("HIST SAME");
            hist_bkg->Draw("HIST SAME");

            cfit->SaveAs(Form("plots/FOM_dx_test_%s_%f.jpg",kin,FOM));

        }
    }

    // 6. Report best
    std::cout<<"Best FOM="<<bestFOM<<" with cuts: \n";
    //std::cout<<" dx ∈ ["<<bestCuts[0]<<","<<bestCuts[1]<<"]\n";
    std::cout<<" dy ∈ ["<<bestCuts[0]<<","<<bestCuts[1]<<"]\n";
    std::cout<<" W2 ∈ ["<<bestCuts[2]<<","<<bestCuts[3]<<"]\n";
    std::cout<<" eHCAL_L : ["<<bestCuts[4]<<"]\n";
    std::cout<<" coin_time ∈ ["<<bestCuts[5]<<","<<bestCuts[6]<<"]\n";
}
