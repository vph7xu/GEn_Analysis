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


struct Store {
    std::vector<double> dx, dy, W2, ct, eHCAL, w;
    std::vector<char>   origin;  // 0=nMC,1=pMC,2=data
};

void OptimizeCuts(const char* dataFile,
                  const char* simFile,
                  int    nbinsDx = 120,
                  double wN_user = -1.0,
                  double wP_user = -1.0)
{
    Store ev;

    std::ofstream outfile("test_FOM.txt");

    gStyle->SetOptFit(0);

    // 1. Load DATA
    {
        TFile fData(dataFile, "READ");
        TTree* tData = (TTree*)fData.Get("Tout");
        double dx, dy, W2, ct, eHCAL;
        tData->SetBranchAddress("dx", &dx);
        tData->SetBranchAddress("dy", &dy);
        tData->SetBranchAddress("W2", &W2);
        tData->SetBranchAddress("coin_time", &ct);
        tData->SetBranchAddress("eHCAL", &eHCAL);
        Long64_t nData = tData->GetEntries();
        ev.dx.reserve(nData);
        ev.dy.reserve(nData);
        ev.W2.reserve(nData);
        ev.ct.reserve(nData);
        ev.eHCAL.reserve(nData);
        ev.w.reserve(nData);
        ev.origin.reserve(nData);
        for(Long64_t i=0; i<0.1*nData; ++i) {
            tData->GetEntry(i);
            ev.dx.push_back(dx);
            ev.dy.push_back(dy);
            ev.W2.push_back(W2);
            ev.ct.push_back(ct);
            ev.eHCAL.push_back(eHCAL);
            ev.w.push_back(1.0);
            ev.origin.push_back(2);
            if (i %1000 == 0 ) std::cout<< (i * 100.0/ nData)<<" data file" << "% \r";
            std::cout.flush();
        }
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
        double dx, dy, W2, ctSim, eHCAL, wt, fnucl;
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
            ev.dx.push_back(dx);
            ev.dy.push_back(dy);
            ev.W2.push_back(W2);
            ev.ct.push_back(bCT? ctSim : 0.0);
            ev.eHCAL.push_back(eHCAL);
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
    const int nSteps=8; // adjust grid resolution 
    std::vector<double> dyL_vals, dyH_vals,
                        W2L_vals, W2H_vals, eL_vals,
                        tL_vals, tH_vals;
    // example: dxL from -1 to 0
    //for(int i=0;i<nSteps;++i) dxL_vals.push_back(-1.0 + i*(1.0/nSteps));
    //for(int i=0;i<nSteps;++i) dxH_vals.push_back( 0.0 + i*(1.0/nSteps));
    for(int i=0;i<nSteps;++i) dyL_vals.push_back(-0.5 + i*(1.0/nSteps));
    for(int i=0;i<nSteps;++i) dyH_vals.push_back( 0.5 - i*(1.0/nSteps));
    for(int i=0;i<nSteps;++i) W2L_vals.push_back(-1.5 + i*(2.0/nSteps));
    for(int i=0;i<nSteps;++i) W2H_vals.push_back( 1.0 + i*(3.0/nSteps));
    for(int i=0;i<nSteps;++i) eL_vals.push_back(0.025 + i*(0.2/nSteps));
    //for(int i=0;i<nSteps;++i) eH_vals.push_back(0.5 + i*(2.0/nSteps));
    for(int i=0;i<nSteps;++i) tL_vals.push_back(120 - i*(8/nSteps));
    for(int i=0;i<nSteps;++i) tH_vals.push_back(120 + i*(8/nSteps));

    // 5. Loop over grid, compute FOM
    double bestFOM = 0.1;
    std::array<double,10> bestCuts{};
    size_t total = dyL_vals.size() * dyH_vals.size() * W2L_vals.size()* W2H_vals.size() * eL_vals.size() * tL_vals.size() * tH_vals.size();
    size_t count = 0;


    //for(double dxL:dxL_vals) for(double dxH:dxH_vals)
    for(double dyL:dyL_vals) for(double dyH:dyH_vals)
    for(double W2L:W2L_vals) for(double W2H:W2H_vals)
    for(double eL:eL_vals)   //for(double eH:eH_vals)
    for(double tL:tL_vals)   for(double tH:tH_vals) {
        // build histograms
        TH1D *h_bkg = new TH1D("b",";dx;",nbinsDx,-6,6);
        TH1D *h_data = new TH1D("d",";dx;",nbinsDx,-6,6);
        TH1D *h_n = new TH1D("n",";dx;",nbinsDx,-6,6);
        TH1D *h_p = new TH1D("p",";dx;",nbinsDx,-6,6);
        // fill
        for(size_t i=0;i<ev.dx.size();++i){
            bool base = //ev.dx[i]>dxL && ev.dx[i]<dxH &&
                        ev.W2[i]>W2L && ev.W2[i]<W2H &&
                        ev.eHCAL[i]>eL;
            bool inct = ev.ct[i]>tL    && ev.ct[i]<tH;
            bool inDy = ev.dy[i]>dyL&&ev.dy[i]<dyH;
            bool antiDy = ev.dy[i]>1.2||ev.dy[i]<-1.5;
            int o=ev.origin[i]; double wt=ev.w[i];
            if(o==2){ if(base && inct && antiDy) h_bkg->Fill(ev.dx[i]);
                      if(base && inct && inDy)  h_data->Fill(ev.dx[i]); }
            if(o==0 && base&&inDy) h_n->Fill(ev.dx[i], wt);
            if(o==1 && base&&inDy) h_p->Fill(ev.dx[i], wt);
        }
        // normalize
        double scale_data = h_data->Integral();

        std::cout<<"scale_data : "<<scale_data<<endl;
        std::cout<<"n : "<<h_n->Integral()<<endl;
        std::cout<<"p : "<<h_p->Integral()<<endl;
        std::cout<<"bkg : "<<h_bkg->Integral()<<endl;


        if(h_bkg->Integral()>0) h_bkg->Scale(1.0/h_bkg->Integral());
        if(h_n->Integral()>0)   h_n->Scale(1.0/h_n->Integral());
        if(h_p->Integral()>0)   h_p->Scale(1.0/h_p->Integral());
        if(h_data->Integral()>0) h_data->Scale(1.0/h_data->Integral());

        hist_bkg = h_bkg;
        hist_p = h_p;
        hist_n = h_n;
        // compose
        //TH1D h_model = h_n; h_model.Add(&h_p); h_model.Add(&h_bkg);
        
        std::cout<<"data : "<<h_data->Integral()<<endl;
        std::cout<<"n : "<<h_n->Integral()<<endl;
        std::cout<<"p : "<<h_p->Integral()<<endl;
        std::cout<<"bkg : "<<h_bkg->Integral()<<endl;

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

        // compute FOM
        double FOM = hist_p->Integral()/sqrt(hist_n->Integral()+hist_bkg->Integral());

        std::cout<<"FOM : "<<FOM<<endl;
        //std::cout<<"dy_H :"<<dyH<<endl;
        std::cout<<"{dyL,dyH,W2L,W2H,eL,tL,tH}" <<"{"<<dyL<<","<<dyH<<","
        <<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<endl;


        if (count % 50000 == 0) {
            std::cout
            << "Scan fit " << count << " / " << total
            << "  (" << (100.0*count/total) << "%)\r"
            << std::flush;
        }

        count++;   

        if(FOM > bestFOM) {
            bestFOM = FOM;
            bestCuts = {dyL,dyH,W2L,W2H,eL,tL,tH};
            std::cout<<"NEW BEST FOM : "<<" {dyL,dyH,W2L,W2H,eL,tL,tH}" <<"{"<<dyL<<","<<dyH<<","
            <<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<endl;

            outfile<<"NEW BEST FOM : "<<bestFOM<<" {dyL,dyH,W2L,W2H,eL,tL,tH}" <<"{"<<dyL<<","<<dyH<<","
            <<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<endl;

        }
    }

    // 6. Report best
    std::cout<<"Best χ²="<<bestFOM<<" with cuts: \n";
    std::cout<<" dx ∈ ["<<bestCuts[0]<<","<<bestCuts[1]<<"]\n";
    std::cout<<" dy ∈ ["<<bestCuts[2]<<","<<bestCuts[3]<<"]\n";
    std::cout<<" W2 ∈ ["<<bestCuts[4]<<","<<bestCuts[5]<<"]\n";
    std::cout<<" eHCAL ∈ ["<<bestCuts[6]<<","<<bestCuts[7]<<"]\n";
    std::cout<<" coin_time ∈ ["<<bestCuts[8]<<","<<bestCuts[9]<<"]\n";
}
