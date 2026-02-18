// OptimizeCuts.C — complete ROOT macro for cut optimization
// ------------------------------------------------------------
// * Loads data & simulation trees via TTree::GetEntry loops
// * Auto-scales simulated neutrons & protons to data (not currently used in fit)
// * For each trial rectangle, builds background, data, sim histograms
// * Composes model = simN + simP + background
// * Fits data with TF1 (fit_sim_n_bkg) to extract relative fractions
// * Computes FOM = Nn / sqrt(Nn + Nbkg) in dx window
// * Scans grid over dy, W2, eHCAL, coin_time
// * Background now taken from a SECOND sim file (bkgSimFile)
// ------------------------------------------------------------

#include "models.h"  // expects fit_sim_n_bkg, hist_bkg, hist_p, hist_n somewhere

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <Math/Functor.h>
#include <Minuit2/Minuit2Minimizer.h>
#include <vector>
#include <array>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <TStyle.h>  // for extern TStyle *gStyle
#include <TROOT.h>
#include <TCanvas.h>

struct Store {
    std::vector<double> vz, ePS, dx, dy, W2, ct, eHCAL, w;
    std::vector<int> helicity, IHWP;
    std::vector<char> origin;  // 0=nMC,1=pMC,2=data
};

struct StoreBkg {
    std::vector<double> vz, ePS, dx, dy, W2, ct, eHCAL, w;
};

void OptimizeCuts(const char* dataFile,
                  const char* simFile,      // neutron/proton sim (fnucl==0,1)
                  const char* bkgSimFile,   // background sim
                  const char* kin,
                  int    nbinsDx = 100,
                  double wN_user = -1.0,
                  double wP_user = -1.0)
{
    Store    ev;
    StoreBkg ev_bkg;

    std::ofstream outfile(Form("test_FOM_%s.txt",kin));
    std::ofstream outfile_all(Form("test_FOM_all_%s.txt",kin));

    gStyle->SetOptFit(0);

    // ------------------------------------------------------------------------
    // 1. Load DATA
    // ------------------------------------------------------------------------
    {
        TFile fData(dataFile, "READ");
        TTree* tData = (TTree*)fData.Get("Tout");
        if (!tData) {
            std::cerr << "ERROR: cannot find Tout in data file " << dataFile << std::endl;
            return;
        }

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
        ev.vz.reserve(nData+100000000);
        ev.ePS.reserve(nData+100000000);
        ev.dx.reserve(nData+100000000);
        ev.dy.reserve(nData+100000000);
        ev.W2.reserve(nData+100000000);
        ev.ct.reserve(nData+100000000);
        ev.eHCAL.reserve(nData+100000000);
        ev.helicity.reserve(nData+100000000);
        ev.IHWP.reserve(nData+100000000);
        ev.w.reserve(nData+100000000);
        ev.origin.reserve(nData+100000000);
        for(Long64_t i=0; i<nData*0.3; ++i) {

            tData->GetEntry(i);
            // NOTE: as written, this condition only rejects events where ALL
            // sub-conditions are true simultaneously. Keeping as-is to match your code.
            if(std::abs(vz)>0.27 || ePS<0.2 ||
               std::abs((eSH+ePS)/trP - 1.0) > 0.2 || eHCAL<0.025 ||
               std::abs(helicity)!=1)
                continue;

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

            if (i % 1000 == 0) {
                std::cout<< (i * 100.0/ nData)<<" data file" << "% \r";
                std::cout.flush();
            }
        }
        std::cout << "\nLoaded " << ev.dx.size() << " total events (data only)" << std::endl;
        std::cout << "Good total events in data: " << good_total_events << std::endl;
    }

    // ------------------------------------------------------------------------
    // 2. Load SIM (neutrons + protons)
    // ------------------------------------------------------------------------
    bool simHasTime = false;
    {
        TFile fSim(simFile, "READ");
        TTree* t = (TTree*)fSim.Get("Tout");
        if (!t) {
            std::cerr << "ERROR: cannot find Tout in sim file " << simFile << std::endl;
            return;
        }
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
            // only accept fnucl 0 or 1 (n or p)
            if(fnucl != 0.0 && fnucl != 1.0) continue;
            if(std::abs(vz)>0.27 || ePS<0.2 ||
               std::abs((eSH+ePS)/trP - 1.0) > 0.2 || eHCAL<0.025)
                continue;

            ev.vz.push_back(vz);
            ev.ePS.push_back(ePS);
            ev.dx.push_back(dx);
            ev.dy.push_back(dy);
            ev.W2.push_back(W2);
            ev.ct.push_back(bCT ? ctSim : 0.0);
            ev.eHCAL.push_back(eHCAL);
            ev.helicity.push_back(0);
            ev.IHWP.push_back(0);
            ev.w.push_back(wt);
            ev.origin.push_back(fnucl==0.0? 0 : 1);

            if (i % 1000 == 0) {
                std::cout<< (i * 100.0/ nSim)<<" sim file" << "% \r";
                std::cout.flush();
            }
        }
        std::cout << "\nLoaded " << ev.dx.size() << " total events (data+np-sim)" << std::endl;
    }

    // ------------------------------------------------------------------------
    // 2b. Load BACKGROUND SIM (for h_bkg instead of antiDy from data)
    // ------------------------------------------------------------------------
    bool bkgHasTime = false;
    {
        TFile fBkg(bkgSimFile, "READ");
        TTree* tB = (TTree*)fBkg.Get("Tout");
        if (!tB) {
            std::cerr << "ERROR: cannot find Tout in background sim file "
                      << bkgSimFile << std::endl;
            return;
        }
        if (tB->GetBranch("coin_time")) bkgHasTime = true;
    }
    {
        TFile fBkg(bkgSimFile, "READ");
        TTree* tBkg = (TTree*)fBkg.Get("Tout");
        if (!tBkg) {
            std::cerr << "ERROR: cannot open background TTree in "
                      << bkgSimFile << std::endl;
            return;
        }

        double vz, ePS, eSH, trP, dx, dy, W2, ctSim, eHCAL, wt;
        tBkg->SetBranchAddress("vz", &vz);
        tBkg->SetBranchAddress("ePS", &ePS);
        tBkg->SetBranchAddress("eSH", &eSH);
        tBkg->SetBranchAddress("trP", &trP);
        tBkg->SetBranchAddress("dx", &dx);
        tBkg->SetBranchAddress("dy", &dy);
        tBkg->SetBranchAddress("W2", &W2);
        TBranch* bCTb = tBkg->GetBranch("coin_time");
        if (bCTb) tBkg->SetBranchAddress("coin_time", &ctSim);
        tBkg->SetBranchAddress("eHCAL", &eHCAL);
        tBkg->SetBranchAddress("weight", &wt);

        Long64_t nBkg = tBkg->GetEntries();
        ev_bkg.dx.reserve(nBkg);

        for (Long64_t i = 0; i < nBkg*0.3; ++i) {
            tBkg->GetEntry(i);

            // same basic preselection form as above (kept identical to your logic)
            if (std::abs(vz)>0.27 || ePS<0.2 ||
                std::abs((eSH+ePS)/trP - 1.0) > 0.2 || eHCAL<0.025)
                continue;

            ev_bkg.vz.push_back(vz);
            ev_bkg.ePS.push_back(ePS);
            ev_bkg.dx.push_back(dx);
            ev_bkg.dy.push_back(dy);
            ev_bkg.W2.push_back(W2);
            ev_bkg.ct.push_back(bkgHasTime ? ctSim : 0.0);
            ev_bkg.eHCAL.push_back(eHCAL);
            ev_bkg.w.push_back(wt);

            if (i % 1000 == 0) {
                std::cout << (i * 100.0 / nBkg) << " bkg sim file" << "% \r";
                std::cout.flush();
            }
        }
        std::cout << "\nLoaded " << ev_bkg.dx.size()
                  << " background-sim events" << std::endl;
    }

    // ------------------------------------------------------------------------
    // 3. Scale factors (still computed, but current fit uses shape-only templates)
    // ------------------------------------------------------------------------
    double sumD=0, sumN=0, sumP=0; long nN=0, nP=0;
    for(size_t i=0;i<ev.dx.size();++i){
        int o=ev.origin[i]; double wt=ev.w[i];
        if(o==2) sumD+=wt;
        else if(o==0){ sumN+=wt; ++nN; }
        else         { sumP+=wt; ++nP; }
    }
    sumN=std::max(sumN,1e-9); sumP=std::max(sumP,1e-9);
    double wN=(wN_user>0? wN_user: sumD/sumN);
    double wP=(wP_user>0? wP_user: sumD/sumP);
    std::cout<<"MC neutrons="<<nN<<" protons="<<nP<<"\n";
    std::cout<<"wN="<<wN<<" wP="<<wP<<"\n";

    // ------------------------------------------------------------------------
    // 4. Define grid ranges and steps
    // ------------------------------------------------------------------------
    const int nSteps=5; // adjust grid resolution
    std::vector<double> dyL_vals, dyH_vals,
                        W2L_vals, W2H_vals, eL_vals,
                        tL_vals, tH_vals;

    for(int i=0;i<nSteps-1;++i) dyL_vals.push_back(-0.5 + i*(0.5/nSteps));
    for(int i=0;i<nSteps-1;++i) dyH_vals.push_back( 0.5 - i*(0.5/nSteps));
    for(int i=0;i<nSteps;++i)   W2L_vals.push_back(-2.0 + i*(1.0/nSteps));
    for(int i=0;i<nSteps;++i)   W2H_vals.push_back( 1.4 + i*(0.5/nSteps));
    for(int i=0;i<nSteps;++i)   eL_vals.push_back(0.025 + i*(0.5/nSteps));
    for(int i=1;i<4;++i)        tL_vals.push_back(184 - i*(6/3)); // 118,116,114
    for(int i=1;i<4;++i)        tH_vals.push_back(184 + i*(6/3)); // 122,124,126

    // ------------------------------------------------------------------------
    // 5. Loop over grid, compute FOM
    // ------------------------------------------------------------------------
    double bestFOM = 0.0001;
    std::array<double,10> bestCuts{};
    size_t total = dyL_vals.size() * dyH_vals.size() *
                   W2L_vals.size()* W2H_vals.size() *
                   eL_vals.size() * tL_vals.size() * tH_vals.size();
    size_t count = 0;

    for(double dyL:dyL_vals) for(double dyH:dyH_vals)
    for(double W2L:W2L_vals) for(double W2H:W2H_vals)
    for(double eL:eL_vals)
    for(double tL:tL_vals)   for(double tH:tH_vals) {

        if(dyL==dyH || W2L==W2H || tL==tH || std::abs(tL-tH)<3){
            continue;
        }

        // build histograms (unique names to avoid ROOT warnings)
        TH1D *h_bkg  = new TH1D(Form("b_%zu",count),";dx;",nbinsDx,-4,3);
        TH1D *h_data = new TH1D(Form("d_%zu",count),";dx;",nbinsDx,-4,3);
        TH1D *h_n    = new TH1D(Form("n_%zu",count),";dx;",nbinsDx,-4,3);
        TH1D *h_p    = new TH1D(Form("p_%zu",count),";dx;",nbinsDx,-4,3);

        double Nplus = 0.0;
        double Nminus = 0.0;

        // fill from data + n/p sim
        for(size_t i=0;i<ev.dx.size();++i){
            bool base = std::abs(ev.vz[i])<0.27 &&
                        ev.ePS[i]>0.2 &&
                        ev.W2[i]>W2L && ev.W2[i]<W2H &&
                        ev.eHCAL[i]>eL;
            bool inct = ev.ct[i]>tL    && ev.ct[i]<tH;
            bool inDy = ev.dy[i]>dyL   && ev.dy[i]<dyH;
            bool inDx = ev.dx[i]>-0.4 && ev.dx[i]<0.4; // NOTE: fixed dx typo

            bool poshel = (ev.IHWP[i]*ev.helicity[i] == 1);
            bool neghel = (ev.IHWP[i]*ev.helicity[i] == -1);

            int o=ev.origin[i]; double wt=ev.w[i];

            if(o==2){
                if(base && inct && inDy)
                    h_data->Fill(ev.dx[i]);
                if(base && inct && inDy && inDx && poshel) ++Nplus;
                if(base && inct && inDy && inDx && neghel) ++Nminus;
            }
            if(o==0 && base && inDy) h_n->Fill(ev.dx[i], wt);
            if(o==1 && base && inDy) h_p->Fill(ev.dx[i], wt);
        }

        // fill background from background sim
        for (size_t j = 0; j < ev_bkg.dx.size(); ++j) {
            bool baseB = std::abs(ev_bkg.vz[j]) < 0.27 &&
                         ev_bkg.ePS[j] > 0.2 &&
                         ev_bkg.W2[j] > W2L && ev_bkg.W2[j] < W2H &&
                         ev_bkg.eHCAL[j] > eL;
            bool inctB = ev_bkg.ct[j] > tL && ev_bkg.ct[j] < tH;
            bool inDyB = ev_bkg.dy[j] > dyL && ev_bkg.dy[j] < dyH;

            if (baseB /*&& inctB*/ && inDyB) {
                h_bkg->Fill(ev_bkg.dx[j], ev_bkg.w[j]);

                //cout<< "Filling bkg with dx: " << ev_bkg.dx[j] << " and weight: " << ev_bkg.w[j] << std::endl;
            }
        }

        // normalize templates to 1 before fitting
        double scale_data = h_data->Integral();

        if(h_bkg->Integral()>0)   h_bkg->Scale(1.0/h_bkg->Integral());
        if(h_n->Integral()>0)     h_n->Scale(1.0/h_n->Integral());
        if(h_p->Integral()>0)     h_p->Scale(1.0/h_p->Integral());
        if(h_data->Integral()>0)  h_data->Scale(1.0/h_data->Integral());

        hist_bkg = h_bkg;
        hist_p   = h_p;
        hist_n   = h_n;

        TF1 *fit_data = new TF1(Form("fit_data_%zu",count),
                                fit_sim_n_bkg,
                                h_data->GetXaxis()->GetXmin(),
                                h_data->GetXaxis()->GetXmax(),3);

        // optional initial guesses (not set right now)
        // double init_param0 = 1;
        // double init_param1 = 1;
        // double init_param2 = 1;
        // fit_data->SetParameters(init_param0, init_param1, init_param2);

        fit_data->SetParLimits(0, 0.1, 100);
        fit_data->SetParLimits(1, 0.1, 100);
        fit_data->SetParLimits(2, 0.0, 100);

        h_data->Fit(fit_data, "RLQ"); // log-likelihood fit

        hist_p->Scale(fit_data->GetParameter(0));
        hist_n->Scale(fit_data->GetParameter(1)*fit_data->GetParameter(0));
        hist_bkg->Scale(fit_data->GetParameter(2)*fit_data->GetParameter(0));

        hist_p->Scale(scale_data);
        hist_n->Scale(scale_data);
        hist_bkg->Scale(scale_data);
        h_data->Scale(scale_data);

        double nN    = hist_n->Integral(hist_n->FindBin(-0.4),hist_n->FindBin(0.4));
        double nData = h_data->Integral(h_data->FindBin(-0.4),h_data->FindBin(0.4));
        double nBkg  = hist_bkg->Integral(hist_bkg->FindBin(-0.4),hist_bkg->FindBin(0.4));

        // compute FOMs
        double FOM   = (nN>0.0 ? nN/std::sqrt(nN+nBkg) : 0.0);
        double FOM_1 = (nData>0.0 ? nN/(nData*std::sqrt(nData)) : 0.0);
        double FOM_2 = (Nplus>0.0 && Nminus>0.0 ?
                        nN/std::sqrt(Nplus*Nminus*(Nplus+Nminus)) : 0.0);

        outfile_all<<"FOM : "<<FOM<<" FOM_1 : "<<FOM_1<<" FOM_2 : "<<FOM_2
                   <<" {dyL,dyH,W2L,W2H,eL,tL,tH} "
                   <<"{"<<dyL<<","<<dyH<<","
                   <<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<std::endl;

        std::cout<< "Scan fit " << count << " / " << total
                 << " : " << (100.0*count/total) << std::endl;
        count++;

        // some fail-safes as in your original code
        bool okRegion =
            (h_data->Integral(h_data->FindBin(-0.4),h_data->FindBin(0.4)) >
             hist_n->Integral(hist_n->FindBin(-0.4),hist_n->FindBin(0.4))) &&
            (h_data->Integral(h_data->FindBin(-1.7),h_data->FindBin(-1.3)) >
             hist_p->Integral(hist_p->FindBin(-1.7),hist_p->FindBin(-1.3)));

            if(FOM > bestFOM && okRegion) {
                bestFOM = FOM;
                bestCuts = {dyL,dyH,W2L,W2H,eL,tL,tH};

                outfile<<"NEW BEST FOM : "<<bestFOM
                    <<" FOM_1 : "<<FOM_1<<" FOM_2 : "<<FOM_2
                    <<" {dyL,dyH,W2L,W2H,eL,tL,tH} "
                    <<"{"<<dyL<<","<<dyH<<","
                    <<W2L<<","<<W2H<<","<<eL<<","<<tL<<","<<tH<<"}"<<std::endl;

                // Style components
                hist_p->SetFillColorAlpha(6,0.5);
                hist_p->SetFillStyle(3004);
                hist_n->SetFillColorAlpha(9,0.5);
                hist_n->SetFillStyle(3005);
                hist_bkg->SetFillColorAlpha(7,0.5);
                hist_bkg->SetFillStyle(3003);
                h_data->SetMarkerStyle(20);
                h_data->SetMarkerColor(kBlack);

                // ---- NEW: build total model histogram = n + p + bkg ----
                TH1D *h_model = (TH1D*)hist_n->Clone(Form("h_model_best_%s", kin));
                h_model->SetTitle("Total model");
                h_model->Add(hist_p);
                h_model->Add(hist_bkg);
                h_model->SetLineColor(kBlue+2);
                h_model->SetLineWidth(3);
                h_model->SetFillStyle(0);    // no fill, just a line

                // Canvas and drawing
                TCanvas* cfit = new TCanvas("cfit","Fit with models",3600,3000);
                cfit->Divide(2,2);

                // Pad 1: data + components + total model
                cfit->cd(1);
                h_data->SetXTitle("HCAL_X(exp)-HCAL_X(act) (m)");
                h_data->Draw("P");
                hist_p->Draw("HIST SAME");
                hist_n->Draw("HIST SAME");
                hist_bkg->Draw("HIST SAME");
                h_model->Draw("HIST SAME");   // <-- total fit histogram

                // Pad 2 (optional): just the model histogram
                cfit->cd(2);
                h_model->SetTitle("Total fit model;HCAL_X(exp)-HCAL_X(act) (m);Counts");
                h_model->Draw("HIST");

                cfit->SaveAs(Form("plots/FOM_dx_test_%s_%f.jpg",kin,FOM));
            }

        // no deletion of hist/fit here so that "best" ones stay around if needed
        // (ROOT will clean up on process exit; you can add explicit delete if desired)
    }

    // ------------------------------------------------------------------------
    // 6. Report best
    // ------------------------------------------------------------------------
    std::cout<<"Best FOM="<<bestFOM<<" with cuts: \n";
    std::cout<<" dy ∈ ["<<bestCuts[0]<<","<<bestCuts[1]<<"]\n";
    std::cout<<" W2 ∈ ["<<bestCuts[2]<<","<<bestCuts[3]<<"]\n";
    std::cout<<" eHCAL_L : ["<<bestCuts[4]<<"]\n";
    std::cout<<" coin_time ∈ ["<<bestCuts[5]<<","<<bestCuts[6]<<"]\n";
}

