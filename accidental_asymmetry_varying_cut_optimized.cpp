#include "models.h"
#include "parse.h"
#include "plotdxdy.h"
#include <cmath>


void accidental_asymmetry_varying_cut_optimized(const char* filename, const char* printfilename, const char* kin){

    std::map<std::string, std::string> config1 = parseConfig(Form("cuts/cut_%s.txt", kin));
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

    // Load the ROOT file and TTree
    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

    // Histogram for asymmetry calculation
    TH1D* hNplus = new TH1D("hNplus", "Nplus", 100, coin_time_L, coin_time_H); // for helicity == 1
    TH1D* hNminus = new TH1D("hNminus", "Nminus", 100, coin_time_L, coin_time_H); // for helicity == -1

    // Define the base selection cut based on dx, dy, W2
    TString baseCut = Form("(W2 > %f && W2 < %f) && (dx > %f && dx < %f) && (dy > %f && dy < %f)", 
                            W2_L, W2_H, dx_L, dx_H, dy_L, dy_H);
    
    // Draw with cuts for helicity == 1 (Nplus)
    tree->Draw(Form("coin_time >> hNplus"), baseCut + " && (-1*helicity*IHWP) == 1", "goff"); // only for GEN2
    
    // Draw with cuts for helicity == -1 (Nminus)
    tree->Draw(Form("coin_time >> hNminus"), baseCut + " && (-1*helicity*IHWP) == -1", "goff"); // only for GEN2

    // Create graph for asymmetry plot
    TGraphErrors *gAsymCut = new TGraphErrors();

    // Loop over bins to calculate asymmetry for each coin time bin
    int nbins = hNplus->GetNbinsX();
    for (int bin = 1; bin <= nbins; bin++) {
        double Nplus = hNplus->GetBinContent(bin);
        double Nminus = hNminus->GetBinContent(bin);
        double coin_time_cut = hNplus->GetBinCenter(bin);

        if (bin %2 == 0 ) std::cout << (bin * 100.0/ nbins) << "% \r"; 
        std::cout.flush();

        if (Nplus + Nminus > 0) {
            double Aexp = (Nplus - Nminus) * 100 / (Nplus + Nminus);  // Asymmetry calculation
            double errAexp = 2 * 100 * sqrt((Nplus * Nminus)) / ((Nplus + Nminus) * (Nplus + Nminus));  // Error on asymmetry

            gAsymCut->SetPoint(bin - 1, coin_time_cut, Aexp);
            gAsymCut->SetPointError(bin - 1, 0, errAexp);

            std::cout << "Coin time: " << coin_time_cut << " Asymmetry: " << Aexp << " ± " << errAexp << std::endl;
        }
    }

    // Plot the asymmetry graph
    TCanvas *cCut = new TCanvas("cCut", "Asymmetry vs Coin Time Cut", 1200, 900);
    gAsymCut->SetTitle("Asymmetry vs Coin Time Cut");
    gAsymCut->GetXaxis()->SetTitle("Coin Time (ns)");
    gAsymCut->GetYaxis()->SetTitle("Asymmetry (%)");
    gAsymCut->SetMarkerStyle(7);
    gAsymCut->SetMarkerSize(1.0);
    gAsymCut->SetMarkerColor(kBlue);
    gAsymCut->Draw("AP");

    cCut->SaveAs(Form("plots/%s_varying_cut_optimized.png", printfilename));
    cCut->Print(Form("plots/%s_varying_cut_optimized.pdf", printfilename));

    // Clean up
    delete hNplus;
    delete hNminus;
    delete gAsymCut;
}