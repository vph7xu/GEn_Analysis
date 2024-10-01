#include <iostream>
#include <cmath>
#include "cuts.h"
#include "plotdxdy.h"

struct TrackingData {
    TTree* tree;
    double dx, dy, W2, Q2, coin_time, ptheta, pphi, ptheta_cal, pphi_cal, pN_expect, trP_sbs;
    TH1D *h_delta_ptheta, *h_delta_pphi, *h_delta_pmag, *h_coin_time, *h_W2;
    TH2D *h_ptheta_corr, *h_pphi_corr, *h_pmag_corr, *h_dxdy;

    TrackingData(TTree* t, const char* suffix) : tree(t), dx(0.0), dy(0.0), W2(0.0), Q2(0.0), coin_time(0.0), 
                                                 ptheta(0.0), pphi(0.0), ptheta_cal(0.0), pphi_cal(0.0),
                                                 pN_expect(0.0), trP_sbs(0.0) {
        // Set branches
        tree->SetBranchAddress("dx", &dx);
        tree->SetBranchAddress("dy", &dy);
        tree->SetBranchAddress("W2", &W2);
        tree->SetBranchAddress("Q2", &Q2);
        tree->SetBranchAddress("coin_time", &coin_time);
        tree->SetBranchAddress("ptheta", &ptheta);
        tree->SetBranchAddress("pphi", &pphi);
        tree->SetBranchAddress("ptheta_cal", &ptheta_cal);
        tree->SetBranchAddress("pphi_cal", &pphi_cal);
        tree->SetBranchAddress("pN_expect", &pN_expect);
        tree->SetBranchAddress("trP_sbs", &trP_sbs);

        // Initialize histograms with suffixes to avoid name clashes
        h_delta_ptheta = new TH1D(Form("h_delta_ptheta_%s", suffix), "h_delta_ptheta", 500, -0.1, 0.1);
        h_delta_pphi = new TH1D(Form("h_delta_pphi_%s", suffix), "h_delta_pphi", 500, -0.2, 0.1);
        h_delta_pmag = new TH1D(Form("h_delta_pmag_%s", suffix), "h_delta_pmag", 500, -1, 1);
        h_ptheta_corr = new TH2D(Form("h_ptheta_corr_%s", suffix), "h_ptheta_corr", 200, 0.5, 0.7, 200, 0.5, 0.7);
        h_pphi_corr = new TH2D(Form("h_pphi_corr_%s", suffix), "h_pphi_corr", 200, 2.8, 4, 200, -0.1, 0.4);
        h_pmag_corr = new TH2D(Form("h_pmag_corr_%s", suffix), "h_pmag_corr", 200, 2, 3, 200, 1.5, 3);
        h_coin_time = new TH1D(Form("h_coin_time_%s", suffix), "h_coin_time", 1000, 75, 120);
        h_W2 = new TH1D(Form("h_W2_%s", suffix), "h_W2", 1000, -2, 4);
        h_dxdy = new TH2D(Form("h_dxdy_%s", suffix), "h_dxdy", 200, -2, 2, 200, -6, 2);
    }
};

void processTree(TrackingData &data, const cuts &cutsobject) {
    bool cutcoin_time = false;
    bool cutW2_coin_time = false;
    bool cutQE = false;

    int nentries = data.tree->GetEntries();
    for (int i = 0; i < nentries; i++) {
        data.tree->GetEntry(i);

        // Fill histograms
        data.h_coin_time->Fill(data.coin_time);

        cutcoin_time = abs(data.coin_time - 100) < 10;
        if (cutcoin_time) data.h_W2->Fill(data.W2);

        cutW2_coin_time = abs(data.W2 - 0.88) < 0.3 && abs(data.coin_time - 101) < 11;
        if (cutW2_coin_time) data.h_dxdy->Fill(data.dy, data.dx);

        cutQE = (abs(data.W2 - 0.88) < 0.3 && abs(data.coin_time - 100) < 10 && 
                 (sqrt(pow((data.dx + 2.7) / 0.4, 2) + pow((data.dy - 0.25) / 0.2, 2))) < 1);
        if (cutQE) {
            data.h_delta_ptheta->Fill(data.ptheta - data.ptheta_cal);
            data.h_delta_pphi->Fill(data.pphi - data.pphi_cal + M_PI);
            data.h_delta_pmag->Fill(data.pN_expect - data.trP_sbs);
            data.h_ptheta_corr->Fill(data.ptheta_cal, data.ptheta);
            data.h_pphi_corr->Fill(data.pphi_cal, data.pphi);
            data.h_pmag_corr->Fill(data.pN_expect, data.trP_sbs);
        }
    }
}

void plotComparisons(TrackingData &data1, TrackingData &data2, const char* printfilename) {
    TCanvas* c = new TCanvas("c", "c", 2400, 1800);
    c->Divide(2, 2); // Two columns, four rows for comparison

    // Create a legend
    TLegend *legend = new TLegend(0.1, 0.8, 0.3, 0.9); // x1, y1, x2, y2 in normalized coordinates
    legend->AddEntry(data1.h_delta_ptheta, "global z - 4.670m (old)", "l");
    legend->AddEntry(data2.h_delta_ptheta, "global z - 4.615m (new)", "l");

    // Position and draw the statboxes
    gPad->Update(); // Ensure the statboxes are created
    TPaveStats *st1 = (TPaveStats*)data1.h_coin_time->FindObject("stats");
    TPaveStats *st2 = (TPaveStats*)data2.h_coin_time->FindObject("stats");

    // Plot ptheta correlation comparison
    //c->cd(1);
    //data1.h_ptheta_corr->SetMarkerColor(kRed);
    //data1.h_ptheta_corr->Draw("COLZ");
    //data2.h_ptheta_corr->SetMarkerColor(kBlue);
    //data2.h_ptheta_corr->Draw("COLZ SAME");

    // Plot pphi correlation comparison
    //c->cd(2);
    //data1.h_pphi_corr->SetMarkerColor(kRed);
    //data1.h_pphi_corr->Draw("COLZ");
    //data2.h_pphi_corr->SetMarkerColor(kBlue);
    //data2.h_pphi_corr->Draw("COLZ SAME");

    // Plot pmag correlation comparison
    //c->cd(3);
    //data1.h_pmag_corr->SetMarkerColor(kRed);
    //data1.h_pmag_corr->Draw("COLZ");
    //data2.h_pmag_corr->SetMarkerColor(kBlue);
    //data2.h_pmag_corr->Draw("COLZ SAME");

    // Plot coin_time comparison
    //c->cd(4);
    //data1.h_coin_time->SetLineColor(kRed);
    //data1.h_coin_time->Draw();
    //data2.h_coin_time->SetLineColor(kBlue);
    //data2.h_coin_time->Draw("SAME");

    // Plot W2 comparison
    //c->cd(5);
    //data1.h_W2->SetLineColor(kRed);
    //data1.h_W2->Draw();
    //data2.h_W2->SetLineColor(kBlue);
    //data2.h_W2->Draw("SAME");

    // Plot dx-dy comparison
    //c->cd(6);
    //data1.h_dxdy->SetMarkerColor(kRed);
    //data1.h_dxdy->Draw("COLZ");
    //data2.h_dxdy->SetMarkerColor(kBlue);
    //data2.h_dxdy->Draw("COLZ SAME");

    // Plot delta_ptheta comparison
    c->cd(1);
    data1.h_delta_ptheta->GetXaxis()->SetRangeUser(-1,1);
    data1.h_delta_ptheta->GetYaxis()->SetRangeUser(0,100);
    data1.h_delta_ptheta->SetLineColor(kRed);
    data1.h_delta_ptheta->SetXTitle("ptheta-ptheta_cal (rad)");
    data1.h_delta_ptheta->Draw();
    data2.h_delta_ptheta->SetLineColor(kBlue);
    data2.h_delta_ptheta->Draw("SAME");
    legend->Draw();

    // Plot delta_pphi comparison
    c->cd(2);
    data1.h_delta_pphi->GetXaxis()->SetRangeUser(-1,1);
    data1.h_delta_pphi->GetYaxis()->SetRangeUser(0,100);
    data1.h_delta_pphi->SetLineColor(kRed);
    data1.h_delta_pphi->SetXTitle("pphi-pphi_cal+pi (rad)");
    data1.h_delta_pphi->Draw();
    data2.h_delta_pphi->SetLineColor(kBlue);
    data2.h_delta_pphi->Draw("SAME");
    legend->Draw();

    // Plot delta_pmag comparison
    c->cd(3);
    data2.h_delta_pmag->GetXaxis()->SetRangeUser(-3,1);    
    data2.h_delta_pmag->SetLineColor(kBlue);
    data2.h_delta_pmag->Draw();
    data2.h_delta_pmag->SetXTitle("pN_expect-trP_sbs (GeV)");
    data1.h_delta_pmag->SetLineColor(kRed);
    data1.h_delta_pmag->Draw("SAME");
    legend->Draw();

    // Save the comparison plots
    c->Print(Form("../plots/%s_sbstracking_comparison.pdf", printfilename));
    c->SaveAs(Form("../plots/%s_sbstracking_comparison.png", printfilename));
}

void compare_sbs_tracking_plots(const char* filename1, const char* printfilename, const char *kin1, const char* filename2, const char* kin2) {
    // Parse configuration for both kinematics
    std::map<std::string, std::string> config1 = parseConfig(Form("cuts/cut_%s.txt", kin1)); // parse cuts for first file
    cuts cutsobject1;    
    cutsobject1.parsecuts(config1);

    std::map<std::string, std::string> config2 = parseConfig(Form("cuts/cut_%s.txt", kin2)); // parse cuts for second file
    cuts cutsobject2;    
    cutsobject2.parsecuts(config2);

    // Load both ROOT files
    TFile* file1 = TFile::Open(filename1);
    TTree* tree1 = (TTree*)file1->Get("Tout");

    TFile* file2 = TFile::Open(filename2);
    TTree* tree2 = (TTree*)file2->Get("Tout");

    // Initialize tracking data for both trees
    TrackingData data1(tree1, "1");
    TrackingData data2(tree2, "2");

    // Process both trees
    processTree(data1, cutsobject1);
    processTree(data2, cutsobject2);

    // Plot and compare the results
    plotComparisons(data1, data2, printfilename);
}
