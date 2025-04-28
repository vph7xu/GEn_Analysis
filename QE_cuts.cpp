#pragma once

//#include "parse.h"
#include "cuts.h"
#include "plotdxdy.h"
//#include "parse.h"

void plotW2(const char* filename, const char* printfilename, const char *kin, bool flag_eHCAL_cut){
	
    std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt", kin)); //parse the cuts

    cuts cutsobject;
    cutsobject.parsecuts(config);

    //double coin_time_L = 90;//getDoubleValue(config,"coin_time_L");
    //double coin_time_H = 115;//getDoubleValue(config,"coin_time_H");

    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

    double dx = 0.0;
    double dy = 0.0;
    double W2 = 0.0;
    double Q2 = 0.0;
    double coin_time = 0.0;
    double xSH = 0.0;
    double ySH = 0.0;
    double eHCAL = 0.0;

    tree->SetBranchAddress("dx", &dx);
    tree->SetBranchAddress("dy", &dy);
    tree->SetBranchAddress("W2", &W2);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("coin_time", &coin_time);
    tree->SetBranchAddress("xSH", &xSH);
    tree->SetBranchAddress("ySH", &ySH);
    tree->SetBranchAddress("eHCAL", &eHCAL);

    // bins for the xSH
    const double binMin = -1.1;
    const double binMax = 1.1;
    const double binWidth = 0.1;
    const int nBins = static_cast<int>((binMax - binMin) / binWidth) + 1;

    double eHCAL_L = cutsobject.eHCAL_L;
    if (flag_eHCAL_cut == false) eHCAL_L = 0.0;

    std::cout << "eHCAL_L: " << eHCAL_L << std::endl;

    std::vector<TH1D*> W2_histograms;
    for (int i = 0; i < nBins; ++i) {
        double binLowEdge = binMin + i * binWidth;
        double binHighEdge = binLowEdge + binWidth;
        W2_histograms.push_back(new TH1D(Form("hist_bin_%d", i),
                                         Form("W2 for Bin[%f,%f)", binLowEdge, binHighEdge),
                                         1000, -4, 8));
    }

    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        if (coin_time < cutsobject.coin_time_L or coin_time > cutsobject.coin_time_H) {
            int binIndex = static_cast<int>((xSH - binMin) / binWidth);
            if (binIndex >= 0 && binIndex < nBins) {
                W2_histograms[binIndex]->Fill(W2);
            }
        }
        if (i % 1000 == 0)
            std::cout << (i * 100.0 / nentries) << "% \r";
        std::cout.flush();
    }

    TFile outFile(Form("rootfiles/%s_histograms.root", printfilename), "RECREATE");
    for (auto hist : W2_histograms) {
        hist->Write();
    }
    outFile.Close();

    std::vector<double> binCenters, means, meanErrors, sigmas, sigmaErrors;
    const int canvasesNeeded = (nBins + 15) / 16;

    for (int canvasIndex = 0; canvasIndex < canvasesNeeded; ++canvasIndex) {
        TCanvas *canvas = new TCanvas(Form("canvas_%d", canvasIndex), "Histograms", 1200, 1200);
        canvas->Divide(4, 4);

        for (int histIndex = 0; histIndex < 16; ++histIndex) {
            int overallIndex = canvasIndex * 16 + histIndex;
            if (overallIndex < nBins) {
                canvas->cd(histIndex + 1);
                W2_histograms[overallIndex]->Draw();

                TF1 *gausFit = new TF1("gausFit", "gaus", -2, 1.2);
                W2_histograms[overallIndex]->Fit(gausFit, "RQ");

                gausFit->Draw("SAME");
                gPad->Update(); // make sure stats box is created

                double mean = gausFit->GetParameter(1);
                double sigma = gausFit->GetParameter(2);
                double meanError = gausFit->GetParError(1);
                double sigmaError = gausFit->GetParError(2);

                double binCenter = (binMin + overallIndex * binWidth) + binWidth / 2;
                binCenters.push_back(binCenter);
                means.push_back(mean);
                meanErrors.push_back(meanError);
                sigmas.push_back(sigma);
                sigmaErrors.push_back(sigmaError);

                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.05);
                latex.SetTextColor(kRed);
                latex.DrawLatex(0.6, 0.7, Form("Mean = %.2f", mean));
                latex.DrawLatex(0.6, 0.6, Form("Sigma = %.2f", sigma));
            }
        }
        if (canvasIndex == 0)
            canvas->Print(Form("%s_W2_bigbite_acceptance.pdf(", printfilename));
        else if (canvasIndex == (canvasesNeeded - 1))
            canvas->Print(Form("%s_W2_bigbite_acceptance.pdf)", printfilename));
        else
            canvas->Print(Form("%s_W2_bigbite_acceptance.pdf", printfilename));

        TGraphErrors *meanGraph = new TGraphErrors(nBins, &binCenters[0], &means[0], nullptr, &meanErrors[0]);
        TGraphErrors *sigmaGraph = new TGraphErrors(nBins, &binCenters[0], &sigmas[0], nullptr, &sigmaErrors[0]);

        TCanvas *graphCanvas = new TCanvas("graphCanvas", "Fit Parameters", 1200, 600);
        graphCanvas->Divide(2, 1);

        graphCanvas->cd(1);
        meanGraph->SetTitle("Mean of Gaussian Fits");
        meanGraph->GetXaxis()->SetTitle("Bin Center");
        meanGraph->GetYaxis()->SetTitle("Mean");
        meanGraph->GetYaxis()->SetRangeUser(0.8, 1.2);
        meanGraph->SetMarkerStyle(7);
        meanGraph->Draw("AP");

        graphCanvas->cd(2);
        sigmaGraph->SetTitle("Sigma of Gaussian Fits");
        sigmaGraph->GetXaxis()->SetTitle("Bin Center");
        sigmaGraph->GetYaxis()->SetTitle("Sigma");
        sigmaGraph->GetYaxis()->SetRangeUser(0.1, 0.4);
        sigmaGraph->SetMarkerStyle(7);
        sigmaGraph->Draw("AP");

        graphCanvas->Print(Form("plots/%s_W2_QE_peak_eHCAL_cut_%s.pdf", printfilename, std::to_string(flag_eHCAL_cut).c_str()));
    }
	
    for (auto hist : W2_histograms) {
        delete hist;
    }
}

void QE_cuts(const char* filename, const char* printfilename, const char* kin, bool flag_eHCAL_cut){
	
    std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt", kin)); //parse the cuts
    cuts cutsobject;	
    cutsobject.parsecuts(config);

    std::cout << "coin_time : " << cutsobject.coin_time_L << std::endl;
	
    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

    double dx = 0.0; 
    double dy = 0.0; 
    double W2 = 0.0;
    double Q2 = 0.0;
    double vz = 0.0;
    double coin_time = 0.0;
    double eHCAL = 0.0;
    int runnum = 0;
    double ntrack_sbs = 0.0;
    double ntrack = 0.0;
    double ePS = 0.0;

    tree->SetBranchAddress("dx", &dx);
    tree->SetBranchAddress("dy", &dy);
    tree->SetBranchAddress("W2", &W2);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("vz", &vz);
    tree->SetBranchAddress("coin_time", &coin_time);
    tree->SetBranchAddress("eHCAL", &eHCAL);
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("ntrack_sbs", &ntrack_sbs);
    tree->SetBranchAddress("ntrack", &ntrack);
    tree->SetBranchAddress("ePS", &ePS);

    TH1D *h_dx = new TH1D("h_dx", "delta-x distribution;delta-x (m)", 200, -6, 3);
    TH1D *h_dx_gcuts = new TH1D("h_dx_gcuts", "delta-x distribution;delta-x (m)", 200, -6, 3);
    TH1D *h_dx_g_eHCAL_cuts = new TH1D("h_dx_g_eHCAL_cuts", "delta-x distribution;delta-x (m)", 200, -6, 3);
    TH1D *h_dx_g_eHCAL_coin_cuts = new TH1D("h_dx_g_eHCAL_coin_cuts", "delta-x distribution;delta-x (m)", 200, -6, 3);
    TH1D *h_dx_g_eHCAL_coin_W2_cuts = new TH1D("h_dx_g_eHCAL_coin_W2_cuts", "delta-x distribution;delta-x (m)", 200, -6, 3);
    TH1D *h_dx_g_eHCAL_coin_W2_dy_cuts = new TH1D("h_dx_g_eHCAL_coin_W2_dy_cuts", "delta-x distribution;delta-x (m)", 200, -6, 3);
    TH1D *h_dx_g_eHCAL_coin_W2_dy_dx_cuts = new TH1D("h_dx_g_eHCAL_coin_W2_dy_dx_cuts", "delta-x distribution;delta-x (m)", 200, -6, 3);

    TH1D *h_dy = new TH1D("h_dy", "delta-y distribution;delta-y (m)", 200, -4, 4);
    TH1D *h_dy_gcuts = new TH1D("h_dy_gcuts", "delta-y distribution;delta-y (m)", 200, -4, 4);
    TH1D *h_dy_g_eHCAL_cuts = new TH1D("h_dy_g_eHCAL_cuts", "delta-y distribution;delta-y (m)", 200, -4, 4);
    TH1D *h_dy_g_eHCAL_coin_cuts = new TH1D("h_dy_g_eHCAL_coin_cuts", "delta-y distribution;delta-y (m)", 200, -4, 4);
    TH1D *h_dy_g_eHCAL_coin_W2_cuts = new TH1D("h_dy_g_eHCAL_coin_W2_cuts", "delta-y distribution;delta-y (m)", 200, -4, 4);
    TH1D *h_dy_g_eHCAL_coin_W2_dy_cuts = new TH1D("h_dy_g_eHCAL_coin_W2_dy_cuts", "delta-y distribution;delta-y (m)", 200, -4, 4);
    TH1D *h_dy_g_eHCAL_coin_W2_dy_dx_cuts = new TH1D("h_dy_g_eHCAL_coin_W2_dy_dx_cuts", "delta-y distribution;delta-y (m)", 200, -4, 4);

    TH1D *h_W2 = new TH1D("h_W2", "W^{2} distribution; W^{2} (GeV^{2})", 100, 0, 4);
    TH1D *h_W2_gcuts = new TH1D("h_W2_gcuts", "W^{2} distribution; W^{2} (GeV^{2})", 100, 0, 4);
    TH1D *h_W2_g_eHCAL_cuts = new TH1D("h_W2_g_eHCAL_cuts", "W^{2} distribution; W^{2} (GeV^{2})", 100, 0, 4);
    TH1D *h_W2_g_eHCAL_coin_cuts = new TH1D("h_W2_g_eHCAL_coin_cuts", "W^{2} distribution; W^{2} (GeV^{2})", 100, 0, 4);
    TH1D *h_W2_g_eHCAL_coin_W2_cuts = new TH1D("h_W2_g_eHCAL_coin_W2_cuts", "W^{2} distribution; W^{2} (GeV^{2})", 100, 0, 4);
    TH1D *h_W2_g_eHCAL_coin_W2_dy_cuts = new TH1D("h_W2_g_eHCAL_coin_W2_dy_cuts", "W^{2} distribution; W^{2} (GeV^{2})", 100, 0, 4);
    TH1D *h_W2_g_eHCAL_coin_W2_dy_dx_cuts = new TH1D("h_W2_g_eHCAL_coin_W2_dy_dx_cuts", "W^{2} distribution; W^{2} (GeV^{2})", 100, 0, 4);

    TH1D *h_coin_time = new TH1D("h_coin_time", "coincidence time distribution; coincidence time (ns)", 200, 0, 200);
    TH1D *h_coin_time_gcuts = new TH1D("h_coin_time_gcuts", "coincidence time distribution; coincidence time (ns)", 200, 0, 200);
    TH1D *h_coin_time_g_eHCAL_cuts = new TH1D("h_coin_time_g_eHCAL_cuts", "coincidence time distribution; coincidence time (ns)", 200, 0, 200);
    TH1D *h_coin_time_g_eHCAL_coin_cuts = new TH1D("h_coin_time_g_eHCAL_coin_cuts", "coincidence time distribution; coincidence time (ns)", 200, 0, 200);
    TH1D *h_coin_time_g_eHCAL_coin_W2_cuts = new TH1D("h_coin_time_g_eHCAL_coin_W2_cuts", "coincidence time distribution; coincidence time (ns)", 200, 0, 200);
    TH1D *h_coin_time_g_eHCAL_coin_W2_dy_cuts = new TH1D("h_coin_time_g_eHCAL_coin_W2_dy_cuts", "coincidence time distribution; coincidence time (ns)", 200, 0, 200);
    TH1D *h_coin_time_g_eHCAL_coin_W2_dy_dx_cuts = new TH1D("h_coin_time_g_eHCAL_coin_W2_dx_cuts", "coincidence time distribution; coincidence time (ns)", 200, 0, 200);

    double eHCAL_L = cutsobject.eHCAL_L;
    if (flag_eHCAL_cut == false) eHCAL_L = 0.0;

    std::cout << "eHCAL_L: " << eHCAL_L << std::endl;
    std::cout << "run_num_L: " << cutsobject.run_num_L << " run_num_H: " << cutsobject.run_num_H << std::endl;
    std::cout << "coin_time_L: " << cutsobject.coin_time_L << " coin_time_H: " << cutsobject.coin_time_H << std::endl;
    std::cout << "W2_L: " << cutsobject.W2_L << " W2_H: " << cutsobject.W2_H << std::endl;
    std::cout << "dx_L: " << cutsobject.dx_L << " dx_H: " << cutsobject.dx_H << std::endl;
    std::cout << "dy_L: " << cutsobject.dy_L << " dy_H: " << cutsobject.dy_H << std::endl;

    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        bool basic_cut = (cutsobject.run_num_L < runnum and runnum < cutsobject.run_num_H and ntrack > 0); 
        bool gcut = (basic_cut and abs(vz) < 0.27 and ePS > 0.2);
        bool g_eHCAL_cut = (gcut and eHCAL > eHCAL_L);
        bool g_eHCAL_coin_cut = (g_eHCAL_cut and cutsobject.coin_time_L < coin_time and coin_time < cutsobject.coin_time_H);
        bool g_eHCAL_coin_W2_cut = (g_eHCAL_coin_cut and cutsobject.W2_L < W2 and W2 < cutsobject.W2_H);
        bool g_eHCAL_coin_W2_dy_cut = (g_eHCAL_coin_W2_cut and cutsobject.dy_L < dy and dy < cutsobject.dy_H);
        bool g_eHCAL_coin_W2_dy_dx_cut = (g_eHCAL_coin_W2_dy_cut and cutsobject.dx_L < dx and dx < cutsobject.dx_H);

        if (basic_cut) {
            h_dx->Fill(dx);
            h_dy->Fill(dy);
            h_W2->Fill(W2);
            h_coin_time->Fill(coin_time);
        }
        if (gcut) {
            h_dx_gcuts->Fill(dx);
            h_dy_gcuts->Fill(dy);
            h_W2_gcuts->Fill(W2);
            h_coin_time_gcuts->Fill(coin_time);
        }
        if (g_eHCAL_cut) {
            h_dx_g_eHCAL_cuts->Fill(dx);
            h_dy_g_eHCAL_cuts->Fill(dy);
            h_W2_g_eHCAL_cuts->Fill(W2);
            h_coin_time_g_eHCAL_cuts->Fill(coin_time);
        }
        if (g_eHCAL_coin_cut) {
            h_dx_g_eHCAL_coin_cuts->Fill(dx);
            h_dy_g_eHCAL_coin_cuts->Fill(dy);
            h_W2_g_eHCAL_coin_cuts->Fill(W2);
            h_coin_time_g_eHCAL_coin_cuts->Fill(coin_time);
        }
        if (g_eHCAL_coin_W2_cut) {
            h_dx_g_eHCAL_coin_W2_cuts->Fill(dx);
            h_dy_g_eHCAL_coin_W2_cuts->Fill(dy);
            h_W2_g_eHCAL_coin_W2_cuts->Fill(W2);
            h_coin_time_g_eHCAL_coin_W2_cuts->Fill(coin_time);
        }
        if (g_eHCAL_coin_W2_dy_cut) {
            h_dx_g_eHCAL_coin_W2_dy_cuts->Fill(dx);
            h_dy_g_eHCAL_coin_W2_dy_cuts->Fill(dy);
            h_W2_g_eHCAL_coin_W2_dy_cuts->Fill(W2);
            h_coin_time_g_eHCAL_coin_W2_dy_cuts->Fill(coin_time);
        }
        if (g_eHCAL_coin_W2_dy_dx_cut) {
            h_dx_g_eHCAL_coin_W2_dy_dx_cuts->Fill(dx);
            h_dy_g_eHCAL_coin_W2_dy_dx_cuts->Fill(dy);
            h_W2_g_eHCAL_coin_W2_dy_dx_cuts->Fill(W2);
            h_coin_time_g_eHCAL_coin_W2_dy_dx_cuts->Fill(coin_time);
        }
        if (i % 1000 == 0)
            std::cout << (i * 100.0 / nentries) << "% \r";
        std::cout.flush();
    }
	
    // Create canvases for drawing the histograms
    TCanvas* c = new TCanvas("c", "c", 3600, 3000);
    TCanvas* clog = new TCanvas("clog", "clog", 5400, 3000);
    gStyle->SetOptStat(0);

    ////////////////////////////////////////////////////////////////////////////
    h_dx->SetLineWidth(2);
    h_dx_gcuts->SetLineWidth(2);
    h_dx_g_eHCAL_cuts->SetLineWidth(2);
    h_dx_g_eHCAL_coin_cuts->SetLineWidth(2);
    h_dx_g_eHCAL_coin_W2_cuts->SetLineWidth(2);
    h_dx_g_eHCAL_coin_W2_dy_cuts->SetLineWidth(2);
    h_dx_g_eHCAL_coin_W2_dy_dx_cuts->SetLineWidth(2);

    h_dy->SetLineWidth(2);
    h_dy_gcuts->SetLineWidth(2);
    h_dy_g_eHCAL_cuts->SetLineWidth(2);
    h_dy_g_eHCAL_coin_cuts->SetLineWidth(2);
    h_dy_g_eHCAL_coin_W2_cuts->SetLineWidth(2);
    h_dy_g_eHCAL_coin_W2_dy_cuts->SetLineWidth(2);
    h_dy_g_eHCAL_coin_W2_dy_dx_cuts->SetLineWidth(2);

    h_W2->SetLineWidth(2);
    h_W2_gcuts->SetLineWidth(2);
    h_W2_g_eHCAL_cuts->SetLineWidth(2);
    h_W2_g_eHCAL_coin_cuts->SetLineWidth(2);
    h_W2_g_eHCAL_coin_W2_cuts->SetLineWidth(2);
    h_W2_g_eHCAL_coin_W2_dy_cuts->SetLineWidth(2);
    h_W2_g_eHCAL_coin_W2_dy_dx_cuts->SetLineWidth(2);

    h_coin_time->SetLineWidth(2);
    h_coin_time_gcuts->SetLineWidth(2);
    h_coin_time_g_eHCAL_cuts->SetLineWidth(2);
    h_coin_time_g_eHCAL_coin_cuts->SetLineWidth(2);
    h_coin_time_g_eHCAL_coin_W2_cuts->SetLineWidth(2);
    h_coin_time_g_eHCAL_coin_W2_dy_cuts->SetLineWidth(2);
    h_coin_time_g_eHCAL_coin_W2_dy_dx_cuts->SetLineWidth(2);

    ////////////////////////////////////////////////////////////////////////////
    h_dx->SetFillStyle(3008);
    h_dx_gcuts->SetFillStyle(3008);
    h_dx_g_eHCAL_cuts->SetFillStyle(3008);
    h_dx_g_eHCAL_coin_cuts->SetFillStyle(3008);
    h_dx_g_eHCAL_coin_W2_cuts->SetFillStyle(3008);
    h_dx_g_eHCAL_coin_W2_dy_cuts->SetFillStyle(3008);
    h_dx_g_eHCAL_coin_W2_dy_dx_cuts->SetFillStyle(3008);

    h_dy->SetFillStyle(3008);
    h_dy_gcuts->SetFillStyle(3008);
    h_dy_g_eHCAL_cuts->SetFillStyle(3008);
    h_dy_g_eHCAL_coin_cuts->SetFillStyle(3008);
    h_dy_g_eHCAL_coin_W2_cuts->SetFillStyle(3008);
    h_dy_g_eHCAL_coin_W2_dy_cuts->SetFillStyle(3008);
    h_dy_g_eHCAL_coin_W2_dy_dx_cuts->SetFillStyle(3008);

    h_W2->SetFillStyle(3008);
    h_W2_gcuts->SetFillStyle(3008);
    h_W2_g_eHCAL_cuts->SetFillStyle(3008);
    h_W2_g_eHCAL_coin_cuts->SetFillStyle(3008);
    h_W2_g_eHCAL_coin_W2_cuts->SetFillStyle(3008);
    h_W2_g_eHCAL_coin_W2_dy_cuts->SetFillStyle(3008);
    h_W2_g_eHCAL_coin_W2_dy_dx_cuts->SetFillStyle(3008);

    h_coin_time->SetFillStyle(3008);
    h_coin_time_gcuts->SetFillStyle(3008);
    h_coin_time_g_eHCAL_cuts->SetFillStyle(3008);
    h_coin_time_g_eHCAL_coin_cuts->SetFillStyle(3008);
    h_coin_time_g_eHCAL_coin_W2_cuts->SetFillStyle(3008);
    h_coin_time_g_eHCAL_coin_W2_dy_cuts->SetFillStyle(3008);
    h_coin_time_g_eHCAL_coin_W2_dy_dx_cuts->SetFillStyle(3008);

    ////////////////////////////////////////////////////////////////////////////
    h_dx->SetFillColor(kGray);
    h_dx_gcuts->SetFillColor(kBlue);
    h_dx_g_eHCAL_cuts->SetFillColor(kGreen);
    h_dx_g_eHCAL_coin_cuts->SetFillColor(kMagenta);
    h_dx_g_eHCAL_coin_W2_cuts->SetFillColor(kCyan);
    h_dx_g_eHCAL_coin_W2_dy_cuts->SetFillColor(kOrange);
    h_dx_g_eHCAL_coin_W2_dy_dx_cuts->SetFillColor(kRed);

    h_dy->SetFillColor(kGray);
    h_dy_gcuts->SetFillColor(kBlue);
    h_dy_g_eHCAL_cuts->SetFillColor(kGreen);
    h_dy_g_eHCAL_coin_cuts->SetFillColor(kMagenta);
    h_dy_g_eHCAL_coin_W2_cuts->SetFillColor(kCyan);
    h_dy_g_eHCAL_coin_W2_dy_cuts->SetFillColor(kOrange);
    h_dy_g_eHCAL_coin_W2_dy_dx_cuts->SetFillColor(kRed);

    h_W2->SetFillColor(kGray);
    h_W2_gcuts->SetFillColor(kBlue);
    h_W2_g_eHCAL_cuts->SetFillColor(kGreen);
    h_W2_g_eHCAL_coin_cuts->SetFillColor(kMagenta);
    h_W2_g_eHCAL_coin_W2_cuts->SetFillColor(kCyan);
    h_W2_g_eHCAL_coin_W2_dy_cuts->SetFillColor(kOrange);
    h_W2_g_eHCAL_coin_W2_dy_dx_cuts->SetFillColor(kRed);

    h_coin_time->SetFillColor(kGray);
    h_coin_time_gcuts->SetFillColor(kBlue);
    h_coin_time_g_eHCAL_cuts->SetFillColor(kGreen);
    h_coin_time_g_eHCAL_coin_cuts->SetFillColor(kMagenta);
    h_coin_time_g_eHCAL_coin_W2_cuts->SetFillColor(kCyan);
    h_coin_time_g_eHCAL_coin_W2_dy_cuts->SetFillColor(kOrange);
    h_coin_time_g_eHCAL_coin_W2_dy_dx_cuts->SetFillColor(kRed);

    ////////////////////////////////////////////////////////////////////////////
    // Draw on canvas "c" (normal scale)
    c->Divide(2, 2);
    c->cd(1);
    h_dx->Draw();
    h_dx_gcuts->Draw("SAME");
    h_dx_g_eHCAL_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");

    c->cd(2);
    h_dy->Draw();
    h_dy_gcuts->Draw("SAME");
    h_dy_g_eHCAL_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");

    c->cd(3);
    h_W2->Draw();
    h_W2_gcuts->Draw("SAME");
    h_W2_g_eHCAL_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");

    c->cd(4);
    h_coin_time->Draw();
    h_coin_time_gcuts->Draw("SAME");
    h_coin_time_g_eHCAL_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");

    ////////////////////////////////////////////////////////////////////////////
    // Draw on canvas "clog" (log scale) with 3*2 layout and one pad reserved solely for the legend
    clog->Divide(3, 2);
    
    // Pad 1: h_dx group
    clog->cd(1);
    gPad->SetLogy();
    h_dx->Draw();
    h_dx_gcuts->Draw("SAME");
    h_dx_g_eHCAL_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_dx_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");
    
    // Pad 2: h_dy group
    clog->cd(2);
    gPad->SetLogy();
    h_dy->Draw();
    h_dy_gcuts->Draw("SAME");
    h_dy_g_eHCAL_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_dy_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");
    
    // Pad 3: h_W2 group
    clog->cd(4);
    gPad->SetLogy();
    h_W2->Draw();
    h_W2_gcuts->Draw("SAME");
    h_W2_g_eHCAL_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_W2_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");
    
    // Pad 4: h_coin_time group
    clog->cd(5);
    gPad->SetLogy();
    h_coin_time->Draw();
    h_coin_time_gcuts->Draw("SAME");
    h_coin_time_g_eHCAL_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_W2_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_W2_dy_cuts->Draw("SAME");
    h_coin_time_g_eHCAL_coin_W2_dy_dx_cuts->Draw("SAME");
    
    // Pad 5: leave blank
    clog->cd(6);
    gPad->SetLogy();
    gPad->Clear();

    // Pad 6: Only the legend
    clog->cd(3);
    gPad->SetLogy(); // Optional for legend pad
    gPad->Clear();
    TLegend* leg = new TLegend(0.1, 0.1, 0.9, 0.9); // Adjust coordinates as needed
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.04);
    leg->AddEntry(h_dx, "tracking vertex cut", "f");
    leg->AddEntry(h_dx_gcuts, "preshower energy > 0.2GeV cut added", "f");
    leg->AddEntry(h_dx_g_eHCAL_cuts, "HCAL energy > 0.025GeV cut added", "f");
    leg->AddEntry(h_dx_g_eHCAL_coin_cuts, "coincidence time cut added", "f");
    leg->AddEntry(h_dx_g_eHCAL_coin_W2_cuts, "W^{2} < 1.6GeV^{2} cut added", "f");
    leg->AddEntry(h_dx_g_eHCAL_coin_W2_dy_cuts, "abs(dy) < 0.5m cut added", "f");
    leg->AddEntry(h_dx_g_eHCAL_coin_W2_dy_dx_cuts, "abs(dx) < 0.5m cut added", "f");
    leg->Draw();

    c->Print(Form("plots/QE_cuts_%s.png", printfilename));
    clog->Print(Form("plots/QE_cuts_%s_logy.png", printfilename));
}
