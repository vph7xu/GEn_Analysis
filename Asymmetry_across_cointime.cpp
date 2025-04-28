#include "models.h"
#include "parse.h"
#include "plotdxdy.h"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <TDatime.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TPad.h>
#include <TGraphErrors.h>
#include <TMath.h>

// -------------------------------------------------------------------
// CalculateAsymmetry(...) returns a TGraphErrors* for the cointime bins
// -------------------------------------------------------------------
TGraphErrors* CalculateAsymmetry(std::vector<TH1D*>& Helicity_histograms,
                                 const char* printfilename,
                                 const char* kin,
                                 bool flag_eHCAL_cut)
{
    // Prepare arrays to store the bin centers, asymmetries, and errors
    int nBins = Helicity_histograms.size();
    std::vector<double> cointime_bin_centers(nBins);
    std::vector<double> asymmetries(nBins);
    std::vector<double> errors(nBins);
    
    // Open a file to save the results
    std::ofstream outfile(Form("txt/%s_asymmetry_graph_cointime_binned_eHCAL_cut_%s.txt",
                               kin,
                               std::to_string(flag_eHCAL_cut).c_str()));
    outfile << "cointime_bin_low_edge, Asymmetry, Error\n";

    for (int i = 0; i < nBins; ++i) {
        TH1D* h = Helicity_histograms[i];
        
        // Get the bin content for helicity = +1 and helicity = -1
        int bin_plus  = h->FindBin( 1.0);  // Helicity +1
        int bin_minus = h->FindBin(-1.0);  // Helicity -1
        
        double N_plus  = h->GetBinContent(bin_plus);
        double N_minus = h->GetBinContent(bin_minus);
        double N_total = N_plus + N_minus;

        double asymmetry = 0.0;
        double error     = 0.0;

        if (N_total > 0.0) {
            asymmetry = 100.0 * (N_plus - N_minus) / N_total;
            // Binomial error (approx)
            error = 100.0 * std::sqrt((4.0 * N_plus * N_minus) / 
                                      (N_total * N_total * N_total));
        }

        // For demonstration, define the bin center
        // (You might want to base this on real cointime bin edges.)
        double cointime_bin_low_edge = 0.0;
        double cointime_bin_width    = 10.0;
        double cointime_bin_center   = i * cointime_bin_width 
                                     + cointime_bin_low_edge 
                                     + 0.5 * cointime_bin_width;

        // Store in vectors
        cointime_bin_centers[i] = cointime_bin_center;
        asymmetries[i]          = asymmetry;
        errors[i]               = error;

        // Write to file
        outfile << cointime_bin_low_edge << "," 
                << asymmetry << "," 
                << error << "\n";
    }
    
    outfile.close();

    // Create the TGraphErrors
    TGraphErrors* graph = new TGraphErrors(nBins,
                                           &cointime_bin_centers[0],
                                           &asymmetries[0],
                                           nullptr,  // no x errors
                                           &errors[0]);
    graph->SetTitle("Helicity Asymmetry vs Cointime; Cointime (ns); Asymmetry (%)");
    graph->SetMarkerStyle(21);
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    // Return the graph to be drawn elsewhere
    return graph;
}

// -------------------------------------------------------------------
// Main function: Asymmetry_across_cointime
// -------------------------------------------------------------------
void Asymmetry_across_cointime(const char* filename,
                               const char* printfilename,
                               const char* kin,
                               bool flag_eHCAL_cut)
{
    // -- Parse your config, read DB, etc. --
    std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin));
    std::map<int, int> HelicityCheck  = readCSVToMap("DB/Helicity_quality.csv");
    std::map<int, int> MollerQuality  = readCSVToMap("DB/Moller_quality.csv");

    double coin_time_L   = getDoubleValue(config,"coin_time_L");
    double coin_time_H   = getDoubleValue(config,"coin_time_H");
    double coin_time_ac_L= getDoubleValue(config,"coin_time_ac_L");
    double coin_time_ac_H= getDoubleValue(config,"coin_time_ac_H");
    double W2_L          = getDoubleValue(config,"W2_L");
    double W2_H          = getDoubleValue(config,"W2_H");
    double dy_L          = getDoubleValue(config,"dy_L");
    double dy_H          = getDoubleValue(config,"dy_H");
    double dx_L          = getDoubleValue(config,"dx_L");
    double dx_H          = getDoubleValue(config,"dx_H");
    double eHCAL_L       = getDoubleValue(config,"eHCAL_L");
    double dy_ac_L       = getDoubleValue(config,"dy_ac_L");
    double dy_ac_H       = getDoubleValue(config,"dy_ac_H");
    double run_num_L     = getDoubleValue(config,"run_num_L");
    double run_num_H     = getDoubleValue(config,"run_num_H");
    double IHWP_flip     = getDoubleValue(config,"IHWP_flip");

    double coin_time_offset_L = coin_time_L + 35.0;
    double coin_time_offset_H = coin_time_H + 35.0;

    // Possibly override eHCAL_L if no cut
    if (!flag_eHCAL_cut) {
        eHCAL_L = 0.0;
    }

    // -- Open ROOT file and set up TTree branches --
    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

    int    runnum    = 0;
    int    helicity  = 0;
    int    IHWP      = 0;
    double dx        = 0.0;
    double dy        = 0.0;
    double W2        = 0.0;
    double Q2        = 0.0;
    double coin_time = 0.0;
    double ntrack    = 0.0;
    double vz        = 0.0;
    double eHCAL     = 0.0;
    double ePS       = 0.0;

    tree->SetBranchAddress("runnum",    &runnum);
    tree->SetBranchAddress("helicity",  &helicity);
    tree->SetBranchAddress("IHWP",      &IHWP);
    tree->SetBranchAddress("dx",        &dx);
    tree->SetBranchAddress("dy",        &dy);
    tree->SetBranchAddress("W2",        &W2);
    tree->SetBranchAddress("Q2",        &Q2);
    tree->SetBranchAddress("coin_time", &coin_time);
    tree->SetBranchAddress("ntrack",    &ntrack);
    tree->SetBranchAddress("vz",        &vz);
    tree->SetBranchAddress("eHCAL",     &eHCAL);
    tree->SetBranchAddress("ePS",       &ePS);

    // Histograms
    TH1D* h_coin_time = new TH1D("h_coin_time","Coincidence Time",125,0,250);
    TH2D* h_dxdy_cut_W2_coin = new TH2D("h_dxdy_cut_W2_coin",
                                        "dx-dy with W2+coin cuts",
                                        200,-4,4,
                                        200,-4,4);

    double N_plus_all      = 0.0;
    double N_minus_all     = 0.0;
    double accidental_events = 0.0;
    double QE_events         = 0.0;

    // Setup bins for cointime
    const double binMin   = 0.0;
    const double binMax   = 250.0;
    const double binWidth = 20.0;
    const int nBins = static_cast<int>((binMax - binMin) / binWidth) + 1;

    // Prepare histograms for helicity in each cointime bin
    std::vector<TH1D*> Helicity_histograms;
    Helicity_histograms.reserve(nBins);
    for (int i = 0; i < nBins; ++i) {
        double binLowEdge  = binMin + i*binWidth;
        double binHighEdge = binLowEdge + binWidth;
        TH1D* hist = new TH1D(Form("hist_bin_%d", i),
                              Form("Helicity for bin [%.1f,%.1f)", 
                                   binLowEdge, binHighEdge),
                              5, -2.5, 2.5);
        Helicity_histograms.push_back(hist);
    }


    std::cout<<"eHCAL_L: "<< eHCAL_L <<endl;
    std::cout<<"run_num_L: "<< run_num_L <<" run_num_H: "<< run_num_H <<endl;
    std::cout<<"coin_time_L: "<< coin_time_L <<" coin_time_H: "<< coin_time_H <<endl;
    std::cout<<"W2_L: "<< W2_L <<" W2_H: "<< W2_H <<endl;
    std::cout<<"dx_L: "<< dx_L<<" dx_H: "<< dx_H<<endl;
    std::cout<<"dy_L: "<< dy_L<<" dy_H: "<< dy_H<<endl;

    // -- Loop over entries, fill histos --
    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; i++){
        tree->GetEntry(i);

        // some conditions from your DB checks:
        bool goodHelicity = (lookupValue(HelicityCheck, runnum) == 1);
        bool goodMoller   = (lookupValue(MollerQuality, runnum) == 1);
        bool goodVz       = (std::abs(vz) < 0.27);
        bool goodPS       = (ePS > 0.2);
        bool goodRunRange = (run_num_L < runnum && runnum < run_num_H);
        bool goodEHCAL    = (eHCAL > 0.025); 
        bool validHel     = (helicity == -1 || helicity == 1);

        if (goodHelicity && goodMoller && goodVz && goodPS && validHel && goodRunRange && goodEHCAL)
        {
            if ( (W2_L < W2 && W2 < W2_H) &&
                 (dx_L < dx && dx < dx_H) &&
                 (dy_L < dy && dy < dy_H) )
            {
                // Fill cointime histogram
                h_coin_time->Fill(coin_time);

                // Count accidentals region
                if ( ((coin_time_ac_L < coin_time) && (coin_time < coin_time_L - 10)) ||
                     ((coin_time_H + 10 < coin_time) && (coin_time < coin_time_ac_H)) )
                {
                    // incorporate sign from IHWP etc.
                    int sign = -1 * IHWP * IHWP_flip * helicity;
                    if (sign == 1) {
                        N_plus_all += 1.0;
                    } else if (sign == -1) {
                        N_minus_all += 1.0;
                    }
                }

                // eHCAL cut if user specified
                if (eHCAL > eHCAL_L) {
                    // offset region
                    if (coin_time_offset_L < coin_time && coin_time < coin_time_offset_H){
                        accidental_events += 1.0;
                    }
                    // QE region
                    if (coin_time_L < coin_time && coin_time < coin_time_H){
                        QE_events += 1.0;
                    }
                }

                // Fill the helicity histogram for cointime bin
                int binIndex = static_cast<int>((coin_time - binMin) / binWidth);
                if (binIndex >= 0 && binIndex < nBins) {
                    int sign = -1 * IHWP * IHWP_flip * helicity;
                    Helicity_histograms[binIndex]->Fill(sign);
                }
            }

            // Fill dx-dy histogram with coin_time cut
            if ( (eHCAL > eHCAL_L) &&
                 (W2_L < W2 && W2 < W2_H) &&
                 (coin_time_L < coin_time && coin_time < coin_time_H) )
            {
                h_dxdy_cut_W2_coin->Fill(dy, dx);
            }
        }

        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries) << "% processed...\r";
            std::cout.flush();
        }
    }

    // -- final accidental fraction & asymmetry --
    double f_acc = accidental_events / QE_events;
    double A_acc = (N_plus_all - N_minus_all) / (N_plus_all + N_minus_all);

    double err_f_acc = std::sqrt(accidental_events) / QE_events; // Poisson approx
    double err_A_acc = std::sqrt( (4.0 * N_plus_all * N_minus_all)
                                  / std::pow(N_plus_all + N_minus_all,3) );

    // Print summary to terminal
    std::cout << "\n";
    std::cout << "N_plus_all       = " << N_plus_all       << "\n";
    std::cout << "N_minus_all      = " << N_minus_all      << "\n";
    std::cout << "A_acc            = " << A_acc            << " +- " << err_A_acc   << "\n";
    std::cout << "accidental_events= " << accidental_events << "\n";
    std::cout << "QE_events        = " << QE_events        << "\n";
    std::cout << "f_acc            = " << f_acc            << " +- " << err_f_acc   << "\n";

    // Write accidental results to file
    {
        std::ofstream outfile(Form("txt/%s_accidental_asymmetry_eHCAL_cut_%s.txt",
                                   kin,
                                   std::to_string(flag_eHCAL_cut).c_str()));
        outfile << "N_plus_all      = " << N_plus_all       << "\n";
        outfile << "N_minus_all     = " << N_minus_all      << "\n";
        outfile << "A_acc           = " << A_acc            << "\n";
        outfile << "err_A_acc       = " << err_A_acc        << "\n";
        outfile << "accidental_events = " << accidental_events << "\n";
        outfile << "QE_events       = " << QE_events        << "\n";
        outfile << "f_acc           = " << f_acc            << "\n";
        outfile << "err_f_acc       = " << err_f_acc        << "\n";
    }

    // ---------------------------------------------------------------------
    // Now create your main TCanvas with 2×2 = 4 main pads
    // ---------------------------------------------------------------------
    TCanvas* ccoin = new TCanvas("ccoin","ccoin",3600,3000);
    ccoin->Divide(2,2);

    // Some boxes for the cointime histogram
    double y1 = -10.0;
    double y2 = h_coin_time->GetMaximum();
    TBox* box_offset = new TBox(coin_time_offset_L, y1, coin_time_offset_H, y2);
    box_offset->SetFillColorAlpha(6, 0.3);
    box_offset->SetLineColor(6);

    TBox* box_anti1 = new TBox(coin_time_ac_L, y1, coin_time_L-10, y2);
    box_anti1->SetFillColorAlpha(kRed, 0.3);
    box_anti1->SetLineColor(kRed);

    TBox* box_anti2 = new TBox(coin_time_H+10, y1, coin_time_ac_H, y2);
    box_anti2->SetFillColorAlpha(kRed, 0.3);
    box_anti2->SetLineColor(kRed);

    // Box for dx-dy
    TBox* box_dxdy = new TBox(dy_L, dx_L, dy_H, dx_H);
    box_dxdy->SetFillStyle(0);
    box_dxdy->SetLineColor(kRed);

    // ------------------------------------------------------------
    // pad (1): Just cointime + offset box
    // ------------------------------------------------------------
    ccoin->cd(1);
    h_coin_time->SetXTitle("cointime (ns)");
    h_coin_time->Draw("hist");
    box_offset->Draw("SAME");

    // ------------------------------------------------------------
    // pad (2): we subdivide it into top and bottom sub‐pads
    // ------------------------------------------------------------
    ccoin->cd(2);

    // top sub‐pad ( ~70% of vertical space )
    TPad* p2_top = new TPad("p2_top", "p2_top", 0.0, 0.3, 1.0, 1.0);
    p2_top->SetBottomMargin(0.02); // small bottom margin
    p2_top->Draw();
    p2_top->cd();

    // Draw h_coin_time + accidental boxes in top sub‐pad
    // Optionally hide the x label here, since we'll put it on bottom
    h_coin_time->GetXaxis()->SetTitle("");
    h_coin_time->Draw("hist");
    box_anti1->Draw("SAME");
    box_anti2->Draw("SAME");
    box_offset->Draw("SAME");

    p2_top->Update();

    // --------------------------------------------------
    // bottom sub‐pad ( ~30% ) for TGraphErrors
    // --------------------------------------------------
    ccoin->cd(2); // go back to main pad #2
    TPad* p2_bottom = new TPad("p2_bottom","p2_bottom", 0.0, 0.0, 1.0, 0.3);
    p2_bottom->SetTopMargin(0.02);
    p2_bottom->SetBottomMargin(0.25); // enough space for axis labels
    p2_bottom->Draw();
    p2_bottom->cd();

    // Build the TGraphErrors from Helicity_histograms
    TGraphErrors* gAsym = CalculateAsymmetry(Helicity_histograms,
                                             printfilename,
                                             kin,
                                             flag_eHCAL_cut);

    // Draw the TGraph in the bottom pad
    gAsym->Draw("AP");
    gAsym->SetMarkerStyle(kFullCircle);
    gAsym->SetMarkerColor(kBlue);
    gAsym->SetLineColor(kBlue);

    // Force x range to match the histogram range (0 → 250)
    gAsym->GetXaxis()->SetRangeUser(0, 250);
    gAsym->SetMinimum(-10);
    gAsym->SetMaximum(10);

    // Now set the bottom pad’s X axis label
    gAsym->GetXaxis()->SetTitle("cointime (ns)");

    p2_bottom->Update();

    // ------------------------------------------------------------
    // pad (3): dx‐dy 2D hist
    // ------------------------------------------------------------
    ccoin->cd(3);
    h_dxdy_cut_W2_coin->SetYTitle("dx (m)");
    h_dxdy_cut_W2_coin->SetXTitle("dy (m)");
    h_dxdy_cut_W2_coin->Draw("COLZ");
    box_dxdy->Draw("SAME");

    // pad (4): could use for something else, or leave blank
    // ccoin->cd(4);
    // ... 

    // ------------------------------------------------------------
    // Finally, save the entire canvas
    // ------------------------------------------------------------
    ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%s.pdf",
                       kin, std::to_string(flag_eHCAL_cut).c_str()));
    ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%s.png",
                       kin, std::to_string(flag_eHCAL_cut).c_str()));
    ccoin->SaveAs(Form("plots/%s_cointime_plots_for_accidentals_eHCAL_cut_%s.jpg",
                       kin, std::to_string(flag_eHCAL_cut).c_str()));

    // Done!
}
