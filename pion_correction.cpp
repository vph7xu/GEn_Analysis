#include "models.h"
#include "parse.h"

// A helper function to perform the fit on a given data histogram (h_data) 
// using the pion (h_pion) and QE (h_QE) simulation histograms.
//
// - Returns a *combined fit* histogram (shape-only by default).
// - The fitted weights (pion_weight, qe_weight) are set by reference.
TH1D* performFit(TH1D* h_data, TH1D* h_pion, TH1D* h_QE,
                 double &pion_weight, double &qe_weight)
{
    // 1) Clone each histogram so we can re-scale them for shape fitting
    //    without affecting the originals.
    TH1D* h_data_norm = (TH1D*)h_data->Clone("h_data_norm");
    TH1D* h_pion_norm = (TH1D*)h_pion->Clone("h_pion_norm");
    TH1D* h_QE_norm   = (TH1D*)h_QE->Clone("h_QE_norm");

    // 2) Normalize each histogram to unity (so we are doing a shape-only fit).
    h_data_norm->Scale(1.0 / h_data_norm->Integral());
    h_pion_norm->Scale(1.0 / h_pion_norm->Integral());
    h_QE_norm->Scale(1.0 / h_QE_norm->Integral());

    // 3) Define a TF1 that returns pion_weight * (pion shape) + qe_weight * (QE shape).
    TF1* fitFunction = new TF1("fitFunction",
        [h_pion_norm, h_QE_norm](double *x, double *par) {
            double pion_w   = par[0];
            double qe_w     = par[1];
            int bin         = h_pion_norm->FindBin(x[0]);
            double pion_val = h_pion_norm->GetBinContent(bin);
            double qe_val   = h_QE_norm->GetBinContent(bin);
            return pion_w * pion_val + qe_w * qe_val;
        },
        0.01, 5, 2);

    fitFunction->SetParNames("Pion Weight", "QE Weight");
    fitFunction->SetParameters(0.5, 0.5); // initial guesses
    fitFunction->SetParLimits(0, 0, 1);
    fitFunction->SetParLimits(1, 0, 1);

    // 4) Fit the normalized data histogram using the shape model.
    h_data_norm->Fit(fitFunction, "RQ"); // R=Range-limit, Q=Quiet

    // 5) Retrieve the fitted weights.
    pion_weight = fitFunction->GetParameter(0);
    qe_weight   = fitFunction->GetParameter(1);

    // 6) Construct a shape-only combined fit histogram:
    //    combined_fit(bin) = pion_weight*pion_norm(bin) + qe_weight*QE_norm(bin).
    TH1D* h_combined_fit = (TH1D*)h_pion_norm->Clone("h_combined_fit");
    h_combined_fit->Reset();
    for (int bin = 1; bin <= h_combined_fit->GetNbinsX(); ++bin) {
        double val = pion_weight * h_pion_norm->GetBinContent(bin)
                   + qe_weight   * h_QE_norm->GetBinContent(bin);
        h_combined_fit->SetBinContent(bin, val);
    }

    // 7) At this point, h_combined_fit is also normalized to area = 1
    //    (just like h_data_norm). We return it for you to scale as needed.
    return h_combined_fit;
}


void pion_correction(const char *pion_filename, const char *QE_filename, 
                     const char *data_filename, const char *kin, bool cutQE)
{
    // ------------------------------------------------------
    // 1. READ CONFIGURATIONS, FILES, AND SET UP BRANCHES
    // ------------------------------------------------------
    std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt", kin));
    std::map<int, int> HelicityCheck = readCSVToMap("DB/Helicity_quality.csv");
    std::map<int, int> MollerQuality = readCSVToMap("DB/Moller_quality.csv");

    TFile* pion_file = TFile::Open(pion_filename);
    TTree* pion_tree = (TTree*)pion_file->Get("Tout");

    TFile* QE_file = TFile::Open(QE_filename);
    TTree* QE_tree = (TTree*)QE_file->Get("Tout");

    TFile* data_file = TFile::Open(data_filename);
    TTree* data_tree = (TTree*)data_file->Get("Tout");

    // Get cut values.
    double eHCAL_L     = getDoubleValue(config,"eHCAL_L");
    double coin_time_L = getDoubleValue(config,"coin_time_L");
    double coin_time_H = getDoubleValue(config,"coin_time_H");
    double W2_L        = getDoubleValue(config,"W2_L");
    double W2_H        = getDoubleValue(config,"W2_H");
    double dx_L        = getDoubleValue(config,"dx_L");
    double dx_H        = getDoubleValue(config,"dx_H");
    double dy_L        = getDoubleValue(config,"dy_L");
    double dy_H        = getDoubleValue(config,"dy_H");
    double IHWP_flip   = getDoubleValue(config,"IHWP_flip"); 

    // Simulation tree variables.
    double PSe_pion=0.0, W2_pion=0.0, eHCAL_pion=0.0, SHe_pion=0.0, trP_pion=0.0;
    double PSe_QE  =0.0, W2_QE  =0.0, eHCAL_QE  =0.0, dx_QE =0.0, dy_QE =0.0, SHe_QE=0.0, trP_QE=0.0;

    // Data tree variables.
    int runnum_data=0, IHWP_data=0, helicity_data=0;
    double dx_data=0.0, dy_data=0.0, ntrack_data=0.0, trP_data=0.0;
    double PSe_data=0.0, W2_data=0.0, eHCAL_data=0.0, SHe_data=0.0;
    double coin_time_data=0.0, vz_data=0.0;

    // Set branch addresses for simulation.
    pion_tree->SetBranchAddress("ePS",    &PSe_pion);
    pion_tree->SetBranchAddress("W2",     &W2_pion);
    pion_tree->SetBranchAddress("eHCAL",  &eHCAL_pion);
    pion_tree->SetBranchAddress("eSH", &SHe_pion);
    pion_tree->SetBranchAddress("trP", &trP_pion);
    QE_tree->SetBranchAddress("ePS",      &PSe_QE);
    QE_tree->SetBranchAddress("W2",       &W2_QE);
    QE_tree->SetBranchAddress("eHCAL",    &eHCAL_QE);
    QE_tree->SetBranchAddress("dx",       &dx_QE);
    QE_tree->SetBranchAddress("dy",    &dy_QE);
    QE_tree->SetBranchAddress("eSH",    &SHe_QE);
    QE_tree->SetBranchAddress("trP", &trP_QE);

    // Set branch addresses for data.
    data_tree->SetBranchAddress("ePS",       &PSe_data);
    data_tree->SetBranchAddress("W2",        &W2_data);
    data_tree->SetBranchAddress("eHCAL",     &eHCAL_data);
    data_tree->SetBranchAddress("helicity",  &helicity_data);
    data_tree->SetBranchAddress("coin_time", &coin_time_data);
    data_tree->SetBranchAddress("vz",        &vz_data);
    data_tree->SetBranchAddress("IHWP",      &IHWP_data);
    data_tree->SetBranchAddress("dx",        &dx_data);
    data_tree->SetBranchAddress("dy",        &dy_data);
    data_tree->SetBranchAddress("ntrack",    &ntrack_data);
    data_tree->SetBranchAddress("runnum",    &runnum_data);
    data_tree->SetBranchAddress("eSH",       &SHe_data);
    data_tree->SetBranchAddress("trP",       &trP_data);

    // ------------------------------------------------------
    // 2. CREATE HISTOGRAMS & FILL THEM
    // ------------------------------------------------------
    TH1D *h_PSe_pion = new TH1D("h_PSe_pion",  "Preshower cluster E (pion sim)",  200, 0.01, 2);
    TH1D *h_PSe_QE   = new TH1D("h_PSe_QE",    "Preshower Energy (QE sim)",      200, 0.01, 2);
    TH1D *h_PSe_data = new TH1D("h_PSe_data",  "Preshower Energy (All data)",    200, 0.01, 2);
    TH1D *h_PSe_data_pos = new TH1D("h_PSe_data_pos","Preshower Energy (Helicity +1)", 200,0.01,2);
    TH1D *h_PSe_data_neg = new TH1D("h_PSe_data_neg","Preshower Energy (Helicity -1)", 200,0.01,2);

    // Fill the pion-simulation histogram.
    int nentries_pion = pion_tree->GetEntries();
    for (int i = 0; i < nentries_pion; i++) {

        pion_tree->GetEntry(i);
        double eoverp_pion = (PSe_pion+SHe_pion)/trP_pion;

        bool pass_W2_pion = (W2_L < W2_pion) && (W2_pion < W2_H);
        bool pass_eoverp_pion = 1; /*abs(eoverp_pion-1)<0.2;*/
        if (pass_W2_pion && pass_eoverp_pion) {
            h_PSe_pion->Fill(PSe_pion);
        }
        std::cout << (i * 100.0 / nentries_pion) << "% of pion simulation processed\r";
        std::cout.flush();
    }

    // Fill the QE-simulation histogram.
    int nentries_QE = QE_tree->GetEntries();
    for (int i = 0; i < nentries_QE; i++) {
        QE_tree->GetEntry(i);
        double eoverp_QE = (PSe_QE+SHe_QE)/trP_QE;

        bool pass_W2_QE   = (W2_L < W2_QE) && (W2_QE < W2_H);
        bool pass_eHCAL_QE = (eHCAL_L < eHCAL_QE);
        bool pass_dx_QE = (dx_L<dx_QE) && (dx_QE<dx_H);
        bool pass_dy_QE = (dy_L<dy_QE) && (dy_QE<dy_H);
        bool pass_eoverp_QE = 1; /*abs(eoverp_QE-1)<0.2;*/

        if (cutQE == 0) { 
            //pass_W2_QE = 1;
            pass_eHCAL_QE = 1;
            pass_dx_QE = 1;
            pass_dy_QE = 1;
            pass_eoverp_QE = 1;
        }

        if (pass_W2_QE && pass_eHCAL_QE && pass_dx_QE && pass_dy_QE && pass_eoverp_QE) {
            h_PSe_QE->Fill(PSe_QE);
        }
        std::cout << (i * 100.0 / nentries_QE) << "% of QE simulation processed\r";
        std::cout.flush();
    }

    // Fill the data histograms (all, pos helicity, neg helicity).
    int nentries_data = data_tree->GetEntries();
    for (int i = 0; i < nentries_data; i++) {
        data_tree->GetEntry(i);
        
        double eoverp_data = (PSe_data+SHe_data)/trP_data;

        // Check quality DB:
        if (lookupValue(HelicityCheck,runnum_data)!=1) continue;
        if (lookupValue(MollerQuality, runnum_data)!=1) continue;

        // Check that helicity is ±1
        if (helicity_data != -1 && helicity_data != 1) continue;

        // Check physics cuts
        bool pass_W2_data   = (W2_L < W2_data) && (W2_data < W2_H);
        bool pass_eHCAL_data = (eHCAL_L < eHCAL_data);
        bool pass_coin_data = (coin_time_L < coin_time_data) && (coin_time_data < coin_time_H);
        bool pass_dx_data = (dx_L<dx_data) && (dx_data<dx_H);
        bool pass_dy_data = (dy_L<dy_data) && (dy_data<dy_H);
        bool pass_eoverp_data = 1; /*abs(eoverp_data-1)<0.2;*/

        if (cutQE == 0) { 
            //pass_W2_data = 1;
            pass_eHCAL_data = 1;
            //pass_coin_data = 1;
            pass_dx_data = 1;
            pass_dy_data = 1;
            pass_eoverp_data = 1;
        }

        if (! (pass_W2_data && pass_eHCAL_data && pass_coin_data && pass_dx_data && pass_dy_data && pass_eoverp_data)) continue;

        // Correct helicity for IHWP
        helicity_data = -1 * IHWP_data * IHWP_flip * helicity_data;
        // Fill histograms
        h_PSe_data->Fill(PSe_data);
        if (helicity_data == 1) {
            h_PSe_data_pos->Fill(PSe_data);
        } else {
            h_PSe_data_neg->Fill(PSe_data);
        }
        std::cout << (i * 100.0 / nentries_data) << "% of data processed\r";
        std::cout.flush();
    }

    // ------------------------------------------------------
    // 3. PERFORM FITS USING THE HELPER FUNCTION
    // ------------------------------------------------------
    double pion_weight_all=0, qe_weight_all=0;
    TH1D* h_combined_all = performFit(h_PSe_data, h_PSe_pion, h_PSe_QE,
                                      pion_weight_all, qe_weight_all);

    double pion_weight_pos=0, qe_weight_pos=0;
    TH1D* h_combined_pos = performFit(h_PSe_data_pos, h_PSe_pion, h_PSe_QE,
                                      pion_weight_pos, qe_weight_pos);

    double pion_weight_neg=0, qe_weight_neg=0;
    TH1D* h_combined_neg = performFit(h_PSe_data_neg, h_PSe_pion, h_PSe_QE,
                                      pion_weight_neg, qe_weight_neg);

    // ------------------------------------------------------
    // 4. SCALE HISTOGRAMS TO MATCH DATA AND PLOT
    // ------------------------------------------------------
    // Below is one approach: 
    //  - Keep the *raw data* in true counts,
    //  - Scale the combined fit and the individual pion/QE shape to match data’s total area
    //    based on the fitted fractions.

    // For "all data" 
    double dataIntegral_all = h_PSe_data->Integral();

    // Create scaled copies of the MC shape so we can superimpose them on data:
    TH1D* h_pion_scaled_all = (TH1D*)h_PSe_pion->Clone("h_pion_scaled_all");
    TH1D* h_QE_scaled_all   = (TH1D*)h_PSe_QE->Clone("h_QE_scaled_all");
    TH1D* h_combined_scaled_all = (TH1D*)h_combined_all->Clone("h_combined_scaled_all");

    // Each one is area=1 originally (from performFit), so multiply by dataIntegral_all:
    h_pion_scaled_all     ->Scale(pion_weight_all * dataIntegral_all/h_pion_scaled_all->Integral());
    h_QE_scaled_all       ->Scale(qe_weight_all   * dataIntegral_all/h_QE_scaled_all->Integral());
    h_combined_scaled_all ->Scale(dataIntegral_all); // The combined sum is 1, times total data.

    // Now do the same for positive helicity 
    double dataIntegral_pos = h_PSe_data_pos->Integral();
    TH1D* h_pion_scaled_pos = (TH1D*)h_PSe_pion->Clone("h_pion_scaled_pos");
    TH1D* h_QE_scaled_pos   = (TH1D*)h_PSe_QE->Clone("h_QE_scaled_pos");
    TH1D* h_combined_scaled_pos = (TH1D*)h_combined_pos->Clone("h_combined_scaled_pos");

    h_pion_scaled_pos     ->Scale(pion_weight_pos * dataIntegral_pos/h_pion_scaled_pos->Integral());
    h_QE_scaled_pos       ->Scale(qe_weight_pos   * dataIntegral_pos/h_QE_scaled_pos->Integral());
    h_combined_scaled_pos ->Scale(dataIntegral_pos);

    // And for negative helicity
    double dataIntegral_neg = h_PSe_data_neg->Integral();
    TH1D* h_pion_scaled_neg = (TH1D*)h_PSe_pion->Clone("h_pion_scaled_neg");
    TH1D* h_QE_scaled_neg   = (TH1D*)h_PSe_QE->Clone("h_QE_scaled_neg");
    TH1D* h_combined_scaled_neg = (TH1D*)h_combined_neg->Clone("h_combined_scaled_neg");

    h_pion_scaled_neg     ->Scale(pion_weight_neg * dataIntegral_neg/h_pion_scaled_neg->Integral());
    h_QE_scaled_neg       ->Scale(qe_weight_neg   * dataIntegral_neg/h_QE_scaled_neg->Integral());
    h_combined_scaled_neg ->Scale(dataIntegral_neg);

    // ------------------------------------------------------
    // 5. MAKE A 2x2 CANVAS, DRAW RESULTS
    // ------------------------------------------------------
    TCanvas* c2 = new TCanvas("c2", "2x2 Canvas", 2400, 1800);
    c2->Divide(2, 2);

    // 5.1 Pad1 -> Overall data
    c2->cd(1);
    h_PSe_data->SetLineColor(kGreen);
    h_PSe_data->Draw("HIST"); // data in true counts
    h_pion_scaled_all->SetLineColor(kRed);
    h_pion_scaled_all->Draw("HIST SAME");
    h_combined_scaled_all->SetLineColor(kBlack);
    h_combined_scaled_all->Draw("HIST SAME");
    h_QE_scaled_all->SetLineColor(kBlue);
    h_QE_scaled_all->Draw("HIST SAME");


    TLegend* leg1 = new TLegend(0.7,0.7,0.9,0.9);
    leg1->AddEntry(h_PSe_data,          "Data", "l");
    leg1->AddEntry(h_pion_scaled_all,   "Pion MC", "l");
    leg1->AddEntry(h_QE_scaled_all,     "QE MC", "l");
    leg1->AddEntry(h_combined_scaled_all, "Combined Fit", "l");
    leg1->Draw();

    // 5.2 Pad2 -> Positive helicity
    c2->cd(2);
    h_PSe_data_pos->SetLineColor(kGreen);
    h_PSe_data_pos->Draw("HIST");
    h_pion_scaled_pos->SetLineColor(kRed);
    h_pion_scaled_pos->Draw("HIST SAME");
    h_QE_scaled_pos->SetLineColor(kBlue);
    h_QE_scaled_pos->Draw("HIST SAME");
    h_combined_scaled_pos->SetLineColor(kBlack);
    h_combined_scaled_pos->Draw("HIST SAME");

    TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
    leg2->AddEntry(h_PSe_data_pos,        "Data (+)", "l");
    leg2->AddEntry(h_pion_scaled_pos,     "Pion MC", "l");
    leg2->AddEntry(h_QE_scaled_pos,       "QE MC", "l");
    leg2->AddEntry(h_combined_scaled_pos, "Combined Fit", "l");
    leg2->Draw();

    // 5.3 Pad3 -> Negative helicity
    c2->cd(3);
    h_PSe_data_neg->SetLineColor(kGreen);
    h_PSe_data_neg->Draw("HIST");
    h_pion_scaled_neg->SetLineColor(kRed);
    h_pion_scaled_neg->Draw("HIST SAME");
    h_QE_scaled_neg->SetLineColor(kBlue);
    h_QE_scaled_neg->Draw("HIST SAME");
    h_combined_scaled_neg->SetLineColor(kBlack);
    h_combined_scaled_neg->Draw("HIST SAME");

    TLegend* leg3 = new TLegend(0.7,0.7,0.9,0.9);
    leg3->AddEntry(h_PSe_data_neg,        "Data (-)", "l");
    leg3->AddEntry(h_pion_scaled_neg,     "Pion MC", "l");
    leg3->AddEntry(h_QE_scaled_neg,       "QE MC", "l");
    leg3->AddEntry(h_combined_scaled_neg, "Combined Fit", "l");
    leg3->Draw();

    // 5.4 Pad4 -> Print Fit Results
    c2->cd(4);
    TPaveText* pt = new TPaveText(0.05, 0.2, 0.95, 0.8, "NDC");
    pt->SetTextSize(0.04);
    pt->AddText(Form("Overall:   Pion Weight = %.3f, QE Weight = %.3f, Total Pions =%.3f ",
                     pion_weight_all, qe_weight_all, h_pion_scaled_all->Integral()));
    pt->AddText(Form("Positive:  Pion Weight = %.3f, QE Weight = %.3f, Pos Helicity Pions =%.3f ",
                     pion_weight_pos, qe_weight_pos, h_pion_scaled_pos->Integral()));
    pt->AddText(Form("Negative:  Pion Weight = %.3f, QE Weight = %.3f, Neg Helicity Pions =%.3f ",
                     pion_weight_neg, qe_weight_neg, h_pion_scaled_neg->Integral()));
    pt->AddText(Form("Overall : Pions (ePS>0.2) =%.1f, QE (ePS>0.2) =%.1f, Fraction =%.6f"
        ,h_pion_scaled_all->Integral(h_pion_scaled_all->FindBin(0.2),h_pion_scaled_all->FindBin(3.0))
        ,h_PSe_data->Integral(h_PSe_data->FindBin(0.2),h_PSe_data->FindBin(3.0))
        ,h_pion_scaled_all->Integral(h_pion_scaled_all->FindBin(0.2),h_pion_scaled_all->FindBin(3.0))
        /h_PSe_data->Integral(h_PSe_data->FindBin(0.2),h_PSe_data->FindBin(3.0))));
    pt->Draw();

    // Save the canvas as a PNG.
    c2->SaveAs(Form("plots/%s_PS_plots_2x2_%s.png", kin, std::to_string(cutQE).c_str()));
}
