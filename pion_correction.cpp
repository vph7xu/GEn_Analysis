#include "models.h"
#include "parse.h"

void pion_correction(const char *pion_filename, const char *QE_filename, const char *data_filename){
    TFile* pion_file = TFile::Open(pion_filename);
    TTree* pion_tree = (TTree*)pion_file->Get("Tout");

    TFile* QE_file = TFile::Open(QE_filename);
    TTree* QE_tree = (TTree*)QE_file->Get("Tout");

    TFile* data_file = TFile::Open(data_filename);
    TTree* data_tree = (TTree*)data_file->Get("Tout");

    double PSe_pion = 0.0;
    double W2_pion = 0.0; 

    double PSe_QE = 0.0;
    double W2_QE = 0.0;

    double PSe_data = 0.0;
    double W2_data = 0.0;
    double eHCAL_data = 0.0;

    pion_tree->SetBranchAddress("ePS", &PSe_pion);
    pion_tree->SetBranchAddress("W2", &W2_pion);

    QE_tree->SetBranchAddress("ePS", &PSe_QE);
    QE_tree->SetBranchAddress("W2", &W2_QE);

    data_tree->SetBranchAddress("ePS", &PSe_data);
    data_tree->SetBranchAddress("W2", &W2_data);
    data_tree->SetBranchAddress("eHCAL", &eHCAL_data);

    TH1D *h_PSe_pion = new TH1D("h_PSe_pion", "Preshower cluster E from pion simulation", 500, 0.01, 5);

    TH1D *h_PSe_QE = new TH1D("h_PSe_QE","Preshower Energy",500,0.01,5);

    TH1D *h_PSe_data = new TH1D("h_PSe_data","Preshower Energy",500,0.01,5);

    int nentries_pion = pion_tree->GetEntries();

    int nentries_QE = QE_tree->GetEntries();

    int nentries_data = data_tree->GetEntries();

    for (int i = 0; i < nentries_pion; i++) {
        pion_tree->GetEntry(i);
        if(abs(W2_pion)<2){
           h_PSe_pion->Fill(PSe_pion);
        }
        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries_pion) << "% of simulation file processed\r";
            std::cout.flush();
        }
    }

    for (int i = 0; i < nentries_QE; i++) {
        QE_tree->GetEntry(i);
        if(abs(W2_QE)<2){
           h_PSe_QE->Fill(PSe_QE);
        }
        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries_QE) << "% of simulation file processed\r";
            std::cout.flush();
        }
    }

    for (int i = 0; i < nentries_data/5; i++) {//only process one fifth of the data
        data_tree->GetEntry(i);
        if(abs(W2_data)<2 /*and eHCAL_data>0.025*/){
           h_PSe_data->Fill(PSe_data);
        }
        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries_data) << "% of simulation file processed\r";
            std::cout.flush();
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);

    h_PSe_QE->Scale(1/h_PSe_QE->Integral());
    h_PSe_pion->Scale(1/h_PSe_pion->Integral());
    h_PSe_data->Scale(1/h_PSe_data->Integral());

    // Create a function for fitting: weighted sum of pion and QE histograms
    auto fitFunction = new TF1("fitFunction", 
        [&](double *x, double *par) {
            double pion_weight = par[0];
            double qe_weight = par[1];
            int bin = h_PSe_pion->FindBin(x[0]);
            double pion_value = h_PSe_pion->GetBinContent(bin);
            double qe_value = h_PSe_QE->GetBinContent(bin);
            return pion_weight * pion_value + qe_weight * qe_value;
        }, 
        0.01, 5, 2);

    fitFunction->SetParNames("Pion Weight", "QE Weight");
    fitFunction->SetParameters(0.5, 0.5); // Initial guesses for weights
    fitFunction->SetParLimits(0, 0, 1);   // Pion weight should be between 0 and 1
    fitFunction->SetParLimits(1, 0, 1);   // QE weight should be between 0 and 1

    // Fit the combined model to the data histogram
    h_PSe_data->Fit(fitFunction, "RL");

    double pion_weight = fitFunction->GetParameter(0);
    double qe_weight = fitFunction->GetParameter(1);

    // Create a new histogram to represent the combined fit
    TH1D* h_combined_fit = (TH1D*)h_PSe_pion->Clone("h_combined_fit");
    h_combined_fit->Reset();
    h_combined_fit->SetTitle("Combined Fit (Pion + QE)");
    h_combined_fit->SetLineColor(kBlack);

    // Fill the combined fit histogram with the weighted sum of pion and QE histograms
    for (int bin = 1; bin <= h_combined_fit->GetNbinsX(); bin++) {
        double pion_value = h_PSe_pion->GetBinContent(bin);
        double qe_value = h_PSe_QE->GetBinContent(bin);
        double combined_value = pion_weight * pion_value + qe_weight * qe_value;
        h_combined_fit->SetBinContent(bin, combined_value);
    }


    h_PSe_QE->SetLineColor(kBlue);
    h_PSe_pion->SetLineColor(kRed);
    h_PSe_data->SetLineColor(kGreen);
    fitFunction->SetLineColor(kBlack);

    h_PSe_QE->Scale(qe_weight);
    h_PSe_pion->Scale(pion_weight);

    h_PSe_QE->Draw("HIST");
    h_PSe_pion->Draw("HIST SAME");
    h_PSe_data->Draw("HIST SAME");
    h_combined_fit->Draw("HIST SAME");

    h_PSe_QE->SetMaximum(0.1);

    // Print fit results
    std::cout << "Fit Results:\n";
    std::cout << "Pion Weight: " << fitFunction->GetParameter(0) << std::endl;
    std::cout << "QE Weight: " << fitFunction->GetParameter(1) << std::endl;
    
    
    // Set up TFractionFitter with data and templates (pion and QE)
    /*TObjArray *mc = new TObjArray(2); // Number of MC histograms
    mc->Add(h_PSe_pion);
    mc->Add(h_PSe_QE);

    TFractionFitter* fit = new TFractionFitter(h_PSe_data, mc);
    fit->Constrain(0, 0.0, 1.0); // Constrain pion weight between 0 and 1
    fit->Constrain(1, 0.0, 1.0); // Constrain QE weight between 0 and 1

    // Perform the fit
    int status = fit->Fit();

    if (status == 0) { // Check if fit is successful
        double pion_weight, qe_weight;
        fit->GetResult(0, pion_weight, qe_weight);
        fit->GetResult(1, qe_weight, qe_weight);

        std::cout << "Fit Results:" << std::endl;
        std::cout << "Pion Weight: " << pion_weight << std::endl;
        std::cout << "QE Weight: " << qe_weight << std::endl;

        // Draw data and fit results
        TCanvas *c1 = new TCanvas("c1", "Preshower Energy Fit", 800, 600);

        h_PSe_QE->Scale(qe_weight);
        h_PSe_pion->Scale(pion_weight);

        h_PSe_data->SetLineColor(kGreen);
        h_PSe_data->Draw("HIST");

        h_PSe_pion->Scale(pion_weight); // Scale pion histogram by its weight
        h_PSe_QE->Scale(qe_weight); // Scale QE histogram by its weight

        h_PSe_pion->SetLineColor(kRed);
        h_PSe_pion->Draw("HIST SAME");

        h_PSe_QE->SetLineColor(kBlue);
        h_PSe_QE->Draw("HIST SAME");

        TH1D *h_combined_fit = (TH1D*)h_PSe_pion->Clone("h_combined_fit");
        h_combined_fit->Add(h_PSe_QE); // Combined fit
        h_combined_fit->SetLineColor(kBlack);
        h_combined_fit->Draw("HIST SAME");
    } else {
        std::cout << "Fit failed with status: " << status << std::endl;
    }
*/

}
