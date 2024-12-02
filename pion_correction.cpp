#include "models.h"
#include "parse.h"

void pion_correction(const char *pion_filename, const char *QE_filename, const char *data_filename,const char *kin){
    
    std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin));

    TFile* pion_file = TFile::Open(pion_filename);
    TTree* pion_tree = (TTree*)pion_file->Get("Tout");

    TFile* QE_file = TFile::Open(QE_filename);
    TTree* QE_tree = (TTree*)QE_file->Get("Tout");

    TFile* data_file = TFile::Open(data_filename);
    TTree* data_tree = (TTree*)data_file->Get("Tout");

    //cuts from cuts file
    double eHCAL_L = getDoubleValue(config,"eHCAL_L");
    double coin_time_L = getDoubleValue(config,"coin_time_L");
    double coin_time_H = getDoubleValue(config,"coin_time_H");

    //varibales from root files
    double PSe_pion = 0.0;
    double W2_pion = 0.0;
    double eHCAL_pion = 0.0; 

    double PSe_QE = 0.0;
    double W2_QE = 0.0;
    double eHCAL_QE = 0.0;

    double PSe_data = 0.0;
    double W2_data = 0.0;
    double eHCAL_data = 0.0;
    int helicity_data = 0;
    double coin_time_data = 0.0;
    double vz_data = 0.0;

    pion_tree->SetBranchAddress("ePS", &PSe_pion);
    pion_tree->SetBranchAddress("W2", &W2_pion);
    pion_tree->SetBranchAddress("eHCAL", &eHCAL_pion);

    QE_tree->SetBranchAddress("ePS", &PSe_QE);
    QE_tree->SetBranchAddress("W2", &W2_QE);
    QE_tree->SetBranchAddress("eHCAL", &eHCAL_QE);

    data_tree->SetBranchAddress("ePS", &PSe_data);
    data_tree->SetBranchAddress("W2", &W2_data);
    data_tree->SetBranchAddress("eHCAL", &eHCAL_data);
    data_tree->SetBranchAddress("helicity", &helicity_data);
    data_tree->SetBranchAddress("coin_time", &coin_time_data);
    data_tree->SetBranchAddress("vz", &vz_data);

    TH1D *h_PSe_pion = new TH1D("h_PSe_pion", "Preshower cluster E from pion simulation", 200, 0.01, 3);

    TH1D *h_PSe_QE = new TH1D("h_PSe_QE","Preshower Energy",200,0.01,3);

    TH1D *h_PSe_data = new TH1D("h_PSe_data","Preshower Energy",200,0.01,3);

    int nentries_pion = pion_tree->GetEntries();

    int nentries_QE = QE_tree->GetEntries();

    int nentries_data = data_tree->GetEntries();

    //bool global_cut_pion_sim = 0;
    //bool global_cut_QE_sim = 0;
    //bool global_cut_data = 0;

    for (int i = 0; i < nentries_pion; i++) {
        pion_tree->GetEntry(i);
        
        //if (eHCAL_pion>eHCAL_L){
            if(abs(W2_pion)<2){
                h_PSe_pion->Fill(PSe_pion);
            }
        //}

        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries_pion) << "% of simulation file processed\r";
            std::cout.flush();
        }
    }

    for (int i = 0; i < nentries_QE; i++) {
        QE_tree->GetEntry(i);
        
        if (eHCAL_QE>eHCAL_L){
            if(abs(W2_QE)<2){
                h_PSe_QE->Fill(PSe_QE);
           }
        }
        
        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries_QE) << "% of simulation file processed\r";
            std::cout.flush();
        }
    }

    for (int i = 0; i < nentries_data; i++) {//only process one fifth of the data
        data_tree->GetEntry(i);
        
        if (eHCAL_data>eHCAL_L and abs(vz_data)<0.27/*and (helicity_data == -1 or helicity_data ==1)*/){
            if(abs(W2_data)<2 /*and eHCAL_data>0.025*/){
                h_PSe_data->Fill(PSe_data);
            }
        }

        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries_data) << "% of simulation file processed\r";
            std::cout.flush();
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Canvas", 2400, 1800);

    double data_scale = h_PSe_data->Integral();

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

    h_PSe_QE->Scale(qe_weight*data_scale);
    h_PSe_pion->Scale(pion_weight*data_scale);
    h_PSe_data->Scale(data_scale);
    h_combined_fit->Scale(data_scale);



    // Create the legend
    TLegend* legend = new TLegend(0.7, 0.7, 0.9, 0.9); // Coordinates for the legend box
    legend->AddEntry(h_PSe_pion, "Pion Simulation", "l");
    legend->AddEntry(h_PSe_QE, "QE Simulation", "l");
    legend->AddEntry(h_PSe_data, "Data", "l");
    legend->AddEntry(h_combined_fit, "Combined Fit", "l");
    legend->SetBorderSize(1);
    legend->SetTextSize(0.03);
    legend->SetFillColor(0);

    c1->Divide(1,1);
    c1->cd(1);
    h_PSe_pion->Draw("HIST");
    h_PSe_QE->Draw("HIST SAME");
    h_PSe_data->Draw("HIST SAME");
    h_combined_fit->Draw("HIST SAME");
    legend->Draw();

    h_PSe_QE->SetMaximum();

    // Print fit results
    std::cout << "Fit Results:\n";
    std::cout << "Pion Weight: " << fitFunction->GetParameter(0) << std::endl;
    std::cout << "QE Weight: " << fitFunction->GetParameter(1) << std::endl;
    

    c1->SaveAs(Form("plots/%s_PS_plots.png",kin));

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
