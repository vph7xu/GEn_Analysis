#include "cuts.h"
#include "plotdxdy.h"
#include <cmath>
#include <fstream>

void Hcal_fsam(const char* filename, const char* printfilename, const char* kin){

	gStyle->SetOptFit(1111);

    std::map<std::string, std::string> config = parseConfig(Form("cuts/calib/cut_%s.txt",kin)); //parse the cuts
    cuts cutsobject;
    cutsobject.parsecuts(config);

    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

    double eHCAL = 0.0;
    double xHCAL = 0.0;
    double yHCAL = 0.0;
    double Q2 = 0.0;
    double dx = 0.0;
    double dy = 0.0;
    double W2 = 0.0;
    double coin_time = 0.0;
    double trP_sbs = 0.0;
    double nblk_HCAL = 0.0;
    double ntrack_sbs = 0.0;
    double vz = 0.0;
    double vz_sbs = 0.0;
    double hcal_clus_id[1000];
    double hcal_clus_mem_id[1000];
    double hcal_clus_mem_e[1000];
    int runnum = 0;
    
    tree->SetBranchAddress("eHCAL",&eHCAL);
    tree->SetBranchAddress("xHCAL",&xHCAL);
    tree->SetBranchAddress("yHCAL",&yHCAL);
    tree->SetBranchAddress("Q2",&Q2);
    tree->SetBranchAddress("dx",&dx);
    tree->SetBranchAddress("dy",&dy);
    tree->SetBranchAddress("W2",&W2);
    tree->SetBranchAddress("coin_time",&coin_time);
    tree->SetBranchAddress("trP_sbs",&trP_sbs);
    tree->SetBranchAddress("hcal_clus_id",&hcal_clus_id);
    tree->SetBranchAddress("hcal_clus_mem_id",&hcal_clus_mem_id);
    tree->SetBranchAddress("hcal_clus_mem_e",&hcal_clus_mem_e);
    tree->SetBranchAddress("nblk_HCAL",&nblk_HCAL);
    tree->SetBranchAddress("ntrack_sbs",&ntrack_sbs);
    tree->SetBranchAddress("vz",&vz);
    tree->SetBranchAddress("vz_sbs",&vz_sbs);
    tree->SetBranchAddress("runnum",&runnum);

    double mp = 0.938272;
    double fs = 0;
    double seedsumratio = 0.0;
    double seedsecondratio = 0.0;
    bool elastic_cut = false;
    int nentries = tree->GetEntries();

    TH1D *hfs = new TH1D("hfs","Sampling fraction distribution",100,0,0.5);
    TH1D *her = new TH1D("her","Actual energy dep/expected energy dep",200,0.0,0.4);
    TH1D *hpN = new TH1D("hpN","Outgoing nucleon momentum",1000,-0.1,4);

    TH2D *hfsample = new TH2D("hfsample","Sampling fraction vs xHCAL distribution",100,-2.7,1.2,100,0,0.4);
    TH2D *hfsampley = new TH2D("hfsampley","Sampling fraction vs yHCAL distribution",100,-0.8,0.8,100,0,0.4);
    TH2D *hdxdy = new TH2D("hdxdy","dx vs dy distribution",400,-2,2,800,-4,4);
    TH2D *heratio = new TH2D("heratio", "Actual energy dep/expected energy dep vs xHCAL",200,-3,1.5,200,0.0,0.4);
    TH2D *heratioy = new TH2D("heratioy", "Actual energy dep/expected energy dep vs xHCAL",200,-0.8,0.8,200,0.0,0.4);

    TH2D *hfsvblkid = new TH2D("hfsvblkid","sampling fraction vs blk id",300,0,300,200,0.025,0.3);

    TH1D *hehcal = new TH1D("hehcal","energy measured by Hcal",100,0,2);
    TH1D *hekin = new TH1D("hekin","kinetic energy (tracking)",500,0,10);

    TH1D *hcointime = new TH1D("hcointime", "cointime distribution",1000,40,200);
    TH1D *hW2 = new TH1D("hW2","W2 distribution",1000,-4,4);

    TH1D *hnblk_hcal = new TH1D("hnblk_hcal","number of blocks in the primary cluster", 100, 0, 15);
    TH1D *hseedsumratio = new TH1D("hseedsumratio","seed block energy / primary cluster energy",200,0,1.2);
    TH1D *hseedsecondratio = new TH1D("hseedsecondratio","second block energy / seed block energy",200,0,1.2);

    TH1D *hxhcal = new TH1D("hxhcal","xhcal",200,-3,1.5);
    TH1D *hyhcal = new TH1D("hyhcal","yhcal",200,-0.8,0.8);

    TH2D *hvertexcorr = new TH2D("hvertexcorr", "vertex correlation", 100,-0.5,0.5,100,-0.5,0.5);
    TH1D *hdeltavertex = new TH1D("hdeltavertex","delta vertex", 100, -0.5, 0.5);

    TGraphErrors *graph = new TGraphErrors(hfsample->GetNbinsX());
    
    for (int i = 0; i<nentries; i++){
        tree->GetEntry(i);
        
        if (ntrack_sbs>0){
            
            elastic_cut = (fabs(coin_time - cutsobject.coin_time_mean) < cutsobject.coin_time_width 
                           && fabs(vz - vz_sbs) < 0.1 /*and 0<W2 and W2<5*/);
        
            //elastic_cut = (((pow((dy-cutsobject.dy_C)/cutsobject.dy_R,2)+pow((dx-cutsobject.dx_C)/cutsobject.dx_R,2))<=2.5)
            //&&(abs(coin_time-cutsobject.coin_time_mean)<cutsobject.coin_time_width)&&(abs(W2-cutsobject.W2_mean)<cutsobject.W2_width));//? true:false;

            // Fill a few diagnostic histograms
            hcointime->Fill(coin_time);
            hvertexcorr->Fill(vz,vz_sbs);
            hdeltavertex->Fill(vz - vz_sbs);

            // Check cointime and vertex cuts for W2 fill
            if (fabs(coin_time - cutsobject.coin_time_mean) < cutsobject.coin_time_width 
                && fabs(vz - vz_sbs) < 0.1 /*and 0<W2 and W2<5*/)
            {
                hW2->Fill(W2);
            }

            // dx/dy fill
            if((fabs(coin_time - cutsobject.coin_time_mean) < cutsobject.coin_time_width) 
               && fabs(vz - vz_sbs) < 0.1)
            {
                hdxdy->Fill(dy,dx);
            }

            // Now the main SF fill if it passes "elastic_cut"
            if(elastic_cut){
                
                /*std::cout << "Original std::array: ";
                for(int i = 0; i < nblk_HCAL; ++i){
                    std::cout << hcal_clus_mem_e[i] << " ";
                }
                std::cout << std::endl;
                */

                std::sort(hcal_clus_mem_e, hcal_clus_mem_e+static_cast<int>(std::floor(nblk_HCAL)), std::greater<double>());
                
                /*std::cout << "Sorted std::array (high to low): ";
                for(int i = 0; i < nblk_HCAL; ++i){
                    std::cout << hcal_clus_mem_e[i] << " ";
                }
                std::cout << std::endl;
                */

                seedsumratio = hcal_clus_mem_e[0]/eHCAL;
                hseedsumratio->Fill(seedsumratio);

                seedsecondratio = hcal_clus_mem_e[1]/hcal_clus_mem_e[0];
                hseedsecondratio->Fill(seedsecondratio);

                hxhcal->Fill(xHCAL);
                hyhcal->Fill(yHCAL);

                fs = eHCAL * 2.0 * mp / Q2;  
                hfs->Fill(fs);
                hfsample->Fill(xHCAL, fs);
                hfsampley->Fill(yHCAL, fs);

                double Eexp = sqrt(mp*mp + trP_sbs*trP_sbs) - mp;
                double ratio = 0;
                
                if(Eexp > 0.001 /*&& abs(Eexp-5.25)<=0.75*/ ) {
                    ratio = eHCAL / Eexp;
                
                    heratio->Fill(xHCAL, ratio);
                    heratioy->Fill(yHCAL, ratio);
                    her->Fill(ratio);

                    hehcal->Fill(eHCAL);
                    hekin->Fill(Eexp);
                }

                

                hpN->Fill(trP_sbs);

                hnblk_hcal->Fill(nblk_HCAL);

                // Fill block ID SF
                for (int iblk = 0; iblk < nblk_HCAL; iblk++){
                    hfsvblkid->Fill(hcal_clus_mem_id[iblk], ratio);
                }

                // For completeness, also fill dx/dy
                hdxdy->Fill(dy,dx);
            }
        }
        if (i %1000 == 0) {
            std::cout << (i * 100.0 / nentries) << "% \r";
            std::cout.flush();
        }
    }

    // Create TProfiles
    TProfile *profX   = hfsample->ProfileX();
    TProfile *profXy  = hfsampley->ProfileX();
    TProfile *profX1  = heratio->ProfileX();
    TProfile *profXy1 = heratioy->ProfileX();
    TProfile *profX2  = hfsvblkid->ProfileX();

    // Write sampling fractions per block to file
    std::ofstream outFile(Form("txt/%s_sampling_fractions_each_blk.txt", kin));
    double sampling_fraction_means[300];

    for (int i = 1; i <= profX2->GetNbinsX(); ++i){
        sampling_fraction_means[i-1] = profX2->GetBinContent(i);
        outFile << i-1 << " " << sampling_fraction_means[i-1] << "\n";
    }
    outFile.close();

    // Create elliptical cut
    TCutG *cut_elipse = CreateOvalCut("cut_elipse",
                                      cutsobject.dy_C, 
                                      cutsobject.dx_C,
                                      cutsobject.dy_R*sqrt(2.5),
                                      cutsobject.dx_R*sqrt(2.5),
                                      100); 

    // =========================
    // Create and fill canvases
    // =========================
    TCanvas *c  = new TCanvas("c", "c", 3200, 2400);
    TCanvas *cc = new TCanvas("cc","cc",3200,2400);
    TCanvas *c1 = new TCanvas("c1","c1",3200,2400);
    TCanvas *c2 = new TCanvas("c2","c2",3200,2400);
    TCanvas *cc2= new TCanvas("cc2","cc2",3200,2400);
    TCanvas *c3 = new TCanvas("c3","c3",3200,2400);
    TCanvas *c4 = new TCanvas("c4","c4",3200,2400);

    // ---- Canvas c ----
    c->Divide(1,2);
    c->cd(1);
    hfsample->SetXTitle("xHCAL (m)");
    hfsample->SetYTitle("eHCAL*2*mp/Q^2");
    hfsample->Draw("COLZ");
    profX->SetMarkerStyle(20);
    profX->SetLineColor(kRed);
    profX->SetMarkerColor(kRed);
    // Increase marker size
    profX->SetMarkerSize(1.5);
    profX->Draw("SAME");

    c->cd(2);
    hfsampley->SetXTitle("yHCAL (m)");
    hfsampley->SetYTitle("eHCAL*2*mp/Q^2");
    hfsampley->Draw("COLZ");
    profXy->SetMarkerStyle(20);
    profXy->SetLineColor(kRed);
    profXy->SetMarkerColor(kRed);
    profXy->SetMarkerSize(1.5);
    profXy->Draw("SAME");

    // ---- Canvas cc ----
    cc->Divide(1,1);
    cc->cd(1);
    hfs->SetLineWidth(3);
    hfs->SetXTitle("eHCAL*2*mp/Q^2");
    hfs->Draw();

    // ---- Canvas c1 ----
    c1->Divide(2,2);
    c1->cd(1);
    hcointime->SetXTitle("cointime (ns)");
    hcointime->Draw();
    TLine *line01 = new TLine(cutsobject.coin_time_mean-cutsobject.coin_time_width, 0,
                              cutsobject.coin_time_mean-cutsobject.coin_time_width, 2750);
    TLine *line02 = new TLine(cutsobject.coin_time_mean+cutsobject.coin_time_width, 0,
                              cutsobject.coin_time_mean+cutsobject.coin_time_width, 2750);
    line01->SetLineColor(kRed);
    line02->SetLineColor(kRed);
    line01->SetLineWidth(2);
    line02->SetLineWidth(2);
    line01->Draw("same");
    line02->Draw("same");

    c1->cd(2);
    TLine *line1 = new TLine(cutsobject.W2_mean-cutsobject.W2_width, 0,
                             cutsobject.W2_mean-cutsobject.W2_width, 500);
    TLine *line2 = new TLine(cutsobject.W2_mean+cutsobject.W2_width, 0,
                             cutsobject.W2_mean+cutsobject.W2_width, 500);
    hW2->Draw();
    hW2->SetXTitle("W^2 (GeV^2)");
    line1->SetLineColor(kRed);
    line2->SetLineColor(kRed);
    line1->SetLineWidth(2);
    line2->SetLineWidth(2);
    line1->Draw("same");
    line2->Draw("same");

    c1->cd(3);
    hdxdy->Draw("COLZ");
    hdxdy->SetXTitle("dy (m)");
    hdxdy->SetYTitle("dx (m)");
    cut_elipse->Draw("L");

    c1->cd(4);
    hpN->Draw();

    // ---- Canvas c2 ----
    c2->Divide(1,2);
    c2->cd(1);
    heratio->SetYTitle("eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp)");
    heratio->SetXTitle("xHCAL (m)");
    heratio->Draw("COLZ");
    profX1->SetMarkerStyle(20);
    profX1->SetLineColor(kRed);
    profX1->SetMarkerColor(kRed);
    // Increase marker size
    profX1->SetMarkerSize(1.5);
    profX1->Draw("SAME");

    c2->cd(2);
    heratioy->SetYTitle("eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp)");
    heratioy->SetXTitle("yHCAL (m)");
    heratioy->Draw("COLZ");
    profXy1->SetMarkerStyle(20);
    profXy1->SetLineColor(kRed);
    profXy1->SetMarkerColor(kRed);
    profXy1->SetMarkerSize(1.5);
    profXy1->Draw("SAME");

    // ---- Canvas cc2 ----
    cc2->Divide(1,1);
    cc2->cd(1);
    her->SetXTitle("eHCAL/(sqrt(mp*mp+trP_sbs*trP_sbs)-mp)");
    her->Draw();

    // ---- Canvas c3 ----
    c3->Divide(2,2);
    c3->cd(1);
    hfsvblkid->Draw("COLZ");
    profX2->SetMarkerStyle(20);
    profX2->SetLineColor(kRed);
    profX2->SetMarkerColor(kRed);
    // Increase marker size
    profX2->SetMarkerSize(1.5);
    profX2->Draw("SAME");

    // ---- Canvas c4 ----
    c4->Divide(2,2);
    c4->cd(1);
    hdeltavertex->SetXTitle("vz_bb-vz_sbs (m)");
    hdeltavertex->Draw();
    c4->cd(2);
    hvertexcorr->SetXTitle("vz_bb (m)");
    hvertexcorr->SetYTitle("vz_sbs (m)");
    hvertexcorr->Draw("COLZ");
    c4->cd(3);
    hcointime->Draw();
    c4->cd(4);
    hW2->Draw();

    // =====================================
    // Canvas for Gaussian Fits on hfs & her
    // =====================================
    TCanvas *c5 = new TCanvas("c5","Gaussian Fits: hfs & her",1200,600);
    c5->Divide(2,1);

    // Left pad: Fit hfs
    c5->cd(1);
    hfs->Draw();
    TF1 *gaus_hfs = new TF1("gaus_hfs","gaus",hfs->GetMean()-1.5*hfs->GetRMS(),hfs->GetMean()+1.5*hfs->GetRMS());
    hfs->Fit(gaus_hfs,"R"); 
    gaus_hfs->SetLineColor(kRed);
    gaus_hfs->SetLineWidth(2);

    // (Optional) Move stats box if overlapping your data
    TPaveStats *st1 = (TPaveStats*)hfs->FindObject("stats");
    if(st1) {
        st1->SetX1NDC(0.15); 
        st1->SetX2NDC(0.45);
        st1->SetY1NDC(0.65);
        st1->SetY2NDC(0.88);
    }

    // Right pad: Fit her
    c5->cd(2);
    her->Draw();
    TF1 *gaus_her = new TF1("gaus_her","gaus",her->GetMean()-0.02-1.2*her->GetRMS(),her->GetMean()-0.02+1.2*her->GetRMS());
    her->Fit(gaus_her,"R");
    gaus_her->SetLineColor(kBlue);
    gaus_her->SetLineWidth(2);

    // (Optional) Move stats box for her
    TPaveStats *st2 = (TPaveStats*)her->FindObject("stats");
    if(st2) {
        st2->SetX1NDC(0.55);
        st2->SetX2NDC(0.85);
        st2->SetY1NDC(0.65);
        st2->SetY2NDC(0.88);
    }

    TCanvas *c6 = new TCanvas("c6","c6",3200,2400);
    c6->Divide(2,2);
    c6->cd(1);
    hnblk_hcal->Draw();
    hnblk_hcal->SetXTitle("number of blocks");
    c6->cd(2);
    hseedsumratio->Draw();
    hseedsumratio->SetXTitle("seed block energy / cluster energy");
    c6->cd(3);
    hseedsecondratio->Draw();
    hseedsecondratio->SetXTitle("secondary block energy / seed block energy");

    TCanvas *c7 = new TCanvas("c7","c7",3200,2400);
    c7->Divide(2,2);
    c7->cd(1);
    hxhcal->Draw();
    hxhcal->SetXTitle("x position of the energy weighted centroid of the primary cluster (m)");
    c7->cd(2);
    hyhcal->Draw();
    hyhcal->SetYTitle("y position of the energy weighted centroid of the primary cluster (m)");

    TCanvas *c8 = new TCanvas("c8","c8",3200,2400);
    c8->Divide(2,2);
    c8->cd(1);
    hehcal->Draw();
    hehcal->SetXTitle("eHCAL(GeV)");
    c8->cd(2);
    hekin->Draw();
    hekin->SetXTitle("sqrt(mp*mp+trP_sbs*trP_sbs)-mp");

    // Save individual canvases as PNG (if desired)
    c->SaveAs(Form("../plots/%s_test_fs.png",printfilename));
    cc->SaveAs(Form("../plots/%s_test_fs_1D.png",printfilename));
    c2->SaveAs(Form("../plots/%s_test_fs_sbstracking.png",printfilename));
    cc2->SaveAs(Form("../plots/%s_test_fs_sbstracking_1D.png",printfilename));
    c4->SaveAs(Form("../plots/%s_test_fs_goodtrackscuts.png",printfilename));
    c5->SaveAs(Form("../plots/%s_test_fs_gausfits.png",printfilename));

    // Multi-page PDF output
    // Open with '(' and close with ')'
    c->Print(Form("../plots/%s_test_fs.pdf(",printfilename));
    cc->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    c1->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    c2->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    cc2->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    c3->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    c5->Print(Form("../plots/%s_test_fs.pdf",printfilename)); // Gauss fits page
    c8->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    c6->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    c7->Print(Form("../plots/%s_test_fs.pdf",printfilename));
    c4->Print(Form("../plots/%s_test_fs.pdf)",printfilename));

    TFile *hcal_file = new TFile(Form("../rootfiles/%s_Hcal.root",printfilename), "RECREATE");
    
    hfsample->Write();
    hfsampley->Write();
    heratio->Write();
    heratioy->Write();

    hcal_file->Close();


    // ============================================================
    // ADDITIONAL CODE: Produce 1D distributions for each X bin
    // but now in 5×5 sub-pads per page
    // ============================================================

    // 1) hfsample -> for each bin along xHCAL, project the sampling fraction distribution
    {
        TCanvas *cx_hfsample = new TCanvas("cx_hfsample","hfsample X-Bin Projections (5x5)",1200,900);
        // Open multi-page PDF with "("
        cx_hfsample->Print(Form("../plots/%s_hfsample_xbins.pdf(",printfilename));

        int nbinsX = hfsample->GetNbinsX();
        int nPage  = 0;
        for(int iBin = 1; iBin <= nbinsX; iBin++){
            // Every 25 bins, start a new canvas page
            if((iBin-1) % 25 == 0){
                // Clear and subdivide
                cx_hfsample->Clear();
                cx_hfsample->Divide(5,5);
                nPage++;
            }
            
            // Select the pad [1..25]
            int iPad = (iBin-1) % 25 + 1;
            cx_hfsample->cd(iPad);

            // Project just this bin along Y
            TH1D *projY = hfsample->ProjectionY(
                Form("hfsample_projY_bin%d_page%d", iBin, nPage),
                iBin, iBin
            );
            projY->SetTitle(Form("hfsample: X-bin %d (%.2f<x<%.2f)",
                                 iBin,
                                 hfsample->GetXaxis()->GetBinLowEdge(iBin),
                                 hfsample->GetXaxis()->GetBinUpEdge(iBin)));
            projY->GetXaxis()->SetTitle("Sampling Fraction");
            projY->Draw();

            // If we've just filled the 25th pad or reached the last bin, print this page
            if(((iBin % 25) == 0) || (iBin == nbinsX)){
                cx_hfsample->Print(Form("../plots/%s_hfsample_xbins.pdf",printfilename));
            }
        }

        // Close multi-page PDF with ")"
        cx_hfsample->Print(Form("../plots/%s_hfsample_xbins.pdf)",printfilename));
        delete cx_hfsample;
    }


    // 2) hfsampley -> for each bin along the x-axis (which is yHCAL), project the sampling fraction
    {
        TCanvas *cx_hfsampley = new TCanvas("cx_hfsampley","hfsampley X-Bin Projections (5x5)",1200,900);
        cx_hfsampley->Print(Form("../plots/%s_hfsampley_xbins.pdf(",printfilename));

        int nbinsX = hfsampley->GetNbinsX();  // x-axis is yHCAL
        int nPage  = 0;
        for(int iBin = 1; iBin <= nbinsX; iBin++){
            if((iBin-1) % 25 == 0){
                cx_hfsampley->Clear();
                cx_hfsampley->Divide(5,5);
                nPage++;
            }
            int iPad = (iBin-1) % 25 + 1;
            cx_hfsampley->cd(iPad);

            TH1D *projY = hfsampley->ProjectionY(
                Form("hfsampley_projY_bin%d_page%d", iBin, nPage),
                iBin, iBin
            );
            projY->SetTitle(Form("hfsampley: Y-bin %d (%.2f<y<%.2f)",
                                 iBin,
                                 hfsampley->GetXaxis()->GetBinLowEdge(iBin),
                                 hfsampley->GetXaxis()->GetBinUpEdge(iBin)));
            projY->GetXaxis()->SetTitle("Sampling Fraction");
            projY->Draw();

            if(((iBin % 25) == 0) || (iBin == nbinsX)){
                cx_hfsampley->Print(Form("../plots/%s_hfsampley_xbins.pdf",printfilename));
            }
        }
        cx_hfsampley->Print(Form("../plots/%s_hfsampley_xbins.pdf)",printfilename));
        delete cx_hfsampley;
    }


    // 3) heratio -> for each bin along xHCAL, project the ratio distribution
    {
        TCanvas *cx_heratio = new TCanvas("cx_heratio","heratio X-Bin Projections (5x5)",1200,900);
        cx_heratio->Print(Form("../plots/%s_heratio_xbins.pdf(",printfilename));

        int nbinsX = heratio->GetNbinsX();
        int nPage  = 0;
        for(int iBin = 1; iBin <= nbinsX; iBin++){
            if((iBin-1) % 25 == 0){
                cx_heratio->Clear();
                cx_heratio->Divide(5,5);
                nPage++;
            }
            int iPad = (iBin-1) % 25 + 1;
            cx_heratio->cd(iPad);

            TH1D *projY = heratio->ProjectionY(
                Form("heratio_projY_bin%d_page%d", iBin, nPage),
                iBin, iBin
            );
            projY->SetTitle(Form("heratio: X-bin %d (%.2f<x<%.2f)",
                                 iBin,
                                 heratio->GetXaxis()->GetBinLowEdge(iBin),
                                 heratio->GetXaxis()->GetBinUpEdge(iBin)));
            projY->GetXaxis()->SetTitle("HCAL energy / Expected KE");
            projY->Draw();

            if(((iBin % 25) == 0) || (iBin == nbinsX)){
                cx_heratio->Print(Form("../plots/%s_heratio_xbins.pdf",printfilename));
            }
        }
        cx_heratio->Print(Form("../plots/%s_heratio_xbins.pdf)",printfilename));
        delete cx_heratio;
    }


    // 4) heratioy -> for each bin along yHCAL, project the ratio distribution
    {
        TCanvas *cx_heratioy = new TCanvas("cx_heratioy","heratioy X-Bin Projections (5x5)",1200,900);
        cx_heratioy->Print(Form("../plots/%s_heratioy_xbins.pdf(",printfilename));

        int nbinsX = heratioy->GetNbinsX(); // x-axis is yHCAL for heratioy
        int nPage  = 0;
        for(int iBin = 1; iBin <= nbinsX; iBin++){
            if((iBin-1) % 25 == 0){
                cx_heratioy->Clear();
                cx_heratioy->Divide(5,5);
                nPage++;
            }
            int iPad = (iBin-1) % 25 + 1;
            cx_heratioy->cd(iPad);

            TH1D *projY = heratioy->ProjectionY(
                Form("heratioy_projY_bin%d_page%d", iBin, nPage),
                iBin, iBin
            );
            projY->SetTitle(Form("heratioy: YHCAL-bin %d (%.2f<y<%.2f)",
                                 iBin,
                                 heratioy->GetXaxis()->GetBinLowEdge(iBin),
                                 heratioy->GetXaxis()->GetBinUpEdge(iBin)));
            projY->GetXaxis()->SetTitle("HCAL energy / Expected KE");
            projY->Draw();

            if(((iBin % 25) == 0) || (iBin == nbinsX)){
                cx_heratioy->Print(Form("../plots/%s_heratioy_xbins.pdf",printfilename));
            }
        }
        cx_heratioy->Print(Form("../plots/%s_heratioy_xbins.pdf)",printfilename));
        delete cx_heratioy;
    }


}
