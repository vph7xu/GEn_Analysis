// plot_thesis_qe.C
// Usage (example):
//   root -l -b -q 'plot_thesis_qe.C("QE_data_GEN2_sbs100p_nucleon_np_model2_sbstrackingon.root","(WCut&&fiduCut&&coinCut&&ntrack>0)")'
//
// If you want no cut, pass "" as the 2nd argument.

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLine.h"
#include "TSystem.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TLegend.h"

static void SetThesisStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetPadRightMargin(0.14);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.06);

  gStyle->SetTitleFont(42, "XYZ");
  gStyle->SetLabelFont(42, "XYZ");
  gStyle->SetTextFont(42);

  gStyle->SetTitleSize(0.050, "XYZ");
  gStyle->SetLabelSize(0.045, "XYZ");
  gStyle->SetTitleOffset(1.10, "X");
  gStyle->SetTitleOffset(1.35, "Y");

  gStyle->SetHistLineWidth(3);
  gStyle->SetLineWidth(2);

  gStyle->SetNumberContours(60);
  gStyle->SetPalette(kRainBow);
}

static void DrawCutLabel(const TString &cut) {
  if (cut.Length() == 0) return;
  TLatex lat;
  lat.SetNDC(true);
  lat.SetTextSize(0.030);
  lat.SetTextAlign(13);
  lat.DrawLatex(0.13, 0.94, Form("Cut: %s", cut.Data()));
}

static void Save1D(TTree *T,
                   const TString &var,
                   const TString &cut,
                   int nbins, double xmin, double xmax,
                   const TString &xtitle,
                   const TString &outpdf,
                   bool logy=false)
{
  TString hname = "h1_" + var;
  TH1D *h = new TH1D(hname, "", nbins, xmin, xmax);
  h->Sumw2();
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle("Counts");

  // Fill
  TString drawcmd = Form("%s>>%s", var.Data(), hname.Data());
  Long64_t n = T->Draw(drawcmd, cut, "goff");
  (void)n;

  // Canvas
  TCanvas *c = new TCanvas("c1_"+var, "c1_"+var, 1000, 750);
  c->SetTicks(1,1);
  if (logy) c->SetLogy(true);

  // For the coin time histogram, fill in gray and overlay a Gaussian fit
  TF1 *fgaus = nullptr;
  if (var.EqualTo("coin_time")) {
    h->SetFillColor(kGray);
    h->SetFillStyle(1001);
    h->SetLineColor(kGray+2);
    h->Draw("hist");

    // Perform Gaussian fit in -5..5 ns (fit range only)
    fgaus = new TF1("fgaus", "gaus", -3.0, 3.0);
    h->Fit(fgaus, "RQ", "", -3.0, 3.0);
    fgaus->SetLineColor(kRed);
    fgaus->SetLineWidth(2);
    fgaus->Draw("same");
  } else {
    h->SetLineColor(kBlack);
    h->Draw("hist");
  }

  // Simple stats box (nice-looking) with entries
  TPaveText *pt = new TPaveText(0.55, 0.80, 0.75, 0.92, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  //pt->AddText(Form("Entries: %.0f", h->GetEntries()));
  if (fgaus) {
    double mean = fgaus->GetParameter(1);
    double mean_err = fgaus->GetParError(1);
    double sigma = fgaus->GetParameter(2);
    double sigma_err = fgaus->GetParError(2);
    pt->AddText(Form("Mean = %.3f #pm %.3f ns", mean, mean_err));
    pt->AddText(Form("#sigma = %.3f #pm %.3f ns", sigma, sigma_err));
  }
  pt->Draw();

  //DrawCutLabel(cut);

  c->SaveAs(outpdf);

  if (fgaus) delete fgaus;
  delete pt;
  delete c;
  delete h;
}

static void SaveCoinTimeComparison(TTree *T,
                                    const TString &cut_tight,
                                    const TString &cut_loose,
                                    int nbins, double xmin, double xmax,
                                    const TString &outpdf)
{
  // Create two histograms: tight cut and loose cut
  TH1D *h_tight = new TH1D("h_coin_tight", "", nbins, xmin, xmax);
  TH1D *h_loose = new TH1D("h_coin_loose", "", nbins, xmin, xmax);
  h_tight->Sumw2();
  h_loose->Sumw2();
  h_tight->GetXaxis()->SetTitle("Coincidence time (ns)");
  h_tight->GetYaxis()->SetTitle("Counts");
  h_loose->GetXaxis()->SetTitle("Coincidence time (ns)");
  h_loose->GetYaxis()->SetTitle("Counts");

  // Fill histograms
  T->Draw("coin_time>>h_coin_tight", cut_tight, "goff");
  T->Draw("coin_time>>h_coin_loose", cut_loose, "goff");

  // Canvas
  TCanvas *c = new TCanvas("c_coin_comparison", "coin_time comparison", 1000, 750);
  c->SetTicks(1,1);
  c->SetLogy(true);

  // Draw loose cut histogram (outline only, different color)
  h_loose->SetFillStyle(0);
  h_loose->SetLineColor(kBlue);
  h_loose->SetLineWidth(2);
  h_loose->Draw("hist");

  // Draw tight cut histogram (gray fill)
  h_tight->SetFillColor(kGray+2);
  h_tight->SetFillStyle(1001);
  h_tight->SetLineColor(kBlack);
  h_tight->Draw("hist same");

  // Fit Gaussian to tight cut in -5..5 ns range
  TF1 *fgaus = new TF1("fgaus", "gaus", -3.0, 3.0);
  h_tight->Fit(fgaus, "RQ", "", -3.0, 3.0);
  fgaus->SetLineColor(kRed);
  fgaus->SetLineWidth(2);
  fgaus->Draw("same");

  // Add legend
  TLegend *leg = new TLegend(0.13, 0.80, 0.35, 0.92);
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.035);
  leg->AddEntry(h_tight, "Tight cut", "f");
  leg->AddEntry(h_loose, "Loose cut", "l");
  leg->AddEntry(fgaus, "Gaussian fit", "l");
  leg->Draw();

  // Add fit parameters
  TPaveText *pt = new TPaveText(0.55, 0.80, 0.75, 0.92, "NDC");
  pt->SetFillStyle(0);
  pt->SetBorderSize(0);
  pt->SetTextAlign(12);
  pt->SetTextSize(0.04);
  if (fgaus) {
    double mean = fgaus->GetParameter(1);
    double mean_err = fgaus->GetParError(1);
    double sigma = fgaus->GetParameter(2);
    double sigma_err = fgaus->GetParError(2);
    pt->AddText(Form("Mean = %.3f #pm %.3f ns", mean, mean_err));
    pt->AddText(Form("#sigma = %.3f #pm %.3f ns", sigma, sigma_err));
  }
  pt->Draw();

  c->SaveAs(outpdf);

  delete leg;
  delete pt;
  delete fgaus;
  delete c;
  delete h_tight;
  delete h_loose;
}

static void Save2D(TTree *T,
                   const TString &yvar,  // plotted on Y
                   const TString &xvar,  // plotted on X
                   const TString &cut,
                   int nxbins, double xmin, double xmax,
                   int nybins, double ymin, double ymax,
                   const TString &xtitle,
                   const TString &ytitle,
                   const TString &outpdf,
                   bool logz=false)
{
  TString hname = "h2_" + yvar + "_vs_" + xvar;
  TH2D *h = new TH2D(hname, "", nxbins, xmin, xmax, nybins, ymin, ymax);
  h->Sumw2();
  h->GetXaxis()->SetTitle(xtitle);
  h->GetYaxis()->SetTitle(ytitle);

  // Fill
  TString drawcmd = Form("%s:%s>>%s", yvar.Data(), xvar.Data(), hname.Data());
  Long64_t n = T->Draw(drawcmd, cut, "goff");
  (void)n;

  // Canvas
  TCanvas *c = new TCanvas("c2_"+hname, "c2_"+hname, 1050, 800);
  c->SetTicks(1,1);
  if (logz) c->SetLogz(true);

  // Give colorbar enough room
  c->SetRightMargin(0.16);

  h->Draw("colz");

  //DrawCutLabel(cut);

  c->SaveAs(outpdf);

  delete c;
  delete h;
}

void plot_thesis_qe(const char *fname,
                    const char *cut_in = "",
                    const char *outdir = "thesis_plots_qe")
{
  SetThesisStyle();

  gSystem->mkdir(outdir, true);

  TFile *f = TFile::Open(fname, "READ");
  if (!f || f->IsZombie()) {
    Error("plot_thesis_qe", "Could not open file: %s", fname);
    return;
  }

  TTree *T = (TTree*)f->Get("Tout");
  if (!T) {
    Error("plot_thesis_qe", "Could not find tree 'Tout' in file.");
    f->Close();
    return;
  }

  TString cut = cut_in;
  cut = cut.Strip(TString::kBoth);

  int total_plots = 6;
  int current_plot = 0;

  // -------------------------
  // 1D examples (edit as you like)
  // -------------------------
  // Compare coin_time with tight vs loose cut
  current_plot++;
  std::cout << "[" << current_plot << "/" << total_plots << "] Generating coin_time comparison..." << std::endl;
  TString loose_cut = "ntrack>0&&abs(vz)<0.27&&ePS>0.2&&eHCAL>0.1&&abs((eSH+ePS)/trP-1)<0.25";
  SaveCoinTimeComparison(T, cut, loose_cut,
                         250, -20, 20,
                         TString::Format("%s/coin_time_comparison.pdf", outdir));
  std::cout << "  ✓ coin_time_comparison.pdf" << std::endl;

  current_plot++;
  std::cout << "[" << current_plot << "/" << total_plots << "] Generating coin_time..." << std::endl;
  Save1D(T, "coin_time", cut,
         250, -20, 20,
         "Coincidence time (ns)",
         TString::Format("%s/coin_time.pdf", outdir),
         true);
  std::cout << "  ✓ coin_time.pdf" << std::endl;

  current_plot++;
  std::cout << "[" << current_plot << "/" << total_plots << "] Generating W2..." << std::endl;
  Save1D(T, "W2", cut,
         250, -1.0, 5.0,
         "W^{2} (GeV^{2})",
         TString::Format("%s/W2.pdf", outdir),
         false);
  std::cout << "  ✓ W2.pdf" << std::endl;

  current_plot++;
  std::cout << "[" << current_plot << "/" << total_plots << "] Generating dx..." << std::endl;
  Save1D(T, "dx", cut,
         200, -4.0, 4.0,
         "dx (m)",
         TString::Format("%s/dx.pdf", outdir),
         false);
  std::cout << "  ✓ dx.pdf" << std::endl;

  // -------------------------
  // 2D examples (edit as you like)
  // ROOT draw syntax "y:x"
  // -------------------------
  current_plot++;
  std::cout << "[" << current_plot << "/" << total_plots << "] Generating coin_time vs W2..." << std::endl;
  Save2D(T, "coin_time", "W2", loose_cut,
         250, -1.0, 5.0,
         250, -20, 20,
         "W^{2} (GeV^{2})",
         "Coincidence time (ns)",
         TString::Format("%s/coin_time_vs_W2.pdf", outdir),
         true);
  std::cout << "  ✓ coin_time_vs_W2.pdf" << std::endl;

  current_plot++;
  std::cout << "[" << current_plot << "/" << total_plots << "] Generating dx vs dy..." << std::endl;
  Save2D(T, "dx", "dy", loose_cut,
         200, -4.0, 4.0,
         200, -3.0, 3.0,
         "dy (m)",
         "dx (m)",
         TString::Format("%s/dx_vs_dy.pdf", outdir),
         true);
  std::cout << "  ✓ dx_vs_dy.pdf" << std::endl;

  std::cout << "\n[✓] All plots completed! Output directory: " << outdir << std::endl;

  f->Close();
}
