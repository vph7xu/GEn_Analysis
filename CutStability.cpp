
// CutStability.cpp — STREAMING single-variable cut stability plots for asymmetry
// ------------------------------------------------------------------------------
// - Streams directly from DATA TTree (no per-event arrays), so memory stays low.
// - Branches used: vz, ePS, eSH, trP, dx, dy, W2, coin_time, eHCAL, helicity, IHWP.
// - Computes raw helicity asymmetry A = (N+ - N-) / (N+ + N-) for each scan point,
//   with σ_A = sqrt((1 - A^2) / (N+ + N-)).
// - Scans one variable at a time: dyL, dyH, W2L, W2H, eHCAL_L, tL, tH, dxL, dxH.
// - Saves multi-panel canvas and per-scan CSVs under plots/.
//
// Build/Run:
//   root -l -q 'CutStability.cpp("data.root","sim.root","bkg.root","kinA")'
// or interactively:
//   root -l
//   root [0] .L CutStability.cpp+
//   root [1] CutStability("data.root","sim.root","bkg.root","kinA")
// ------------------------------------------------------------------------------

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>
#include <TLatex.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <string>

// -----------------------------
// Result container per scan
// -----------------------------
struct ScanResult {
  std::vector<double> x;      // cut values
  std::vector<double> A;      // asymmetry
  std::vector<double> Aerr;   // uncertainty
  std::vector<long long> Np;  // Nplus
  std::vector<long long> Nm;  // Nminus
};

// ------------------------------------------------------
// Helper: build a TGraphErrors with consistent styling, values in %
// ------------------------------------------------------
static TGraphErrors* make_graph(const char* title,
                                const std::vector<double>& xs,
                                const std::vector<double>& As,
                                const std::vector<double>& Aerrs)
{
  auto g = new TGraphErrors(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) {
    g->SetPoint(i, xs[i], As[i] * 100.0);  // Convert to percentage
    g->SetPointError(i, 0.0, Aerrs[i] * 100.0);
  }
  g->SetTitle(title);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlue+2);
  g->SetLineWidth(2);
  return g;
}


// Print simple percent progress on one line (carriage return).
inline void print_percent_progress(Long64_t i, Long64_t N, const char* label = nullptr) {
  if (N == 0) return;
  // Update every ~1M entries or 0.5% (tune these to your dataset size)
  const Long64_t step = std::max<Long64_t>(1, N / 200); // ~0.5% steps
  if ((i % step) == 0 || i + 1 == N) {
    const double pct = 100.0 * double(i + 1) / double(N);
    if (label) {
      std::cout << label << " " << std::fixed << std::setprecision(1)
                << pct << "%\r";
    } else {
      std::cout << std::fixed << std::setprecision(1)
                << pct << "%\r";
    }
    std::cout.flush();
  }
}


// ------------------------------------------------------
// Helper: build CSV filename safely
// ------------------------------------------------------
static std::string make_csv_filename(const char* kin, const char* name) {
  return std::string(Form("plots/%s_%s.csv", name, kin));
}

// ------------------------------------------------------
// STREAMING asymmetry scan: loops over TTree without storing events
// varType: "dyL","dyH","W2L","W2H","eL","tL","tH","dxL","dxH"
// ------------------------------------------------------
ScanResult StreamScanA(TTree* tData,
                       const std::vector<double>& xvals,
                       const char* varType,
                       // Baseline (fixed) cuts
                       double dyL0, double dyH0,
                       double W2L0, double W2H0,
                       double eL0,
                       double tL0, double tH0,
                       double dxL0, double dxH0)
{
  // Branch variables (float for speed/memory)
  double vz, ePS, eSH, trP, dx, dy, W2, ct, eHCAL;
  int helicity, IHWP;

  // Set addresses
  tData->SetBranchAddress("vz",       &vz);
  tData->SetBranchAddress("ePS",      &ePS);
  tData->SetBranchAddress("eSH",      &eSH);
  tData->SetBranchAddress("trP",      &trP);
  tData->SetBranchAddress("dx",       &dx);
  tData->SetBranchAddress("dy",       &dy);
  tData->SetBranchAddress("W2",       &W2);
  tData->SetBranchAddress("coin_time",&ct);
  tData->SetBranchAddress("eHCAL",    &eHCAL);
  tData->SetBranchAddress("helicity", &helicity);
  tData->SetBranchAddress("IHWP",     &IHWP);

  const Long64_t nData = 1.00 * tData->GetEntries();

  // Counters per scan point
  std::vector<long long> Nplus(xvals.size(), 0), Nminus(xvals.size(), 0);

  // Stream over all events once
  for (Long64_t i = 0; i < nData; ++i) {
    tData->GetEntry(i);


    // Update in-place percent
    print_percent_progress(i, nData, "Streaming events:");


    // Preserve the original "reject only if ALL sub-conditions are true" preselection
    if (std::abs(vz) > 0.27 && ePS < 0.2 &&
        std::abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL < 0.025 &&
        std::abs(helicity) != 1)
      continue;

    // Fixed baseline cuts (except the scanned variable)
    bool baseCommon = std::abs(vz) < 0.27 &&
                      ePS > 0.2 &&
                      W2 > W2L0 && W2 < W2H0 &&
                      eHCAL > eL0;
    bool inCT = ct > tL0 && ct < tH0;
    bool inDX = dx > dxL0 && dx < dxH0;

    if (!baseCommon || !inCT || !inDX) continue;

    // Helicity sign with IHWP
    bool poshel = (IHWP * helicity == 1);
    bool neghel = (IHWP * helicity == -1);

    const std::string vtype(varType);

    if (vtype == "dyL") {
      // dy lower edge scan: dy > x (ensure dy < dyH0)
      if (dy < dyH0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (dy > xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    } else if (vtype == "dyH") {
      // dy upper edge scan: dy < x (ensure dy > dyL0)
      if (dy > dyL0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (dy < xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    } else if (vtype == "W2L") {
      // W2 lower edge scan: W2 > x (ensure W2 < W2H0)
      if (W2 < W2H0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (W2 > xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    } else if (vtype == "W2H") {
      // W2 upper edge scan: W2 < x (ensure W2 > W2L0)
      if (W2 > W2L0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (W2 < xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    } else if (vtype == "eL") {
      // eHCAL low threshold scan: eHCAL > x
      for (size_t k = 0; k < xvals.size(); ++k)
        if (eHCAL > xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
    } else if (vtype == "tL") {
      // coin_time lower edge scan: ct > x (ensure ct < tH0)
      if (ct < tH0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (ct > xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    } else if (vtype == "tH") {
      // coin_time upper edge scan: ct < x (ensure ct > tL0)
      if (ct > tL0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (ct < xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    } else if (vtype == "dxL") {
      // dx lower edge scan: dx > x (ensure dx < dxH0)
      if (dx < dxH0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (dx > xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    } else if (vtype == "dxH") {
      // dx upper edge scan: dx < x (ensure dx > dxL0)
      if (dx > dxL0) {
        for (size_t k = 0; k < xvals.size(); ++k)
          if (dx < xvals[k]) { if (poshel) ++Nplus[k]; else if (neghel) ++Nminus[k]; }
      }
    }
  } // end streaming loop

  std::cout << std::endl; // finish the line

  // Compose results
  ScanResult out;
  out.x = xvals;
  out.Np = Nplus;
  out.Nm = Nminus;
  out.A.resize(xvals.size());
  out.Aerr.resize(xvals.size());

  for (size_t k = 0; k < xvals.size(); ++k) {
    const double Np = (double)Nplus[k];
    const double Nm = (double)Nminus[k];
    const double N  = Np + Nm;
    if (N > 0.0) {
      out.A[k]    = (Np - Nm) / N;
      out.Aerr[k] = std::sqrt(std::max(0.0, (1.0 - out.A[k]*out.A[k]) / N));
    } else {
      out.A[k]    = 0.0;
      out.Aerr[k] = 1.0;
    }
  }
  return out;
}

// ------------------------------------------------------
// Helper: write CSV file (no lambda to avoid HTML escapes)
// ------------------------------------------------------
static void write_csv_file(const char* kin, const char* base, const ScanResult& r) {
  const std::string fname = make_csv_filename(kin, base);
  std::ofstream csv(fname.c_str());
  csv << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (size_t i=0;i<r.x.size();++i)
    csv << r.x[i] << "," << r.A[i] << "," << r.Aerr[i] << ","
        << r.Np[i] << "," << r.Nm[i] << "\n";
  csv.close();
}

// ------------------------------------------------------
// Main entry: loads data and runs streaming stability scans
// ------------------------------------------------------
void CutStability(const char* dataFile,
                  const char* kin)
{
  gStyle->SetOptFit(0);
  gSystem->mkdir("plots", kTRUE);

  // Open data
  TFile fData(dataFile, "READ");
  TTree* tData = (TTree*)fData.Get("Tout");
  if (!tData) {
    std::cerr << "ERROR: cannot find Tout in data file " << dataFile << std::endl;
    return;
  }
  std::cout << "Streaming from data TTree: " << dataFile << std::endl;

  // -----------------------------
  // 1) Define scan ranges
  // -----------------------------
  const int nSteps = 10;

  std::vector<double> dyL_vals, dyH_vals,
                      W2L_vals, W2H_vals, eL_vals,
                      tL_vals, tH_vals,
                      dxL_vals, dxH_vals;

  for (int i = 0; i <= nSteps; ++i) dyL_vals.push_back(-0.8 + i*(0.8/nSteps));
  for (int i = 0; i <= nSteps; ++i) dyH_vals.push_back( 0.8 - i*(0.8/nSteps));
  for (int i = 0; i <= nSteps; ++i) W2L_vals.push_back(-2.0 + i*(2.0/nSteps));
  for (int i = 0; i <= nSteps; ++i) W2H_vals.push_back( 1.0 + i*(1.0/nSteps));
  for (int i = 0; i <= nSteps; ++i) eL_vals.push_back( 0.025 + i*(0.5/nSteps));
  for (int i = 0; i <= nSteps; ++i) tL_vals.push_back(130.0 - i*(6.0/nSteps));
  for (int i = 0; i <= nSteps; ++i) tH_vals.push_back(130.0 + i*(6.0/nSteps));
  for (int i = 0; i <= nSteps; ++i) dxL_vals.push_back(-0.8 + i*(0.8/nSteps)); // -0.8 to 0.0
  for (int i = 0; i <= nSteps; ++i) dxH_vals.push_back( 0.8 - i*(0.8/nSteps)); //  0.8 to 0.0

  // -----------------------------
  // 2) Baseline cuts (working point)
  // -----------------------------
  double dyL0 = -0.8;
  double dyH0 =  0.8;
  double W2L0 = -2.0;
  double W2H0 =  2.0;
  double eL0  =  0.025;
  double tL0  = 124.0;
  double tH0  = 136.0;
  double dxL0 = -0.8;
  double dxH0 =  0.8;

  std::cout << "Baseline cuts:\n"
            << " dy in [" << dyL0 << ", " << dyH0 << "]\n"
            << " W2 in [" << W2L0 << ", " << W2H0 << "]\n"
            << " eHCAL_L = " << eL0 << "\n"
            << " coin_time in [" << tL0 << ", " << tH0 << "]\n"
            << " dx in [" << dxL0 << ", " << dxH0 << "]\n";

  // -----------------------------
  // 3) Run streaming scans and produce two-panel side-by-side plots
  // -----------------------------

  // small helper to set bold, larger axis fonts on a graph
  auto style_graph = [](TGraphErrors* g){
    g->GetXaxis()->SetTitleFont(62);
    g->GetXaxis()->SetLabelFont(62);
    g->GetXaxis()->SetTitleSize(0.05);
    g->GetXaxis()->SetLabelSize(0.045);
    g->GetYaxis()->SetTitleFont(62);
    g->GetYaxis()->SetLabelFont(62);
    g->GetYaxis()->SetTitleSize(0.05);
    g->GetYaxis()->SetLabelSize(0.045);
    g->GetYaxis()->SetTitleOffset(1.2);
  };

  // Helper to add horizontal Y-axis title on left pad at top
  auto add_y_label = [](const char* label) {
    TLatex tex;
    tex.SetNDC(true);
    tex.SetTextFont(62);
    tex.SetTextSize(0.045);
    tex.SetTextAlign(11); // left-aligned, centered
    tex.DrawLatex(0.08, 0.95, label);
  };

  // --- dyL scan
  // std::cout << "[Scan] dyL (" << dyL_vals.size() << " points) ..." << std::endl;
  // auto r_dyL = StreamScanA(tData, dyL_vals, "dyL",
  //                          dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  // auto g_dyL = make_graph("", r_dyL.x, r_dyL.A, r_dyL.Aerr);
  // g_dyL->GetXaxis()->SetTitle("dy_{L}");
  // g_dyL->GetYaxis()->SetTitle("");
  // style_graph(g_dyL);
  // TF1 fconst_dyL("fconst_dyL","pol0",
  //   *std::min_element(r_dyL.x.begin(), r_dyL.x.end()),
  //   *std::max_element(r_dyL.x.begin(), r_dyL.x.end()));
  // g_dyL->Fit(&fconst_dyL, "Q");

  // // --- dyH scan
  // std::cout << "[Scan] dyH (" << dyH_vals.size() << " points) ..." << std::endl;
  // auto r_dyH = StreamScanA(tData, dyH_vals, "dyH",
  //                          dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  // auto g_dyH = make_graph("", r_dyH.x, r_dyH.A, r_dyH.Aerr);
  // g_dyH->GetXaxis()->SetTitle("dy_{H}");
  // g_dyH->GetYaxis()->SetTitle("");
  // style_graph(g_dyH);
  // TF1 fconst_dyH("fconst_dyH","pol0",
  //   *std::min_element(r_dyH.x.begin(), r_dyH.x.end()),
  //   *std::max_element(r_dyH.x.begin(), r_dyH.x.end()));
  // g_dyH->Fit(&fconst_dyH, "Q");

  // --- W2L scan
  std::cout << "[Scan] W2L (" << W2L_vals.size() << " points) ..." << std::endl;
  auto r_W2L = StreamScanA(tData, W2L_vals, "W2L",
                           dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  auto g_W2L = make_graph("", r_W2L.x, r_W2L.A, r_W2L.Aerr);
  g_W2L->GetXaxis()->SetTitle("W^{2}_{L}");
  g_W2L->GetYaxis()->SetTitle("");
  style_graph(g_W2L);
  // TF1 fconst_W2L("fconst_W2L","pol0",
  //   *std::min_element(r_W2L.x.begin(), r_W2L.x.end()),
  //   *std::max_element(r_W2L.x.begin(), r_W2L.x.end()));
  // g_W2L->Fit(&fconst_W2L, "Q");

  // --- W2H scan
  std::cout << "[Scan] W2H (" << W2H_vals.size() << " points) ..." << std::endl;
  auto r_W2H = StreamScanA(tData, W2H_vals, "W2H",
                           dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  auto g_W2H = make_graph("", r_W2H.x, r_W2H.A, r_W2H.Aerr);
  g_W2H->GetXaxis()->SetTitle("W^{2}_{H}");
  g_W2H->GetYaxis()->SetTitle("");
  style_graph(g_W2H);
  // TF1 fconst_W2H("fconst_W2H","pol0",
  //   *std::min_element(r_W2H.x.begin(), r_W2H.x.end()),
  //   *std::max_element(r_W2H.x.begin(), r_W2H.x.end()));
  // g_W2H->Fit(&fconst_W2H, "Q");

  // --- eHCAL low threshold (eL) scan
  std::cout << "[Scan] eL (" << eL_vals.size() << " points) ..." << std::endl;
  auto r_eL = StreamScanA(tData, eL_vals, "eL",
                          dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  auto g_eL = make_graph("", r_eL.x, r_eL.A, r_eL.Aerr);
  g_eL->GetXaxis()->SetTitle("eHCAL low cut");
  g_eL->GetYaxis()->SetTitle("");
  style_graph(g_eL);
  // TF1 fconst_eL("fconst_eL","pol0",
  //   *std::min_element(r_eL.x.begin(), r_eL.x.end()),
  //   *std::max_element(r_eL.x.begin(), r_eL.x.end()));
  // g_eL->Fit(&fconst_eL, "Q");

  // --- tL scan (coincidence time lower edge)
  // std::cout << "[Scan] tL (" << tL_vals.size() << " points) ..." << std::endl;
  // auto r_tL = StreamScanA(tData, tL_vals, "tL",
  //                         dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  // auto g_tL = make_graph("", r_tL.x, r_tL.A, r_tL.Aerr);
  // g_tL->GetXaxis()->SetTitle("coincidence time lower edge");
  // g_tL->GetYaxis()->SetTitle("");
  // style_graph(g_tL);
  // TF1 fconst_tL("fconst_tL","pol0",
  //   *std::min_element(r_tL.x.begin(), r_tL.x.end()),
  //   *std::max_element(r_tL.x.begin(), r_tL.x.end()));
  // g_tL->Fit(&fconst_tL, "Q");

  // // --- tH scan (coincidence time upper edge)
  // std::cout << "[Scan] tH (" << tH_vals.size() << " points) ..." << std::endl;
  // auto r_tH = StreamScanA(tData, tH_vals, "tH",
  //                         dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  // auto g_tH = make_graph("", r_tH.x, r_tH.A, r_tH.Aerr);
  // g_tH->GetXaxis()->SetTitle("coincidence time upper edge");
  // g_tH->GetYaxis()->SetTitle("");
  // style_graph(g_tH);
  // TF1 fconst_tH("fconst_tH","pol0",
  //   *std::min_element(r_tH.x.begin(), r_tH.x.end()),
  //   *std::max_element(r_tH.x.begin(), r_tH.x.end()));
  // g_tH->Fit(&fconst_tH, "Q");

  // --- dxL scan (lower edge of dx window)
  // std::cout << "[Scan] dxL (" << dxL_vals.size() << " points) ..." << std::endl;
  // auto r_dxL = StreamScanA(tData, dxL_vals, "dxL",
  //                          dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  // auto g_dxL = make_graph("", r_dxL.x, r_dxL.A, r_dxL.Aerr);
  // g_dxL->GetXaxis()->SetTitle("dx_{L}");
  // g_dxL->GetYaxis()->SetTitle("");
  // style_graph(g_dxL);
  // TF1 fconst_dxL("fconst_dxL","pol0",
  //   *std::min_element(r_dxL.x.begin(), r_dxL.x.end()),
  //   *std::max_element(r_dxL.x.begin(), r_dxL.x.end()));
  // g_dxL->Fit(&fconst_dxL, "Q");

  // // --- dxH scan (upper edge of dx window)
  // std::cout << "[Scan] dxH (" << dxH_vals.size() << " points) ..." << std::endl;
  // auto r_dxH = StreamScanA(tData, dxH_vals, "dxH",
  //                          dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
  // auto g_dxH = make_graph("", r_dxH.x, r_dxH.A, r_dxH.Aerr);
  // g_dxH->GetXaxis()->SetTitle("dx_{H}");
  // g_dxH->GetYaxis()->SetTitle("");
  // style_graph(g_dxH);
  // TF1 fconst_dxH("fconst_dxH","pol0",
  //   *std::min_element(r_dxH.x.begin(), r_dxH.x.end()),
  //   *std::max_element(r_dxH.x.begin(), r_dxH.x.end()));
  // g_dxH->Fit(&fconst_dxH, "Q");

  // create and save two-panel canvases: left=low, right=high (right empty for eHCAL_L)

  // dy canvas
  // TCanvas* c_dy = new TCanvas(Form("c_dy_%s", kin), "dy scans", 1400, 700);
  // c_dy->Divide(2,1);
  // c_dy->SetLeftMargin(0.14);
  // c_dy->SetRightMargin(0.05);
  // c_dy->SetTopMargin(0.08);
  // c_dy->SetBottomMargin(0.12);
  // c_dy->cd(1);
  // gPad->SetLeftMargin(0.22);
  // gPad->SetRightMargin(0.03);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_dyL->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_dy->cd(2);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.06);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_dyH->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_dy->SaveAs(Form("plots/asym_dy_%s.png", kin));
  // c_dy->SaveAs(Form("plots/asym_dy_%s.jpg", kin));
  

  // W2 canvas
  TCanvas* c_W2 = new TCanvas(Form("c_W2_%s", kin), "W2 scans", 1400, 700);
  c_W2->Divide(2,1);
  c_W2->SetLeftMargin(0.14);
  c_W2->SetRightMargin(0.05);
  c_W2->SetTopMargin(0.08);
  c_W2->SetBottomMargin(0.12);
  c_W2->cd(1);
  gPad->SetLeftMargin(0.22);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  g_W2L->Draw("AP");
  add_y_label("Asymmetry A (%)");
  c_W2->cd(2);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  g_W2H->Draw("AP");
  add_y_label("Asymmetry A (%)");
  c_W2->SaveAs(Form("plots/asym_W2_%s.png", kin));
  c_W2->SaveAs(Form("plots/asym_W2_%s.jpg", kin));
  

  // eHCAL canvas (keep right pad empty)
  TCanvas* c_eL = new TCanvas(Form("c_eL_%s", kin), "eHCAL low scan", 1400, 700);
  c_eL->Divide(2,1);
  c_eL->SetLeftMargin(0.14);
  c_eL->SetRightMargin(0.05);
  c_eL->SetTopMargin(0.08);
  c_eL->SetBottomMargin(0.12);
  c_eL->cd(1);
  gPad->SetLeftMargin(0.22);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  g_eL->Draw("AP");
  add_y_label("Asymmetry A (%)");
  // leave pad 2 empty
  c_eL->SaveAs(Form("plots/asym_eHCAL_L_%s.png", kin));
  c_eL->SaveAs(Form("plots/asym_eHCAL_L_%s.jpg", kin));
  

  // // t canvas
  // TCanvas* c_t = new TCanvas(Form("c_t_%s", kin), "time scans", 1400, 700);
  // c_t->Divide(2,1);
  // c_t->SetLeftMargin(0.14);
  // c_t->SetRightMargin(0.05);
  // c_t->SetTopMargin(0.08);
  // c_t->SetBottomMargin(0.12);
  // c_t->cd(1);
  // gPad->SetLeftMargin(0.22);  // // t canvas
  // TCanvas* c_t = new TCanvas(Form("c_t_%s", kin), "time scans", 1400, 700);
  // c_t->Divide(2,1);
  // c_t->SetLeftMargin(0.14);
  // c_t->SetRightMargin(0.05);
  // c_t->SetTopMargin(0.08);
  // c_t->SetBottomMargin(0.12);
  // c_t->cd(1);
  // gPad->SetLeftMargin(0.22);
  // gPad->SetRightMargin(0.03);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_tL->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_t->cd(2);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.06);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_tH->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_t->SaveAs(Form("plots/asym_time_%s.png", kin));
  // c_t->SaveAs(Form("plots/asym_time_%s.jpg", kin));
  

  // // dx canvas
  // TCanvas* c_dx = new TCanvas(Form("c_dx_%s", kin), "dx scans", 1400, 700);
  // c_dx->Divide(2,1);
  // c_dx->SetLeftMargin(0.14);
  // c_dx->SetRightMargin(0.05);
  // c_dx->SetTopMargin(0.08);
  // c_dx->SetBottomMargin(0.12);
  // c_dx->cd(1);
  // gPad->SetLeftMargin(0.22);
  // gPad->SetRightMargin(0.03);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_dxL->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_dx->cd(2);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.06);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_dxH->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_dx->SaveAs(Form("plots/asym_dx_%s.png", kin));
  // c_dx->SaveAs(Form("plots/asym_dx_%s.jpg", kin));
  // gPad->SetRightMargin(0.03);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_tL->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_t->cd(2);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.06);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_tH->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_t->SaveAs(Form("plots/asym_time_%s.png", kin));
  // c_t->SaveAs(Form("plots/asym_time_%s.jpg", kin));
  

  // // dx canvas
  // TCanvas* c_dx = new TCanvas(Form("c_dx_%s", kin), "dx scans", 1400, 700);
  // c_dx->Divide(2,1);
  // c_dx->SetLeftMargin(0.14);
  // c_dx->SetRightMargin(0.05);
  // c_dx->SetTopMargin(0.08);
  // c_dx->SetBottomMargin(0.12);
  // c_dx->cd(1);
  // gPad->SetLeftMargin(0.22);
  // gPad->SetRightMargin(0.03);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_dxL->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_dx->cd(2);
  // gPad->SetLeftMargin(0.08);
  // gPad->SetRightMargin(0.06);
  // gPad->SetTopMargin(0.08);
  // gPad->SetBottomMargin(0.12);
  // g_dxH->Draw("AP");
  // add_y_label("Asymmetry A (%)");
  // c_dx->SaveAs(Form("plots/asym_dx_%s.png", kin));
  // c_dx->SaveAs(Form("plots/asym_dx_%s.jpg", kin));
  
  // Also create a single combined PDF containing all canvases
  const TString combined_pdf = Form("plots/asym_all_%s.pdf", kin);
  //c_dy->SaveAs(Form("%s(", combined_pdf.Data()));
  c_W2->SaveAs(Form("%s(", combined_pdf.Data()));
  c_eL->SaveAs(Form("%s)", combined_pdf.Data()));
  //c_t->SaveAs(combined_pdf);
  //c_dx->SaveAs(Form("%s)", combined_pdf.Data()));

  // write_csv_file(kin, "asym_scan_dyL",      r_dyL);
  // write_csv_file(kin, "asym_scan_dyH",      r_dyH);
  write_csv_file(kin, "asym_scan_W2L",      r_W2L);
  write_csv_file(kin, "asym_scan_W2H",      r_W2H);
  write_csv_file(kin, "asym_scan_eHCAL_L",  r_eL);
  // write_csv_file(kin, "asym_scan_tL",       r_tL);
  // write_csv_file(kin, "asym_scan_tH",       r_tH);
  // write_csv_file(kin, "asym_scan_dxL",      r_dxL);
  // write_csv_file(kin, "asym_scan_dxH",      r_dxH);

  std::cout << "Saved plots/asym_stability_" << kin
            << ".png/.jpg and per-scan CSVs in plots/." << std::endl;
}
