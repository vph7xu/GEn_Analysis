
// CutStability.cpp — single-variable cut stability plots for asymmetry
// --------------------------------------------------------------------
// - Loads DATA (and SIM/BKG for continuity; SIM/BKG not used for A)
// - Reuses your variables/branches: vz, ePS, eSH, trP, dx, dy, W2, coin_time, eHCAL, helicity, IHWP
// - Uses your selection booleans & fixed signal dx window [-0.4, 0.4]
// - Scans one cut at a time (dyL, dyH, W2L, W2H, eHCAL_L, tL, tH) and plots A vs cut
// - Saves multi-panel stability canvas to plots/asym_stability_<kin>.png/.jpg
//
// Build/Run (as a ROOT macro):
//   root -l -q 'CutStability.cpp("data.root","sim.root","bkg.root","kinA")'
// or interactively:
//   root -l
//   root [0] .L CutStability.cpp+
//   root [1] CutStability("data.root","sim.root","bkg.root","kinA")
// --------------------------------------------------------------------

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TH1D.h>
#include <TGraphErrors.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TString.h>

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <string>

// -----------------------------
// Event containers (same names)
// -----------------------------
struct Store {
  std::vector<double> vz, ePS, dx, dy, W2, ct, eHCAL, w;
  std::vector<int>    helicity, IHWP;
  std::vector<char>   origin; // 0=nMC,1=pMC,2=data
};

struct StoreBkg {
  std::vector<double> vz, ePS, dx, dy, W2, ct, eHCAL, w;
};

// -----------------------------
// Asymmetry result container
// -----------------------------
struct AsymmetryResult {
  double cut_value = 0.0;
  double Nplus     = 0.0;
  double Nminus    = 0.0;
  double A         = 0.0;
  double Aerr      = 0.0;
};

// ------------------------------------------------------
// Helper: compute A and sigma_A for a given set of cuts
// ------------------------------------------------------
static AsymmetryResult compute_asymmetry_single(
  const Store& ev,
  double dyL, double dyH,
  double W2L, double W2H,
  double eL,
  double tL, double tH,
  double dxL = -0.4, double dxH = 0.4)
{
  double Nplus = 0.0;
  double Nminus = 0.0;

  // Reuse selection logic and fixed signal dx window
  for (size_t i = 0; i < ev.dx.size(); ++i) {
    if (ev.origin[i] != 2) continue; // data only
    bool base   = std::abs(ev.vz[i]) < 0.27 &&
                  ev.ePS[i] > 0.2 &&
                  ev.W2[i] > W2L && ev.W2[i] < W2H &&
                  ev.eHCAL[i] > eL;
    bool inct   = ev.ct[i] > tL && ev.ct[i] < tH;
    bool inDy   = ev.dy[i] > dyL && ev.dy[i] < dyH;
    bool inDx   = ev.dx[i] > dxL && ev.dx[i] < dxH;      // signal window
    bool poshel = (ev.IHWP[i] * ev.helicity[i] == 1);
    bool neghel = (ev.IHWP[i] * ev.helicity[i] == -1);

    if (base && inct && inDy && inDx && poshel) ++Nplus;
    if (base && inct && inDy && inDx && neghel) ++Nminus;
  }

  AsymmetryResult r;
  r.Nplus = Nplus;
  r.Nminus = Nminus;
  double Nsum = Nplus + Nminus;

  if (Nsum > 0.0) {
    r.A = (Nplus - Nminus) / Nsum;
    // Binomial-like asymmetry uncertainty
    r.Aerr = std::sqrt(std::max(0.0, (1.0 - r.A * r.A) / Nsum));
  } else {
    r.A = 0.0;
    r.Aerr = 1.0;
  }
  return r;
}

// ------------------------------------------------------
// Helper: build a TGraphErrors with consistent styling
// ------------------------------------------------------
static TGraphErrors* make_graph(const char* title,
                                const std::vector<double>& xs,
                                const std::vector<double>& As,
                                const std::vector<double>& Aerrs)
{
  auto g = new TGraphErrors(xs.size());
  for (size_t i = 0; i < xs.size(); ++i) {
    g->SetPoint(i, xs[i], As[i]);
    g->SetPointError(i, 0.0, Aerrs[i]);
  }
  g->SetTitle(title);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlack);
  g->SetLineColor(kBlue+2);
  g->SetLineWidth(2);
  return g;
}

// ------------------------------------------------------
// Helper: build CSV filename safely (ROOT's Form -> std::string)
// ------------------------------------------------------
static std::string make_csv_filename(const char* kin, const char* name) {
  return std::string(Form("plots/%s_%s.csv", name, kin));
}

// ------------------------------------------------------
// Produce single-variable cut stability plots for A
// ------------------------------------------------------
void MakeCutStabilityPlots(
  const Store& ev,
  const char* kin,
  // Baseline cuts (pick your working point or "best" from prior study)
  double dyL0, double dyH0,
  double W2L0, double W2H0,
  double eL0,
  double tL0, double tH0,
  double dxL0, double dxH0,
  // Scan vectors (reused from your definitions)
  const std::vector<double>& dyL_vals,
  const std::vector<double>& dyH_vals,
  const std::vector<double>& W2L_vals,
  const std::vector<double>& W2H_vals,
  const std::vector<double>& eL_vals,
  const std::vector<double>& tL_vals,
  const std::vector<double>& tH_vals,
  const std::vector<double>& dxL_vals,
  const std::vector<double>& dxH_vals)
{
  gStyle->SetOptFit(0);

  std::vector<double> xs, As, Aerrs;

  // Canvas (4x2 pads)
  TCanvas* c = new TCanvas(Form("c_asym_stab_%s", kin),
                           "Asymmetry cut stability (single-variable scans)",
                           4500, 3000);
  c->Divide(5,2);

  // Ensure output folder exists
  gSystem->mkdir("plots", kTRUE);

  // 1) dyL scan
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_dyL = make_csv_filename(kin, "asym_scan_dyL");
  std::ofstream csv_dyL; csv_dyL.open(fn_dyL.c_str());
  csv_dyL << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : dyL_vals) {
    auto r = compute_asymmetry_single(ev, v, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_dyL << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_dyL.close();
  c->cd(1);
  auto g_dyL = make_graph("A vs dyL", xs, As, Aerrs);
  g_dyL->GetXaxis()->SetTitle("dy_{L}");
  g_dyL->GetYaxis()->SetTitle("Asymmetry A");
  g_dyL->Draw("AP");
  TF1 fconst_dyL("fconst_dyL","pol0", *std::min_element(xs.begin(),xs.end()),
                                   *std::max_element(xs.begin(),xs.end()));
  g_dyL->Fit(&fconst_dyL, "Q");

  // 2) dyH scan
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_dyH = make_csv_filename(kin, "asym_scan_dyH");
  std::ofstream csv_dyH; csv_dyH.open(fn_dyH.c_str());
  csv_dyH << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : dyH_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, v, W2L0, W2H0, eL0, tL0, tH0, dxL0, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_dyH << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_dyH.close();
  c->cd(2);
  auto g_dyH = make_graph("A vs dyH", xs, As, Aerrs);
  g_dyH->GetXaxis()->SetTitle("dy_{H}");
  g_dyH->GetYaxis()->SetTitle("Asymmetry A");
  g_dyH->Draw("AP");
  TF1 fconst_dyH("fconst_dyH","pol0", *std::min_element(xs.begin(),xs.end()),
                                   *std::max_element(xs.begin(),xs.end()));
  g_dyH->Fit(&fconst_dyH, "Q");

  // 3) W2L scan
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_W2L = make_csv_filename(kin, "asym_scan_W2L");
  std::ofstream csv_W2L; csv_W2L.open(fn_W2L.c_str());
  csv_W2L << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : W2L_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, dyH0, v, W2H0, eL0, tL0, tH0, dxL0, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_W2L << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_W2L.close();
  c->cd(3);
  auto g_W2L = make_graph("A vs W2L", xs, As, Aerrs);
  g_W2L->GetXaxis()->SetTitle("W^{2}_{L}");
  g_W2L->GetYaxis()->SetTitle("Asymmetry A");
  g_W2L->Draw("AP");
  TF1 fconst_W2L("fconst_W2L","pol0", *std::min_element(xs.begin(),xs.end()),
                                     *std::max_element(xs.begin(),xs.end()));
  g_W2L->Fit(&fconst_W2L, "Q");

  // 4) W2H scan
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_W2H = make_csv_filename(kin, "asym_scan_W2H");
  std::ofstream csv_W2H; csv_W2H.open(fn_W2H.c_str());
  csv_W2H << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : W2H_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, dyH0, W2L0, v, eL0, tL0, tH0, dxL0, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_W2H << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_W2H.close();
  c->cd(4);
  auto g_W2H = make_graph("A vs W2H", xs, As, Aerrs);
  g_W2H->GetXaxis()->SetTitle("W^{2}_{H}");
  g_W2H->GetYaxis()->SetTitle("Asymmetry A");
  g_W2H->Draw("AP");
  TF1 fconst_W2H("fconst_W2H","pol0", *std::min_element(xs.begin(),xs.end()),
                                     *std::max_element(xs.begin(),xs.end()));
  g_W2H->Fit(&fconst_W2H, "Q");

  // 5) eHCAL low threshold (eL) scan
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_eL = make_csv_filename(kin, "asym_scan_eHCAL_L");
  std::ofstream csv_eL; csv_eL.open(fn_eL.c_str());
  csv_eL << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : eL_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, dyH0, W2L0, W2H0, v, tL0, tH0, dxL0, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_eL << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_eL.close();
  c->cd(5);
  auto g_eL = make_graph("A vs eHCAL_{L}", xs, As, Aerrs);
  g_eL->GetXaxis()->SetTitle("eHCAL low cut");
  g_eL->GetYaxis()->SetTitle("Asymmetry A");
  g_eL->Draw("AP");
  TF1 fconst_eL("fconst_eL","pol0", *std::min_element(xs.begin(),xs.end()),
                                   *std::max_element(xs.begin(),xs.end()));
  g_eL->Fit(&fconst_eL, "Q");

  // 6) tL scan (coincidence time lower edge)
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_tL = make_csv_filename(kin, "asym_scan_tL");
  std::ofstream csv_tL; csv_tL.open(fn_tL.c_str());
  csv_tL << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : tL_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, dyH0, W2L0, W2H0, eL0, v, tH0, dxL0, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_tL << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_tL.close();
  c->cd(6);
  auto g_tL = make_graph("A vs coin_time_{L}", xs, As, Aerrs);
  g_tL->GetXaxis()->SetTitle("coincidence time lower edge");
  g_tL->GetYaxis()->SetTitle("Asymmetry A");
  g_tL->Draw("AP");
  TF1 fconst_tL("fconst_tL","pol0", *std::min_element(xs.begin(),xs.end()),
                                   *std::max_element(xs.begin(),xs.end()));
  g_tL->Fit(&fconst_tL, "Q");

  // 7) tH scan (coincidence time upper edge)
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_tH = make_csv_filename(kin, "asym_scan_tH");
  std::ofstream csv_tH; csv_tH.open(fn_tH.c_str());
  csv_tH << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : tH_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, dyH0, W2L0, W2H0, eL0, tL0, v, dxL0, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_tH << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_tH.close();
  c->cd(7);
  auto g_tH = make_graph("A vs coin_time_{H}", xs, As, Aerrs);
  g_tH->GetXaxis()->SetTitle("coincidence time upper edge");
  g_tH->GetYaxis()->SetTitle("Asymmetry A");
  g_tH->Draw("AP");
  TF1 fconst_tH("fconst_tH","pol0", *std::min_element(xs.begin(),xs.end()),
                                   *std::max_element(xs.begin(),xs.end()));
  g_tH->Fit(&fconst_tH, "Q");

  // 8) dxL scan (lower edge of dx window)
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_dxL = make_csv_filename(kin, "asym_scan_dxL");
  std::ofstream csv_dxL; csv_dxL.open(fn_dxL.c_str());
  csv_dxL << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : dxL_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, v, dxH0);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_dxL << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_dxL.close();
  c->cd(8);
  auto g_dxL = make_graph("A vs dxL", xs, As, Aerrs);
  g_dxL->GetXaxis()->SetTitle("dx_{L}");
  g_dxL->GetYaxis()->SetTitle("Asymmetry A");
  g_dxL->Draw("AP");
  TF1 fconst_dxL("fconst_dxL","pol0", *std::min_element(xs.begin(),xs.end()),
                                    *std::max_element(xs.begin(),xs.end()));
  g_dxL->Fit(&fconst_dxL, "Q");

  // 9) dxH scan (upper edge of dx window)
  xs.clear(); As.clear(); Aerrs.clear();
  std::string fn_dxH = make_csv_filename(kin, "asym_scan_dxH");
  std::ofstream csv_dxH; csv_dxH.open(fn_dxH.c_str());
  csv_dxH << "cut_value,A,Aerr,Nplus,Nminus\n";
  for (double v : dxH_vals) {
    auto r = compute_asymmetry_single(ev, dyL0, dyH0, W2L0, W2H0, eL0, tL0, tH0, dxL0, v);
    xs.push_back(v); As.push_back(r.A); Aerrs.push_back(r.Aerr);
    csv_dxH << v << "," << r.A << "," << r.Aerr << "," << r.Nplus << "," << r.Nminus << "\n";
  }
  csv_dxH.close();
  c->cd(9);
  auto g_dxH = make_graph("A vs dxH", xs, As, Aerrs);
  g_dxH->GetXaxis()->SetTitle("dx_{H}");
  g_dxH->GetYaxis()->SetTitle("Asymmetry A");
  g_dxH->Draw("AP");
  TF1 fconst_dxH("fconst_dxH","pol0", *std::min_element(xs.begin(),xs.end()),
                                    *std::max_element(xs.begin(),xs.end()));
  g_dxH->Fit(&fconst_dxH, "Q");

  // Save
  c->SaveAs(Form("plots/asym_stability_%s.png", kin));
  c->SaveAs(Form("plots/asym_stability_%s.jpg", kin));
}

// ------------------------------------------------------
// Main entry: loads data/sim/bkg and runs stability scans
// ------------------------------------------------------
void CutStability(const char* dataFile,
                  const char* simFile,
                  const char* bkgSimFile,
                  const char* kin,
                  int    nbinsDx   = 100,
                  double wN_user   = -1.0,
                  double wP_user   = -1.0)
{
  gStyle->SetOptFit(0);

  Store    ev;
  StoreBkg ev_bkg;

  // ---------------------------
  // 1) Load DATA
  // ---------------------------
  {
    TFile fData(dataFile, "READ");
    TTree* tData = (TTree*)fData.Get("Tout");
    if (!tData) {
      std::cerr << "ERROR: cannot find Tout in data file " << dataFile << std::endl;
      return;
    }

    double vz, ePS, eSH, trP, dx, dy, W2, ct, eHCAL;
    int helicity, IHWP;

    tData->SetBranchAddress("vz",      &vz);
    tData->SetBranchAddress("ePS",     &ePS);
    tData->SetBranchAddress("eSH",     &eSH);
    tData->SetBranchAddress("trP",     &trP);
    tData->SetBranchAddress("dx",      &dx);
    tData->SetBranchAddress("dy",      &dy);
    tData->SetBranchAddress("W2",      &W2);
    tData->SetBranchAddress("coin_time",&ct);
    tData->SetBranchAddress("eHCAL",   &eHCAL);
    tData->SetBranchAddress("helicity",&helicity);
    tData->SetBranchAddress("IHWP",    &IHWP);

    Long64_t nData = tData->GetEntries();

    ev.vz.reserve(nData);
    ev.ePS.reserve(nData);
    ev.dx.reserve(nData);
    ev.dy.reserve(nData);
    ev.W2.reserve(nData);
    ev.ct.reserve(nData);
    ev.eHCAL.reserve(nData);
    ev.helicity.reserve(nData);
    ev.IHWP.reserve(nData);
    ev.w.reserve(nData);
    ev.origin.reserve(nData);

    double good_total_events = 0.0;

    for (Long64_t i = 0; i < nData*0.5; ++i) {
      tData->GetEntry(i);

      // Preserved preselection structure (rejects only if ALL sub-conditions are true)
      if (std::abs(vz) > 0.27 && ePS < 0.2 &&
          std::abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL < 0.025 &&
          std::abs(helicity) != 1)
        continue;

      ++good_total_events;
      ev.vz.push_back(vz);
      ev.ePS.push_back(ePS);
      ev.dx.push_back(dx);
      ev.dy.push_back(dy);
      ev.W2.push_back(W2);
      ev.ct.push_back(ct);
      ev.eHCAL.push_back(eHCAL);
      ev.helicity.push_back(helicity);
      ev.IHWP.push_back(IHWP);
      ev.w.push_back(1.0);
      ev.origin.push_back(2);

      if (i % 500000 == 0) {
        std::cout << std::fixed << std::setprecision(1)
                  << (i * 100.0 / nData) << "% data file\r";
        std::cout.flush();
      }
    }
    std::cout << "\nLoaded " << ev.dx.size() << " total events (data only)\n";
    std::cout << "Good total events in data: " << good_total_events << std::endl;
  }

  // -----------------------------------------
  // 2) (Optional) Load SIM and BACKGROUND SIM
  // -----------------------------------------
  // Not used for asymmetry plots; kept for continuity
  {
    TFile fSim(simFile, "READ");
    TTree* tSim = (TTree*)fSim.Get("Tout");
    if (tSim) {
      double vz, ePS, eSH, trP, dx, dy, W2, ctSim, eHCAL, wt, fnucl;
      tSim->SetBranchAddress("vz",     &vz);
      tSim->SetBranchAddress("ePS",    &ePS);
      tSim->SetBranchAddress("eSH",    &eSH);
      tSim->SetBranchAddress("trP",    &trP);
      tSim->SetBranchAddress("dx",     &dx);
      tSim->SetBranchAddress("dy",     &dy);
      tSim->SetBranchAddress("W2",     &W2);
      if (tSim->GetBranch("coin_time")) tSim->SetBranchAddress("coin_time", &ctSim);
      tSim->SetBranchAddress("eHCAL",  &eHCAL);
      tSim->SetBranchAddress("weight", &wt);
      tSim->SetBranchAddress("fnucl",  &fnucl);

      Long64_t nSim = tSim->GetEntries();
      for (Long64_t i = 0; i < nSim; ++i) {
        tSim->GetEntry(i);
        if (fnucl != 0.0 && fnucl != 1.0) continue;
        if (std::abs(vz) > 0.27 && ePS < 0.2 &&
            std::abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL < 0.025)
          continue;
        ev.vz.push_back(vz);
        ev.ePS.push_back(ePS);
        ev.dx.push_back(dx);
        ev.dy.push_back(dy);
        ev.W2.push_back(W2);
        ev.ct.push_back(tSim->GetBranch("coin_time") ? ctSim : 0.0);
        ev.eHCAL.push_back(eHCAL);
        ev.helicity.push_back(0);
        ev.IHWP.push_back(0);
        ev.w.push_back(wt);
        ev.origin.push_back(fnucl==0.0 ? 0 : 1);
      }
      std::cout << "Loaded " << nSim << " sim entries (added; not used for A).\n";
    } else {
      std::cout << "SIM file " << simFile << " has no Tout (skipping).\n";
    }
  }

  {
    TFile fBkg(bkgSimFile, "READ");
    TTree* tBkg = (TTree*)fBkg.Get("Tout");
    if (tBkg) {
      double vz, ePS, eSH, trP, dx, dy, W2, ctSim, eHCAL, wt;
      tBkg->SetBranchAddress("vz",     &vz);
      tBkg->SetBranchAddress("ePS",    &ePS);
      tBkg->SetBranchAddress("eSH",    &eSH);
      tBkg->SetBranchAddress("trP",    &trP);
      tBkg->SetBranchAddress("dx",     &dx);
      tBkg->SetBranchAddress("dy",     &dy);
      tBkg->SetBranchAddress("W2",     &W2);
      if (tBkg->GetBranch("coin_time")) tBkg->SetBranchAddress("coin_time", &ctSim);
      tBkg->SetBranchAddress("eHCAL",  &eHCAL);
      tBkg->SetBranchAddress("weight", &wt);

      Long64_t nBkg = tBkg->GetEntries();
      ev_bkg.dx.reserve(nBkg);
      for (Long64_t i = 0; i < nBkg*0.2; ++i) {
        tBkg->GetEntry(i);
        if (std::abs(vz) > 0.27 && ePS < 0.2 &&
            std::abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL < 0.025)
          continue;
        ev_bkg.vz.push_back(vz);
        ev_bkg.ePS.push_back(ePS);
        ev_bkg.dx.push_back(dx);
        ev_bkg.dy.push_back(dy);
        ev_bkg.W2.push_back(W2);
        ev_bkg.ct.push_back(tBkg->GetBranch("coin_time") ? ctSim : 0.0);
        ev_bkg.eHCAL.push_back(eHCAL);
        ev_bkg.w.push_back(wt);
      }
      std::cout << "Loaded " << ev_bkg.dx.size() << " background-sim events (not used for A).\n";
    } else {
      std::cout << "BKG SIM file " << bkgSimFile << " has no Tout (skipping).\n";
    }
  }

  // ----------------------------------------------
  // 3) Define scan ranges
  // ----------------------------------------------
  const int nSteps = 10;

  std::vector<double> dyL_vals, dyH_vals,
                      W2L_vals, W2H_vals, eL_vals,
                      tL_vals, tH_vals,
                      dxL_vals, dxH_vals;


  //these cut values are hardcoded for now, but can be modified to user input later
  for (int i = 0; i < nSteps+1; ++i) dyL_vals.push_back(-0.8 + i*(0.8/nSteps));
  for (int i = 0; i < nSteps+1; ++i) dyH_vals.push_back( 0.8 - i*(0.8/nSteps));
  for (int i = 0; i < nSteps+1;   ++i) W2L_vals.push_back(-2.0 + i*(2.0/nSteps));
  for (int i = 0; i < nSteps+1;   ++i) W2H_vals.push_back( 1.0 + i*(1.0/nSteps));
  for (int i = 0; i < nSteps+1;   ++i) eL_vals.push_back( 0.025 + i*(0.5/nSteps));
  for (int i = 1; i < nSteps+1;        ++i) tL_vals.push_back(120.0 - i*(6.0/nSteps));
  for (int i = 1; i < nSteps+1;        ++i) tH_vals.push_back(120.0 + i*(6.0/nSteps));
  // dx scans: vary lower edge from -0.8 to 0.0, upper edge from 0.0 to 0.8
  for (int i = 0; i < nSteps+1; ++i) dxL_vals.push_back( -0.8 + i*(0.8/(nSteps)));
  for (int i = 0; i < nSteps+1; ++i) dxH_vals.push_back( 0.8 - i*(0.8/(nSteps)));

  // ----------------------------------------------
  // 4) Choose baseline cuts (working point)
  // ----------------------------------------------
  double dyL0 = -0.8;
  double dyH0 =  0.8;
  double W2L0 = -2.0;
  double W2H0 =  2.0;
  double eL0  =  0.025;
  double tL0  = 115.0;
  double tH0  = 125.0;
  double dxL0 = -0.8;
  double dxH0 =  0.8;

  std::cout << "Baseline cuts:\n"
            << " dy in [" << dyL0 << ", " << dyH0 << "]\n"
            << " W2 in [" << W2L0 << ", " << W2H0 << "]\n"
            << " eHCAL_L = " << eL0 << "\n"
            << " coin_time in [" << tL0 << ", " << tH0 << "]\n";

  // ----------------------------------------------
  // 5) Run single-variable scans and make plots
  // ----------------------------------------------
  MakeCutStabilityPlots(ev, kin,
                        dyL0, dyH0,
                        W2L0, W2H0,
                        eL0,
                        tL0, tH0,
                        dxL0, dxH0,
                        dyL_vals, dyH_vals,
                        W2L_vals, W2H_vals,
                        eL_vals, tL_vals, tH_vals,
                        dxL_vals, dxH_vals);

  std::cout << "Saved plots/asym_stability_" << kin << ".png/.jpg and per-scan CSVs in plots/.\n";
}

