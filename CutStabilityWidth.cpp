
// CutStabilityWidth.cpp — Streaming cut stability for dx, dy, and coin_time
// -------------------------------------------------------------------------
// - Streams directly from DATA TTree (no large arrays).
// - For each variable (dx, dy, coin_time), midpoint is fixed, width varies.
// - Computes A = (N+ - N-) / (N+ + N-) and σ_A for each width.
// - Adds progress bars and constant fits for plateau checks.
// - Saves plots and CSVs under plots/. Also writes a single consolidated PDF.
//
// Run:
//   root -l
//   root [0] .L CutStabilityWidth.cpp+
//   root [1] CutStabilityWidth("data.root","kinA")
// -------------------------------------------------------------------------

#include <TF1.h>
#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TStopwatch.h>
#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sstream>

// Result container
struct ScanResult {
  std::vector<double> x;      // width values
  std::vector<double> A;      // asymmetry
  std::vector<double> Aerr;   // uncertainty
  std::vector<long long> Np;  // N+
  std::vector<long long> Nm;  // N-
};

// Progress bar helper
inline void progress_bar(Long64_t i, Long64_t N,
                         const char* label = nullptr,
                         int width = 50) {
  if (N == 0) return;
  double frac = double(i + 1) / double(N);
  int filled = int(frac * width);
  std::ostringstream oss;
  if (label) oss << label << " ";
  oss << "[";
  for (int j = 0; j < width; ++j) oss << (j < filled ? '#' : '-');
  oss << "] " << std::fixed << std::setprecision(1) << (frac * 100.0) << "%\r";
  std::cout << oss.str();
  std::cout.flush();
}

// Make graph
static TGraphErrors* make_graph(const char* title,
                                const std::vector<double>& xs,
                                const std::vector<double>& As,
                                const std::vector<double>& Aerrs) {
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

// Write CSV
static void write_csv(const char* kin, const char* base, const ScanResult& r) {
  std::string fname = Form("plots/%s_%s.csv", base, kin);
  std::ofstream csv(fname.c_str());
  csv << "width,A,Aerr,Nplus,Nminus\n";
  for (size_t i=0;i<r.x.size();++i)
    csv << r.x[i] << "," << r.A[i] << "," << r.Aerr[i] << ","
        << r.Np[i] << "," << r.Nm[i] << "\n";
  csv.close();
}

// Streaming width scan with progress bar
ScanResult StreamWidthScan(TTree* tData,
                           const std::vector<double>& widths,
                           const char* varType,
                           double dx_mid, double dy_mid, double ct_mid,
                           // Baseline cuts
                           double W2L0, double W2H0, double eL0,
                           double dyL0, double dyH0, double dxL0, double dxH0,
                           double ctL0, double ctH0,
                           const char* kin) {
  double vz, ePS, eSH, trP, dx, dy, W2, ct, eHCAL;
  int helicity, IHWP;

  tData->SetBranchAddress("vz", &vz);
  tData->SetBranchAddress("ePS", &ePS);
  tData->SetBranchAddress("eSH", &eSH);
  tData->SetBranchAddress("trP", &trP);
  tData->SetBranchAddress("dx", &dx);
  tData->SetBranchAddress("dy", &dy);
  tData->SetBranchAddress("W2", &W2);
  tData->SetBranchAddress("coin_time", &ct);
  tData->SetBranchAddress("eHCAL", &eHCAL);
  tData->SetBranchAddress("helicity", &helicity);
  tData->SetBranchAddress("IHWP", &IHWP);

  Long64_t nData = 1.0 * tData->GetEntries();
  std::vector<long long> Nplus(widths.size(), 0), Nminus(widths.size(), 0);

  int p_kin = 1;

  if (std::strcmp(kin, "GEN2_He3") == 0){
    p_kin = -1;
  }
  else if (std::strcmp(kin, "GEN3_He3") == 0){
    p_kin = 1;
  }
  else if(std::strcmp(kin, "GEN4_He3") == 0){
    p_kin = 1;
  }
  else if(std::strcmp(kin, "GEN4b_He3") == 0){
    p_kin = 1;
  }

  p_kin=-1*p_kin;

  std::cout<<"p_kin : "<<p_kin<<endl;

  Long64_t step = std::max<Long64_t>(1, nData / 200); // update ~0.5%
  TStopwatch sw; sw.Start();

  for (Long64_t i = 0; i < nData; ++i) {
    tData->GetEntry(i);

    if (std::abs(vz) > 0.27 && ePS < 0.2 &&
        std::abs((eSH+ePS)/trP - 1.0) > 0.2 && eHCAL < 0.025 &&
        std::abs(helicity) != 1)
      continue;

    bool baseCommon = std::abs(vz) < 0.27 &&
                      ePS > 0.2 &&
                      W2 > W2L0 && W2 < W2H0 &&
                      eHCAL > eL0;
    if (!baseCommon) continue;

    bool poshel = (p_kin * IHWP * helicity == 1);
    bool neghel = (p_kin * IHWP * helicity == -1);

    for (size_t k = 0; k < widths.size(); ++k) {
      double w = widths[k];
      bool pass = false;
      if (std::string(varType) == "dx") {
        pass = (dx > dx_mid - w && dx < dx_mid + w && dyL0 < dy && dy < dyH0 && ctL0< ct && ct< ctH0);
      } else if (std::string(varType) == "dy") {
        pass = (dy > dy_mid - w && dy < dy_mid + w && dxL0 < dx && dx < dxH0 && ctL0< ct && ct< ctH0);
      } else if (std::string(varType) == "ct") {
        pass = (ct > ct_mid - w && ct < ct_mid + w && dyL0 < dy && dy < dyH0 && dxL0< dx && dx< dxH0);
      }
      if (pass) {
        if (poshel) ++Nplus[k];
        else if (neghel) ++Nminus[k];
      }
    }

    if ((i % step) == 0 || i + 1 == nData) {
      progress_bar(i, nData, varType);
    }
  }
  std::cout << "\nDone in " << sw.RealTime() << " s.\n";

  ScanResult out;
  out.x = widths;
  out.Np = Nplus;
  out.Nm = Nminus;
  out.A.resize(widths.size());
  out.Aerr.resize(widths.size());
  for (size_t k = 0; k < widths.size(); ++k) {
    double Np = (double)Nplus[k], Nm = (double)Nminus[k], N = Np + Nm;
    if (N > 0) {
      out.A[k] = (Np - Nm) / N;
      out.Aerr[k] = std::sqrt(std::max(0.0, (1.0 - out.A[k]*out.A[k]) / N));
    } else {
      out.A[k] = 0.0;
      out.Aerr[k] = 1.0;
    }
  }
  return out;
}

// Main entry
void CutStabilityWidth(const char* dataFile, const char* kin) {
  gStyle->SetOptFit(0);
  gSystem->mkdir("plots", kTRUE);

  TFile fData(dataFile, "READ");
  TTree* tData = (TTree*)fData.Get("Tout");
  if (!tData) {
    std::cerr << "ERROR: cannot find Tout in " << dataFile << std::endl;
    return;
  }

  // Baseline midpoints
  double dx_mid = 0.0;
  double dy_mid = 0.0;
  double ct_mid = 120.0; // ns
  double W2L0 = -2, W2H0 = 1.6, eL0 = 0.025;
  double dyL0 = -0.8, dyH0=0.8;
  double dxL0 = -0.8, dxH0=0.8;
  double ctL0 = 126, ctH0=134;

  if(std::strcmp(kin, "GEN2_He3") == 0){
    ct_mid = 130.0; W2L0 = -2; W2H0 = 1.3; eL0 = 0.025;
    dyL0 = -0.8; dyH0 = 0.8; dxL0 = -0.8; dxH0 = 0.8; ctL0 = 126; ctH0 = 134;
  }
  else if(std::strcmp(kin, "GEN3_He3") == 0){
    ct_mid = 120.0; W2L0 = -2; W2H0 = 1.5; eL0 = 0.085;
    dyL0 = -0.5; dyH0 = 0.5; dxL0 = -0.5; dxH0 = 0.5; ctL0 = 117; ctH0 = 123;
  }
  else if(std::strcmp(kin, "GEN4_He3") == 0){
    ct_mid = 121.0; W2L0 = -2; W2H0 = 1.4; eL0 = 0.225;
    dyL0 = -0.4; dyH0 = 0.4; dxL0 = -0.4; dxH0 = 0.4; ctL0 = 118; ctH0 = 124;
  }
  else if(std::strcmp(kin, "GEN4b_He3") == 0){
    ct_mid = 185.0; W2L0 = -2; W2H0 = 1.4; eL0 = 0.325;
    dyL0 = -0.4; dyH0 = 0.4; dxL0 = -0.4; dxH0 = 0.4; ctL0 = 181; ctH0 = 188;
  }
  else{
    ct_mid = 120.0; W2L0 = -2; W2H0 = 1.6; eL0 = 0.025;
    dyL0 = -0.8; dyH0 = 0.8; dxL0 = -0.8; dxH0 = 0.8; ctL0 = 126; ctH0 = 134;
  }

  // Baseline cuts for other vars
  //double W2L0 = -2, W2H0 = 1.6, eL0 = 0.025;

  // Width scan ranges
  std::vector<double> dx_widths = {0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.5,0.55,0.60,0.65,0.70,0.75,0.80};
  std::vector<double> dy_widths = {0.05,0.10,0.15,0.20,0.25,0.30,0.35,0.40,0.45,0.5,0.55,0.60,0.65,0.70,0.75,0.80};
  std::vector<double> ct_widths = {1.0,2.0,3.0,4.0,5.0,6.0,7.0};

  std::cout << "Running width-based stability scans...\n";

  auto r_dx = StreamWidthScan(tData, dx_widths, "dx", dx_mid, dy_mid, ct_mid, W2L0, W2H0, eL0, dyL0, dyH0, dxL0, dxH0, ctL0, ctH0, kin);
  auto r_dy = StreamWidthScan(tData, dy_widths, "dy", dx_mid, dy_mid, ct_mid, W2L0, W2H0, eL0, dyL0, dyH0, dxL0, dxH0, ctL0, ctH0, kin);
  auto r_ct = StreamWidthScan(tData, ct_widths, "ct", dx_mid, dy_mid, ct_mid, W2L0, W2H0, eL0, dyL0, dyH0, dxL0, dxH0, ctL0, ctH0, kin);

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
    //g->GetYaxis()->SetNdivisions(505);
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


  // Plot
  TCanvas* c = new TCanvas(Form("c_width_stab_%s", kin),
                           "Width-based cut stability", 2400, 800);
  c->Divide(3,1);
  c->SetLeftMargin(0.14);
  c->SetRightMargin(0.05);
  c->SetTopMargin(0.08);
  c->SetBottomMargin(0.12);

  c->cd(1);
  gPad->SetLeftMargin(0.22);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  auto g_dx = make_graph("", r_dx.x, r_dx.A, r_dx.Aerr);
  g_dx->GetXaxis()->SetTitle("dx half-width (m)");
  g_dx->GetYaxis()->SetTitle("");
  style_graph(g_dx);
  g_dx->Draw("AP");
  add_y_label("Asymmetry A (%)");
  
  // TF1 fconst_dx("fconst_dx","pol0",
  //               *std::min_element(r_dx.x.begin(),r_dx.x.end()),
  //               *std::max_element(r_dx.x.begin(),r_dx.x.end()));
  // g_dx->Fit(&fconst_dx, "Q");

  c->cd(2);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.03);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  auto g_dy = make_graph("", r_dy.x, r_dy.A, r_dy.Aerr);
  g_dy->GetXaxis()->SetTitle("dy half-width (m)");
  g_dy->GetYaxis()->SetTitle("");
  style_graph(g_dy);
  g_dy->Draw("AP");
  add_y_label("Asymmetry A (%)");
  
  // TF1 fconst_dy("fconst_dy","pol0",
  //               *std::min_element(r_dy.x.begin(),r_dy.x.end()),
  //               *std::max_element(r_dy.x.begin(),r_dy.x.end()));
  // g_dy->Fit(&fconst_dy, "Q");

  c->cd(3);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.06);
  gPad->SetTopMargin(0.08);
  gPad->SetBottomMargin(0.12);
  auto g_ct = make_graph("", r_ct.x, r_ct.A, r_ct.Aerr);
  g_ct->GetXaxis()->SetTitle("coincidence time half-width (ns)");
  g_ct->GetYaxis()->SetTitle("");
  style_graph(g_ct);
  g_ct->Draw("AP");
  add_y_label("Asymmetry A (%)");

  // TF1 fconst_ct("fconst_ct","pol0",
  //               *std::min_element(r_ct.x.begin(),r_ct.x.end()),
  //               *std::max_element(r_ct.x.begin(),r_ct.x.end()));
  // g_ct->Fit(&fconst_ct, "Q");

  // Save raster images (unchanged)
  c->SaveAs(Form("plots/width_stability_%s.png", kin));
  c->SaveAs(Form("plots/width_stability_%s.jpg", kin));

  // NEW: Save to a single consolidated PDF
  // If you only have one canvas, this will produce a single-page PDF.
  // If later you add more canvases, use the "open / append / close" syntax shown below.
  const char* outpdf = Form("plots/width_stability_%s_all.pdf", kin);

  // Single page (current canvas):
  c->Print(outpdf);  // a single-page PDF

  // --- If you add more canvases, do it like this ---
  // c->Print(Form("plots/width_stability_%s_all.pdf(", kin));  // OPEN multi-page
  // c->Print(Form("plots/width_stability_%s_all.pdf",  kin));  // APPEND middle pages
  // otherCanvas->Print(Form("plots/width_stability_%s_all.pdf",  kin));  // more pages
  // lastCanvas->Print(Form("plots/width_stability_%s_all.pdf)", kin));  // CLOSE multi-page
  // --------------------------------------------------

  // Write CSVs
  write_csv(kin, "width_scan_dx", r_dx);
  write_csv(kin, "width_scan_dy", r_dy);
  write_csv(kin, "width_scan_ct", r_ct);

  std::cout << "Saved plots and CSVs in plots/ for " << kin << std::endl;
}
