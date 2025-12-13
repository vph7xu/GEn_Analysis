// overlay_W2_3.C
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <memory>
#include <algorithm>

static TH1D* makeHist(TFile* f,
                      const char* treename,
                      const char* var,
                      const char* cut,
                      int nbins, double xmin, double xmax,
                      const char* hname)
{
  if (!f || f->IsZombie()) { ::Error("makeHist","Bad TFile"); return nullptr; }

  // If a histogram named `var` already exists in the file, clone and return it
  if (TObject* obj = f->Get(var)) {
    if (obj->InheritsFrom(TH1::Class())) {
      auto* h = dynamic_cast<TH1D*>(obj);
      if (h) return (TH1D*)h->Clone(hname);
    }
  }

  // Otherwise, draw from a TTree branch/expression
  auto* t = dynamic_cast<TTree*>(f->Get(treename));
  if (!t) { ::Error("makeHist","Tree '%s' not found", treename); return nullptr; }

  gROOT->cd(); // keep histogram in memory, not inside the file
  auto* h = new TH1D(hname, Form("%s;%s;Events", var, var), nbins, xmin, xmax);
  h->Sumw2();
  t->Draw(Form("%s>>%s", var, hname), cut, "goff");
  return h;
}

// Overlay three W2 plots (different files and cuts)
void overlay_W2_3(const char* file1,
                  const char* file2,
                  const char* file3,
                  const char* treename="Tout",
                  const char* var="W2",
                  const char* cut1="1",
                  const char* cut2="1",
                  const char* cut3="1",
                  int nbins=120, double xmin=0.0, double xmax=3.0,
                  bool normalize=true,
                  const char* out="W2_overlay.pdf")
{
  gStyle->SetOptStat(0);

  std::unique_ptr<TFile> f1(TFile::Open(file1,"READ"));
  std::unique_ptr<TFile> f2(TFile::Open(file2,"READ"));
  std::unique_ptr<TFile> f3(TFile::Open(file3,"READ"));
  if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie() || !f3 || f3->IsZombie()) {
    ::Error("overlay_W2_3","Could not open one of the files.");
    return;
  }

  const char *cut1a = "ePS>0.2&&abs(vz)<0.27&&abs((ePS+eSH)/trP-1)<0.2&&eHCAL>0.025&&abs(dy)<0.4&&abs(dx)<0.4&&abs(coin_time-121)<3";
  const char *cut2a = "ePS>0.2&&abs(vz)<0.27&&abs((ePS+eSH)/trP-1)<0.2&&eHCAL>0.125&&abs(dy)<0.4&&abs(dx)<0.4&&abs(coin_time-121)<3";
  const char *cut3a = "ePS>0.2&&abs(vz)<0.27&&abs((ePS+eSH)/trP-1)<0.2&&eHCAL>0.325&&abs(dy)<0.4&&abs(dx)<0.4&&abs(coin_time-184)<2";  

  auto* h1 = makeHist(f1.get(), treename, var, cut1a, nbins, xmin, xmax, "h1");
  auto* h2 = makeHist(f2.get(), treename, var, cut2a, nbins, xmin, xmax, "h2");
  auto* h3 = makeHist(f3.get(), treename, var, cut3a, nbins, xmin, xmax, "h3");
  if (!h1 || !h2 || !h3) { ::Error("overlay_W2_3","Failed to build histograms."); return; }

  if (normalize) {
    double i1=h1->Integral(), i2=h2->Integral(), i3=h3->Integral();
    if (i1>0) h1->Scale(1.0/i1);
    if (i2>0) h2->Scale(1.0/i2);
    if (i3>0) h3->Scale(1.0/i3);
    h1->GetYaxis()->SetTitle("Normalized events");
  }

  // Styling
  h1->SetLineColor(kAzure+1);
  h2->SetLineColor(kOrange+7);
  h3->SetLineColor(kGreen+2);
  h1->SetLineWidth(3); h2->SetLineWidth(3); h3->SetLineWidth(3);

  double ymax = std::max({h1->GetMaximum(), h2->GetMaximum(), h3->GetMaximum()});
  h1->SetMaximum(1.15*ymax);

  TCanvas c("c","W2 overlay (3 files)",900,700);
  h1->Draw("HIST");
  h2->Draw("HIST SAME");
  h3->Draw("HIST SAME");

  TLegend leg(0.58,0.72,0.88,0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(h1, gSystem->BaseName(file1), "l");
  leg.AddEntry(h2, gSystem->BaseName(file2), "l");
  leg.AddEntry(h3, gSystem->BaseName(file3), "l");
  leg.Draw();

  if (out && strlen(out)) c.SaveAs(out);
}

