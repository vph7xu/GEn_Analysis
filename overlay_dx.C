// overlay_dx.C
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
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

  // Otherwise, try to draw from a TTree branch/expression
  auto* t = dynamic_cast<TTree*>(f->Get(treename));
  if (!t) { ::Error("makeHist","Tree '%s' not found", treename); return nullptr; }

  gROOT->cd(); // keep histogram in memory, not inside the file
  auto* h = new TH1D(hname, Form("%s;%s;Events", var, var), nbins, xmin, xmax);
  h->Sumw2();
  t->Draw(Form("%s>>%s", var, hname), cut, "goff");
  return h;
}

void overlay_dx(const char* file1,
                const char* file2,
                const char* treename="Tout",
                const char* var="dx",
                const char* cut="1",
                int nbins=120, double xmin=-4, double xmax=3,
                bool normalize=true,
                const char* out="dx_overlay.pdf")
{
  gStyle->SetOptStat(0);

  std::unique_ptr<TFile> f1(TFile::Open(file1,"READ"));
  std::unique_ptr<TFile> f2(TFile::Open(file2,"READ"));
  if (!f1 || f1->IsZombie() || !f2 || f2->IsZombie()) {
    ::Error("overlay_dx","Could not open one of the files.");
    return;
  }

  const char *cut1 = "ePS>0.2&&abs(vz)<0.27&&abs((ePS+eSH)/trP-1)<0.2&&eHCAL>0.225&&W2<1.5&&abs(dy)<0.4&&abs(coin_time-121)<3";
  const char *cut2 = "ePS>0.2&&abs(vz)<0.27&&abs((ePS+eSH)/trP-1)<0.2&&eHCAL>0.225&&W2<1.5&&abs(dy)<0.4&&abs(coin_time-184)<2";  

  auto* h1 = makeHist(f1.get(), treename, var, cut1, nbins, xmin, xmax, "h1");
  auto* h2 = makeHist(f2.get(), treename, var, cut2, nbins, xmin, xmax, "h2");
  if (!h1 || !h2) { ::Error("overlay_dx","Failed to build histograms."); return; }

  if (normalize) {
    double i1 = h1->Integral();
    double i2 = h2->Integral();
    if (i1>0) h1->Scale(1.0/i1);
    if (i2>0) h2->Scale(1.0/i2);
    h1->GetYaxis()->SetTitle("Normalized events");
  }

  h1->SetLineColor(kAzure+1);
  h2->SetLineColor(kOrange+7);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);

  double ymax = std::max(h1->GetMaximum(), h2->GetMaximum());
  h1->SetMaximum(1.15*ymax);

  TCanvas c("c","dx overlay",900,700);
  h1->Draw("HIST");
  h2->Draw("HIST SAME");

  TLegend leg(0.60,0.75,0.88,0.88);
  leg.SetBorderSize(0);
  leg.SetFillStyle(0);
  leg.AddEntry(h1, gSystem->BaseName(file1), "l");
  leg.AddEntry(h2, gSystem->BaseName(file2), "l");
  leg.Draw();

  if (out && strlen(out)) c.SaveAs(out);
}

