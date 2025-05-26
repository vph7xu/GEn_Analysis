#ifndef MODELS_H
#define MODELS_H

#include <TH1D.h>    // bring in the definition of TH1D
#include <TF1.h>     // if you ever need TF1* in this header
#include <Rtypes.h>  // for Double_t, etc.  (often pulled in transitively)


TH1D *hist_p = nullptr;
TH1D *hist_n = nullptr;
TH1D *hist_bkg = nullptr;

TH1D *hist_sim_bkg_p = nullptr;
TH1D *hist_sim_bkg_n = nullptr;
TH1D *hist_sim_bkg = nullptr;

Double_t fit_sim_n_bkg(Double_t *x, Double_t *par){
	
	int bin1 = hist_p->FindBin(x[0]);
	int bin2 = hist_n->FindBin(x[0]);
	int bin3 = hist_bkg->FindBin(x[0]);

	double value =	par[0] * (hist_p->GetBinContent(bin1) +
			par[1] * hist_n->GetBinContent(bin2) +
			par[2] * hist_bkg->GetBinContent(bin3));

	return value;
}

Double_t fit_sim_n_bkg_sim(Double_t *x, Double_t *par){
	
	int bin1 = hist_p->FindBin(x[0]);
	int bin2 = hist_n->FindBin(x[0]);
	int bin3 = hist_sim_bkg->FindBin(x[0]);

	double value =	par[0] * (hist_p->GetBinContent(bin1) +
			par[1] * hist_n->GetBinContent(bin2) +
			par[2] * hist_sim_bkg->GetBinContent(bin3));

	return value;
}

#endif // MODELS_H
