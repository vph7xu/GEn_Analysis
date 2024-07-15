TH1D *hist_p = nullptr;
TH1D *hist_n = nullptr;
TH1D *hist_bkg = nullptr;

Double_t fit_sim_n_bkg(Double_t *x, Double_t *par){
	
	int bin1 = hist_p->FindBin(x[0]);
	int bin2 = hist_n->FindBin(x[0]);
	int bin3 = hist_bkg->FindBin(x[0]);

	double value =	par[0] * hist_p->GetBinContent(bin1) +
			par[1] * hist_n->GetBinContent(bin2) +
			par[2] * hist_bkg->GetBinContent(bin3);

	return value;
}
