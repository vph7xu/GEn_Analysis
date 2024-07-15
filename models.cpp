//this is to return models from the simulation and for the bkg

#include "models.h"
#include "parse.h"

TH1D* gExistingHist = nullptr;

Double_t customFunction(Double_t *x, Double_t *par) {
    if (gExistingHist == nullptr) return 0.0;
    Int_t bin = gExistingHist->FindBin(x[0]);
    return gExistingHist->GetBinContent(bin);
}

TH1D* bkg_model(TH2D *hist2D, double xmin, double xmax, const char* kin){

	//std::map<std::string, std::string> config = parseConfig(Form("cut_%s.txt",kin));

	TH2D *h2_anticut = (TH2D*)hist2D->Clone("h2_anticut");
   	h2_anticut->Reset();
	
	int nx = hist2D->GetNbinsX();
	int ny = hist2D->GetNbinsY();
    	
	for (int ix = 1; ix <= nx; ++ix) {
        	for (int iy = 1; iy <= ny; ++iy) {
            		double x = hist2D->GetXaxis()->GetBinCenter(ix);
            		double y = hist2D->GetYaxis()->GetBinCenter(iy);

            		if (!(x > xmin && x < xmax)) {
                		double content = hist2D->GetBinContent(ix, iy);
              	  		h2_anticut->SetBinContent(ix, iy, content);
         	   	}
       		}
	}

	TH1D *projY = h2_anticut->ProjectionX();
        TH1D *projX = h2_anticut->ProjectionY();

	TCanvas *c = new TCanvas("c","c",1600,1200);
	c->Divide(2,2);
	c->cd(1);
	projY->SetTitle("dy projection");
	projY->Draw("HIST");
	c->cd(2);
	projX->SetTitle("dx projection");
	projX->Draw("HIST");
	c->cd(3);
	h2_anticut->Draw();
	c->cd(4);
	hist2D->GetXaxis()->SetRangeUser(-4,4);
        hist2D->GetYaxis()->SetRangeUser(-4,4);
	hist2D->Draw();

	c->Print(Form("../plots/bkg_%s.pdf",kin));

	gExistingHist = projY;

	TF1 *model = new TF1("model", customFunction,projY->GetXaxis()->GetXmin(),projY->GetXaxis()->GetXmax(),0);

	return projX;

}

std::pair<TH1D*, TH1D*> sim_hist(const char* sim_filename, double W2_L_sim, double W2_H_sim){
	TFile* sim_file = TFile::Open(sim_filename);

	TTree* sim_tree = (TTree*)sim_file->Get("Tout");

	double dx = 0.0;
	double dy = 0.0;
	double W2 = 0.0;
	double Q2 = 0.0;
	double weight = 0.0;
	double fnucl = 0.0;

	sim_tree->SetBranchAddress("dx",&dx);
	sim_tree->SetBranchAddress("dy",&dy);
	sim_tree->SetBranchAddress("W2",&W2);
	sim_tree->SetBranchAddress("Q2",&Q2);
	sim_tree->SetBranchAddress("weight",&weight);
	sim_tree->SetBranchAddress("fnucl",&fnucl);

        TH1D *h_dx_n = new TH1D("h_dx_n","dx neutrons",200,-10,10);
        TH1D *h_dx_p = new TH1D("h_dx_p","dx_protons",200,-10,10);

	int nentries = sim_tree->GetEntries();

	for (int i = 0; i<nentries; i++){
                sim_tree->GetEntry(i);

		if (W2_L_sim<W2 and W2<W2_H_sim){
                        if (fnucl == 0.0){
                                h_dx_n->Fill(dx,weight);
                        }
                        else if (fnucl == 1.0){
                                h_dx_p->Fill(dx,weight);
                        }
                }

                if (i %1000 == 0 ) std::cout<< (i * 100.0/ nentries)<<"sim file" << "% \r";
                std::cout.flush();
	}

	//TCanvas* ctest = new TCanvas("ctest","ctest",800,1200);
	//h_dx_n->Draw();

	return std::make_pair(h_dx_n,h_dx_p);
}


void models(const char* filename,const char* sim_filename,const char* printfilename, const char* kin){
        
	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
	
	double coin_time_L = getDoubleValue(config,"coin_time_L");
	double coin_time_H = getDoubleValue(config,"coin_time_H");

	double W2_L = getDoubleValue(config,"W2_L");
	double W2_H = getDoubleValue(config,"W2_H");

	double dy_L = getDoubleValue(config,"dy_L");
	double dy_H = getDoubleValue(config,"dy_H");

	//end of parsing cuts
	
	TFile* file = TFile::Open(filename);
	//TFile* sim_file = TFile::Open(sim_filename);

        TTree* tree = (TTree*)file->Get("Tout");
	//TTree* sim_tree = (TTree*)sim_file->Get("Tout");

	//p and n dist
	auto sim_histograms = sim_hist(sim_filename,W2_L,W2_H);
	hist_p = sim_histograms.second;
	hist_n = sim_histograms.first;

        double dx = 0.0;
        double dy = 0.0;
        double W2 = 0.0;
        double Q2 = 0.0;
	double coin_time = 0.0;

        tree->SetBranchAddress("dx",&dx);
        tree->SetBranchAddress("dy",&dy);
        tree->SetBranchAddress("W2",&W2);
        tree->SetBranchAddress("Q2",&Q2);
	tree->SetBranchAddress("coin_time",&coin_time);

        TH1D *h_dx = new TH1D("h_dx","dx",200,-10,10);
        TH1D *h_dy = new TH1D("h_dy","dy",200,-10,10);
        TH1D *h_W2 = new TH1D("h_W2","W2",250,-4,8);
        TH2D *h_dxdy = new TH2D("h_dxdy","dxdy",250,-2,2,250,-4,4);

        TH1D *h_dx_W2_cut = new TH1D("h_dx_W2_cut","dx",200,-10,10);
        TH1D *h_dy_W2_cut = new TH1D("h_dy_W2_cut","dy",200,-10,10);
        //TH1D *h_W2 = new TH1D("h_W2","W2",1000,-4,8);
        TH2D *h_dxdy_W2_cut = new TH2D("h_dxdy_W2_cut","dx v dy",200,-4,4,200,-10,10);

	TH1D *h_dx_sim_n_bkg = new TH1D("h_dx_sim_n_bkg","dx from simulation",200,-10,10);

        int nentries = tree->GetEntries();

	std::cout<<"coin_L and H : "<<coin_time_L<<"  "<<coin_time_H<<endl;
	std::cout<<"dx_L and H : "<<dy_L <<" "<<dy_H <<endl;

	//get the bkg from data
	for (int i = 0; i<nentries; i++){
                tree->GetEntry(i);

                //before adding a cut on W2
                h_dx->Fill(dx);
                h_dy->Fill(dy);
                h_dxdy->Fill(dy,dx);
                h_W2->Fill(W2);
	        
		//add a cut on W2
                if ((W2_L<W2 and W2<W2_H)and(coin_time_L<coin_time and coin_time<coin_time_H)){
                        h_dx_W2_cut->Fill(dx);
                        h_dy_W2_cut->Fill(dy);
                        h_dxdy_W2_cut->Fill(dy,dx);
                }	

                if (i %1000 == 0 ) std::cout << (i * 100.0/ nentries) << "% \r";
                std::cout.flush();

        }
        //add lines to W2 plot to show the cut
	//
	
	hist_bkg=bkg_model(h_dxdy_W2_cut, dy_L, dy_H, kin);//background dist
	
	//scale everything
	double scale_data = h_dx_W2_cut->Integral();
	h_dx_W2_cut->Scale(1.0/h_dx_W2_cut->Integral());
	hist_p->Scale(1.0/hist_p->Integral());
	hist_n->Scale(1.0/hist_n->Integral());
	hist_bkg->Scale(1.0/hist_bkg->Integral());

	//bkg_model(h_dxdy_W2_cut,-1.5,1.5, kin);

	TCanvas *c1 = new TCanvas("c1","c1",800,1200);
	h_dxdy_W2_cut->Draw();

	TF1 *fit_data = new TF1("fit_data", fit_sim_n_bkg, h_dx_W2_cut->GetXaxis()->GetXmin(), h_dx_W2_cut->GetXaxis()->GetXmax(),3);

	fit_data->SetParameters(1.0,1.0,1.0);//starting parameters

	h_dx_W2_cut->Fit(fit_data,"R");//chi2 minimization happen here

	hist_p->Scale(fit_data->GetParameter(0));
	hist_n->Scale(fit_data->GetParameter(1));
	hist_bkg->Scale(fit_data->GetParameter(2));

	hist_p->SetLineColor(46);
	hist_n->SetLineColor(38);
	hist_bkg->SetLineColor(16);
	h_dx_sim_n_bkg->SetLineColor(30);

	hist_p->Scale(scale_data);
	hist_n->Scale(scale_data);
	hist_bkg->Scale(scale_data);
	h_dx_W2_cut->Scale(scale_data);
	
	h_dx_W2_cut->GetXaxis()->SetRangeUser(-6,6);

	hist_p->SetFillColorAlpha(46,0.5);
	hist_p->SetFillStyle(3004);
	hist_n->SetFillColorAlpha(38,0.5);
	hist_n->SetFillStyle(3005);
	hist_bkg->SetFillColorAlpha(16,0.5);
	hist_bkg->SetFillStyle(3003);
	h_dx_sim_n_bkg->SetFillColorAlpha(30,0.1);
	h_dx_sim_n_bkg->SetFillStyle(3009);

	//fill up the total model histogram
	
	int nBins = h_dx_sim_n_bkg->GetNbinsX();

	for (int i = 1; i<=nBins; ++i){
		h_dx_sim_n_bkg->SetBinContent(i,(hist_p->GetBinContent(i)+hist_n->GetBinContent(i)+hist_bkg->GetBinContent(i)));	
	}

        TLegend *legend = new TLegend(0.6,0.6,0.9,0.9);
        legend->AddEntry(h_dx_W2_cut,"Data","lf");
        legend->AddEntry(hist_p,"sim proton","lf");
        legend->AddEntry(hist_n,"sim neutron","lf");
        legend->AddEntry(hist_bkg,"bkg data","lf");
        legend->AddEntry(h_dx_sim_n_bkg,"full model","lf");
        //legend->Draw();

	TCanvas* cfit = new TCanvas("cfit","Fit with models",800,1200);
	cfit->Divide(1,2);
	cfit->cd(1);
	h_dx_W2_cut->Draw("HIST");
	hist_p->Draw("HIST SAME");
	hist_n->Draw("HIST SAME");
	hist_bkg->Draw("HIST SAME");
	h_dx_sim_n_bkg->Draw("HIST SAME");
	legend->Draw();
	//h_dx_sim_n_bkg->Add(hist_p);
	//h_dx_sim_n_bkg->Add(hist_n);
	//h_dx_sim_n_bkg->Add(hist_bkg);

	cfit->cd(2);
	h_dx_sim_n_bkg->Draw();

	std::cout<<" Neutrons : "<<hist_n->Integral()<<endl;

	//fit_data->Draw();

	cfit->Print(Form("models_%s.pdf",kin));

} 

