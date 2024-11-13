#include "cuts.h"
#include "plotdxdy.h"


void Hcal_time_energy_v1(const char* filename,const char* printfilename,const char* kin){

	std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt",kin)); //parse the cuts
        cuts cutsobject;
        cutsobject.parsecuts(config);

        TFile* file = TFile::Open(filename);
        TTree* tree = (TTree*)file->Get("Tout");

	double hcal_clus_mem_e[1000];
	double hcal_clus_mem_atime[1000];
	double hcal_clus_mem_tdctime[1000];
	double hcal_time = 0.0;
	double hcal_clus_mem_id[1000];
	double nblk_HCAL = 0.0;
	//double hcal_clus_nblk[1000];

	tree->SetBranchAddress("hcal_clus_mem_e",&hcal_clus_mem_e);
	tree->SetBranchAddress("hcal_clus_mem_atime",&hcal_clus_mem_atime);
	tree->SetBranchAddress("hcal_clus_mem_tdctime",&hcal_clus_mem_tdctime);
	tree->SetBranchAddress("hcal_time",&hcal_time);
	//tree->SetBranchAddress("hcal_clus_nblk",&hcal_clus_nblk);
	tree->SetBranchAddress("hcal_clus_mem_id",&hcal_clus_mem_id);
	tree->SetBranchAddress("nblk_HCAL",&nblk_HCAL);

	TH2D *hcal_clus_avg_time_offset = new TH2D("hcal_clus_avg_time_offset","hcal clus average time offset",300,0,300,400,-100,100);

	//TH1D *h_cal_clus_e = new TH1D("h_cal_clus_e","hcal cluster energy",2000,0,2);

	//TGraph *hcal_clus_atime_g = new TGraph();
	//TGraph *hcal_clus_avg_time_g = new TGraph();

	//vector of ids
	std::vector<double> ids;
	for(int i=0; i<289;i++){
		ids.push_back(i);
	}

	//vector to store TGraphs
	std::vector<TGraph*> graphs;

	int nentries = tree->GetEntries();
	int nentriesMAX = 100000;	//limiting the number of events
	
	//double clus_avg_atime = 0.0;

	double clus_avg_times[nentriesMAX];

	for (double id : ids){
		TGraph* graph = new TGraph();
		graph->SetTitle(Form("Graph for blk ID %f",id));
		for (int i = 0; i<nentriesMAX; i++){
			tree->GetEntry(i);
			//std::cout<<"number of blks in the primary cluster : "<<nblk_HCAL<<endl;
			double clus_atime_avg = 0.0;
			for (int j = 0; j<nblk_HCAL;j++){
				if (static_cast<int>(id)==0){
					clus_atime_avg+=(hcal_clus_mem_atime[j]/(nblk_HCAL));
					//std::cout<<"i : "<< i<<" j : "<<j<<" hcal_clus_atime_avg : "<<clus_atime_avg<<endl;
					//std::cout<<"j : "<<j<<" i : "<<i<<" hcal_adc_time : "<<hcal_clus_mem_atime[j]<<" blk_id : "<<id<<endl;
					clus_avg_times[i]=clus_atime_avg;
				}
				//std::cout<<"hcal clus id : "<<hcal_clus_id[j]<<endl;
				if (static_cast<int>(hcal_clus_mem_id[j])==static_cast<int>(id) and static_cast<int>(id)>0){
					//std::cout<<"hcal_adc_time_avg : "<<clus_atime_avg<<endl;
					//std::cout<<"****hcal_time : "<<hcal_time<<"hcal_time_blk_0 : "<<hcal_clus_mem_atime[0]<<endl;
					//std::cout<<"hit, atime : "<<clus_avg_times[i]<<" entry : "<<i<<" blk_id : "<<id<<endl;
					graph->SetPoint(i,i,(hcal_clus_mem_atime[j]-clus_avg_times[i]));
					hcal_clus_avg_time_offset->Fill(static_cast<int>(id),(hcal_clus_mem_atime[j]-clus_avg_times[i]));
					//count++;
				}
				else graph->SetPoint(1,-10,-10);
			}
			//reset the array
			//for (int i =0; i<10){

			//}

		}
		graphs.push_back(graph);

		std::cout << (id * 100.0/288) << "% \r";
                std::cout.flush();

	}


   	const int graphsPerCanvas = 25;
   	int canvasCount = 0;

   	for (size_t i = 0; i < graphs.size(); i += graphsPerCanvas) {
   	    // Create a new canvas
   	    TCanvas* canvas = new TCanvas(Form("canvas_%d", canvasCount), Form("Graphs Canvas %d", canvasCount), 800, 600);
   	    canvas->Divide(5, 5); // Divide the canvas into 25 pads

   	    for (size_t j = 0; j < graphsPerCanvas && (i + j) < graphs.size(); ++j) {
   	        canvas->cd(j + 1);
   	        if (graphs[i + j]) {
   	            graphs[i + j]->Draw("AP");
   	        }
   	    }
  	    // Save the canvas
	    if (i==0) canvas->Print("../plots/hcal_primary_adc_time_per_blk.pdf(");
	    else if (i==275) canvas->Print("../plots/hcal_primary_adc_time_per_blk.pdf)");
	    else canvas->Print("../plots/hcal_primary_adc_time_per_blk.pdf");
	    
	    
	    canvasCount++;

  	    // Clean up the canvas
  	    delete canvas;
  	}

  	// Clean up the graphs
  	for (TGraph* graph : graphs) {
  	    delete graph;
  	}

  	gStyle->SetPalette(kRainBow);
	hcal_clus_avg_time_offset->SetMinimum(0);
	hcal_clus_avg_time_offset->SetMaximum(500);

  	TCanvas* c1 = new TCanvas("c1","2d distribution of timing offset",800,600);
	hcal_clus_avg_time_offset->Draw("COLZ");
	c1->Print("../plots/timimg_offset_distribution.pdf");
	
	
	/*
	int count = 0;

	for (int i = 0; i<nentriesMAX; i++){
		tree->GetEntry(i);
		std::cout<<"hcal e : "<<hcal_clus_e[9]<<endl;
		std::cout<<"hcal time : "<<hcal_time<<endl;
		std::cout<<"hcal adc time : "<<hcal_clus_atime[0]<<endl;
		std::cout<<"hcal tdc time : "<<hcal_clus_tdctime[0]<<endl;
		std::cout<<"hcal nblk : "<<hcal_clus_nblk[0]<<" "<<hcal_clus_nblk[1]<<endl;
		std::cout<<"hcal blk id : "<<hcal_clus_nblk[0]<<endl;

		for (int j = 0; j<9; j++){
			std::cout<<"adc time : "<<hcal_clus_atime[j]<<endl;
			hcal_clus_atime_g->SetPoint(count,count,hcal_clus_atime[j]);
			hcal_clus_avg_time_g->SetPoint(count,count,hcal_time);
			count++;
		}
	}
	
	TCanvas* c = new TCanvas("c","c",1200,1600);
	hcal_clus_atime_g->SetMarkerStyle(8);
	hcal_clus_avg_time_g->SetMarkerStyle(7);

	hcal_clus_atime_g->SetMarkerSize(1);
	hcal_clus_avg_time_g->SetMarkerSize(1);

	hcal_clus_avg_time_g->SetMarkerColor(kRed);

	hcal_clus_atime_g->Draw("AP");
	hcal_clus_avg_time_g->Draw("P");
	c->Print(Form("%s_Hcal_atime.pdf",printfilename));
	*/
}
