#include "models.h"
#include "parse.h"
#include "include/KinematicVar.h"

void pion_correction(const char *pion_filename){
    TChain *TC_pion = new TChain("T");
    TC_pion->Add(pion_filename);

    double PSe_pion = 0.0;

    if (!TC_pion->GetBranch("bb.ps.clus.e")) {
        std::cerr << "Branch bb.ps.clus.e does not exist!" << std::endl;
        return;
    }
    TC_pion->SetBranchAddress("bb.ps.clus.e", &PSe_pion);

    TH1D *h_PSe_pion = new TH1D("h_PSe_pion", "Preshower cluster E from pion simulation", 100, 0, 1);

    int nentries_pion = TC_pion->GetEntries();
    if (nentries_pion == 0) {
        std::cerr << "No entries in TChain!" << std::endl;
        return;
    }

    for (int i = 0; i < nentries_pion; i++) {
        TC_pion->GetEntry(i);
        h_PSe_pion->Fill(PSe_pion);

        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries_pion) << "% of simulation file processed\r";
            std::cout.flush();
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Canvas", 800, 600);
    h_PSe_pion->Draw();
}
