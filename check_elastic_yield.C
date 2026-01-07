// check_elastic_yield.C
// Usage (single file):
//   root -l -q 'check_elastic_yield.C("g4sbs_elastic.root")'
// Usage (many files with wildcard):
//   root -l -q 'check_elastic_yield.C("g4sbs_elastic_*.root")'

#include <TChain.h>
#include <TTree.h>
#include <TString.h>
#include <TMath.h>
#include <iostream>
#include <TTreeFormula.h>
using namespace std;

// ---- USER: fill from your run card / notes ----
static const double kGenerationVolume_sr = 0.21377;     // [sr]
static const double kLuminosity_sinv_cm2 = 1.7754e37;//1.1543e37;//1.7754e37;   // [s^-1 cm^-2]
static const double kNtries              = 2.0e4;       // total trials thrown
static const double kNormalization       = (kGenerationVolume_sr * kLuminosity_sinv_cm2) / kNtries; // [s^-1 cm^-2 * sr] -> [s^-1] after *sigma
static const double kBeamCurrent_Cps     = 60e-6;       // 60 µA = 60e-6 C/s

// ---- USER: put your *exact* coincidence selection here ----
// (Pick one "best pair" per event in your own analysis; here we assume
// a simple, single-candidate selection using cluster-level branches.)
static const char* kCoincCut = "";//"Earm.BBGEM.Track.ntracks>0&&ev.pmperp<0.1&&abs(ev.pmpar)<0.1&&ev.W2<1.6&&Earm.BBPSTF1.det.esum>0.2&&Earm.BBSHTF1.det.esum>0.01&&Harm.HCalScint.det.esum>0.325&&abs(((Earm.BBPSTF1.det.esum+Earm.BBSHTF1.det.esum)/Earm.BBGEM.Track.P)-1)<0.2";
 /* "Earm.BBGEM.Track.ntracks>0"
 &&Earm.BBGEM.Track.ntracks>0
 &&Earm.BBSHTF1.det.esum>0
 ev.earmaccept==1&&
  " && Earm.BBCal.clus.nclus>0 && Earm.BBCal.clus.e>0.5"         // e-arm calorimeter threshold (example)
  " && Harm.HCalScint.clus.nclus>0 && Harm.HCalScint.clus.e>0.03"// HCAL threshold (example)
  // timing: replace 'TOFcorr' with your expression or leave as a loose window
  " && abs(Harm.HCalScint.clus.t - Earm.BBCal.clus.t) < 15";
*/
void check_elastic_yield(const char* files)
{
  TChain C("T");
  C.Add(files);


  int nentries = C.GetEntries();

  C.SetBranchStatus("ev*",1);           // enable the whole ev branch (and its leaves)
  TLeaf* lSigma = C.GetLeaf("ev.sigma");
  TLeaf* lRate  = C.GetLeaf("ev.rate");

  // Create a TTreeFormula to evaluate the coincidence cut per-entry
  TTreeFormula *fCoinc = new TTreeFormula("fCoinc", kCoincCut, &C);

  double total_sigma=0, total_rate=0;
  Long64_t n = C.GetEntries();
  for (Long64_t i=0;i<n;++i) {
    C.GetEntry(i);
    // Only add to total_sigma if this entry passes the coincidence cut
    if (lSigma) {
      if (fCoinc->EvalInstance()) total_sigma += lSigma->GetValue();
    }
    // total_rate remains the full sum over all entries (leave unchanged)
    if (fCoinc->EvalInstance()) total_rate  += lRate ->GetValue();
  }

  delete fCoinc;

  // --- Sum ev.sigma over your coincidence selection ---
  C.Draw("ev.sigma", kCoincCut, "goff");
  Long64_t nSelSigma = C.GetSelectedRows();
  if (nSelSigma<=0) {
    cout << "[WARN] No events selected for sigma with your cut.\n";
  }
  double sum_sigma = 0.0;
  for (Long64_t i=0;i<nSelSigma;++i) sum_sigma += C.GetV1()[i]; // ev.sigma is cm^2/sr for elastic

  // --- Sum ev.rate (elastic-only shortcut; already has GenVol*Lumi/Ntries baked in) ---
  C.Draw("ev.rate", kCoincCut, "goff");
  Long64_t nSelRate = C.GetSelectedRows();
  if (nSelRate<=0) {
    cout << "[WARN] No events selected for rate with your cut.\n";
  }
  double sum_rate = 0.0;
  for (Long64_t i=0;i<nSelRate;++i) sum_rate += C.GetV1()[i]; // Hz (s^-1)

  // --- Yields per Coulomb ---
  const double Y_perC_sigma = (kNormalization / kBeamCurrent_Cps) * sum_sigma; // [events/C]
  const double Y_perC_rate  = (1.0 / kBeamCurrent_Cps) * sum_rate;             // [events/C]

  const double Y_perC_sigma_total = (kNormalization / kBeamCurrent_Cps) * total_sigma; // [events/C]
  const double Y_perC_rate_total  = (1.0 / kBeamCurrent_Cps) * total_rate;             // [events/C]

  // --- Consistency ratio (should be ~1 for elastic) ---
  const double expect_rate_from_sigma = sum_sigma * kNormalization; // should equal sum_rate
  const double ratio = (sum_rate>0.0) ? ( sum_rate / expect_rate_from_sigma ) : 0.0;

  const double expect_rate_from_sigma_total = total_sigma * kNormalization; // should equal sum_rate
  const double ratio_total = (total_rate>0.0) ? ( total_rate / expect_rate_from_sigma ) : 0.0;

  cout.setf(std::ios::scientific); cout.precision(6);
  cout << "\n=== Elastic coincidence yield check ===\n";
  cout << "  Files           : " << files << "\n";
  cout << "  Cut             : " << kCoincCut << "\n";
  cout << "  N_selected (σ)  : " << nSelSigma << "\n";
  cout << "  N_selected (rate): " << nSelRate  << "\n";
  cout << "  Sum(ev.sigma)   : " << sum_sigma << "  [cm^2/sr]\n";
  cout << "  Sum(ev.rate)    : " << sum_rate  << "  [s^-1]\n";
  cout << "  Total(ev.sigma) : " << total_sigma<< " [cm^2/sr]\n";
  cout << "  Total(ev.rate)  : " << total_rate<< " [s^-1]\n";
  
  cout << "  Normalization   : " << kNormalization << "  [s^-1]\n";
  cout << "---------------------------------------------\n";
  cout << "  Yield per C (σ) : " << Y_perC_sigma << "  [events/C]\n";
  cout << "  Yield per C (rate): " << Y_perC_rate  << "  [events/C]\n";
  cout << "  Yield per C (σ) total: " << Y_perC_sigma_total << "  [events/C]\n";
  cout << "  Yield per C (rate) total: " << Y_perC_rate_total  << "  [events/C]\n";
  cout << "  Consistency ratio = sum(rate) / (sum(sigma)*Norm) = " << ratio << "\n";
  cout << "  (Expect ~1.0 for elastic if cuts/units match and no double-counting)\n\n";
}

