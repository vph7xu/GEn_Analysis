#include "cuts.h"
#include "plotdxdy.h"

void targetpol(const char* filename, const char* printfilename, const char *kin) {

    std::map<std::string, std::string> config = parseConfig(Form("cuts/cut_%s.txt", kin)); //parse the cuts

    cuts cutsobject;
    cutsobject.parsecuts(config);

    TFile* file = TFile::Open(filename);
    TTree* tree = (TTree*)file->Get("Tout");

    double He3Pol = 0.0;
    int runnum = 0;
    double ntrack = 0;
    double dx = 0.0; 
    double dy = 0.0; 
    double W2 = 0.0;
    double Q2 = 0.0;
    double coin_time = 0.0;
    double eHCAL = 0.0;

    tree->SetBranchAddress("He3Pol", &He3Pol);
    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("ntrack", &ntrack);
    tree->SetBranchAddress("dx", &dx);
    tree->SetBranchAddress("dy", &dy);
    tree->SetBranchAddress("W2", &W2);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("coin_time", &coin_time);
    tree->SetBranchAddress("eHCAL", &eHCAL);

    // -- Axes limits for histograms by Kinematic
    double upper_limit = 0;
    double lower_limit = 0;

    std::string kin_str = kin;

    if (kin_str == "GEN2_He3") {
        lower_limit = 2000;
        upper_limit = 2400;
    }
    else if (kin_str == "GEN3_He3") {
        lower_limit = 2500;
        upper_limit = 3300;
        std::cout << "HERE" << std::endl;
    }
    else if (kin_str == "GEN4_He3") {
        lower_limit = 3500;
        upper_limit = 4600;
    }
    else if (kin_str == "GEN4b_He3") {
        lower_limit = 5000;
        upper_limit = 6100;
    }

    // Create histograms (numeric axis from lower_limit to upper_limit)
    TH1I* h_runnum = new TH1I("h_runnum", "Run Number Distribution;Run Number;Counts",
                              upper_limit - lower_limit, lower_limit, upper_limit);
    TH1I* h_runnum_no_norm = new TH1I("h_runnum_no_norm", "Run Number Distribution;Run Number;Counts",
                                      upper_limit - lower_limit, lower_limit, upper_limit);
    TH1D* h_He3Pol = new TH1D("h_He3Pol", "He3 Polarization", 100, 0, 100);
    TH2D* h_He3Polvrunnum = new TH2D("h_He3Polvrunnum", "He3Pol vs Run Number;Run Number;He3Pol",
                                     upper_limit - lower_limit, lower_limit, upper_limit,
                                     100, 0, 100);

    // TGraphs
    TGraph* He3Polgraph = new TGraph();
    TGraph* numberofeventsgraph = new TGraph();

    int nentries = tree->GetEntries();
    int runx = 0;
    bool QE_cut = false;

    // We'll also collect all actual run numbers here for the bar-graph approach:
    std::set<int> runSet;

    // -------------------
    //  MAIN EVENT LOOP
    // -------------------
    for (int i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // Kinematic cut
        QE_cut = ( eHCAL > cutsobject.eHCAL_L &&
                   (cutsobject.W2_L < W2 && W2 < cutsobject.W2_H) &&
                   (cutsobject.coin_time_L < coin_time && coin_time < cutsobject.coin_time_H) &&
                   (dx < cutsobject.dx_H && cutsobject.dx_L < dx) &&
                   (dy < cutsobject.dy_H && cutsobject.dy_L < dy) );

        if (ntrack > 0 && QE_cut) {
            h_runnum->Fill(runnum);
            h_runnum_no_norm->Fill(runnum);
            h_He3Polvrunnum->Fill(runnum, He3Pol);
            h_He3Pol->Fill(He3Pol);

            // Collect actual run numbers
            runSet.insert(runnum);
        }

        // Build "He3Polgraph" for the first TCanvas
        if (i == 0) {
            runx = runnum;
            He3Polgraph->SetPoint(He3Polgraph->GetN(), runx, He3Pol);    
        }
        if (runx != runnum) {
            He3Polgraph->SetPoint(He3Polgraph->GetN(), runx, He3Pol);
            runx = runnum;
        }

        if (i % 1000 == 0) {
            std::cout << (i * 100.0 / nentries) << "% \r";
            std::cout.flush();
        }
    } // end event loop

    // Fill "numberofeventsgraph" from the histogram "h_runnum"
    int nBins = h_runnum->GetNbinsX();
    numberofeventsgraph->Set(0); // Clear existing points in the graph
    for (int i = 1; i <= nBins; ++i) {
        double xCenter = h_runnum->GetBinCenter(i);
        double yContent = h_runnum->GetBinContent(i);
        if (yContent > 0) {
            numberofeventsgraph->SetPoint(numberofeventsgraph->GetN(), xCenter, yContent);
        }
    }

    // Optionally set the metadata for consistency
    numberofeventsgraph->SetTitle(h_runnum->GetTitle());
    numberofeventsgraph->GetXaxis()->SetTitle(h_runnum->GetXaxis()->GetTitle());
    numberofeventsgraph->GetYaxis()->SetTitle(h_runnum->GetYaxis()->GetTitle());

    // -------------------
    //  Canvas 1: Histos
    // -------------------
    TCanvas* c1 = new TCanvas("c1", "Run Number and He3Pol", 800, 600);
    c1->Divide(2, 2);

    c1->cd(1);
    h_He3Polvrunnum->Draw();
    c1->cd(2);
    h_He3Pol->Draw();
    c1->cd(3);
    h_runnum_no_norm->Draw();
    c1->cd(4);
    // (Unused pad, or add something else here if you want)

    c1->Update();
    c1->SaveAs(Form("plots/targetpol_%s.pdf", kin));

    // -------------------------------
    //  Canvas 2: Dual TGraph Overlay
    // -------------------------------
    TCanvas* c2 = new TCanvas("c2", "Overlay TGraphs with Dual Y-Axes", 800, 600);

    // Draw the first graph (He3Polgraph)
    He3Polgraph->SetMarkerStyle(20);
    He3Polgraph->SetMarkerColor(kBlue);
    He3Polgraph->SetLineColor(kBlue);  // line color won't be visible if no "L" in draw option
    He3Polgraph->Draw("AP");

    // Adjust the primary y-axis
    He3Polgraph->GetYaxis()->SetTitle(Form("%s Polarization percentage",kin));
    He3Polgraph->GetYaxis()->SetRangeUser(0, 100);
    He3Polgraph->GetXaxis()->SetTitle("Run Number");
    He3Polgraph->GetXaxis()->SetLimits(lower_limit, upper_limit);

    // Create the secondary y-axis
    double yMinPrimary = 0, yMaxPrimary = 100;  // range of primary y-axis
    double yMinSecondary = 0, yMaxSecondary = h_runnum->GetMaximum(); // range for # events

    // Scale the second graph (red graph)
    TGraph* scaledGraph = new TGraph();
    for (int i = 0; i < numberofeventsgraph->GetN(); ++i) {
        double x, y;
        numberofeventsgraph->GetPoint(i, x, y);
        double scaledY = yMinPrimary + (y - yMinSecondary) *
                         (yMaxPrimary - yMinPrimary) /
                         (yMaxSecondary - yMinSecondary + 1e-9);
        scaledGraph->SetPoint(i, x, scaledY);
    }

    // Draw the scaled second graph
    scaledGraph->SetMarkerStyle(21);
    scaledGraph->SetMarkerColor(kRed);
    scaledGraph->SetLineColor(kRed);
    scaledGraph->Draw("P SAME");

    // Add the secondary Y-axis on the right
    double leftX  = gPad->GetUxmin();  
    double rightX = gPad->GetUxmax();

    TGaxis* rightAxis = new TGaxis(upper_limit - 50, yMinPrimary, upper_limit - 50, yMaxPrimary,
                                   yMinSecondary, yMaxSecondary, 510, "+L");
    rightAxis->SetLineColor(kRed);
    rightAxis->SetLabelColor(kRed);
    rightAxis->SetTitle("Number of QE Events");
    rightAxis->SetTitleColor(kRed);
    rightAxis->Draw();

    c2->Update();
    c2->SaveAs(Form("plots/dual_axis_%s.pdf", kin));

    // ------------------------------------------------------------------------
// ADDITIONAL CODE to generate BAR GRAPHS with a COMPRESSED X-axis
// (No missing run numbers in between), skipping some labels
// ------------------------------------------------------------------------

// 1) We have all real runs in runSet (unique, sorted).
//    Let's create a "bin index" for each run.
int nRuns = runSet.size();
if (nRuns == 0) {
    std::cerr << "No valid runs found. Skipping bar graph.\n";
    return;
}

// 2) Create new histograms with exactly nRuns bins, from 0.5 to nRuns+0.5
TH1D* h_bar_He3Pol = new TH1D("h_bar_He3Pol",
                              Form("%s ;Run Number;He3Pol percentage",kin),
                              nRuns, 0.5, nRuns + 0.5);
TH1D* h_bar_nEvents = new TH1D("h_bar_nEvents",
                               Form("%s ;Run Number;Events",kin),
                               nRuns, 0.5, nRuns + 0.5);

// 3) Map each run number => integer bin index
//    Then label that bin with the run number
std::map<int, int> runToBin;
{
    int binIndex = 1;
    for (auto r : runSet) {
        runToBin[r] = binIndex;

        // Label the bin with the actual run number
        // We'll skip most labels so they don't overlap: only label every 5th bin
        // Adjust "5" if you want a different spacing
        if ((binIndex % 25) == 1) {
            h_bar_He3Pol->GetXaxis()->SetBinLabel(binIndex, Form("%d", r));
            h_bar_nEvents->GetXaxis()->SetBinLabel(binIndex, Form("%d", r));
        }
        else {
            // Remove the label for bins we are skipping
            h_bar_He3Pol->GetXaxis()->SetBinLabel(binIndex, "");
            h_bar_nEvents->GetXaxis()->SetBinLabel(binIndex, "");
        }

        binIndex++;
    }
}

// 4) Fill the bar histograms from the TGraphs (He3Polgraph, numberofeventsgraph)
//    Each TGraph point has (x=run, y=value).
//    Convert run => bin index => fill
for (int i = 0; i < He3Polgraph->GetN(); ++i) {
    double x, y;
    He3Polgraph->GetPoint(i, x, y);
    int run = (int)x;  // run number
    if (runToBin.find(run) != runToBin.end()) {
        int bin = runToBin[run];
        h_bar_He3Pol->SetBinContent(bin, y);
    }
}
for (int i = 0; i < numberofeventsgraph->GetN(); ++i) {
    double x, y;
    numberofeventsgraph->GetPoint(i, x, y);
    int run = (int)x;  // run number
    if (runToBin.find(run) != runToBin.end()) {
        int bin = runToBin[run];
        h_bar_nEvents->SetBinContent(bin, y);
    }
}

// 5) Apply semi-transparent fill colors so both bar graphs remain visible
h_bar_He3Pol->SetLineColor(kBlue);
h_bar_He3Pol->SetFillColorAlpha(kBlue, 0.40);  // 40% opacity
h_bar_nEvents->SetLineColor(kRed);
h_bar_nEvents->SetFillColorAlpha(kRed, 0.40);  // 40% opacity

h_bar_He3Pol->SetStats(false);
h_bar_nEvents->SetStats(false);

// 6) Optional: zero-suppress leading/trailing empty bins
int firstBinHe3    = h_bar_He3Pol->FindFirstBinAbove(0);
int lastBinHe3     = h_bar_He3Pol->FindLastBinAbove(0);
int firstBinEvents = h_bar_nEvents->FindFirstBinAbove(0);
int lastBinEvents  = h_bar_nEvents->FindLastBinAbove(0);

if (firstBinHe3 <= 0)     firstBinHe3 = 1;
if (firstBinEvents <= 0)  firstBinEvents = 1;
if (lastBinHe3 <= 0)      lastBinHe3 = 1;
if (lastBinEvents <= 0)   lastBinEvents = 1;

int firstBin = TMath::Min(firstBinHe3, firstBinEvents);
int lastBin  = TMath::Max(lastBinHe3,  lastBinEvents);
if (lastBin > nRuns) lastBin = nRuns;

h_bar_He3Pol->GetXaxis()->SetRange(firstBin, lastBin);
h_bar_nEvents->GetXaxis()->SetRange(firstBin, lastBin);

// 7) Increase bar width (and optionally offset them)
h_bar_He3Pol->SetBarWidth(0.7);
h_bar_nEvents->SetBarWidth(0.7);

// 8) Create a new canvas for the bar graphs
TCanvas* c3 = new TCanvas("c3", "Bar Graphs with Dual Y-Axes (Compressed X)", 900, 600);

// (A) Increase bottom margin for large labels
gPad->SetBottomMargin(0.25);

// (B) Reduce label font size
h_bar_He3Pol->GetXaxis()->SetLabelSize(0.03);
h_bar_nEvents->GetXaxis()->SetLabelSize(0.03);

// (C) Rotate labels (vertical: "v", or angled: ">" or "d", etc.)
h_bar_He3Pol->LabelsOption("v");
h_bar_nEvents->LabelsOption("v");

// 9) Draw the first bar graph (He3Pol) on the left axis
double he3polMax = 1.2 * h_bar_He3Pol->GetMaximum();
if (he3polMax < 100) he3polMax = 100;
h_bar_He3Pol->SetMaximum(he3polMax);

h_bar_He3Pol->Draw("bar");

// 10) Scale the second histogram (h_bar_nEvents) to overlay
double primaryMin   = 0.0;
double primaryMax   = he3polMax;
double secondaryMin = 0.0;
double secondaryMax = h_bar_nEvents->GetMaximum();

TH1D* h_bar_nEvents_scaled = (TH1D*)h_bar_nEvents->Clone("h_bar_nEvents_scaled");
for (int b = 1; b <= h_bar_nEvents_scaled->GetNbinsX(); ++b) {
    double val = h_bar_nEvents_scaled->GetBinContent(b);
    double valScaled = primaryMin + (val - secondaryMin) *
                       (primaryMax - primaryMin) /
                       (secondaryMax - secondaryMin + 1e-9);
    h_bar_nEvents_scaled->SetBinContent(b, valScaled);
}

// Draw the scaled second bar histogram on top of the first
h_bar_nEvents_scaled->Draw("bar same");

// 11) Add the secondary Y-axis on the right
gPad->Update();
double axisX   = gPad->GetUxmax();
double axisYmin = gPad->GetUymin();
double axisYmax = gPad->GetUymax();

TGaxis* rightAxis2 = new TGaxis(axisX, axisYmin, axisX, axisYmax,
                                secondaryMin, secondaryMax, 510, "+L");
rightAxis2->SetTitle("Number of QE Events");
rightAxis2->SetLineColor(kRed);
rightAxis2->SetTitleColor(kRed);
rightAxis2->SetLabelColor(kRed);
rightAxis2->Draw();

// 12) Create and draw a legend
TLegend* legend = new TLegend(0.65, 0.75, 0.88, 0.88);
legend->SetBorderSize(0);
legend->SetFillStyle(0);
legend->AddEntry(h_bar_He3Pol,         "He3Pol (%)",       "f");
legend->AddEntry(h_bar_nEvents_scaled, "Number of Events", "f");
legend->Draw();

// Finally update and save
c3->Update();
c3->SaveAs(Form("plots/bar_graphs_with_dual_yaxes_%s.pdf", kin));

}

