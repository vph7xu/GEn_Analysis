#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TStyle.h>

void plot_f_bkg_with_error() {
    // Data
    const int nPoints = 10;
    double eHCAL_cut[nPoints] = {0.000, 0.100, 0.125, 0.150, 0.175, 0.200, 0.225, 0.250, 0.275, 0.300};
    double f_bkg[nPoints] = {0.643748, 0.611040, 0.596170, 0.571912, 0.544291, 0.521407, 0.502419, 0.475389, 0.458343, 0.426655};
    double error_f_bkg[nPoints] = {0.004971, 0.005033, 0.005139, 0.005235, 0.005336, 0.005480, 0.005685, 0.005824, 0.006093, 0.006221};
    double eHCAL_error[nPoints] = {0}; // No error for x-axis values

    // Create a TGraphErrors
    TGraphErrors *graph = new TGraphErrors(nPoints, eHCAL_cut, f_bkg, eHCAL_error, error_f_bkg);

    // Create a canvas
    TCanvas *c1 = new TCanvas("c1", "f_bkg vs. eHCAL cut with error bars", 1600, 1200);
    c1->SetGrid();

    // Set graph style
    graph->SetMarkerStyle(21); // Circle markers
    graph->SetMarkerColor(kBlue);
    graph->SetLineColor(kBlue);

    // Set axis titles
    graph->GetXaxis()->SetTitle("eHCAL cut");
    graph->GetYaxis()->SetTitle("f_bkg");
    graph->SetTitle("f_bkg vs. eHCAL cut");

    // Draw the graph
    graph->Draw("AP"); // "A" for axis, "P" for points

    // Save the canvas as an image
    c1->SaveAs("plots/f_bkg_vs_eHCAL_cut.png");
}