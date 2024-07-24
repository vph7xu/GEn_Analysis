// Function to create an oval cut
TCutG* CreateOvalCut(const char* name, double centerX, double centerY, double a, double b, int nPoints) {
    double *x = new double[nPoints];
    double *y = new double[nPoints];
    
    for (int i = 0; i < nPoints; i++) {
        double angle = 2 * TMath::Pi() * i / nPoints;
        x[i] = centerX + a * cos(angle);
        y[i] = centerY + b * sin(angle);
    }
    
    TCutG *cutg = new TCutG(name, nPoints, x, y);
    cutg->SetLineColor(kRed);
    cutg->SetLineWidth(2);
    
    delete[] x;
    delete[] y;
    
    return cutg;
}
