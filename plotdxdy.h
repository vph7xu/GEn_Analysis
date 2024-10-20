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

// Function to create and return a square cut
inline void CreateSquareCut(double x1, double y1, double x2, double y2) {
    // Define the vertices of the square cut
    double x[5] = {x1, x2, x2, x1, x1};
    double y[5] = {y1, y1, y2, y2, y1};

    // Create a TCutG object for the square cut
    //TCutG *cut = new TCutG("squareCut", 5, x, y);
    TBox *box = new TBox(x1,y1,x2,y2);
    box->SetLineColor(kRed); // Set the line color to red
    box->SetLineWidth(2);    // Set the line width
    box->SetFillStyle(0);

    box->Draw();
}
