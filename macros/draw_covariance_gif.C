#include <TLatex.h>
#include <TRandom3.h>
#include <TText.h> 
#include <TCanvas.h> 
#include <TStyle.h> 
#include <TH2D.h> 
// Eigen 
#include <eigen3/Eigen/Core>
#include <eigen3/Eigen/Dense>
// stdlib
#include <memory> 


void draw_covariance_gif(double cov=0.)
{   
    TRandom3 mrand; 

    unsigned stats = 1e4; 

    auto canv = new TCanvas("c", "canv", 500, 500); 
    gStyle->SetOptStat(0); 

    double min_cov =0.; 
    double max_cov =1.; 
    double d_cov = 0.05; 
    
    int n_steps=21; 
    
    double ymax = 5.5; 

    auto hist = new TH2D("h_test", "Covariance of gaussian variables", 100, -ymax,+ymax, 100, -ymax,+ymax); 

    for (unsigned j=0; j<stats; j++) {

        double v1 = std::sqrt(1. + cov)*mrand.Gaus();
        double v2 = std::sqrt(1. - cov)*mrand.Gaus();   

        if (v2 != v2) v2 =0.; 
        
        hist->Fill( v1+v2, v1-v2 ); 
    }
    hist->SetTitle(Form("Cov(y_{1}, y_{2}) = %.3f;y_{1};y_{2}",cov)); 

    hist->Draw("col"); 

    return; 
}