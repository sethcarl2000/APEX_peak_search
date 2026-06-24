
//project headers
#include <log_likelihood.hpp> 

//ROOT headers
#include <TCanvas.h>
#include <TH1D.h>s
#include <TRandom3.h> 
#include <Math/ProbFunc.h>
#include <TF1.h> 
#include <TAxis.h>
#include <TSystem.h> 
#include <TError.h> 

//stdlib headers
#include <cmath> 
#include <cstdio> 
#include <functional> 

/// @brief normalized PDF
/// @param x input var
/// @param x0
/// @param sigma 
inline double normal_dist(double x, double x0, double sigma)
{
    return 0.398942280401*std::exp( -0.5*std::pow((x-x0)/sigma, 2) ) / sigma; 
} 

int draw_toy_plot(int n_events_b, int n_events_s)
{
    if (gSystem->Load("build/libApexPeakSearch.so") < 0) {
        Error(__func__, "Something went wrong trying to load the ApexPeakSearch library");
        return -1; 
    };

    using namespace std; 

    using ROOT::Math::normal_cdf; 

    TRandom3 rand;

    auto hist = new TH1D("h", "toy data", 50, 3, 5);

    auto xax = hist->GetXaxis(); 
    const double dx = (xax->GetXmax()-xax->GetXmin())/((double)xax->GetNbins());
    const double sum_b = (double)n_events_b; 
    const double sum_s = (double)n_events_s; 

    const double sigma_b = 4;   
    const double sigma_s = 0.05; 

    const double x0_b = 3.1;
    const double x0_s = 4.2; 

    const double norm_s 
        = ROOT::Math::normal_cdf(xax->GetXmax(), sigma_s, x0_s) 
        - ROOT::Math::normal_cdf(xax->GetXmin(), sigma_s, x0_s); 
        
    const double norm_b 
        = ROOT::Math::normal_cdf(xax->GetXmax(), sigma_b, x0_b) 
        - ROOT::Math::normal_cdf(xax->GetXmin(), sigma_b, x0_b); 

    double xmin = xax->GetXmin();
    double xmax = xax->GetXmax(); 

    const double x_range = xmax - xmin; 

    std::function<double(double,double,double,double)> bin_expect = [x_range, dx, xax, xmin,xmax](double x, double total_events, double sigma, double x0)
    {
        //fraction of total events we expect to be in this bin
        double frac = normal_cdf(x+dx/2., sigma, x0) - normal_cdf(x-dx/2., sigma, x0); 
        frac = frac / (normal_cdf(xmax, sigma, x0) - normal_cdf(xmin, sigma, x0)); 
        //expectation value for number of events
        
        return frac * total_events; 
    };

    /*std::printf("result: %f\n", bin_expect(3.01, 1e4, 4, 3.1)); 
    return 0; */ 

   
    auto f_b = new TF1(
        "fcn_b", [bin_expect, sum_b](double *x, double *par){ 
            return bin_expect(x[0], sum_b, par[0], par[1]); }, 
        xax->GetXmin(), xax->GetXmax(), 2
    ); 
    f_b->SetParameters(sigma_b, x0_b);

    auto f_b_plus_s = new TF1(
        "fcn_bs", [bin_expect, sum_b, sum_s](double *x, double *par){ 
            return bin_expect(x[0], sum_b, par[0], par[1]) + bin_expect(x[0], sum_s, par[2], par[3]); }, 
        xax->GetXmin(), xax->GetXmax(), 4
    );

    

    f_b_plus_s->SetParameters(sigma_b,x0_b, sigma_s,x0_s); 

    for (int i=0; i<n_events_b; i++) { 
        double x; 
        do { x = rand.Gaus()*sigma_b + x0_b; } while (x < xax->GetXmin() || x > xax->GetXmax() ); 
        hist->Fill(x);
    }

    for (int i=0; i<n_events_s; i++) {
        double x; 
        do { x = rand.Gaus()*sigma_s + x0_s; } while (x < xax->GetXmin() || x > xax->GetXmax() ); 
        hist->Fill(x);
    }

    new TCanvas; 

    hist->SetMaximum( hist->GetMaximum()*1.2 );
    hist->SetMinimum(0.);
    hist->Draw("E"); 

    f_b_plus_s->SetLineColor(kRed);
    f_b_plus_s->Draw("SAME"); 
    
    f_b->SetLineColor(kBlack); 
    f_b->SetLineStyle(kDashed);
    f_b->Draw("SAME"); 
    

    return 1; 
}
