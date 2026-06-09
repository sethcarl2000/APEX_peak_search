
#include "fit_exponential_poly.hpp"

#include <TRandom3.h> 
#include <Math/ProbFunc.h>
#include <Math/Functor.h>
#include <Math/Factory.h>
#include <TF1.h> 
#include <TCanvas.h>
#include <TSystem.h> 
#include <TH1D.h> 
#include <TError.h> 

#include <cmath>


//this macro is designed to test different Profile Likelihood Ratio (PLR) tests on varieties of TOY (made-up, unrealistic) data.

namespace {

    //snr = the likelihood that any event is 'accidental' or not.
    //x is our 'primary' search variable
    constexpr double x_min = 100.;
    constexpr double x_max = 400.; 
    
    constexpr double pi = 3.1415926536;

    //standerd toy 'event'
    struct event_t { 
        double x; //x-value
        double p_coinc; // probability that this event represents a true 'coincidence'   
    };

    //this will be the background event generator 
    double generate_bg_accidental(TRandom3& rand) {

        static constexpr double x0_bg = 250.; 
        static constexpr double sigma_bg = 45.; 

        double x; 
        do {
            x = rand.Gaus()*sigma_bg + x0_bg;
    
            //now, bias it a bit 
            x += 50.*std::sin( pi*(x-x_min)/((x_max-x_min)) ); 

        } while (x < x_min || x > x_max); 

        return x; 
    }

}

int toy_peak_search(const long n_events)
{
    using peak_search::ExponentialPoly; 

    if (gSystem->Load("build/libApexPeakSearch.so") < 0) {
        Error(__func__, "Error loading ApexPeakSearch lib"); 
        return -1; 
    };  
    
    auto hist = new TH1D("h", "background;x;n. events", 200, x_min, x_max);

    TRandom3 rand; 

    for (long i=0; i<n_events; i++) hist->Fill( generate_bg_accidental(rand) );

    auto fit_result = peak_search::fit_exponential_poly(hist, 5); 
    
    
    if (!(bool)fit_result) {
        Error(__func__, "Something went wrong doing the background fit"); 
        return -1; 
    }
    auto poly = (ExponentialPoly)fit_result; 

    auto xax = hist->GetXaxis(); 
    double xmin = xax->GetXmin();
    double xmax = xax->GetXmax(); 

    double dx = (xmax - xmin)/((double)xax->GetNbins());

    auto tf1 = new TF1("exp_poly", [poly, dx](double *x, double *par){ return dx*poly(x[0]); }, xmin, xmax, 0);

    new TCanvas; 
    hist->SetMaximum(hist->GetMaximum()*1.2);
    hist->SetMinimum(0.);
    hist->Draw("E"); 

    tf1->Draw("SAME"); 

    return 0; 
}