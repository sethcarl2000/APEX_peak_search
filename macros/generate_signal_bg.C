#ifndef generate_signal_bg_C
#define generate_signal_bg_C

#include <bininfo.hpp>
#include <fit_exponential_poly.hpp>
//ROOT headers
#include <TRandom3.h> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TAxis.h> 
#include <functional> 
//stdlib headers
#include <cstdio> 
#include <iostream> 
#include <cmath> 

namespace {

    const double xmin{-7}, xmax{+7}; 

    const double sigma_bg = 2.5;
    const double bg_skew = 0.5; 

    const double sigma_signal = 0.05; 
    const double x0_signal = 2.; 

    double gen_bg(TRandom3& rand) {
        double ret; 
        do { 
            ret = rand.Gaus()*sigma_bg; 
            ret += bg_skew * std::sin( ret*3.1415925636/(sigma_bg*3.) ); 
        } while (ret > xmax || ret < xmin); 
        return ret; 
    }

    double gen_signal(TRandom3& rand) {
        double ret; 
        do { ret=rand.Gaus()*sigma_signal + x0_signal; } while (ret > xmax || ret < xmin); 
        return ret; 
    }

}

void generate_signal_bg(unsigned int n_bg, unsigned int n_signal)
{   

    TRandom3 myrand; 

    auto hist = new TH1D("h_combined", "Signal + background;X;", 100, xmin, xmax); 

    std::cout << 
        "generating bg..." << std::flush;
    
    for (unsigned int i=0; i<n_bg; i++) { hist->Fill(gen_bg(myrand)); }


    std::cout << 
        "done.\n"
        "generating signal..." << std::flush; 
    
    for (unsigned int i=0; i<n_signal; i++) { hist->Fill(gen_signal(myrand)); }
    
    std::cout << 
        "done.\n"; 


    new TCanvas; 
    hist->Draw(); 

    hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum()*1.1); 

    auto data = peak_search::copy_1D_hist(hist); 

    auto exp_poly = (peak_search::ExponentialPoly)peak_search::fit_exponential_poly(hist, 6); 

    int i=0; 
    std::cout << "coeffs:\n"; 
    for (auto& coeff : exp_poly.GetParams()) {
        std::printf("   %3i : %f\n", i++, coeff); 
    }

    const double dx = (xmax - xmin)/((double)hist->GetXaxis()->GetNbins()); 

    auto fcn = (std::function<double(double*,double*)>)[exp_poly,dx](double *x, double *par){ 
        //std::cout << "<tf1>: dx = " << dx << "\n"; 
        return dx * exp_poly(x[0]);
    }; 

    auto tf1 = new TF1("exp_poly", fcn, xmin, xmax, 0); 

    tf1->Draw("SAME"); 

}

#endif