#ifndef generate_signal_bg_C
#define generate_signal_bg_C

#include "ToyEventGenerator.h"

#include <bininfo.hpp>
#include <fit_exponential_poly.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <newton_optimizer.hpp> 
#include <fit_parameter.hpp> 
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


void fit_toy_peak(unsigned int n_events)
{   
    double xmin{-7.}, xmax{+7.};

    double x0 = 2.; 
    double sigma = 0.05; 

    ToyEventGenerator generator;
    generator.Range(xmin, xmax);
    generator.Signal_sigma(sigma); 
    generator.Signal_center(x0);

    generator.Signal_fraction(1e-3);
    
    auto hist = new TH1D("h_combined", "Signal + background;X;", 100, xmin, xmax); 

    std::cout << 
        "generating events..." << std::flush;
    
    for (unsigned int i=0; i<n_events; i++) { hist->Fill(generator()); }


    
    std::cout << 
        "done.\n"; 


    new TCanvas; 
    hist->Draw("E"); 

    hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum()*1.1); 

    auto data = peak_search::copy_1D_hist(hist); 

    auto exp_poly = (peak_search::ExponentialPoly)peak_search::fit_exponential_poly(hist, 6); 

    const double dx = (xmax - xmin)/((double)hist->GetXaxis()->GetNbins()); 

    auto fcn_poly = (std::function<double(double*,double*)>)[exp_poly,dx](double *x, double *par){ 
        //std::cout << "<tf1>: dx = " << dx << "\n"; 
        return dx * exp_poly(x[0]);
    }; 

    auto tf1 = new TF1("exp_poly", fcn_poly, xmin, xmax, 0); 

    int i=0; 
    std::cout << "coeffs:\n"; 
    for (auto& coeff : exp_poly.GetParams()) {
        std::printf("   %3i : %f\n", i++, coeff); 
    }

    auto gaussian_fcn = peak_search::Gauss(0, x0, sigma); 


    std::cout << "exp poly dof: " << exp_poly.GetDoF() << " size: " << exp_poly.GetParams().size() << "\n"; 

    auto fcn_s_plus_b = peak_search::FcnSum(&gaussian_fcn, &exp_poly); 

    tf1->SetLineStyle(kDashed); 
    tf1->Draw("SAME"); 


    auto coeffs_vec = fcn_s_plus_b.GetParams(); 
    std::vector<fit_parameter_t> coeffs; coeffs.reserve(coeffs_vec.size()); 
    
    coeffs.emplace_back(fit_parameter_t{ .val = coeffs_vec[0], .name = "mu", .is_fixed = false });
    for (size_t i=1; i<coeffs_vec.size(); i++) {
        coeffs.emplace_back(fit_parameter_t{ .val = coeffs_vec[i], .name = Form("c%zi",i), .is_fixed = false }); 
    }
    peak_search::newton_optimizer(data, fcn_s_plus_b, coeffs); 


    
    auto fcn = (std::function<double(double*,double*)>)[exp_poly, gaussian_fcn, dx](double *x, double *par){ 
        //std::cout << "<tf1>: dx = " << dx << "\n"; 
        return dx * ( exp_poly(x[0]) + gaussian_fcn(x[0]) );
    }; 
    auto tf1_s_b = new TF1("f_s_b", fcn, xmin, xmax, 0); 

    tf1_s_b->SetLineColor(kBlue); 
    tf1_s_b->SetLineStyle(kSolid); 

    tf1_s_b->Draw("SAME"); 

    return; 

}

#endif