
#include "ToyEventGenerator.h"

#include <bininfo.hpp>
#include <fit_exponential_poly.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <newton_optimizer.hpp> 
#include <fit_parameter.hpp> 
#include <compute_Q0.hpp>
//ROOT headers
#include <TRandom3.h> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TAxis.h> 
#include <TGraph.h> 
#include <TLine.h> 
//stdlib headers
#include <cstdio> 
#include <functional> 
#include <iostream> 
#include <cmath> 


void scan_toy_spectrum(unsigned int n_events)
{   
    int n_steps=200; 
    double m_min{-5.5}, m_max{+5.5};

    double xmin{-7.}, xmax{+7.};

    double x0 = 2.; 
    double sigma = 0.05; 

    ToyEventGenerator generator;
    generator.Range(xmin, xmax);
    generator.Signal_sigma(sigma); 
    generator.Signal_center(x0);

    generator.Signal_fraction(8e-4);
    
    auto hist = new TH1D("h_combined", "Signal + background;X;", 200, xmin, xmax); 

    std::cout << 
        "generating events..." << std::flush;
    
    for (unsigned int i=0; i<n_events; i++) { hist->Fill(generator()); } 

    new TCanvas; 

    hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum()*1.1); 
    hist->Draw("E"); 


    auto data = peak_search::copy_1D_hist(hist); 

    auto exp_poly = (peak_search::ExponentialPoly)peak_search::fit_exponential_poly(hist, 8); 

    const double dx = (xmax - xmin)/((double)hist->GetXaxis()->GetNbins()); 

    auto fcn_poly = (std::function<double(double*,double*)>)[exp_poly,dx](double *x, double *par){ 
        //std::cout << "<tf1>: dx = " << dx << "\n"; 
        return dx * exp_poly(x[0]);
    };  
    
    auto tf1_poly = new TF1("tf1_poly", fcn_poly, xmin, xmax, 0); 

    tf1_poly->SetLineColor(kRed);
    tf1_poly->SetLineStyle(kDashed); 
    tf1_poly->Draw("SAME"); 


    auto poly_coeffs = exp_poly.GetParams(); 
    
    std::vector<fit_parameter_t> coeffs; coeffs.reserve(poly_coeffs.size() + 1); 
    
    coeffs.emplace_back(fit_parameter_t{ .val = 0., .name = "mu", .is_fixed = false });

    for (size_t i=0; i<poly_coeffs.size(); i++) {
        coeffs.emplace_back(fit_parameter_t{ .val = poly_coeffs[i], .name = Form("p_coeff_%zi",i), .is_fixed = false }); 
    }
    
    int i=0; 
    std::cout << "coeffs:\n"; 
    for (auto& coeff : exp_poly.GetParams()) {
        std::printf("   %3i : %f\n", i++, coeff); 
    }

    std::vector<double> pts_x0(n_steps,0.), pts_mu(n_steps,0.), pts_sqrtQ0(n_steps,0.), pts_pQ0(n_steps,0.);

    x0 = m_min;
    for (int i_step=0; i_step < n_steps; i_step++) {

        auto gaussian_fcn = peak_search::Gauss(0, x0, sigma); 

        auto fcn_s_plus_b = peak_search::FcnSum(&gaussian_fcn, &exp_poly); 

        double Q0 = peak_search::compute_Q0(data, fcn_s_plus_b); 
        double pQ0 = peak_search::compute_Q0_p(Q0);

        pts_x0[i_step]      = x0; 
        pts_mu[i_step]      = fcn_s_plus_b.GetParams()[0]; 
        pts_sqrtQ0[i_step]  = (Q0>0.?+1.:-1.) * std::sqrt( std::fabs(Q0) ); 
        pts_pQ0[i_step]     = std::log10(pQ0); 

        x0 += (m_max - m_min)/((double)n_steps-1); 

        std::printf(" x0=%.2f,  mu=%.1f, sqrt(Q0)=%.2f\n", x0, pts_mu[i_step], pts_sqrtQ0[i_step]); 
    }

    TGraph *graph; 
    TLine* line; 

    // 'mu' graph
    new TCanvas;
    graph = new TGraph(n_steps, pts_x0.data(), pts_mu.data()); 
    graph->SetTitle("Best-fit mu vs. x0;x0;best-fit mu"); 
    graph->Draw(); 

    line = new TLine(pts_x0.front(),0., pts_x0.back(),0.); 
    line->SetLineColor(kRed); 
    line->SetLineStyle(kDashed); 
    line->Draw(); 

    // 'Q0' graph
    new TCanvas; 
    graph = new TGraph(n_steps, pts_x0.data(), pts_sqrtQ0.data()); 
    graph->SetTitle("Signficance (Z = #sqrt{|Q0|}) vs. x0;x0;#sqrt{|Q0|}"); 
    graph->Draw(); 

    line = new TLine(pts_x0.front(),0., pts_x0.back(),0.); 
    line->SetLineColor(kRed); 
    line->SetLineStyle(kDashed); 
    line->Draw(); 

    // 'pQ0' graph
    new TCanvas; 
    graph = new TGraph(n_steps, pts_x0.data(), pts_pQ0.data()); 
    graph->SetTitle("log_{10} p(Q0) vs. x0;x0;log_{10} p(Q0)"); 
    graph->Draw(); 

    return; 
}
