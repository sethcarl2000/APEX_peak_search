
#include "ToyEventGenerator.h"

#include <bininfo.hpp>
#include <copy_subhist.hpp>
#include <fit_exponential_poly.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <newton_optimizer.hpp> 
#include <fit_parameter.hpp> 
#include <compute_Q0.hpp>
#include <numbers.hpp>
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

#define VERBOSE 1

void scan_toy_spectrum_window(unsigned int n_events, std::string output_gif="")
{   
    //if the user provided a path, then make a gif. 
    const bool make_gif=(output_gif != ""); 

    int n_steps=200; 
    double m_min{-5.5}, m_max{+5.5};

    double xmin{-7.}, xmax{+7.};

    double x0 = 2.; 
    double sigma = 0.05; 

    double window_size = sigma*10.; 

    ToyEventGenerator generator;
    generator.Range(xmin, xmax);
    generator.Signal_sigma(sigma); 
    generator.Signal_center(x0);

    generator.Signal_fraction(8e-4);
    
    auto hist = new TH1D("h_combined", "Signal + background;X;", 200, xmin, xmax); 

    std::cout << 
        "generating events..." << std::flush;
    
    for (unsigned int i=0; i<n_events; i++) { hist->Fill(generator()); } 

    auto c_fit = new TCanvas; 

    hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum()*1.1); 
    hist->Draw("E"); 

    const double dx = (xmax - xmin)/((double)hist->GetXaxis()->GetNbins()); 
    
    std::vector<double> pts_x0, pts_mu, pts_sqrtQ0, pts_pQ0;

    pts_x0.reserve(n_steps);
    pts_mu.reserve(n_steps);
    pts_sqrtQ0.reserve(n_steps);
    pts_pQ0.reserve(n_steps);

    auto gaussian_fcn   = peak_search::Gauss(0, x0, sigma); 

    x0 = m_min;
    gaussian_fcn.Set_mu(0.); 
    gaussian_fcn.Set_x0(0.);

    for (int i_step=0; i_step < n_steps; i_step++) {

        auto sub_data       = peak_search::copy_subhist(hist, x0-window_size, x0+window_size); 

        if (VERBOSE >= 2) {
            std::cout << "sub histogram xmax=" << x0-window_size << ", " << x0+window_size << "\n"; 
            std::cout << "bins:\n"; 
            for (auto& bin : sub_data.bins) {
                std::cout << "  " << bin.x << " | " << bin.N << "\n"; 
            }
        }

        auto exp_poly_result = peak_search::fit_exponential_poly(sub_data, 3);

        if (!exp_poly_result) {
            Warning(__func__, "Fit for x0=%.3f failed.", x0); 
            continue; 
        }

        auto exp_poly = exp_poly_result.data; 

        auto fcn_s_plus_b   = peak_search::FcnSum(&gaussian_fcn, &exp_poly); 

        double Q0 = peak_search::compute_Q0(sub_data, fcn_s_plus_b); 

        if (peak_search::numbers::is_nan(Q0)) {
            Warning(__func__, "computed Q0 is nan"); 
            continue; 
        }

        double pQ0 = peak_search::compute_Q0_p(Q0);

        if (make_gif) {
            auto fcn_b = (std::function<double(double*,double*)>)[exp_poly,dx](double *x, double *par){ 
                //std::cout << "<tf1>: dx = " << dx << "\n"; 
                return dx * exp_poly(x[0]);
            };  
            auto tf1_b = new TF1("tf1_poly", fcn_b, x0-window_size, x0+window_size, 0); 
            tf1_b->SetBit(kMustCleanup);
            tf1_b->ResetBit(kCanDelete); 

            tf1_b->SetLineColor(kBlue); 
            tf1_b->SetLineStyle(kDashed); 
            tf1_b->Draw("SAME"); 

            auto fcn_sb = (std::function<double(double*,double*)>)[fcn_s_plus_b,dx](double *x, double *par){ 
                //std::cout << "<tf1>: dx = " << dx << "\n"; 
                return dx * fcn_s_plus_b(x[0]);
            }; 
            auto tf1_sb = new TF1("tf1_poly", fcn_sb, x0-window_size, x0+window_size, 0); 
            tf1_sb->SetBit(kMustCleanup);
            tf1_sb->ResetBit(kCanDelete); 

            tf1_sb->SetLineColor(kRed); 
            tf1_sb->SetLineStyle(kSolid); 
            tf1_sb->Draw("SAME"); 

            c_fit->Modified();
            c_fit->Update(); 
            c_fit->SaveAs(std::string{output_gif + "+20"}.c_str()); 

            delete tf1_b; 
            delete tf1_sb; 
        }

        pts_x0      .emplace_back( x0 ); 
        pts_mu      .emplace_back( fcn_s_plus_b.GetParams()[0] ); 
        pts_sqrtQ0  .emplace_back( (Q0>0.?+1.:-1.) * std::sqrt( std::fabs(Q0) ) ); 
        pts_pQ0     .emplace_back( std::log10(pQ0) ); 

        std::printf(" x0=%.2f,  mu=%.1f, sqrt(Q0)=%.2f\n", x0, pts_mu[i_step], pts_sqrtQ0[i_step]); 

        x0 += (m_max - m_min)/((double)n_steps-1); 
        gaussian_fcn.Set_x0(x0);
        gaussian_fcn.Set_mu(0.);


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
