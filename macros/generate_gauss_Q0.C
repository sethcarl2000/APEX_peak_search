
#include "../chisquare.hpp"
#include "gauss_integrate.hpp"
#include "compute_Q0.hpp"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h> 
#include <TRandom3.h> 
#include <Math/ProbFunc.h>
#include <Math/QuantFuncMathCore.h>
#include <Math/SpecFuncMathCore.h>
#include <Math/PdfFuncMathCore.h>
#include <TF1.h> 
#include <TAxis.h> 

#include <cmath> 
#include <cstdio> 
#include <functional>
#include <vector>  

namespace 
{
    constexpr int n_bins = 100;
    constexpr double range = 5.; 

    constexpr double sigma_background = 2.5; 

    constexpr double signal_sigma = 0.20; 
    constexpr double signal_x0 = -4.7; 

    double generate_background(TRandom3& rand) {
        double ret; do { ret = rand.Gaus()*sigma_background; } while (std::fabs(ret) > range); 
        return ret;  
    }

    double gauss_pdf(double Q0, double sigma) {
        double xx{ Q0/sigma }; 
        return std::exp( - xx*xx/2. ) / (2.50662827463 * sigma); 
    }

};

/// generate a toy Q0 distribution, and measure its negative log-likelihood, to see how its distributed. 
int generate_gauss_Q0(const int n_experiments)
{
    const long int n_pts = 4e4; 

    const int n_bins = 200; 
    
    auto hist_Q0 = new TH1D("h_exp", "Q0 dist for gaussian;Q0;", 
        60, -5., 5.
    );

    auto hist_p = new TH1D("h_p", "p-value Q0 dist for gaussian;p(Q0);", 
        60, -0.25, 1.25
    );

    TRandom3 rand; 
    
    for (int i_step=0; i_step<n_experiments; i_step++) {

        //do an experiment    
        TH1D hist_1d("hist_exp", "", n_bins, -range, +range); 
            
        //double dx_chi = (hist_1d.GetXaxis()->GetXmax() - hist_1d.GetXaxis()->GetXmin()) / ((double)hist_1d.GetXaxis()->GetNbins());

        for (long i=0; i<n_pts; i++) { hist_1d.Fill(generate_background(rand)); }

        double n_pts_hist = hist_1d.Integral() / (ROOT::Math::normal_cdf(+range, sigma_background) - ROOT::Math::normal_cdf(-range, sigma_background));
        auto fcn_expect = [n_pts_hist](double Q0){ return gauss_pdf(Q0,sigma_background)*n_pts_hist; };

        if (i_step==0) {
            new TCanvas; 
            hist_1d.DrawCopy();
            auto xax = hist_1d.GetXaxis(); 
            double dx_chi = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()); 

            auto fcn_dist = (std::function<double(double*,double*)>)[dx_chi, fcn_expect](double *x, double *par){
                return dx_chi*fcn_expect(x[0]);
            };

            auto tf11 = new TF1("hhh", fcn_dist, -range, +range, 0); 
            tf11->DrawCopy("SAME"); 
        }

        auto fcn_expect_s_plus_b = (std::function<double(double,double,const double*)>)[&fcn_expect](double x, double mu, const double* par){
            return 
                fcn_expect(x)*par[0] +
                gauss_pdf(x-signal_x0, signal_sigma)*mu; 
        };
        std::vector<double> nuissance_params{ 1. };

        double Q0 = peak_search::compute_Q0(&hist_1d, fcn_expect_s_plus_b, nuissance_params);

        double p = peak_search::compute_Q0_p(Q0); 

        hist_p->Fill( p );
        hist_Q0->Fill( Q0 );
    }

    new TCanvas; 
    hist_Q0->Draw("E"); 

    auto xax = hist_Q0->GetXaxis(); 
    double dx = (xax->GetXmax() - xax->GetXmin())/(xax->GetNbins()); 

    auto fcn_draw = (std::function<double(double*,double*)>)[dx, n_bins, n_experiments](double *x, double *par) {
        //return ROOT::Math::gamma_pdf(x[0], ((double)n_bins)/2., 1.) * dx * (double)experiments_per_step; 
        double absQ0 = std::fabs(x[0]);
        double p_val = std::exp( -absQ0/2. ) / (2. * 2.50662827463 * std::sqrt(absQ0));
        return p_val * dx * (double)n_experiments; 
    };

    auto tf1 = new TF1("hist_f1", fcn_draw, xax->GetXmin(), xax->GetXmax(), 0);
    tf1->Draw("SAME"); 

    new TCanvas; 
    hist_p->Draw("HIST, E"); 

    return 0; 
}