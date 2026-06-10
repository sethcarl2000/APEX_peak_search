
#include "../chisquare.hpp"
#include "gauss_integrate.hpp"

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

namespace 
{
    constexpr int n_bins = 50;
    constexpr double range = 5.; 

    constexpr double sigma = 2.5; 

    double generate_Q0(TRandom3& rand) {
        double ret; do { ret = rand.Gaus()*sigma; } while (std::fabs(ret) > range); 
        return ret;  
    }

    double Q0_pdf(double Q0) {
        double xx{ Q0/sigma }; 
        return std::exp( - xx*xx/2. ) / (2.50662827463 * sigma); 
    }

};

/// generate a toy Q0 distribution, and measure its negative log-likelihood, to see how its distributed. 
int generate_gauss_chi2()
{
    const long int n_pts = 4e4; 

    const int n_bins = 200; 

    int experiments_per_step = 1e5; 
    
    auto hist_exp = new TH1D("h_exp", Form("#chi^2 dist for Q0;#chi^{2}/2 (DoF=%i)",n_bins), 
        200, 0, 250
    );

    auto hist_p = new TH1D("h_p", Form("p-value #chi^2 dist for Q0;p(#chi^{2}/2) (DoF=%i)",n_bins), 
        200, -0.5, 1.5
    );

    TRandom3 rand; 
    
    for (int i_step=0; i_step<experiments_per_step; i_step++) {

        //do an experiment    
        TH1D hist_1d("hist_exp", "", n_bins, -range, +range); 
            
        //double dx_chi = (hist_1d.GetXaxis()->GetXmax() - hist_1d.GetXaxis()->GetXmin()) / ((double)hist_1d.GetXaxis()->GetNbins());

        for (long i=0; i<n_pts; i++) { hist_1d.Fill(generate_Q0(rand)); }

        double n_pts_hist = hist_1d.Integral() / (ROOT::Math::normal_cdf(+range, sigma) - ROOT::Math::normal_cdf(-range, sigma));
        auto fcn_expect = [n_pts_hist](double Q0){ return Q0_pdf(Q0)*n_pts_hist; };

        if (i_step==0) {
            new TCanvas; 
            hist_1d.DrawCopy();
            auto xax = hist_1d.GetXaxis(); 
            double dx_chi = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()); 

            auto fcn_dist = (std::function<double(double*,double*)>)[n_pts, dx_chi, fcn_expect](double *x, double *par){
                return dx_chi*fcn_expect(x[0]);
            };

            auto tf11 = new TF1("hhh", fcn_dist, -range, +range, 0); 
            tf11->DrawCopy("SAME"); 
        }

        double chi2 = peak_search::chisquare(&hist_1d, fcn_expect);

        double p = peak_search::chisquare_p(chi2, hist_1d.GetXaxis()->GetNbins()); 

        hist_p->Fill( p );
        hist_exp->Fill( chi2/2. );
    }

    new TCanvas; 
    hist_exp->Draw("E"); 

    auto xax = hist_exp->GetXaxis(); 
    double dx = (xax->GetXmax() - xax->GetXmin())/(xax->GetNbins()); 

    auto fcn_draw = (std::function<double(double*,double*)>)[dx, n_bins, experiments_per_step](double *x, double *par) {
        //return ROOT::Math::gamma_pdf(x[0], ((double)n_bins)/2., 1.) * dx * (double)experiments_per_step; 
        double p_val = std::exp( -x[0] ) * std::pow( x[0], ((double)n_bins)/2. - 1.) / ROOT::Math::tgamma( ((double)n_bins)/2. );
        return p_val * dx * (double)experiments_per_step; 
    };

    auto tf1 = new TF1("hist_f1", fcn_draw, xax->GetXmin(), xax->GetXmax(), 0);
    tf1->Draw("SAME"); 


    new TCanvas; 
    hist_p->Draw(); 

    return 0; 
}