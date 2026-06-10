
#include "log_likelihood.hpp"

#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h> 
#include <TRandom3.h> 
#include <Math/ProbFunc.h>
#include <Math/QuantFuncMathCore.h>

#include <cmath> 
#include <cstdio> 

namespace 
{
    constexpr int n_bins = 50;
    constexpr double range = 5.; 

    double generate_Q0(TRandom3& rand)
    {
        double p = rand.Rndm(); 
        
        double q0 = ROOT::Math::normal_quantile(p, 1.);
        
        return (p < 0.5 ? -1 : +1)*q0*q0; 
    }

    double Q0_pdf(double Q0) {
        Q0 = std::fabs(Q0);
        return std::exp( -Q0/2. ) / ( 2. * 2.50662827463 * std::sqrt(Q0) ); 
    }

};

/// generate a toy Q0 distribution, and measure its negative log-likelihood, to see how its distributed. 
int generate_toy_Q0_dist()
{
    const double min_pts = 10; 
    const double max_pts = 1e6; 
    int steps = 100; 

    double pts_log_0 = std::log10(min_pts); 
    double pts_log_1 = std::log10(max_pts);

    double dx = (pts_log_1 - pts_log_0)/((double)steps);

    int experiments_per_step = 100; 
    
    auto hist_exp = new TH2D("h_exp", "specific NLL vs N.samples;log_{10}(N.pts); log_{10}(NLL)", 
        steps, pts_log_0-dx/2., pts_log_1+dx/2.,
        200, -10, +10
    );

    TRandom3 rand; 

    double pts_log = pts_log_0 + dx/2.;
    for (int i_step=0; i_step<steps; i_step++) {
        //do an experiment
        double n_pts_d = std::exp(pts_log*std::log(10.));

        long n_pts = (long)n_pts_d;
        
        std::printf("pts: %li\n", n_pts);

        auto fcn_expect = [n_pts_d](double Q0){ return Q0_pdf(Q0)*n_pts_d; }; 

        for (int i_exp=0; i_exp<experiments_per_step; i_exp++) {

            TH1D hist_1d("hist_exp", "", n_bins, -range, +range); 
            
            for (long i=0; i<n_pts; i++) { hist_1d.Fill(generate_Q0(rand)); }
            
            double eta = -peak_search::log_likelihood(&hist_1d, fcn_expect);// / n_pts_d;

            hist_exp->Fill(std::log10(n_pts_d), std::log10(eta)); 
        }

        pts_log += dx; 
    }

    new TCanvas; 
    hist_exp->Draw("col"); 

    return 0; 
}