
#include <compute_statistics.hpp>
#include <gauss_integrate.hpp> 
#include "numbers.hpp"
#include "fit_parameter.hpp"
#include <newton_optimizer.hpp>
#include <make_histogram_copy.hpp>
#include <Fcn1D/FcnSum.hpp>
//Eigen
#include <eigen3/Eigen/Core> 
#include <eigen3/Eigen/Dense>
// ROOT
#include <TError.h>
#include <TString.h> 
#include <Math/ProbFuncMathCore.h>
#include <Math/QuantFuncMathCore.h>
// stdlib
#include <cmath> 
#include <cstdio> 

//#define DEBUG
namespace peak_search
{

#ifdef DEBUG
#define DEBUG_STATS
#endif 

constexpr double muon_mass = 105.66; // MeV; 
constexpr double electron_mass = 0.501; // MeV; 

/// @return the approximate branch ratio for A'-> e+e- (as opposed to A'->mu+mu-)
double get_branch_ratio(double mass)
{
    if (mass < 2.*muon_mass) return 1.; 

    double amp_electron = std::sqrt( 1. - 4.*(electron_mass*electron_mass)/(mass*mass) ) ; 
    double amp_muon     = std::sqrt( 1. - 4.*(muon_mass*muon_mass)/(mass*mass) ) ; 
    
    return amp_electron / (amp_electron + amp_muon); 
}
 
/// @return a _very_ rough estimate of the fraction of our background events that are the 'radiative' type of gamma->e+e- production. 
double estimate_radiative_fraction(double mass)
{
    return 0.205 + (mass - 120)* ((0.155 - 0.205)/(220 - 120)); 
}


FitResult<FitStats> compute_statistics(const Histo1D& data, Fcn1D& fcn_s, Fcn1D& fcn_b, double mass, double CLs, double mass_window_size)
{
    using Eigen::MatrixXd, Eigen::VectorXd;

    FcnSum fcn(&fcn_s, &fcn_b); 

    //first, compute the best fit
    auto& coeffs = fcn.GetParams(); 
    std::vector<fit_parameter_t> fit_params; fit_params.reserve(coeffs.size());

    fit_params.push_back({ .val=coeffs[0], .name="mu", .is_fixed=false});
    for (size_t i=1; i<fcn.GetDoF(); i++) {
        fit_params.push_back({ .val=coeffs[i], .name=Form("param_%zi",i), .is_fixed=false});
    }

#ifdef DEBUG_STATS    
    double total_stats =0.; 
    for (const auto& bin : data.bins) total_stats += bin.N; 
    std::printf("<%s> in body. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "   params:\n"
        "       CL: %.5f\n"
        "       reference mass-window size: %.5f\n"
        "   data: \n"
        "       bins: %zi\n"
        "       x-range: [%.1f, %.1f]\n"
        "       total stats: %.3e\n"
        "   fcn: \n"
        "       n. params: %zi (mu + %zi nuissance params)\n",
    __func__,
        
        CLs, 
        mass_window_size,
        
        data.bins.size(), 
        data.bins.front().xmin, data.bins.back().xmax, 
        total_stats, 

        coeffs.size(), coeffs.size()-1
    );
#endif

    auto& param_mu = fit_params[0]; 
    
    // now, get the value of mu which is the Maximum-likelihood-estimate (MLE)
    param_mu.is_fixed = false; 
    double NLL_mle   = newton_optimizer(data, fcn, fit_params); 

    if (numbers::is_nan(NLL_mle)) {
        Error(__func__, "NLL for mu=best-fit test case is nan (newton-optimizer step failed)."); 
        return FitResult<FitStats>::Fail(); 
    }
    const double mu_mle = param_mu.val; 

    // let's estimate the variance of each parameter around its MLE value. 
    // this formula comes from: https://arxiv.org/pdf/1007.172
    // 'Asymptotic formulae for profile likelihood tests of new physics' Glen Cowan, Kyle Cranmer, Eilam Gross, Ofer Vitells.
    //
    // in particular, see section 3 and equation (28) therein. 
    // 
    
    // this will be the inverse of the covariance matrix for each of our function's paramters (centered around their best-fit values)
    const size_t n_params = fit_params.size(); 
    
    // call this V^{-1}_{ij}
    MatrixXd V_inv = MatrixXd::Zero(n_params, n_params); 

    std::vector<double> dExpect_di(n_params,0.); 

    for (const auto& bin : data.bins) {

        double expect = gauss_integrate(fcn, bin.xmin, bin.xmax); 

        for (size_t i=0; i<n_params; i++) 
            dExpect_di[i] = gauss_integrate([i,&fcn](double x){ return fcn.Di(x,i); }, bin.xmin, bin.xmax); 

        for (size_t i=0; i<n_params; i++) {
            for (size_t j=i; j<n_params; j++) { // <== Note that j starts counting at 'i'. this is symmetric matrix, so we only need to compute the upper-triangle. 

                V_inv(i,j) += dExpect_di[i] * dExpect_di[j] / expect; 
            }
        }
    }
    
    //now, we can fill out the rest of the lower-triangle by copying the symmetric elements from the upper-triangle.  
    for (int i=1; i<n_params; i++) 
        for (int j=0; j<i; j++)
            V_inv(i,j) = V_inv(j,i); 

    //and now we use a sneaky trick to extrat sigma^2 = V_{00}. 
    VectorXd unit_vec = VectorXd::Zero(n_params); 
    unit_vec(0) = 1.; 

    //when we solve the lin. equation w/r/t the above unit vector, and take the first-column of the resultant vector, we have the element V_{00}
    VectorXd X = V_inv.llt().solve(unit_vec); 
    double sigma = std::sqrt( X[0] );  


    //now, fix mu=0, and find the NLL
    param_mu.val =0.; 
    param_mu.is_fixed = true; 
    double NLL_0      = newton_optimizer(data, fcn, fit_params);

    if (numbers::is_nan(NLL_0)) {
        Error(__func__, "NLL for mu=0 test case is nan (newton-optimizer step failed)."); 
        return FitResult<FitStats>::Fail(); 
    }

    double mu_mle_sign = (mu_mle < 0.) ? -1 : +1; 

    // now, we can compute Q0. 
    double Z_b = mu_mle_sign * std::sqrt(2.*(NLL_0 - NLL_mle)); 

    double Q0  = mu_mle_sign * (Z_b*Z_b); 

    // next, we will estimate the upper-limit on the signal-strength parameter, mu_CL, using the CLs technique. 
    // for an outline of this technique, see:  
    //
    //      Presentation of search results: The CL(s) technique
    //      Alexander L. Read (Oslo Univ.), 2002
    //
    // this represents: 
    //
    //      p( mu_mle_obs < mu_mle ; mu_true = 0 ) 
    //
    //  in words: if the true signal rate (mu_true) is 0, then we would expect any other identical measurement to find a vale of best-fit 'mu' less than or equal to mu_mle with probability: p
    //
    using ROOT::Math::normal_cdf, ROOT::Math::normal_quantile; 

    double CL_b = normal_cdf(Z_b, 1.); 

    double mu_CL = mu_mle   -   sigma*normal_quantile(CLs * CL_b, 1.); 

    // get the (expected) background rate for this fcn.
    // so we set the signal rate to 0: 
    fcn.GetParams()[0] = 0.; 

    double expect_b = gauss_integrate(fcn_b, mass - mass_window_size/2., mass + mass_window_size/2.);
    
    fcn_s.GetParams()[0] = mu_CL; 
    double expect_s = gauss_integrate(fcn_s, mass - mass_window_size/2., mass + mass_window_size/2.);
    
    //now, use mu_CL to compute the upper-limit on the coupling, epsilon^2 
    double epsilon2_CL  = (expect_s / expect_b) * (mass_window_size/mass) * (2./137.) / (3.*3.1415926536); 

    //now, apply some corrections. 
    epsilon2_CL = epsilon2_CL / ( get_branch_ratio(mass) * estimate_radiative_fraction(mass) ); 

    
    FitResult<FitStats> result;
    
    result.status = Status::kSuccess; 
    result.data   = FitStats{ .mass=mass, .mu_MLE=mu_mle, .Q0=Q0, .mu_CL=mu_CL, .epsilon2_CL=epsilon2_CL, .mu_sigma=sigma };

    return result; 
}

};