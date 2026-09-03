
#include <solve_for_CLs.hpp>
#include <compute_Q0.hpp>
#include "numbers.hpp"
#include "fit_parameter.hpp"
#include <newton_optimizer.hpp>
#include <make_histogram_copy.hpp>
// ROOT
#include <TError.h>
#include <TString.h> 
#include <Math/ProbFuncMathCore.h>
// stdlib
#include <cmath> 
#include <cstdio> 


namespace peak_search
{

//________________________________________________________________________________________________________________________________________
double solve_for_CLs(const Histo1D& data, Fcn1D& fcn, double CL, double rel_tolerance, double grid_search_size)
{
    //first, compute the best fit
    auto& coeffs = fcn.GetParams(); 
    std::vector<fit_parameter_t> fit_params; fit_params.reserve(coeffs.size());

    fit_params.push_back({ .val=coeffs[0], .name="mu", .is_fixed=false});
    for (size_t i=1; i<fcn.GetDoF(); i++) {
        fit_params.push_back({ .val=coeffs[i], .name=Form("param_%zi",i), .is_fixed=false});
    }

#ifdef DEBUG    
    double total_stats =0.; 
    for (const auto& bin : data.bins) total_stats += bin.N; 
    std::printf("<%s> in body. ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
        "   params:\n"
        "       CL: %.5f\n"
        "       rel_tolerance: %.3e\n"
        "       grid_search_size: %.3e\n"
        "   data: \n"
        "       bins: %zi\n"
        "       x-range: [%.1f, %.1f]\n"
        "       total stats: %.3e\n"
        "   fcn: \n"
        "       n. params: %zi (mu + %zi nuissance params)\n",
    __func__,
        
        CL, 
        rel_tolerance,
        grid_search_size,

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
        return numbers::nan; 
    }
    const double mu_mle = param_mu.val; 

    //now, fix mu=0, and find the NLL
    param_mu.val =0.; 
    param_mu.is_fixed = true; 
    double NLL_0      = newton_optimizer(data, fcn, fit_params);

    if (numbers::is_nan(NLL_0)) {
        Error(__func__, "NLL for mu=0 test case is null (newton-optimizer step failed)."); 
        return numbers::nan; 
    }
    
    // the background-only confidence level. 
    double Q0_b = NLL_0 - NLL_mle; 
    double CL_b = 1. - compute_Q0_p( (mu_mle < 0.) ? -2.*Q0_b : +2.*Q0_b ); 

    double NLL_lo, NLL_hi; 

    // the CL_{s+b} (signal + background) confidence level. this is different for each different 'mu'. 
    double CLs_lo, CLs_hi; 

    // the value of 'mu' which satisfies the given CL
    // if mu_MLE is biggen than 0, that's our starting place, otherwise, we'll start at 0. 
    double mu_CL;
    if (mu_mle > 0.) { 
        NLL_lo = NLL_mle;  mu_CL = mu_mle; 
    } else {
        NLL_lo = NLL_0;    mu_CL = 0; 
    }


    param_mu.is_fixed = true; 

    //return CLs for given value of 'mu' 
    auto get_CLs = [&data, &fcn, &fit_params, &param_mu, CL_b, NLL_mle, mu_mle](double mu) {

        param_mu.val = mu;
        double NLL = newton_optimizer(data, fcn, fit_params);         

        return (1. - compute_Q0_p(2.*(NLL - NLL_mle)))/CL_b;  
    };

    // first, do a basic grid-search until we pass the CL 
    double mu_low{mu_CL}, mu_high; 

    CLs_lo = get_CLs(mu_low);

#ifdef DEBUG
    std::printf("<%s>: preliminary fits:\n"
        "   NLL:\n"
        "       mle:        %+.8e\n"
        "       mu=0:       %+.8e\n"
        "   mu_mle:         %+.5e\n"
        "   mu_CL (start)   %+.5e\n"
        "\n"
        "   CL_sb(mu_CL):   %.5f\n"
        "   CL_b:           %.5f\n",
        __func__, 
        NLL_mle, 
        NLL_0, 
        mu_mle, 
        mu_CL,
        CLs_lo,
        CL_b
    );
#endif

    //first, keep increasing our upper bound until we find a value of 'mu_CL' that's above our CL 
    while (mu_CL < grid_search_size*1e6) {
        
        //set 'mu' to the given value, and find the MLE for all other parameters. 
        mu_high = mu_low + grid_search_size; 

        CLs_hi = get_CLs(mu_high); 

        if (CLs_hi > CL) {
#ifdef DEBUG
            std::printf("   found upper limit:\n"
                "   mu_CL bounds:   [%+.5e, %+.5e]\n"
                "   CLs bounds:     [%+.5e, %+.5e]\n",
                mu_low, mu_high, 
                CLs_lo, CLs_hi
            );
#endif
            break; 
        }
    
        //if we haven't yet found an upper bound, go one step up. 
        mu_low = mu_high; 
        CLs_lo = CLs_hi; 
    }
#ifdef DEBUG
    int n_iterations =0; 
    std::printf("starting bisection iterations...\n"); 
#endif

    double rel_error_bound = 2.*(mu_high - mu_low)/(mu_high + mu_low); 

    while (2.*(mu_high - mu_low)/(mu_high + mu_low) > rel_tolerance) {

#ifdef DEBUG
        std::printf("bisection it: mu [%+.5e, %+.5e]      CL: [%+.5e, %+.5e]    rel. error: %.4e\n",
            mu_low, mu_high, CLs_lo, CLs_hi, rel_error_bound
        );
#endif 
        //do bisection alg. 
        double mu = (mu_low + mu_high)/2.; 

        double residual = get_CLs(mu) - CL; 

        if (residual > 0.) { 
            mu_high = mu; 
            CLs_hi = residual + CL; 
        } else { 
            mu_low = mu; 
            CLs_lo = residual + CL; 
        }

        rel_error_bound = 2.*(mu_high - mu_low)/(mu_high + mu_low); 
    }

    return (mu_low + mu_high)/2.; 
}

};