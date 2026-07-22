
#include "compute_Q0.hpp"
#include "numbers.hpp"
#include "fit_parameter.hpp"
#include "bininfo.hpp"
#include <newton_optimizer.hpp>

// ROOT
#include <TError.h>
#include <TString.h> 
#include <Math/ProbFuncMathCore.h>
// stdlib
#include <cmath> 


namespace peak_search
{

//________________________________________________________________________________________________________________________________________
double compute_Q0(TH1D* hist, Fcn1D& fcn, std::vector<double>& nuissance_params)
{
    //extract the histogram values
    if (!hist) {
        Error(__func__, "hist passed is null");
        return numbers::nan;
    }

    auto data = copy_1D_hist(hist);
    
    return compute_Q0(data, fcn, nuissance_params);
}
//________________________________________________________________________________________________________________________________________
double compute_Q0(histo_1D_t data, Fcn1D& fcn, std::vector<double>& nuissance_params)
{
    //check to make sure that the fcn and params match
    if (fcn.GetDoF() != (int)nuissance_params.size()+1) {
        std::ostringstream oss; 
        oss << "in <"<<__func__<<">: number of fcn parameters ("<<fcn.GetDoF()<<") "
            "does not match size of parameter list (1+"<<nuissance_params.size()<<")"; 
        
        throw std::logic_error(oss.str()); 
        return numbers::nan; 
    }

    //first, compute the best fit
    std::vector<fit_parameter_t> fit_params; fit_params.reserve(nuissance_params.size() + 1);

    fit_params.push_back({ .val=0., .name="mu" }); 

    auto& param_mu = fit_params[0]; 

    for (size_t i=0; i<nuissance_params.size(); i++) {
        fit_params.push_back({ .val=nuissance_params[i], .name=Form("param_%zi",i), .is_fixed=false});
    }
    
    //now, fix mu=0, and find the NLL
    param_mu.val =0.; 
    param_mu.is_fixed = true; 
    double NLL_0      = newton_optimizer(data, fcn, fit_params);
    
    param_mu.is_fixed = false; 
    double NLL_best   = newton_optimizer(data, fcn, fit_params); 

    //return the computed value of Q0 
    double Q0 = NLL_0 - NLL_best; 
    if (param_mu.val < 0) { Q0 *= -2.; } else { Q0 *= +2.; }
    return Q0; 
}
//________________________________________________________________________________________________________________________________________
double compute_Q0_p(double Q0)
{
    double p = ROOT::Math::normal_cdf_c( std::sqrt(std::fabs(Q0)), 1. );
    return (Q0 > 0) ? p : 1 - p; 
}

};