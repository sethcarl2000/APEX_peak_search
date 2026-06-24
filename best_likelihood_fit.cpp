#include "best_likelihood_fit.hpp"
#include "log_likelihood.hpp"
#include "numbers.hpp"
#include "gauss_integrate.hpp"
#include "bininfo.hpp"

#include <TString.h> 
#include <TAxis.h> 
#include <Math/Minimizer.h>
#include <Math/Factory.h>
#include <Math/Functor.h>

#include <stdexcept>
#include <vector>
#include <memory> 

namespace peak_search
{

double best_likelihood_fit(const histo_1D_t& data, const std::function<double(double,const double*)>& fcn, std::vector<fit_parameter_t>& params)
{

    if (params.empty()) {
        throw std::logic_error(Form("in <%s>: 'params' vec (third arg.) is empty",__func__));
        return numbers::nan; 
    }
    auto minimizer = ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad");

    if (!minimizer) {
        throw std::logic_error(Form("in <%s>: minimizer failed to be initialized.",__func__)); 
        return numbers::nan; 
    }

    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(0.1);
    minimizer->SetPrintLevel(0);

    const double dx = (data.xmax - data.xmin)/((double)data.bins.size());

    //__________________________________________________________________________________________________
    //compute the negative log-likelihood (NLL)
    auto objective_fcn = [&fcn, &data, dx](const double* par){

        double NLL =0.;

        auto fcn_wrapper = [&fcn,par](double x){ return fcn(x,par); };
        for (const auto& bin : data.bins) {
            
            double x = bin.x; 
            
            double expect = gauss_integrate(fcn_wrapper, x-dx/2., x+dx/2.);

            NLL += expect - (bin.N * std::log(expect));  
        }
        return NLL; 
    };
    //__________________________________________________________________________________________________

    auto f_minimizer = ROOT::Math::Functor(objective_fcn, params.size());
    
    minimizer->SetFunction(f_minimizer);

    //set the list of variables
    for (int i_var=0; i_var<params.size(); i_var++) { 
        const auto& par = params[i_var];
        minimizer->SetVariable(i_var, Form("x_pow_%i",i_var), par.val, std::pow(1.25, i_var));
        if (par.is_fixed) { minimizer->FixVariable(i_var); }
    }
    
    bool fit_status = minimizer->Minimize();

    if (!fit_status) { 
        Error(__func__, "fit failed");
        return numbers::nan;
    }

    const double* params_result = minimizer->X(); 
    const double* errors = minimizer->Errors(); 

    //copy the result to our vector
    for (int i=0; i<params.size(); i++) params[i].val = params_result[i];
    
    return negative_log_likelihood(data, [params_result,&fcn](double x){ return fcn(x,params_result); });
}
    
};