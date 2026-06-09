#include "best_likelihood_fit.hpp"
#include "log_likelihood.hpp"
#include "numbers.hpp"
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

double best_likelihood_fit(TH1D* hist, double (*fcn)(double,const double*), std::vector<double>& params)
{

    if (!hist) {
        throw std::logic_error(Form("in <%s>: hist passed is null",__func__));
        return numbers::nan; 
    }

    if (params.empty()) {
        throw std::logic_error(Form("in <%s>: 'params' vec (third arg.) is empty",__func__));
        return numbers::nan; 
    }

    //get the histogram axis
    auto xax = hist->GetXaxis(); 

    //make a copy of all bin contents
    auto bins = copy_1D_hist(hist);

    auto minimizer = std::unique_ptr<ROOT::Math::Minimizer>(
        ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad")
    ); 

    if (!minimizer) {
        throw std::logic_error(Form("in <%s>: minimizer failed to be initialized.",__func__)); 
        return numbers::nan; 
    }

    minimizer->SetMaxFunctionCalls(1e7); 
    minimizer->SetMaxIterations(1e6); 
    minimizer->SetTolerance(1e-3);
    minimizer->SetPrintLevel(2);

    auto negative_log_likelihood = [fcn, &bins](const double* par){
        
        return -1.* log_likelihood(bins, fcn, par); 
    };

    auto f_minimizer = ROOT::Math::Functor(negative_log_likelihood, params.size());
    
    minimizer->SetFunction(f_minimizer);

    //set the list of variables
    int i_var=0;
    for (int i_var=0; i_var<params.size(); i_var++) { 
        minimizer->SetVariable(i_var, Form("param_%i",i_var), params[i_var], 1e-4);
    }

    bool fit_status = minimizer->Minimize();

    if (!fit_status) return numbers::nan;

    const double* params_result = minimizer->X(); 
    const double* errors = minimizer->Errors(); 

    //copy the result to our vector
    for (int i=0; i<params.size(); i++) params[i] = params_result[i];
    
    return negative_log_likelihood(const_cast<const double*>(params.data()));
}
    
};