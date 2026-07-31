#include "fit_exponential_legendre.hpp"
#include "best_likelihood_fit.hpp"
#include "log_likelihood.hpp"
#include "fit_parameter.hpp"
#include "bininfo.hpp"
#include <newton_optimizer.hpp>
#include <legendre_polynomial.hpp>

//ROOT header
#include <TAxis.h> 
#include <TError.h> 
#include <TString.h> 
//eigen header
#include <eigen3/Eigen/Core> 
#include <eigen3/Eigen/Dense>
//stdlib headers
#include <vector>
#include <iostream>  
#include <cstdio> 

namespace {
    struct binval_t { double x, N; }; 
};

namespace peak_search 
{
//______________________________________________________________________________________________________________________________
FitResult<ExponentialLegendre> fit_exponential_legendre(TH1D* hist, int degree)
{
    using Eigen::MatrixXd, Eigen::VectorXd; 

    //extract the histogram values
    if (!hist) {
        Error(__func__, "hist passed is null");
        return FitResult<ExponentialLegendre>::Null(); 
    }

    auto data = copy_1D_hist(hist);

    return fit_exponential_legendre(data, degree);
}
//______________________________________________________________________________________________________________________________
FitResult<ExponentialLegendre> fit_exponential_legendre(histo_1D_t data, int degree)
{
    using Eigen::MatrixXd, Eigen::VectorXd; 
    
    double dx = (data.xmax - data.xmin)/((double)data.bins.size()); 

    double x_scale  = (data.xmax - data.xmin)/2.;
    double x_center = (data.xmax + data.xmin)/2.;  
    
    //first, try to fit the polynomial using a chi-square fit. 
    MatrixXd A = MatrixXd::Zero(degree, degree);
    VectorXd B = VectorXd::Zero(degree);

    for (auto& bin : data.bins) {

        //don't want to take the logarithm of 0! 
        if (bin.N <= 1e-6) continue; 

        double Xmu[degree]; 

        double x = (bin.x - x_center)/x_scale; 

         
        for (int i=0; i<degree; i++) Xmu[i] = legendre_polynomial(x, i); 

        for (int i=0; i<degree; i++) {

            B(i) += std::log( bin.N/dx ) * Xmu[i];

            for (int j=0; j<degree; j++) {

                A(i,j) += Xmu[i]*Xmu[j];
            }
        }
    }
    
    //now, solve this system.
    //we call these coeffs 'scaled', cause the x-values were normalized to be on the interval x=[-1,+1] 
    VectorXd coeffs_vec = A.llt().solve(B); 
    
    std::vector<fit_parameter_t> coeffs; coeffs.reserve(coeffs_vec.size()); 
    
    for (size_t i=0; i<coeffs_vec.size(); i++) {
        coeffs.emplace_back(fit_parameter_t{ .val = coeffs_vec[i], .name = Form("c%zi",i), .is_fixed = false }); 
    }
    
    //now, let's 
    ExponentialLegendre poly({}, data.xmin, data.xmax); 
    poly.SetParams(coeffs_vec); 

    //check for NaN
    if (coeffs.size() != (size_t)degree) { return FitResult<ExponentialLegendre>::Fail(); }

    poly.SetParams(coeffs); 

    //now, de-scale each exponent. (! update, for ExponentialLegendre this is not necessary, as it is already done by ExponentialLegendre)
    /*
    std::vector<double> coeffs_descale(degree, 0.);

    for (int n=0; n<degree; n++) {

        double prefactor = coeffs[n].val / numbers::int_pow(x_scale, n);

        for (int k=0; k<=n; k++) { 
            coeffs_descale[k] += prefactor * numbers::n_choose_k(n, k) * numbers::int_pow(-x_center, n-k); 
        }
    }
    
    poly.SetParams(coeffs_descale);     
    */ 
    
    double eta = newton_optimizer(data, poly, coeffs); 

#ifdef DEBUG
    int i=0; 
    std::cout <<"<"<<__func__<">: final coeffs:\n"; 
    for (auto& coeff : exp_poly.GetParams()) {
        std::printf("   %3i : %f\n", i++, coeff); 
    }
#endif 

    //check for NaN vals. 
    if (numbers::is_nan(eta) || numbers::contains_nan(poly.GetParams())) {
        Error(__func__, "Nan value encountered in newton opimization step.");
        return FitResult<ExponentialLegendre>::Fail(); 
    }

    return { poly, Status::kSuccess }; 
}
//______________________________________________________________________________________________________________________________

};