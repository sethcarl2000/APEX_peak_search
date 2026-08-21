#include "fit_exponential_poly.hpp"
#include "best_likelihood_fit.hpp"
#include "log_likelihood.hpp"
#include "fit_parameter.hpp"
#include <newton_optimizer.hpp>
#include <make_histogram_copy.hpp>
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
FitResult<ExponentialPoly> fit_exponential_poly(TH1D* hist, int degree)
{
    using Eigen::MatrixXd, Eigen::VectorXd; 

    //extract the histogram values
    if (!hist) {
        Error(__func__, "hist passed is null");
        return FitResult<ExponentialPoly>::Null(); 
    }

    Histo1D data = make_histogram_copy(hist);

    return fit_exponential_poly(data, degree);
}
//______________________________________________________________________________________________________________________________
FitResult<ExponentialPoly> fit_exponential_poly(const Histo1D& data, int degree)
{
    using Eigen::MatrixXd, Eigen::VectorXd; 
    
    //make a copy of the histogram 
    Histo1D hist = data;

    double dx = (hist.GetXmax() - hist.GetXmin())/((double)hist.GetNbins()); 

    double x_scale  = (hist.GetXmax() - hist.GetXmin())/2.;
    double x_center = (hist.GetXmax() + hist.GetXmin())/2.;  
    
    //first, try to fit the polynomial using a chi-square fit. 
    MatrixXd A = MatrixXd::Zero(degree, degree);
    VectorXd B = VectorXd::Zero(degree);

    for (auto& bin : hist.bins) {

        //don't want to take the logarithm of 0! 
        if (bin.N <= 1e-6) continue; 

        double Xmu[degree]; 

        double x = (0.5*(bin.xmax + bin.xmin) - x_center)/x_scale; 

        Xmu[0] = 1.; 
        for (int i=1; i<degree; i++) Xmu[i] = Xmu[i-1]*x; 

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
    ExponentialPoly poly({}, hist.GetXmin(), hist.GetXmax()); 
    poly.SetParams(coeffs_vec); 

    //check for NaN
    if (coeffs.size() != (size_t)degree) { return FitResult<ExponentialPoly>::Fail(); }

    poly.SetParams(coeffs); 

    //now, de-scale each exponent. (! update, for ExponentialPoly this is not necessary, as it is already done by ExponentialPoly)
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
    
    double eta = newton_optimizer(hist, poly, coeffs); 

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
        return FitResult<ExponentialPoly>::Fail(); 
    }

    return { poly, Status::kSuccess }; 
}
//______________________________________________________________________________________________________________________________

};