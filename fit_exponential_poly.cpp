#include "fit_exponential_poly.hpp"
#include "best_likelihood_fit.hpp"
#include "log_likelihood.hpp"
#include "bininfo.hpp"

//ROOT header
#include <TAxis.h> 
#include <TError.h> 
//eigen header
#include <eigen3/Eigen/Core> 
#include <eigen3/Eigen/Dense>
//stdlib headers
#include <vector>
#include <iostream>  

namespace {
    struct binval_t { double x, N; }; 
};

namespace peak_search 
{
    
FitResult<ExponentialPoly> fit_exponential_poly(TH1D* hist, int degree)
{
    using Eigen::MatrixXd, Eigen::VectorXd; 

    //extract the histogram values
    if (!hist) {
        Error(__func__, "hist passed is null");
        return FitResult<ExponentialPoly>::Null(); 
    }
    
    auto xax = hist->GetXaxis(); 
    
    const double dx = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()); 

    //download all the histogram data locally
    auto data = copy_1D_hist(hist);

    double scale  = (xax->GetXmax() - xax->GetXmin())/2.;    
    double center = (xax->GetXmax() + xax->GetXmin())/2.;

    //re-normalize the points (to make fit go a little easier)
    for (auto& bin : data.bins) { 
        bin.x = (bin.x - center)/scale;
    }
    data.xmax = +1.;
    data.xmin = -1.; 

    //first, try to fit the polynomial using a chi-square fit. 
    MatrixXd A = MatrixXd::Zero(degree, degree);
    VectorXd B = VectorXd::Zero(degree);

    for (const auto& bin : data.bins) {

        //don't want to take the logarithm of 0! 
        if (bin.N <= 1e-6) continue; 

        double Xmu[degree]; 

        Xmu[0] = 1.; 
        for (int i=1; i<degree; i++) Xmu[i] = Xmu[i-1]*bin.x; 

        for (int i=0; i<degree; i++) {

            B(i) += std::log( bin.N ) * Xmu[i];

            for (int j=0; j<degree; j++) {

                A(i,j) += Xmu[i]*Xmu[j];
            }
        }
    }
    
    //now, solve this system 
    VectorXd coeffs_vec = A.llt().solve(B); 
    
    std::vector<double> coeffs( coeffs_vec.data(), coeffs_vec.data() + coeffs_vec.size() );

    //check for NaN
    if (coeffs.size() != (size_t)degree || numbers::contains_nan(coeffs)) { return FitResult<ExponentialPoly>::Fail(); }

    //scale this coefficient to that we can *integrate* over each bin
    double dx_scaled = 2./((double)data.bins.size());
    coeffs[0] += std::log( 1./dx_scaled );

    
#ifdef DEBUG
    std::cout << "coeffs: ";
    for (auto x : coeffs) std::cout << x << " "; 
    std::cout << "\n";
#endif 


#ifdef DEBUG
    std::cout << "coeffs_descale: ";
    for (auto x : coeffs_descale) std::cout << x << " "; 
    std::cout << "\n";
#endif
    
    //now, we will find the max-likelihood estimator 
    auto wrapper = [degree](double x, const double* par){ return ExponentialPoly::Eval(x,par,degree); };
 
    ExponentialPoly poly{coeffs};
    std::printf("Negative log likelihood for least-squares fit: %f\n", negative_log_likelihood(data, poly)); 

    double nll = best_likelihood_fit(data, wrapper, coeffs); 


    if (numbers::is_nan(nll)) { return FitResult<ExponentialPoly>::Fail(); } //*/ 

    //remove this scaling for a moment
    coeffs[0] += -std::log( 1./dx_scaled );

    //now, de-scale each exponent.
    std::vector<double> coeffs_descale(degree, 0.);

    for (int i=0; i<degree; i++) {

        double prefactor = coeffs[i] / numbers::int_pow(scale, i);

        for (int k=0; k<=i; k++) { 
            coeffs_descale[k] 
                += prefactor * numbers::n_choose_k(i, k) * numbers::int_pow(-center, i-k); 
        }
    }
    coeffs_descale[0] += std::log( 1./dx );


    return { ExponentialPoly{coeffs_descale}, Status::kSuccess }; 
}

};