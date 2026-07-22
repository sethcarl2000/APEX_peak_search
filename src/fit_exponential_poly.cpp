#include "fit_exponential_poly.hpp"
#include "best_likelihood_fit.hpp"
#include "log_likelihood.hpp"
#include "fit_parameter.hpp"
#include "bininfo.hpp"
#include <newton_optimizer.hpp>

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

    auto data = copy_1D_hist(hist);

    return fit_exponential_poly(data, degree);
}
//______________________________________________________________________________________________________________________________
FitResult<ExponentialPoly> fit_exponential_poly(histo_1D_t data, int degree)
{
    using Eigen::MatrixXd, Eigen::VectorXd; 
    
    double dx = (data.xmax - data.xmin)/((double)data.bins.size()); 

    double x_scale  = (data.xmax - data.xmin)/2.;
    double x_center = (data.xmax + data.xmin)/2.;  
    
    //first, try to fit the polynomial using a chi-square fit. 
    MatrixXd A = MatrixXd::Zero(degree, degree);
    VectorXd B = VectorXd::Zero(degree);

    for (const auto& bin : data.bins) {

        //don't want to take the logarithm of 0! 
        if (bin.N <= 1e-6) continue; 

        double Xmu[degree]; 

        double x = (bin.x - x_center)/x_scale; 
        Xmu[0] = 1.; 
        for (int i=1; i<degree; i++) Xmu[i] = Xmu[i-1]*x; 

        for (int i=0; i<degree; i++) {

            B(i) += std::log( bin.N ) * Xmu[i];

            for (int j=0; j<degree; j++) {

                A(i,j) += Xmu[i]*Xmu[j];
            }
        }
    }
    
    //now, solve this system 
    VectorXd coeffs_vec = A.llt().solve(B); 
    
    std::vector<fit_parameter_t> coeffs; coeffs.reserve(coeffs_vec.size()); 
    
    for (size_t i=0; i<coeffs_vec.size(); i++) {
        coeffs.emplace_back(fit_parameter_t{ .val = coeffs_vec[i], .name = Form("c%zi",i), .is_fixed = false }); 
    }
    //std::vector<double> coeffs( coeffs_vec.data(), coeffs_vec.data() + coeffs_vec.size() );
    //now, let's 
    ExponentialPoly poly({}, data.xmin, data.xmax); 
    poly.SetParams(coeffs_vec); 

    //check for NaN
    if (coeffs.size() != (size_t)degree) { return FitResult<ExponentialPoly>::Fail(); }

    //scale this coefficient to that we can *integrate* over each bin. 
    coeffs[0].val += std::log( ((double)data.bins.size())/2. );

    newton_optimizer(data, poly, coeffs); 
    
#ifdef DEBUG
    std::cout << "coeffs: ";
    for (auto x : coeffs) std::cout << x << " "; 
    std::cout << "\n";
#endif 
    
    return { poly, Status::kSuccess }; 
}
//______________________________________________________________________________________________________________________________

};