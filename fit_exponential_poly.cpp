#include "fit_exponential_poly.hpp"

//ROOT header
#include <TAxis.h> 
#include <TError.h> 
//eigen header
#include <eigen3/Eigen/Core> 
#include <eigen3/Eigen/Dense>
//stdlib headers
#include <vector> 

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
    
    //download all the histogram data locally
    std::vector<binval_t> bins; 
    bins.reserve(xax->GetNbins());
    for (int b=1; b<=xax->GetNbins(); b++) bins.emplace_back( xax->GetBinCenter(b), hist->GetBinContent(b) );

    //first, try to fit the polynomial using a chi-square fit. 
    MatrixXd A = MatrixXd::Zero(degree, degree);
    VectorXd B = VectorXd::Zero(degree);

    for (const auto& bin : bins) {

        //don't want to take the logarithm of 0! 
        if (bin.N < 1) continue; 

        double Xmu[degree]; 

        Xmu[0] = 1.; 
        for (int i=1; i<degree; i++) Xmu[i] = Xmu[i-1]*bin.x; 

        for (int i=0; i<degree; i++) {

            B(i) += std::log(bin.N) * Xmu[i];

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
    
    FitResult<ExponentialPoly> result{ ExponentialPoly{coeffs}, Status::kSuccess };  
    return result; 
}

};