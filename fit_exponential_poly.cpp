#include "fit_exponential_poly.hpp"

#include <TAxis.h> 

#include <eigen3/Eigen/Core> 

namespace peak_search 
{
    
FitResult<ExponentialPoly> fit_exponential_poly(TH1D* hist, int degree)
{

    return FitResult<ExponentialPoly>{ .status = Status::kNull };
}

};