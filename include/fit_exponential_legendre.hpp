#ifndef peak_search_fit_exponential_legendre_hpp
#define peak_search_fit_exponential_legendre_hpp

#include <Fcn1D/ExponentialLegendre.hpp>
#include "FitResult.hpp"
#include <Histo1D.hpp> 

//ROOT headers
#include <TH1D.h> 


namespace peak_search
{

/// @brief Fit an exponentiated legendre-polynomial to the given histogram 
/// @param hist histogram to fit 
/// @param degree degree of exponentiated polynomial 
/// @return Exponential polynomial fit-result 
FitResult<ExponentialLegendre> fit_exponential_legendre(TH1D* hist, int degree);

/// @brief Fit an exponentiated legendre-polynomial to the given histogram 
/// @param hist histogram to fit 
/// @param degree degree of exponentiated polynomial 
/// @return Exponential polynomial fit-result 
FitResult<ExponentialLegendre> fit_exponential_legendre(const Histo1D& data, int degree);

};

#endif