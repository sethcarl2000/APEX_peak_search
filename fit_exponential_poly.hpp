#ifndef peak_search_fit_exponential_poly_hpp
#define peak_search_fit_exponential_poly_hpp

#include "ExponentialPoly.hpp"

//ROOT headers
#include <TH1D.h> 

namespace peak_search
{
/// @brief Fit an exponentiated-polynomial to the given histogram 
/// @param hist histogram to fit 
/// @param degree degree of exponentiated polynomial 
/// @return Exponential polynomial fit-result 
ExponentialPoly fit_exponential_poly(TH1D* hist, int degree);

};

#endif