#ifndef gauss_integrate_hpp
#define gauss_integrate_hpp

#include <Fcn1D/Fcn1D.hpp> 

#include <functional> 

namespace peak_search
{

/// @brief Do gauss integration (n=4) for region given
/// @param fcn input fcn 
/// @param params relevant params
/// @return numerical gauss integral over region
double gauss_integrate(const std::function<double(double)>& fcn, double xmin,double xmax);

/// @brief Do gauss integration (n=4) for region given
/// @param fcn input fcn 
/// @param params relevant params
/// @return numerical gauss integral over region
double gauss_integrate(const Fcn1D& fcn, double xmin,double xmax);

};

#endif