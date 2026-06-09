#ifndef gauss_integrate_hpp
#define gauss_integrate_hpp

#include <pdf_fcn.hpp>

namespace peak_search
{

/// @brief Do gauss integration (n=4) for region given
/// @param fcn input fcn 
/// @param params relevant params
/// @return numerical gauss integral over region
double gauss_integrate(double (*fcn)(double), double xmin,double xmax);

/// @brief Do gauss integration (n=4) for region given
/// @param fcn input fcn 
/// @return numerical gauss integral over region
double gauss_integrate(double (*fcn)(double,double), double xmin,double xmax, double ymin,double ymax);

};

#endif