#ifndef peak_search_solve_for_CLs_hpp
#define peak_search_solve_for_CLs_hpp

#include <Fcn1D/Fcn1D.hpp>
#include <Histo1D.hpp>

#include <TH1D.h> 
#include <TH2D.h> 

#include <vector> 

namespace peak_search 
{


/// @brief Solve for an upper-limit at the given confidence level, given the 
/// @param data data to fit 
/// @param fcn function to perform MLE fits with. first param is assumed to be signal paramter 'mu'. 
/// @param CL confidence level to test (i.e., 0.95 => 95% frequentist confidence that true parameter value is below this value)
/// @param rel_tolerance relative confidence with which to report parameter value
/// @param grid_search_size initial grid-search size to scan 'mu' with. Once the bounding interval for mu_CL is found
/// @return mu_CL, the upper limit on 'mu', computed for the given CL (and within the requested rel. tolerance)
double solve_for_CLs(const Histo1D& data, Fcn1D& fcn, double CL=0.95, double rel_tolerance=1e-2, double grid_search_size=100); 

};

#endif