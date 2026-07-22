#ifndef peak_search_compute_Q0_hpp
#define peak_search_compute_Q0_hpp

#include <Fcn1D/Fcn1D.hpp>
#include "bininfo.hpp"

#include <TH1D.h> 
#include <TH2D.h> 

#include <vector> 

namespace peak_search 
{

/// @brief Compute the 'Q0' statistic 
/// @param hist histogram representing data 
/// @param fcn fcn representing the expectation value of each bin 
/// @return Q0 statistic (representing degree of agreement with no-signal hypothesis)
double compute_Q0(TH1D* hist, Fcn1D& fcn, std::vector<double>& nuissance_params); 

/// @brief Compute the 'Q0' statistic 
/// @param hist histogram representing data 
/// @param pdf PDF of expecation value for each bin 
/// @param nuissance_params first-guess of nuissance parameters to fit. 
/// @return Q0 statistic (representing degree of agreement with no-signal hypothesis)
double compute_Q0(histo_1D_t data, Fcn1D& fcn, std::vector<double>& nuissance_params); 


/// @brief returns the p-value for a computed Q0
/// @param Q0 a computed value of 'Q0' 
/// @return p-value (prob. that any observation of Q0 is as great or greater than passed value)
double compute_Q0_p(double Q0); 

};

#endif