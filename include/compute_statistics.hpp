#ifndef peak_search_compute_statistics_hpp
#define peak_search_compute_statistics_hpp

#include <Fcn1D/Fcn1D.hpp>
#include <FitResult.hpp>
#include <Histo1D.hpp>

namespace peak_search 
{

struct FitStats {
    double mass, mu_MLE, Q0, mu_CL, epsilon2_CL, mu_sigma; 
};

/// @brief Compute the fit-statistics for the given mass hypothesis. 
/// @param data data to fit
/// @param fcn functional form. assumed that parameter[0] is mu (signal-strength parameter)
/// @param mass mass hypothesis
/// @param CLs confidence level (given epsilon^2_CL is the **true value** of epsilon^2 there is a p = CL% chance that any identical experiment will measure a best-fit 'mu' value _less than or equal to_ our observed best-fit 'mu')
/// @param mass_window_size size of 'mass window' used in ratio technique (This is distict from the 'mass window' over which the fit is performed!)
/// @return List of Fit stats. Null 'FitResult' object if fit fails. 
FitResult<FitStats> compute_statistics(const Histo1D& data, Fcn1D& fcn_s, Fcn1D& fcn_b, double mass, double CLs=0.05, double mass_window_size=1.); 

};

#endif