#ifndef peak_search_compute_epsilon2_hpp
#define peak_search_compute_epsilon2_hpp

#include <Fcn1D/Fcn1D.hpp> 

namespace peak_search
{

/// @brief Compute epsilon^2
/// @param fcn_b background function
/// @param n_signal_events number of signal events 
/// @param mass current mass hypothesis  
/// @param mass_window_size mass window over which to integrate
/// @return estimate for epsilon^2 
double compute_epsilon2(const Fcn1D& fcn_b, double n_signal_events, double mass, double mass_window_size=1.);

/// @brief Given epsilon^2, return the expected number of signal events
/// @param fcn_b background-only fcn 
/// @param epsilon2 epsilon^2 (coupling strength)
/// @param mass mass hypothesis
/// @param mass_window_size mass window size 
/// @return expected number of signal events 
double get_n_signal_events(const Fcn1D& fcn_b, double epsilon2, double mass, double mass_window_size=1.); 


}; 

#endif