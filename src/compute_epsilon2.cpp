
#include <compute_epsilon2.hpp>
#include <gauss_integrate.hpp>

namespace peak_search
{

constexpr double muon_mass = 105.66; // MeV; 
constexpr double electron_mass = 0.501; // MeV; 
constexpr double alpha = 1./137.; 


/// @return the approximate branch ratio for A'-> e+e- (as opposed to A'->mu+mu-)
double get_branch_fraction(double mass)
{
    if (mass < 2.*muon_mass) return 1.; 

    double amp_electron = std::sqrt( 1. - 4.*(electron_mass*electron_mass)/(mass*mass) ) ; 
    double amp_muon     = std::sqrt( 1. - 4.*(muon_mass*muon_mass)/(mass*mass) ) ; 
    
    return amp_electron / (amp_electron + amp_muon); 
}
//________________________________________________________________________________________________________________________________
/// @return a _very_ rough estimate of the fraction of our background events that are the 'radiative' type of gamma->e+e- production. 
double estimate_radiative_fraction(double mass)
{
    return 0.205 + (mass - 120)* ((0.155 - 0.205)/(220 - 120)); 
}
//________________________________________________________________________________________________________________________________
double compute_epsilon2(const Fcn1D& fcn_b, double N_signal, double mass, double mass_window_size)
{
    double N_background = gauss_integrate(
        fcn_b, 
        mass - mass_window_size/2.,
        mass + mass_window_size/2.
    );

    double epsilon2 = (N_signal / N_background) * (mass_window_size/mass) * (2.*alpha) / (3.*3.1415926536); 

    double branch_fraction = get_branch_fraction(mass); 
    double radiative_frac  = estimate_radiative_fraction(mass);

    // apply a few rough corrections
    return epsilon2 / (branch_fraction * radiative_frac); 
}   
//________________________________________________________________________________________________________________________________
double get_n_signal_events(const Fcn1D& fcn_b, double epsilon2, double mass, double mass_window_size)
{
    double N_background = gauss_integrate(
        fcn_b, 
        mass - mass_window_size/2.,
        mass + mass_window_size/2.
    );

    //double epsilon2 = (N_signal / N_background) * (mass_window_size/mass) * (2.*alpha) / (3.*3.1415926536); 

    double branch_fraction = get_branch_fraction(mass); 
    double radiative_frac  = estimate_radiative_fraction(mass);

    epsilon2 *= (branch_fraction * radiative_frac);
    
    return epsilon2 * N_background * (mass/mass_window_size) * (3.*3.1415926536) / (2.*alpha);
}

};
