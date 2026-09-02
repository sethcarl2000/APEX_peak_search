
#include <FitTestFunction.hpp>
#include <FitTestKernel.hpp>
#include <FitTestThreadManager.hpp>
#include <solve_for_CLs.hpp>
#include <SignalFit.hpp>
#include <Histo1D.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <fit_exponential_poly.hpp>
#include <compute_Q0.hpp>

peak_search::SignalFit fit_window(peak_search::FitTestThreadManager* mgr);
///________________________________________________________________________________________________________
void test_scan()
{

    peak_search::FitTestKernel kernel(155., 265., 10); 

    peak_search::FitTestFunction fit_window_f{fit_window}; 

    kernel.RunTest(1, fit_window_f); 

    kernel.DrawResults(); 
}
///________________________________________________________________________________________________________
peak_search::SignalFit fit_window(peak_search::FitTestThreadManager* mgr)
{
    double window_size = 7.; // MeV 

    auto mass = mgr->get_mass(); 

    peak_search::Histo1D spectrum = mgr->get_spectrum(mass - window_size, mass + window_size);

    auto gaussian_fcn = peak_search::Gauss(0, mass, 1.); 

    //fit the background
    auto background_poly = peak_search::fit_exponential_poly(spectrum, 6).data; 

    peak_search::FcnSum f_sum(&gaussian_fcn, &background_poly); 

    double Q0 = peak_search::compute_Q0(spectrum, f_sum); 

    double mu = f_sum.GetParams()[0];

    //get the middle(-ish)bin. this gives us an order-of-magnitude estimate for the natural variance of the signal paramter, mu.
    double N_middle = spectrum.bins.at( spectrum.GetNbins()/2 ).N; 

    //also, find the CLs value
    double mu_CL95 = peak_search::solve_for_CLs(spectrum, f_sum, 0.95, 0.01, std::sqrt(N_middle) );

    std::printf("found upper limit on mu: %+.5e\n", mu_CL95); 

    return peak_search::SignalFit{ .mass=mass, .mu=mu, .Q0=Q0 };
}


