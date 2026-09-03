
#include <FitTestFunction.hpp>
#include <FitTestKernel.hpp>
#include <FitTestThreadManager.hpp>
#include <solve_for_CLs.hpp>
#include <newton_optimizer.hpp>
#include <log_likelihood.hpp>
#include <fit_parameter.hpp>
#include <SignalFit.hpp>
#include <Histo1D.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <fit_exponential_poly.hpp>
#include <compute_Q0.hpp>
#include <gauss_integrate.hpp>
// ROOT headers
#include <TCanvas.h> 
#include <TH2D.h> 
#include <TStyle.h>
#include <TAxis.h> 
// stdlib
#include <cstdio> 

constexpr int poly_order = 6; 

//void fit_window(peak_search::FitTestThreadManager* mgr);
///________________________________________________________________________________________________________
void test_compute_sigma()
{
    using namespace peak_search; 

    auto& kernel = FitTestKernel::Instance(); 

    const double min_mass = 200.; 
    const double max_mass = 200.; 

    int n_steps = 1; 
    int n_bins  = 100; 


    const double stat_scaling = 1./10.;

    const double sigma_est = 6.7e3 * std::sqrt(stat_scaling); 
    const double mu_range  = 3*sigma_est; 

    const double max_y_range = 1.5*std::pow( mu_range/sigma_est, 2 ); 

    auto h_sigma = new TH2D("h", "mu vs sigma;#mu - #hat{#mu};2( NLL(#mu) - NLL(#hat{#mu})", 
        200, -mu_range, +mu_range,
        200, 0, max_y_range
    ); 

    auto ptr_sigma = kernel.AddTH2D(*h_sigma); 

    kernel.SetMassRange(min_mass, max_mass); 
    kernel.SetNSteps(n_steps); 
    kernel.SetMassBinSize(0.5); 

    kernel.SetTotalStats(250e6 * stat_scaling); 

    auto fit_window_fcn = static_cast<FitTestFunction>([ptr_sigma](peak_search::FitTestThreadManager* mgr)
    {
        double window_size = 7.; // MeV 

        auto mass = mgr->get_mass(); 

        peak_search::Histo1D spectrum = mgr->get_spectrum(mass - window_size, mass + window_size);

        auto gaussian_fcn = peak_search::Gauss(0, mass, 1.); 

        //fit the background
        auto background_poly = peak_search::fit_exponential_poly(spectrum, poly_order).data; 

        peak_search::FcnSum f_sum(&gaussian_fcn, &background_poly); 

        //get the best-fit 
        std::vector<fit_parameter_t> params; params.reserve(poly_order+1); 
        for (const auto& val : f_sum.GetParams()) { params.emplace_back(val, "", false); }

        auto& mu_param = params[0]; 

        double NLL_MLE = peak_search::newton_optimizer(spectrum, f_sum, params); 

        double mu_MLE = mu_param.val;
        mu_param.is_fixed = true; 

        std::printf("\nmu_MLE: %+.8e\n", mu_param.val); 

        auto h = mgr->GetUserTH2D(ptr_sigma); 

        double n_events_2MeV_window = peak_search::gauss_integrate(background_poly, mass-1., mass+1); 

        double x_guess = 2.7 * n_events_2MeV_window; 

        std::printf(
            "n. eventsi n 2 MeV window: %.6e\n"
            "guess for sigma: %.6e\n", 
            x_guess,    
            n_events_2MeV_window
        ); 

        mu_param.val = x_guess + mu_MLE; 
        double NLL = peak_search::newton_optimizer(spectrum, f_sum, params); 

        double sigma = x_guess / std::sqrt(2.*(NLL - NLL_MLE)); 
        std::printf("sigma: %.8e:\n", sigma); 
        
        auto x_ax = h->GetXaxis(); 
        const int n_bins = h->GetXaxis()->GetNbins(); 

        for (int i=1; i<=n_bins; i++) {
            
            mu_param.val = x_ax->GetBinCenter(i) + mu_MLE; 

            double NLL = peak_search::newton_optimizer(spectrum, f_sum, params); 


            h->Fill( x_ax->GetBinCenter(i), 2.*(NLL - NLL_MLE) );
        }

        return; 
    });

    kernel.RunTest(1, fit_window_fcn); 

    new TCanvas;
    h_sigma->Draw("col2");
}
///________________________________________________________________________________________________________



