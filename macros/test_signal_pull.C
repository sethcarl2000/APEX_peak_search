

#include "make_brazil_flag_plot.h"
// peak_search headers
#include <FitTest/Function.hpp>
#include <FitTest/Run.hpp>
#include <FitTest/ThreadManager.hpp>
#include <FitTest/Configuration.hpp>
#include <FitTest/Outputs.hpp>
#include <FitTest/ParameterList.hpp>
#include <solve_for_CLs.hpp>
#include <SignalFit.hpp>
#include <Histo1D.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <fit_exponential_poly.hpp>
#include <gauss_integrate.hpp>
#include <compute_Q0.hpp>
#include <compute_statistics.hpp> 
// ROOT headers
#include <TCanvas.h> 
#include <TH2D.h> 
#include <TStyle.h>
#include <TAxis.h> 
#include <TGraphErrors.h>
#include <TLegend.h> 
#include <TLine.h> 
// stdlib
#include <vector> 

/// @brief Returns estimate of mass resolution for the given mass hypothesis (this is a slightly conservative over-estimate)
/// @param mass_hypothesis mass hypothesis (MeV)
/// @return estimate of mass resolution (MeV)
double mass_resolution(double mass_hypothesis)
{
    return 1. + (mass_hypothesis - 140.) * ((0.8 - 1.0)/(270 - 140)); 
}

///________________________________________________________________________________________________________
void test_signal_pull()
{
    using namespace peak_search; 

    const double min_mass = 145.; 
    const double max_mass = 275.; 

    const double mass = 260.; 

    int n_steps = 400; 
    int n_bins  = n_steps/4; 

    //pick a reasonable number of bins
    FitTest::Configuration config; 

    config.total_stats = 76e6; 

    config.n_steps_per_task = 200; 

    const double max_injected_signal = 50e3; 

    auto h_mu_vs_MLE = new TH2D("h_mu_MLE", Form("#mu (true) vs #hat{#mu} (mass = %.0f MeV);#mu (truth);#hat{#mu} (best-fit #mu)",mass),
        n_bins, 0, max_injected_signal, 
        n_bins, -20e3, max_injected_signal + 20e3
    );

    //config.n_threads = 1; 

    //auto p_mass = config.params.Append(400, min_mass, max_mass);   

    // amount of toy-signal to inject. 
    auto p_signal = config.params.Append(800, 0., max_injected_signal); 

    FitTest::Outputs outputs; 

    auto p_mu_vs_MLE = outputs.Add(h_mu_vs_MLE);

    auto fit_window_fcn = static_cast<FitTest::Function>([&](FitTest::ThreadManager* mgr)
    {
        double window_size = 7.; // MeV 

        const auto params = mgr->GetParamList(); 

        double injected_signal = params[p_signal]; 


        double resolution = mass_resolution(mass); 

        double m_min = mass - window_size*resolution; 
        double m_max = mass + window_size*resolution; 

        int n_bins = (m_max - m_min)/(0.5); 

        auto spectrum = mgr->GetSpectrum(n_bins, mass - window_size*resolution, mass + window_size*resolution);

        //add injected signal
        auto injected_signal_fcn = peak_search::Gauss(injected_signal, mass, resolution); 
        
        auto my_rand = mgr->GetRand(); 
        for (auto& bin : spectrum.bins) {
            double S_expect = peak_search::gauss_integrate(injected_signal_fcn, bin.xmin, bin.xmax); 
            bin.N += my_rand->PoissonD(S_expect); 
        }

        /*double N{0.}; 
        for (const auto& bin : spectrum.bins) N += bin.N; 
        std::printf("mass range: [%5.1f, %5.1f], stats: %+.3e\n", m_min,m_max, N); */ 

        auto gaussian_fcn = peak_search::Gauss(0, mass, resolution); 

        //fit the background
        auto background_poly = peak_search::fit_exponential_poly(spectrum, 6).data; 
        
        auto stat_result = peak_search::compute_statistics(spectrum, gaussian_fcn, background_poly, mass, 0.05, 1.); 

        if (stat_result.status != peak_search::Status::kSuccess) {
            Warning("fit_window_function", "Fit failed for mass: %.1f", mass); 
            return; 
        }

        auto stats = stat_result.data; 
        double Q0 = stats.Q0; 
        double mu = stats.mu_MLE; 
        double mu_cl95 = stats.mu_CL; 
        double epsilon2_CL = stats.epsilon2_CL; 

        //double mu_sigma = stats.mu_sigma; 

        double Z  = (Q0<0.?-1:+1) * std::sqrt(std::fabs(Q0)); 
 
        double pQ0 = compute_Q0_p(Q0); 

        //get the middle(-ish)bin. this gives us an order-of-magnitude estimate for the natural variance of the signal paramter, mu.
        double N_middle = spectrum.bins.at( spectrum.GetNbins()/2 ).N; 

        mgr->GetOutput<TH2D>(p_mu_vs_MLE)->Fill( injected_signal, mu ); 

        return;
    });

    FitTest::Run(200, config, outputs, fit_window_fcn); 

    //now, we're going to do one **real** scan (still on the accidental spectrum)

    new TCanvas;
    gStyle->SetOptStat(0); 
    h_mu_vs_MLE->Draw("col"); 

    auto line = new TLine(0,0, max_injected_signal,max_injected_signal); 
    line->SetLineStyle(kDashed); 
    line->SetLineColor(kRed); 
    line->SetLineWidth(2); 
    line->Draw("SAME"); 
}
///________________________________________________________________________________________________________



