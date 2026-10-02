
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
#include <compute_Q0.hpp>
#include <compute_statistics.hpp> 
// ROOT headers
#include <TCanvas.h> 
#include <TH2D.h> 
#include <TStyle.h>
#include <TAxis.h> 
#include <TGraphErrors.h>
#include <TLegend.h> 
// stdlib
#include <vector> 
#include <iostream> 

/// @brief Returns estimate of mass resolution for the given mass hypothesis (this is a slightly conservative over-estimate)
/// @param mass_hypothesis mass hypothesis (MeV)
/// @return estimate of mass resolution (MeV)
double mass_resolution(double mass_hypothesis)
{
    return 1. + (mass_hypothesis - 140.) * ((0.8 - 1.0)/(270 - 140)); 
}

/// @brief Return empty vector of given type, with reserved size of 'n'
template <typename T> std::vector<T> make_resd_vec(const std::size_t n) {
    std::vector<T> vec; vec.reserve(n);
    return vec;  
}

///________________________________________________________________________________________________________
void test_scan()
{
    using namespace peak_search; 

    const double min_mass = 145.; 
    const double max_mass = 275.; 

    int n_steps = 400; 
    int n_bins  = n_steps/4; 

    //pick a reasonable number of bins
    FitTest::Configuration config; 

    config.total_stats = 100e6; 

    config.n_steps_per_task = 200; 

    auto h_m_vs_mu = new TH2D(
        "h_signal", "Best-fit signal parameter '#mu' vs m;signal mass hypothesis (MeV);best-fit #mu", 
        n_bins, min_mass, max_mass,
        100, -40e3, 40e3
    ); 

    auto h_m_vs_uCL = new TH2D(
        "h_uCL", "Signal parameter upper-limit '#mu_{>0.95}' vs m;signal mass hypothesis (MeV);log_{10} #mu_{>0.95}", 
        n_bins, min_mass, max_mass,
        100, -2, 6
    ); 
    
    auto h_m_vs_e2CL = new TH2D(
        "h_e2CL", "Coupling CL_{s} upper limit 0.95;;signal mass hypothesis (MeV);#epsilon^{2}, CL=0.95", 
        n_bins, min_mass, max_mass,
        100, -9, -5
    ); 

    auto h_m_vs_Z = new TH2D(
        "h_Z", "Significance Z ~ #sqrt{Q0} vs m;signal mass hypothesis (MeV);Significance Z (n. #sigma)",
        n_bins, min_mass, max_mass,
        100, -7, 7
    ); 

    auto h_pQ0 = new TH1D(
        "h_pZ", "p(Q0) vs m;signal mass hypothesis (MeV);p(Q0)",  
        50, 0., 1.
    );  
    //config.n_threads = 1;
    
    const int n_mass_tests = 400; 
    auto p_mass = config.params.Append(n_mass_tests, min_mass, max_mass);   

    FitTest::Outputs outputs; 

    auto p_m_vs_mu   = outputs.Add(h_m_vs_mu);
    auto p_m_vs_Z    = outputs.Add(h_m_vs_Z);
    auto p_m_vs_uCL  = outputs.Add(h_m_vs_uCL);
    auto p_m_vs_e2CL = outputs.Add(h_m_vs_e2CL);
    auto p_pQ0       = outputs.Add(h_pQ0);

    const int n_scans = 200; 

    bool fill_pts = false; 
    auto pts_m    = make_resd_vec<double>(n_mass_tests); 
    auto pts_eps2 = make_resd_vec<double>(n_mass_tests); 

    auto fit_window_fcn = static_cast<FitTest::Function>([&](FitTest::ThreadManager* mgr)
    {
        double window_size = 7.; // MeV 

        const auto params = mgr->GetParamList(); 

        double mass = params[p_mass]; 
        double resolution = mass_resolution(mass); 

        double m_min = mass - window_size*resolution; 
        double m_max = mass + window_size*resolution; 

        int n_bins = (m_max - m_min)/(0.5); 

        const auto spectrum = mgr->GetSpectrum(n_bins, mass - window_size*resolution, mass + window_size*resolution);

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

        if (fill_pts) {
            pts_m.emplace_back(mass); 
            pts_eps2.emplace_back(std::log10(epsilon2_CL)); 
        } else {
            mgr->GetOutput<TH2D>(p_m_vs_mu)  ->Fill(mass, mu); 
            mgr->GetOutput<TH2D>(p_m_vs_Z)   ->Fill(mass, Z); 
            mgr->GetOutput<TH2D>(p_m_vs_uCL) ->Fill(mass, std::log10(mu_cl95)); 
            mgr->GetOutput<TH2D>(p_m_vs_e2CL)->Fill(mass, std::log10(epsilon2_CL)); 

            mgr->GetOutput<TH1D>(p_pQ0)->Fill(pQ0);
        } 
        return;
    });


    FitTest::Run(n_scans, config, outputs, fit_window_fcn); 
    //now, we're going to do one **real** scan (still on the accidental spectrum)
    TCanvas *canv; 

    //now, run once and fill test points
    fill_pts = true; 
    config.n_threads = 1;
    FitTest::Run(1, config, outputs, fit_window_fcn, 0); 
    std::cout << "size: " << pts_m.size() << "\n"; 

    canv = new TCanvas;
    gStyle->SetOptStat(0); 

    h_m_vs_mu->Draw("col"); 

    canv = new TCanvas;
    h_m_vs_Z->Draw("col"); 

    canv = new TCanvas;
    h_m_vs_uCL->Draw("col"); 

    canv = new TCanvas;
    canv->SetTopMargin(0.15);
    h_m_vs_e2CL->SetTitle(Form("CL=0.95 upper limits on #varepsilon^{2}, %.1f x 10^{6} events;signal mass hypothesis (MeV);#epsilon^{2}, CL=0.95", config.total_stats/1e6));
    make_brazil_flag_plot(h_m_vs_e2CL, Form("Avg. of %i pseudo-spectra",n_scans)); 
    auto g = new TGraph(pts_m.size(), pts_m.data(), pts_eps2.data());  
    g->Draw("SAME"); 

    canv = new TCanvas; 
    h_pQ0->SetMaximum( h_pQ0->GetMaximum()*1.5 );
    h_pQ0->SetMinimum( 0. );  
    h_pQ0->Draw("HIST"); 
}



