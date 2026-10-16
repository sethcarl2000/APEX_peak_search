
#include "make_brazil_flag_plot.h"
#include "accidental_spectra_generator.h"
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
// analysis_utils headers 
#include <analysis_utils.hpp> 
// ROOT headers
#include <TCanvas.h> 
#include <TH2D.h> 
#include <TStyle.h>
#include <TAxis.h> 
#include <TLine.h> 
#include <TGraphErrors.h>
#include <TLegend.h> 
// stdlib
#include <vector> 
#include <thread> 
#include <iostream> 
#include <memory> 

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
void test_stats_level()
{
    using namespace peak_search; 

    const double min_mass = 150.; 
    const double max_mass = 270.; 


    //let's test at our most sensitive level
    const double test_mass0 = 220.;

    const double min_stats = 20e6; 
    const double max_stats = 275e6; 

    int n_steps = 400; 
    int n_bins  = n_steps/2; 


    using namespace analysis_utils; 
    using ASG = peak_search::AccidentalSpectraGenerator; 


    auto t_stats_vs_mu = root::threadlocal_hist(new TH2D(
        "h_signal", "Best-fit signal parameter '#mu' vs N. prompt stats;total prompt stats;best-fit #mu", 
        n_bins, min_stats, max_stats,
        100, -40e3, 40e3
    )); 

    auto t_stats_vs_uCL = root::threadlocal_hist(new TH2D(
        "h_uCL", "Signal parameter upper-limit '#mu_{>0.90}' vs N. prompt stats;total prompt stats;log_{10} #mu_{>0.90}", 
        n_bins, min_stats, max_stats,
        100, -2, 6
    )); 
    
    auto t_stats_vs_e2CL = root::threadlocal_hist(new TH2D(
        "h_e2CL", "Coupling CL_{s} upper limit 0.90 vs N. prompt stats;total prompt stats;#epsilon^{2}, CL=0.90", 
        n_bins, min_stats, max_stats,
        100, -7.75, -5.50
    )); 

    
    const size_t n_threads = std::thread::hardware_concurrency(); 
    
    //create the mass-spectrum generators
    std::vector<std::unique_ptr<ASG>> gens; gens.reserve(n_threads); 

    for (size_t t=0; t<n_threads; t++) {
        gens.emplace_back(std::make_unique<ASG>(min_stats, t+1)); 
    }

    const double bin_size = 0.5; 
    const size_t total_bins = (ASG::max_mass() - ASG::min_mass())/bin_size; 


    const double window_size = 8.; // MeV

    int n_bins_subhist = (2.*window_size)/bin_size; 
    Histo1D subhist_template; 
    subhist_template.bins.reserve(n_bins_subhist); 
    double m=test_mass0 - window_size; 
    for (int i=0; i<n_bins_subhist; i++) { subhist_template.bins.emplace_back(m, m+bin_size, 0.); m += bin_size; }

    auto t_subhist = thread_local_obj(subhist_template); 

    const int max_scans = 750; 

    task::config_t cfg; 
    cfg.n_scans = max_scans; 
    cfg.verbosity = 1; 
    cfg.chunk_size = 128; 
    //cfg.n_threads = 1; 

    struct stats_t { double total_n_prompt; }; 
    task::param_list<stats_t> plist; 

    plist.add(&stats_t::total_n_prompt, min_stats, max_stats, n_steps); 

    int n_fitpts = n_steps; 


    auto ffcn = [&](const stats_t *__restrict str, std::size_t thread_id){

        auto& subhist = t_subhist(thread_id); 


        double total_stats = str->total_n_prompt; 

        //make a new spectrum 
        auto& bgen = gens[thread_id]; 
        bgen->SetStats(total_stats); 
        bgen->FillSpectra(subhist); 

        double resolution = mass_resolution(test_mass0); 

        auto gaussian_fcn = peak_search::Gauss(0, test_mass0, resolution); 

        //fit the background
        auto background_poly = peak_search::fit_exponential_poly(subhist, 6).data; 
        
        auto stat_result = peak_search::compute_statistics(subhist, gaussian_fcn, background_poly, test_mass0, 0.10, nullptr, 1.); 

        if (stat_result.status != peak_search::Status::kSuccess) {
            Warning("fit_window_function", "Fit failed for mass: %.1f", test_mass0); 
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

        t_stats_vs_mu(thread_id)  ->Fill(total_stats, mu); 
        t_stats_vs_uCL(thread_id) ->Fill(total_stats, std::log10(mu_cl95)); 
        t_stats_vs_e2CL(thread_id)->Fill(total_stats, std::log10(epsilon2_CL)); 

    };

    task::scan_params(plist, ffcn, cfg);

    TCanvas* canv; 
    TGraph* g; 

    canv = new TCanvas;
    gStyle->SetOptStat(0); 
    auto h_stats_vs_uCL = t_stats_vs_uCL.aggregate(); 
    make_brazil_flag_plot(h_stats_vs_uCL, Form("Avg. of %i pseudo-spectra",max_scans));  
    

    canv = new TCanvas;
    auto h_stats_vs_e2CL = t_stats_vs_e2CL.aggregate(); 
    auto g_median_e2CL = make_brazil_flag_plot(h_stats_vs_e2CL, Form("Avg. of %i pseudo-spectra",max_scans));  

    const double current_stats = 111e6;
    const double best_possible = 141e6; 

    double ymin{g_median_e2CL->GetYaxis()->GetXmin()}, ymax{g_median_e2CL->GetYaxis()->GetXmax()}; 
    auto line = new TLine;
    line->SetLineStyle(kDashed); 
    
    line->DrawLine(current_stats,ymin, current_stats,ymax); 
    line->DrawLine(best_possible,ymin, best_possible,ymax); 

    std::printf(
        "\\varepsilon^2 0.90 CL limit achieveable:\n"
        "   %+.5e    111 x 10^6 prompt stats (current)\n"
        "   %+.5e    141 x 10^6 prompt stats (est. stats without VDC hardware issue)\n",
        std::pow( 10., g_median_e2CL->Eval(current_stats) ),
        std::pow( 10., g_median_e2CL->Eval(best_possible) )
    ); 

    canv = new TCanvas; 
    t_stats_vs_mu.aggregate()->Draw("col"); 

}



