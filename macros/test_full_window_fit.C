
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
void test_full_window_fit()
{
    using namespace peak_search; 

    const double min_mass = 150.; 
    const double max_mass = 270.; 

    int n_steps = 400; 
    int n_bins  = n_steps/4; 
    const double dm_fit = (max_mass - min_mass)/((double)n_steps-1); 

    //pick a reasonable number of bins
    FitTest::Configuration config; 

    config.total_stats = 246e6; 

    config.n_steps_per_task = 200; 


    using namespace analysis_utils; 
    using ASG = peak_search::AccidentalSpectraGenerator; 

    auto t_Z_corr = root::threadlocal_hist(new TH2D(
        "h_Z", "Covariance matrix of Z (significance): #braket{Z(m_{1})Z(m_{2})};mass 1;mass 2",
        n_bins, min_mass-dm_fit/2., max_mass+dm_fit/2.,
        n_bins, min_mass-dm_fit/2., max_mass+dm_fit/2.
    )); 

    auto t_m_vs_mu = root::threadlocal_hist(new TH2D(
        "h_signal", "Best-fit signal parameter '#mu' vs m;signal mass hypothesis (MeV);best-fit #mu", 
        n_bins, min_mass, max_mass,
        100, -40e3, 40e3
    )); 

    auto t_m_vs_uCL = root::threadlocal_hist(new TH2D(
        "h_uCL", "Signal parameter upper-limit '#mu_{>0.90}' vs m;signal mass hypothesis (MeV);log_{10} #mu_{>0.90}", 
        n_bins, min_mass, max_mass,
        100, -2, 6
    )); 
    
    auto t_m_vs_e2CL = root::threadlocal_hist(new TH2D(
        "h_e2CL", "Coupling CL_{s} upper limit 0.90;;signal mass hypothesis (MeV);#epsilon^{2}, CL=0.90", 
        n_bins, min_mass, max_mass,
        100, -9, -5
    )); 

    auto t_m_vs_Z = root::threadlocal_hist(new TH2D(
        "h_Z", "Significance Z ~ #sqrt{Q0} vs m;signal mass hypothesis (MeV);Significance Z (n. #sigma)",
        n_bins, min_mass, max_mass,
        100, -7, 7
    )); 


    auto t_pQ0 = root::threadlocal_hist(new TH1D(
        "h_pZ", "p(Q0) vs m;signal mass hypothesis (MeV);p(Q0)",  
        50, 0., 1.
    ));  
    //config.n_threads = 1;

    
    const size_t n_threads = std::thread::hardware_concurrency(); 
    
    //create the mass-spectrum generators
    std::vector<std::unique_ptr<ASG>> gens; gens.reserve(n_threads); 

    for (size_t t=0; t<n_threads; t++) {
        gens.emplace_back(std::make_unique<ASG>(config.total_stats, t+1)); 
    }

    const double bin_size = 0.5; 
    const size_t total_bins = (ASG::max_mass() - ASG::min_mass())/bin_size; 

     


    auto hist1d_template = gens.front()->MakeNewSpectra(total_bins, ASG::min_mass(), ASG::max_mass()); 

    const double window_size = 8.; // MeV

    int n_bins_subhist = (2.*window_size)/bin_size; 
    Histo1D subhist_template; 
    subhist_template.bins.reserve(n_bins_subhist); 
    double m=0.; 
    for (int i=0; i<n_bins_subhist; i++) { subhist_template.bins.emplace_back(m, m+bin_size, 0.); m += bin_size; }

    auto t_hist1d = thread_local_obj(hist1d_template); 
    auto t_subhist = thread_local_obj(subhist_template); 

    //assumes evenly-spaced bins!
    auto copy_subhist = [bin_size](const Histo1D& src, Histo1D& dest, double center) {
        double bin_span = dest.GetXmax() - dest.GetXmin(); 

        int first_bin = (center - bin_span/2. - src.GetXmin())/bin_size; 
        int last_bin  = first_bin + dest.bins.size(); 

        first_bin = std::max<int>( 0, first_bin );
        last_bin  = std::min<int>( src.bins.size()-1, last_bin ); 

        //std::cout << "center: " << center << "  min/max (approx): " << (center-bin_size/2.) << "/" << (center+bin_size/2.) << "\n"; 
        //std::cout << "copying " << (last_bin-first_bin+1) << " bins, in range: ["<<dest.GetXmin()<<","<<dest.GetXmax()<<"]\n"; 
        //std::cout << "bin range: ["<<first_bin<<","<<last_bin<<"]\n"; 
        int i_bin=first_bin; 
        for (auto& bin : dest.bins) {
            bin = src.bins[i_bin]; 
            ++i_bin; 
        }
    };

    const int max_scans = 320; 

    task::config_t cfg; 
    cfg.n_scans = max_scans; 
    cfg.verbosity = 1; 
    cfg.chunk_size = 1; 
    //cfg.n_threads = 1; 

    struct DummyStr {}; 
    task::param_list<DummyStr> plist; 

    int n_fitpts = n_steps; 

    bool fill_pts=false; 
    auto pts_m    = make_resd_vec<double>(n_fitpts); 
    auto pts_eps2 = make_resd_vec<double>(n_fitpts);  
    auto pts_Z    = make_resd_vec<double>(n_fitpts);  
    auto pts_uCL   = make_resd_vec<double>(n_fitpts); 

    auto t_pts_mass = thread_local_obj(std::vector<double>(n_fitpts)); 
    auto t_pts_Z    = thread_local_obj(std::vector<double>(n_fitpts)); 

    auto ffcn = [&](const DummyStr *__restrict str, std::size_t thread_id){

        auto& hist = t_hist1d(thread_id); 
        auto& subhist = t_subhist(thread_id); 

        //make a new spectrum 
        gens[thread_id]->FillSpectra(hist); 

        for (int i=0; i<n_fitpts; i++) {

            double mass = min_mass + ((double)i)*dm_fit; 
            double resolution = mass_resolution(mass); 

            //std::cout << "trial " << i << ", mass: " << mass << "\n"; 

            double m_min = mass - window_size*resolution; 
            double m_max = mass + window_size*resolution; 

            copy_subhist(hist, subhist, mass);              

            int n_bins = (m_max - m_min)/(0.5); 

            auto gaussian_fcn = peak_search::Gauss(0, mass, resolution); 

            //fit the background
            auto background_poly = peak_search::fit_exponential_poly(subhist, 6).data; 
            
            auto stat_result = peak_search::compute_statistics(subhist, gaussian_fcn, background_poly, mass, 0.10, 1.); 

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

            if (fill_pts) {
            
                pts_m   .emplace_back(mass); 
                pts_eps2.emplace_back(std::log10(epsilon2_CL)); 
                pts_Z   .emplace_back(Z); 
                pts_uCL .emplace_back(std::log10(mu_cl95)); 

            } else {

                t_pts_mass(thread_id)[i] = mass; 
                t_pts_Z(thread_id)[i]    = Z; 

                t_m_vs_mu(thread_id)  ->Fill(mass, mu); 
                t_m_vs_Z(thread_id)   ->Fill(mass, Z); 
                t_m_vs_uCL(thread_id) ->Fill(mass, std::log10(mu_cl95)); 
                t_m_vs_e2CL(thread_id)->Fill(mass, std::log10(epsilon2_CL)); 

                t_pQ0(thread_id)      ->Fill(pQ0);
            } 
        }
        if (!fill_pts) {
            const auto& pts_mm = t_pts_mass(thread_id); 
            const auto& pts_Z = t_pts_Z(thread_id); 
            for (int i=0; i<n_fitpts; i++) {
                for (int j=0; j<n_fitpts; j++) {
                    t_Z_corr(thread_id)->Fill( pts_mm[i], pts_mm[j], pts_Z[i]*pts_Z[j] ); 
                }
            }
        }

    };

    task::scan_params(plist, ffcn, cfg);

    //now, run once and fill test points
    fill_pts = true; 
    cfg.n_threads = 1;
    cfg.n_scans = 1; 
    task::scan_params(plist, ffcn, cfg); 


    TCanvas* canv; 
    TGraph* g; 

    canv = new TCanvas;
    gStyle->SetOptStat(0); 
    t_m_vs_mu.aggregate()->Draw("col"); 

    canv = new TCanvas;
    auto h_m_vs_Z = t_m_vs_Z.aggregate(); 
    make_brazil_flag_plot(h_m_vs_Z, Form("Avg. of %i pseudo-spectra",max_scans)); 
    g = new TGraph(pts_m.size(), pts_m.data(), pts_Z.data());  
    g->Draw("SAME"); 

    canv = new TCanvas;
    auto h_m_vs_uCL = t_m_vs_uCL.aggregate(); 
    make_brazil_flag_plot(h_m_vs_uCL, Form("Avg. of %i pseudo-spectra",max_scans)); 
    g = new TGraph(pts_m.size(), pts_m.data(), pts_uCL.data());  
    g->Draw("SAME"); 

    canv = new TCanvas;
    canv->SetTopMargin(0.15);
    auto h_m_vs_e2CL = t_m_vs_e2CL.aggregate();
    h_m_vs_e2CL->SetTitle(Form("CL=0.95 upper limits on #varepsilon^{2}, %.1f x 10^{6} events;signal mass hypothesis (MeV);#epsilon^{2}, CL=0.95", config.total_stats/1e6));
    make_brazil_flag_plot(h_m_vs_e2CL, Form("Avg. of %i pseudo-spectra",max_scans)); 
    g = new TGraph(pts_m.size(), pts_m.data(), pts_eps2.data());  
    g->Draw("SAME"); 

    canv = new TCanvas; 
    auto h_pQ0 = t_pQ0.aggregate(); 
    h_pQ0->SetMaximum( h_pQ0->GetMaximum()*1.5 );
    h_pQ0->SetMinimum( 0. );  
    h_pQ0->Draw("HIST"); 

    canv = new TCanvas; 
    canv->SetRightMargin(0.15); 
    auto h_Z_corr = t_Z_corr.aggregate(); 
    h_Z_corr->Scale(1./((double)(n_fitpts)*std::pow(n_fitpts/n_bins,2))); 
    h_Z_corr->Draw("colz"); 

}



