
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
#include <compute_epsilon2.hpp>
#include <make_histogram_copy.hpp>
#include <compute_Q0.hpp>
#include <compute_statistics.hpp> 
// ROOT headers
#include <TCanvas.h> 
#include <TH2D.h> 
#include <TStyle.h>
#include <TVector3.h> 
#include <TAxis.h> 
#include <TGraphErrors.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <Math/QuantFuncMathCore.h> 
#include <Math/ProbFuncMathCore.h> 
#include <TLegend.h> 
#include <TLine.h> 
// stdlib
#include <vector> 
#include <cstdio> 
#include <memory> 

/// @brief Returns estimate of mass resolution for the given mass hypothesis (this is a slightly conservative over-estimate)
/// @param mass_hypothesis mass hypothesis (MeV)
/// @return estimate of mass resolution (MeV)
double mass_resolution(double mass_hypothesis)
{
    return 1. + (mass_hypothesis - 140.) * ((0.8 - 1.0)/(270 - 140)); 
}

#define MAX_EVENTS 0

///________________________________________________________________________________________________________
void do_peak_search(const std::string& path_input="replay-jul-9.root", ULong64_t max_events=0)
{
    using ROOT::Math::normal_quantile; 
    using ROOT::Math::normal_cdf; 

    using namespace peak_search; 

    const double min_plot_mass = 140.; 
    const double max_plot_mass = 280.; 

    const double dm = 0.25; 
    const int n_bins = (max_plot_mass - min_plot_mass)/dm; 

    const double min_fit_mass = 155.; 
    const double max_fit_mass = 270.; 

    //number of bins on either edge to fit to 
    const double mass_window_size = 7.; 

    const int n_steps = ((max_fit_mass - min_fit_mass)/dm) +1; 

    if (MAX_EVENTS > 0) { 
        if (ROOT::IsImplicitMTEnabled()) ROOT::DisableImplicitMT(); 
    } else {
        ROOT::EnableImplicitMT(); 
    }

    double Z_CL = normal_quantile(0.95, 1.); 
    double Z_power_constrained_limit = 1.; 
    
    //approximate correction for the 'look elsewhere' effect
    double n_regions = (max_fit_mass - min_fit_mass)/0.8; 

    ROOT::RDataFrame df("track_data", path_input); 

    auto df_prompt = df
#if MAX_EVENTS > 0
        .Range(0, MAX_EVENTS)
#endif
        .Define("p_prompt", [](double p){ return p; }, {"p_coinc"})
        .Define("mass", [](const TVector3& Pp, const TVector3& Pe){
  
            double E = std::sqrt( 0.501*0.501 + Pe.Mag2() ) + std::sqrt( 0.501*0.501 + Pp.Mag2() ); 
            double P2 = ( Pe + Pp ).Mag2(); 

            return std::sqrt( E*E - P2 ); 
        }, {"P_p", "P_e"});  


    auto h_spectrum = df_prompt.Histo1D({"h_prompt", "X #rightarrow e^{#pm} Invariant Mass Spectrum;e^{#pm} inv. mass (MeV);counts",
        n_bins, min_plot_mass, max_plot_mass
    }, "mass", "p_prompt"); 

    ROOT::RDF::Experimental::AddProgressBar(df_prompt); 

    std::vector<double> 
        vec_m(n_bins, 0.), 
        vec_0(n_bins, 0.),
        vec_mu(n_bins, 0.),
        vec_PCUL(n_bins, 0.), 
        vec_1s(n_bins, 0.),
        vec_2s(n_bins, 0.), 

        vec_pQ0(n_bins, 0.), 

        vec_epsCL(n_bins, 0.),
        vec_epsCL_1s(n_bins, 0.), 
        vec_epsCL_2s(n_bins, 0.), 
        vec_eps_mu(n_bins, 0.); 

    auto h_m = dynamic_cast<TH1D*>(h_spectrum->Clone("h_spec_cpy")); 
    auto xax = h_m->GetXaxis(); 

    std::cout << "Executing peak-search...\n"; 

    double mu_max =0.; 

    double eps_min =1e6;
    double eps_max =0.;
    double min_p =1.; 

    for (int i=0; i<n_steps; i++) {

        double mass = min_fit_mass + dm*((double)i); 

        vec_m[i] = mass; 

        int center_bin = xax->FindBin(mass); 

        double sigma_m = mass_resolution(mass); 

        Histo1D window = make_histogram_copy(h_m, 
            mass - ((double)mass_window_size)*sigma_m, 
            mass + ((double)mass_window_size)*sigma_m
        ); 

        /*std::printf("mass: %.1f m-range: [%.1f, %.1f]\n", mass, 
            window.GetXmin(), 
            window.GetXmin()
        );*/ 

        auto gaussian_fcn = peak_search::Gauss(0, mass, sigma_m); 

        //fit the background
        auto background_poly = peak_search::fit_exponential_poly(window, 6).data; 
        
        auto stat_result = peak_search::compute_statistics(window, gaussian_fcn, background_poly, mass, 0.05, 1.); 

        if (stat_result.status != peak_search::Status::kSuccess) {
            Warning("fit_window_function", "Fit failed for mass: %.1f", mass); 
            continue; 
        }

        auto stats = stat_result.data; 
        double Q0 = stats.Q0; 
        double mu = stats.mu_MLE; 
        double mu_cl95 = stats.mu_CL; 
        double epsilon2_CL = stats.epsilon2_CL; 

        //double mu_sigma = stats.mu_sigma; 
        vec_mu[i] = mu; 

        double mu_PCUL = std::max(mu + Z_CL*stats.mu_sigma, stats.mu_sigma*Z_power_constrained_limit); 


        // PCUL with power of 1-sigma 
        vec_PCUL[i] = mu_PCUL; 

        mu_max = std::max( mu_max, mu_PCUL ); 
        mu_max = std::max( mu_max, std::fabs(mu) ); 
        mu_max = std::max( mu_max, 2.*stats.mu_sigma ); 
        
        vec_1s[i] = 1.*stats.mu_sigma; 
        vec_2s[i] = 2.*stats.mu_sigma; 

        double Z  = (Q0<0.?-1:+1) * std::sqrt(std::fabs(Q0)); 

        double pQ0 = compute_Q0_p(Q0); 

        vec_pQ0[i] = pQ0; 
        
        min_p = std::min( min_p, pQ0 ); 

        vec_eps_mu[i] = mu > 0. ? compute_epsilon2(background_poly, mu, mass, 1.) : 1e-3; 
        vec_epsCL[i] = compute_epsilon2(background_poly, mu_PCUL, mass, 1.); 
        vec_epsCL_1s[i] = compute_epsilon2(background_poly, 1.*stats.mu_sigma + Z_CL*stats.mu_sigma, mass, 1.); 
        vec_epsCL_2s[i] = compute_epsilon2(background_poly, 2.*stats.mu_sigma + Z_CL*stats.mu_sigma, mass, 1.); 
    
        eps_min = std::min( eps_min, vec_epsCL[i] ); 
        eps_max = std::max( eps_max, vec_epsCL[i] ); 
        
    }
       
    std::vector<double> vec_eps_min(n_bins, eps_min*10e-4); 

    new TCanvas; 
    h_m->DrawCopy(); 
    std::printf("total prompt stats: %.0f\n", h_m->Integral()); 

    new TCanvas;    
    gStyle->SetOptStat(0); 
    //gPad->DrawFrame(min_fit_mass-10.,-1.2*mu_max,  max_fit_mass+10.,+1.2*mu_max); 
    auto h_frame = new TH2D("frame_mu", "X #rightarrow e^{#pm} spectrum fit;inv. mass (MeV); #varepsilon^{2}", 
        100, min_fit_mass-10., max_fit_mass+10.,
        100, -1.2*mu_max, +1.2*mu_max
    );
    h_frame->Draw(); 

    auto g_2s = new TGraphErrors(n_steps, vec_m.data(), vec_0.data(), nullptr, vec_2s.data()); 
    g_2s->SetTitle("X #rightarrow e^{#pm} spectrum fit;inv. mass (MeV); #varepsilon^{2}"); 
    g_2s->SetFillColor(kYellow); 
    g_2s->SetLineColor(kYellow); 
    g_2s->Draw("3");

    auto g_1s = new TGraphErrors(n_steps, vec_m.data(), vec_0.data(), nullptr, vec_1s.data()); 
    g_1s->SetFillColor(kGreen); 
    g_2s->SetLineColor(kGreen); 
    g_1s->Draw("3"); 

    auto g_mu = new TGraph(n_steps, vec_m.data(), vec_mu.data()); 
    g_mu->Draw("SAME");

    auto g_PCUL = new TGraph(n_steps, vec_m.data(), vec_PCUL.data()); 
    g_PCUL->SetLineWidth(2);
    g_PCUL->Draw("SAME"); 

    auto line = new TLine(vec_m.front(),0., vec_m.back(),0.); 
    line->SetLineStyle(kDashed); 
    line->Draw(); 

    new TCanvas;
    h_frame = new TH2D("frame_eps", "X #rightarrow e^{#pm} #varepsilon^{2} CL=0.95 upper limit;inv. mass (MeV); #varepsilon^{2} (best-fit)", 
        100, min_fit_mass-10., max_fit_mass+10.,
        100, eps_min/3., eps_max*3.
    );
    gPad->SetLogy(1);
    h_frame->Draw(); 

    auto ge_mu = new TGraph(n_steps, vec_m.data(), vec_epsCL.data()); 
    ge_mu->Draw("SAME");

    auto ge_PCUL = new TGraph(n_steps, vec_m.data(), vec_epsCL.data()); 
    ge_PCUL->SetLineWidth(2);
    ge_PCUL->Draw("SAME"); 

    new TCanvas;
    auto gp = new TGraph(n_steps, vec_m.data(), vec_pQ0.data()); 
    gPad->SetLogy(1);
    gp->SetTitle("Local p(Q_{0}) significance;inv. mass (MeV);p(Q0)"); 
    gp->Draw(); 
    return; 

}



