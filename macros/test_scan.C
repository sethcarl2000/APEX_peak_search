
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
// ROOT headers
#include <TCanvas.h> 
#include <TH2D.h> 
#include <TStyle.h>

//void fit_window(peak_search::FitTestThreadManager* mgr);
///________________________________________________________________________________________________________
void test_scan()
{
    using namespace peak_search; 

    auto& kernel = FitTestKernel::Instance(); 

    const double min_mass = 155.; 
    const double max_mass = 265.; 

    int n_steps = 10; 
    int n_bins  = 100; 

    auto h_m_vs_mu = new TH2D(
        "h_signal", "Best-fit signal parameter '#mu' vs m;signal mass hypothesis (MeV);best-fit #mu", 
        n_steps, min_mass, max_mass,
        100, -40e3, 40e3
    ); 

    auto h_m_vs_uCL = new TH2D(
        "h_uCL", "Signal parameter upper-limit '#mu_{>0.95}' vs m;signal mass hypothesis (MeV);log_{10} #mu_{>0.95}", 
        n_steps, min_mass, max_mass,
        100, -2, 6
    ); 

    auto h_m_vs_Z = new TH2D(
        "h_Z", "Significance Z ~ #sqrt{Q0} vs m;signal mass hypothesis (MeV);Significance Z (n. #sigma)",
        n_steps, min_mass, max_mass,
        100, -7, 7
    ); 

    auto h_pQ0 = new TH1D(
        "h_pZ", "p(Q0) vs m;signal mass hypothesis (MeV);p(Q0)",  
        50, 0., 1.
    ); 

    kernel.SetMassRange(min_mass, max_mass); 
    kernel.SetNSteps(n_steps); 
    kernel.SetMassBinSize(0.5); 

    auto ptr_m_vs_mu    = kernel.AddTH2D(*h_m_vs_mu); 
    auto ptr_m_vs_Z     = kernel.AddTH2D(*h_m_vs_Z); 
    auto ptr_m_vs_uCL   = kernel.AddTH2D(*h_m_vs_uCL); 

    auto ptr_pQ0        = kernel.AddTH1D(*h_pQ0); 

    auto fit_window_fcn = static_cast<FitTestFunction>([ptr_m_vs_mu, ptr_m_vs_Z, ptr_pQ0, ptr_m_vs_uCL](peak_search::FitTestThreadManager* mgr)
    {
        double window_size = 7.; // MeV 

        auto mass = mgr->get_mass(); 

        peak_search::Histo1D spectrum = mgr->get_spectrum(mass - window_size, mass + window_size);

        auto gaussian_fcn = peak_search::Gauss(0, mass, 1.); 

        //fit the background
        auto background_poly = peak_search::fit_exponential_poly(spectrum, 6).data; 

        peak_search::FcnSum f_sum(&gaussian_fcn, &background_poly); 

        double Q0 = peak_search::compute_Q0(spectrum, f_sum); 

        double Z  = (Q0<0.?-1:+1) * std::sqrt(std::fabs(Q0)); 
 
        double mu = f_sum.GetParams()[0];

        double pQ0 = compute_Q0_p(Q0); 

        //get the middle(-ish)bin. this gives us an order-of-magnitude estimate for the natural variance of the signal paramter, mu.
        double N_middle = spectrum.bins.at( spectrum.GetNbins()/2 ).N; 

        //also, find the CLs value
        // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%M<>       
        //      - muon's comment (2 sep 2026)
        double mu_CL95 = peak_search::solve_for_CLs(spectrum, f_sum, 0.95, 0.01, 5e3);

        mgr->GetUserTH2D(ptr_m_vs_mu) ->Fill(mass, mu); 
        mgr->GetUserTH2D(ptr_m_vs_Z)  ->Fill(mass, Z); 
        mgr->GetUserTH2D(ptr_m_vs_uCL)->Fill(mass, std::log10(mu_CL95)); 

        mgr->GetUserTH1D(ptr_pQ0)    ->Fill(pQ0); 
        return;
    });

    kernel.RunTest(200, fit_window_fcn); 

    new TCanvas;
    gStyle->SetOptStat(0); 

    h_m_vs_mu->Draw("col"); 

    new TCanvas;
    h_m_vs_Z->Draw("col"); 

    new TCanvas;
    h_m_vs_uCL->Draw("col"); 

    new TCanvas;
    h_pQ0->SetMaximum( h_pQ0->GetMaximum()*1.5 );
    h_pQ0->SetMinimum( 0. );  
    h_pQ0->Draw("E"); 

}
///________________________________________________________________________________________________________



