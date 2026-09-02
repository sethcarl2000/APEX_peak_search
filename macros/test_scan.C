
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

//void fit_window(peak_search::FitTestThreadManager* mgr);
///________________________________________________________________________________________________________
void test_scan()
{
    using namespace peak_search; 

    auto& kernel = FitTestKernel::Instance(); 

    const double min_mass = 155.; 
    const double max_mass = 265.; 

    int n_steps = 100; 
    int n_bins  = 100; 

    auto h_m_vs_mu = new TH2D(
        "h_signal", "Best-fit signal parameter '#mu' vs m;signal mass hypothesis (MeV);best-fit #mu", 
        n_bins, min_mass, max_mass,
        150, -40e-3, 40e3
    ); 

    /*
    const double max_signal_significance = 6.; 
    fModel_m_vs_Z = ROOT::RDF::TH2DModel{"h_Z", "Significance Z ~ #sqrt{Q0} vs m;signal mass hypothesis (MeV);Significance Z (n. #sigma)", 
        (int)fN_steps/4, fMinFitMass, fMaxFitMass, 
        150, -max_signal_significance, +max_signal_significance
    }; 

    fModel_m_vs_pQ0 = ROOT::RDF::TH2DModel{"h_pZ", "p(Q0) vs m;signal mass hypothesis (MeV);p(Q0)",  
        (int)fN_steps/4, fMinFitMass, fMaxFitMass, 
        50, 0., 1.
    }; 

    fModel_pQ0 = ROOT::RDF::TH1DModel{"h_pQ0", "p(Q0);p(Q0);", 50, 0., 1.}; 
*/ 
    kernel.SetMassRange(min_mass, max_mass); 
    kernel.SetNSteps(n_steps); 

    kernel.AddTH2D(*h_m_vs_mu); 

    auto fit_window_fcn = static_cast<FitTestFunction>([h_m_vs_mu](peak_search::FitTestThreadManager* mgr)
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

        auto h_m_vs_mu_t = mgr->GetUserTH2D(h_m_vs_mu); 

        h_m_vs_mu->Fill(mass, mu); 
        return;
    });

    kernel.RunTest(20, fit_window_fcn); 

    new TCanvas;
    h_m_vs_mu->Draw("col"); 

}
///________________________________________________________________________________________________________



