
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

/// @brief Returns estimate of mass resolution for the given mass hypothesis (this is a slightly conservative over-estimate)
/// @param mass_hypothesis mass hypothesis (MeV)
/// @return estimate of mass resolution (MeV)
double mass_resolution(double mass_hypothesis)
{
    return 1. + (mass_hypothesis - 140.) * ((0.8 - 1.0)/(270 - 140)); 
}

void make_brazil_flag_plot(TH2D* hist, double cl_1=0.34134475, double cl_2=0.47724987)
{
    //go through each bin, and find the cumulative stats corresponding to each cl given. 
    const std::vector<double> levels{ 0.5-cl_2, 0.5-cl_1, 0.5, 0.5+cl_1, 0.5+cl_2 };
    
    auto x_ax = hist->GetXaxis(); const int n_bins_x = x_ax->GetNbins(); 
    auto y_ax = hist->GetYaxis(); const int n_bins_y = y_ax->GetNbins(); 

    const int n_bins = x_ax->GetNbins(); 

    std::vector<double> x; x.reserve(n_bins_x); 
    std::vector<double> y_cl1, y_err_cl1; y_cl1.reserve(n_bins_x); y_err_cl1.reserve(n_bins_x); 
    std::vector<double> y_cl2, y_err_cl2; y_cl2.reserve(n_bins_x); y_err_cl2.reserve(n_bins_x); 
    std::vector<double> y_med; y_med.reserve(n_bins_x); 

    double minval{+1e30}, maxval{-1e30}; 

    for (int bx=1; bx<=n_bins_x; bx++) {

        double integral =0.; 
        for (int by=1; by<=n_bins_y; by++) integral += hist->GetBinContent(bx,by); 
        
        auto find_cumulant = [hist,bx,y_ax,n_bins_y,integral](double p) {

            double cum=0.; 


            double bin_val; 
            for (int by=1; by<=n_bins_y; by++) {

                bin_val = hist->GetBinContent(bx,by) / integral; 
                
                if (cum + bin_val >= p) {
                    
                    double val = y_ax->GetBinCenter(by-1);
                    
                    double remainder = p - cum; 
                    val += ((remainder/bin_val) * y_ax->GetBinWidth(by)); 
                    
                    return val; 

                } else {
                    cum += bin_val; 
                } 
            }
            
            return -1.; 
        };

        double y1_lo = find_cumulant(0.5-cl_1);
        double y1_hi = find_cumulant(0.5+cl_1);
        
        double y2_lo = find_cumulant(0.5-cl_2);
        double y2_hi = find_cumulant(0.5+cl_2);

        minval = std::min(y2_lo, minval); 
        maxval = std::max(y2_hi, maxval); 

        double y_median = find_cumulant(0.5); 

        x.push_back(x_ax->GetBinCenter(bx)); 

        y_cl1    .emplace_back((y1_hi + y1_lo)/2.);
        y_err_cl1.emplace_back((y1_hi - y1_lo));

        y_cl2    .emplace_back((y2_hi + y2_lo)/2.);
        y_err_cl2.emplace_back((y2_hi - y2_lo));

        y_med.emplace_back(y_median); 
    }

    /*if (!gPad) new TCanvas;
    
    double x_span = x_ax->GetXmax() - x_ax->GetXmin(); 
    gPad->DrawFrame(
        x_ax->GetXmin() - 0.1*x_span,
        minval - 0.1*(maxval-minval), 
        x_ax->GetXmax() + 0.1*x_span,
        maxval + 0.1*(maxval-minval)
    );*/ 

    auto g2 = new TGraphErrors(n_bins_x, x.data(), y_cl2.data(), nullptr, y_err_cl2.data()); 
    
    g2->SetTitle("CL_{s} 0.95 Upper limit on Dark Photon Coupling #varepsilon^{2} (76 M coinc events);mass hypothesis (MeV);#varepsilon^{2} (0.95 CL_{s} upper-limit)"); 
    
    g2->SetFillColor(kYellow);
    g2->Draw("A3"); 

    auto g1 = new TGraphErrors(n_bins_x, x.data(), y_cl1.data(), nullptr, y_err_cl1.data()); 
    g1->SetFillColor(kGreen);
    g1->Draw("3"); 

    auto gmed = new TGraph(n_bins_x, x.data(), y_med.data()); 
    gmed->Draw("SAME"); 

    auto legend = new TLegend; 
    legend->AddEntry(g1, "#pm 1 #sigma");
    legend->AddEntry(g2, "#pm 2 #sigma");
    legend->AddEntry(gmed, "median");
    legend->Draw(); 

    return; 
}


//void fit_window(peak_search::FitTestThreadManager* mgr);
///________________________________________________________________________________________________________
void test_scan()
{
    using namespace peak_search; 

    auto& kernel = FitTestKernel::Instance(); 

    const double min_mass = 150.; 
    const double max_mass = 270.; 

    int n_steps = 400; 
    int n_bins  = n_steps/4; 

    //pick a reasonable number of bins

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

    kernel.SetMassRange(min_mass, max_mass); 
    kernel.SetNSteps(n_steps); 
    kernel.SetMassBinSize(0.5); 

    auto ptr_m_vs_mu    = kernel.AddTH2D(*h_m_vs_mu); 
    auto ptr_m_vs_Z     = kernel.AddTH2D(*h_m_vs_Z); 
    auto ptr_m_vs_uCL   = kernel.AddTH2D(*h_m_vs_uCL); 
    auto ptr_m_vs_e2CL  = kernel.AddTH2D(*h_m_vs_e2CL); 

    auto ptr_pQ0        = kernel.AddTH1D(*h_pQ0); 

    auto fit_window_fcn = static_cast<FitTestFunction>([ptr_m_vs_mu, ptr_m_vs_Z, ptr_pQ0, ptr_m_vs_uCL, ptr_m_vs_e2CL](peak_search::FitTestThreadManager* mgr)
    {
        double window_size = 7.; // MeV 

        auto mass = mgr->get_mass(); 
        double resolution = mass_resolution(mass); 

        peak_search::Histo1D spectrum = mgr->get_spectrum(mass - window_size*resolution, mass + window_size*resolution);

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

        mgr->GetUserTH2D(ptr_m_vs_mu)  ->Fill(mass, mu); 
        mgr->GetUserTH2D(ptr_m_vs_Z)   ->Fill(mass, Z); 
        mgr->GetUserTH2D(ptr_m_vs_uCL) ->Fill(mass, std::log10(mu_cl95)); 
        mgr->GetUserTH2D(ptr_m_vs_e2CL)->Fill(mass, std::log10(epsilon2_CL)); 

        mgr->GetUserTH1D(ptr_pQ0)    ->Fill(pQ0); 
        return;
    });

    kernel.SetTotalStats(100e6);
    kernel.RunTest(200, fit_window_fcn); 

    new TCanvas; 
    make_brazil_flag_plot(h_m_vs_e2CL); 
    return; 

    new TCanvas;
    gStyle->SetOptStat(0); 

    h_m_vs_mu->Draw("col"); 

    new TCanvas;
    h_m_vs_Z->Draw("col"); 

    new TCanvas;
    h_m_vs_uCL->Draw("col"); 

    new TCanvas;
    h_m_vs_e2CL->Draw("col"); 

    new TCanvas;
    h_pQ0->SetMaximum( h_pQ0->GetMaximum()*1.5 );
    h_pQ0->SetMinimum( 0. );  
    h_pQ0->Draw("E"); 
}
///________________________________________________________________________________________________________



