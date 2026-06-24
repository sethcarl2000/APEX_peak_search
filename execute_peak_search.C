#include "bootstrap_sample_hist.C"
#include "macros/scan_masses.C"
#include "bininfo.hpp"
#include "numbers.hpp"

#include <TH2D.h>
#include <TH1D.h>
#include <TFile.h> 
#include <TRandom3.h> 
#include <TString.h> 
#include <TGraph.h> 
#include <TPad.h> 
#include <Math/QuantFuncMathCore.h>
#include <TLine.h> 

constexpr bool yes=true;

/// @brief execute computation of Q0 statistic for APEX m^2 spectrum 
/// @param ARMED if TRUE, the FULL UNBLINDED STATISTICS WILL BE USED
int execute_peak_search(bool ARMED=false)
{
    double frac=0.10; 

    auto file = new TFile("histo_all.root", "READ");

    auto h_m_p_FULL = (TH2D*)file->Get("h_m_p"); 

    TRandom3 rand; 
    
    int bg_poly_degree = 3;

    double m_min = 140.;
    double m_max = 265.; 
    int n_tests = 250; 
    
    //the *global* significance corresponding to 1-sigma
    double p_sigma[5]; 
    for (int i=0; i<5; i++) {
        p_sigma[i] = ROOT::Math::normal_cdf_c((double)i+1, 1., 0.); 
    }
    
    double dm = (m_max - m_min)/((double)n_tests-1);

    auto hist_p_chi2 = new TH1D("h_p_chi2", Form("p(#chi^2) for %i-degree polynomial bg (blind sample size: %.1f%%);p(#chi^2);",bg_poly_degree, frac*100.),
        50, -0.25, 1.25 
    );

    auto hist_p_Q0 = new TH1D("h_p_Q0", Form("p(q_{0}) distribution for %i-degree polynomial bg (blind sample size: %.1f%%);p(#chi^2);",bg_poly_degree, frac*100.),
        50, -0.25, 1.25 
    );



    //first, trim the histogram. 
    auto xax = h_m_p_FULL->GetXaxis();
    auto yax = h_m_p_FULL->GetYaxis();
    
    int max_bin = xax->FindBin(280.); 
    int min_bin = xax->FindBin(140.);
    double dx = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()); 

    double xmin = xax->GetBinCenter(min_bin)-dx/2.;
    double xmax = xax->GetBinCenter(max_bin)+dx/2.;
        

    auto h_m_p_trim = new TH2D("h_m_p_trim", h_m_p_FULL->GetTitle(), 
        max_bin-min_bin+1, xmin, xmax, 
        yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
    );

    std::vector<m_and_Q0_t> ARMED_result; ARMED_result.reserve(n_tests); 

    int n_experiments = ARMED ? 1 : 500; 

    for (int i=0; i<n_experiments; i++) {

        auto h_m_p = ARMED ? bootstrap_sample_hist(h_m_p_FULL, frac, rand) : bootstrap_sample_hist(h_m_p_FULL, frac, rand); 
        h_m_p->SetDirectory(0);
        
        h_m_p_trim->Reset(); 

        //h_m_p_trim.SetDirectory(0);
        for (int ix=min_bin; ix<=max_bin; ix++) {
            for (int iy=1; iy<=yax->GetNbins(); iy++) {
                h_m_p_trim->Fill( 
                    xax->GetBinCenter(ix), 
                    yax->GetBinCenter(iy), 
                    h_m_p->GetBinContent(ix,iy) 
                );
            }
        }

        std::printf(" test %i total statis: %f\n", i, h_m_p_trim->Integral());

        auto results = scan_masses(h_m_p_trim, 300, 140., 265., 0.3, bg_poly_degree); 

        if (ARMED==false) {
            for (auto& result : results) hist_p_chi2->Fill( result.p_chi2 ); 
            for (auto& result : results) hist_p_Q0  ->Fill( result.p_Q0 );
        } else {
            //ARMED_result = results;  
            for (auto result : results) { 
                if (!peak_search::numbers::is_nan(result.p_Q0)) ARMED_result.push_back(result);
                std::printf("my result: %f \n", result.p_Q0 ); 
            }
        }

        delete h_m_p; 
    }

    std::vector<double> m, p_Q0; 
    for (const auto& result : ARMED_result) { m.push_back(result.m); p_Q0.push_back(result.p_Q0); }

    if (ARMED==false) {
        new TCanvas; 
        hist_p_chi2->Draw("col"); 

        new TCanvas; 
        hist_p_Q0->Draw("col"); 
    } else {

        new TCanvas; 

        gPad->SetLogy(1);

        auto g = new TGraph(n_tests, m.data(), p_Q0.data()); 
        g->SetTitle(Form("p-value for mass hypotheses (%.1f%% statistical sample);mass (MeV);p(Q0);",frac*100.)); 
        g->SetMarkerStyle(kOpenCircle);
        g->SetMarkerSize(0.2); 
        g->Draw();
        g->Draw("SAME, P");
        
        //draw lines corresponding to different sigma-values
        auto line = new TLine; 
        line->SetLineStyle(kDashDotted);
        line->SetLineColor(kRed); 
        for (int i=0; i<5; i++) {
            line->DrawLine(m.front(),p_sigma[i], m.back(),p_sigma[i]); 
        }

        //draw global significance lines
        line->SetLineStyle(kSolid);
        line->SetLineColor(kRed); 
        for (int i=0; i<5; i++) {
            line->DrawLine(m.front(),p_sigma[i]/139., m.back(),p_sigma[i]/139.); 
        }
    }

    return 0; 
}