//ROOT
#include <TH1D.h>
#include <ROOT/RDataFrame.hxx>
#include <TGaxis.h> 
#include <TStyle.h>
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TAxis.h> 
#include <TGraph.h> 
#include <TVector3.h>
#include <Math/Vector4D.h> 
#include <TString.h> 
//stdlib
#include <string>
#include <vector> 
#include <cstdio> 
#include <cmath> 

namespace {
    constexpr double me2 = 0.511*0.511; 

    constexpr double m_min = 140.; //MeV
    constexpr double m_max = 280.; //MeV 

    constexpr double dm = 0.500; //MeV 
}

/// @brief Draw cumulative stats for full replay, starting at p_coinc = 100%, and moving down towards p_coinc = 0%.
int draw_coinc_and_accidental(std::string path_infile, double threshold_p_accidental=0.25, double threshold_p_coinc=0.95)
{
    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("track_data", path_infile);

    using FourVec = ROOT::Math::XYZTVector; 

    auto df_invmass = df    

        .Define("invariant_mass", [](const TVector3& Pp, const TVector3& Pe)
        { 
            
            double energy = std::sqrt( Pp.Mag2() + me2 ) + std::sqrt( Pe.Mag2() + me2 ); 

            TVector3 P = Pp + Pe; 

            //return (time-like) mass of Pe + PP, using (+---) metric 
            return std::sqrt( energy*energy - P.Mag2() ); 
        }, {"P_p", "P_e"}); 

    const int n_bins = std::floor((m_max - m_min)/dm); 

    auto hist_coinc = (TH1D*)df_invmass

        .Filter([threshold_p_coinc](double p){ return p > threshold_p_coinc; }, {"p_coinc"})

        .Histo1D<double>({"hh_coinc", "Invariant mass;m_{#pm} (MeV)", n_bins, m_min, m_max}, "invariant_mass")->Clone("h_coinc"); 

    
    auto hist_accidental = (TH1D*)df_invmass

        .Filter([threshold_p_accidental](double p){ return p < threshold_p_accidental; }, {"p_coinc"})

        .Histo1D<double>({"hh_accidental", "Invariant mass;m_{#pm} (MeV)", n_bins, m_min, m_max}, "invariant_mass")->Clone("h_accidental"); 

    
    //get *actual* purity of sample (for both)
    auto h_purity_coinc =  
        df.Filter([threshold_p_coinc](double p){ return p > threshold_p_coinc; }, {"p_coinc"})
        .Histo1D({"hhh_c", "", 200, 0., 1.}, "p_coinc");
    
    auto get_coinc_fraction = [](ROOT::RDF::RNode node)
    {
        auto hist_frac = (TH1D*)node.Histo1D<double>({"hhh", "", 200, 0., 1.}, "p_coinc")->Clone("hhh_cpy"); 
        hist_frac->SetDirectory(0); 

        double sum = 0.; 
        double coinc = 0.; 
        for (int ib=1; ib<=200; ib++) {
            sum += hist_frac->GetBinContent(ib);
            coinc += hist_frac->GetBinContent(ib) * hist_frac->GetXaxis()->GetBinCenter(ib);
        }
        delete hist_frac; 
        
        return coinc / sum;
    };

    double coinc_frac_coinc = get_coinc_fraction(
        df.Filter([threshold_p_coinc](double p){ return p > threshold_p_coinc; }, {"p_coinc"})
    );

    double coinc_frac_accidental = get_coinc_fraction(
        df.Filter([threshold_p_accidental](double p){ return p < threshold_p_accidental; }, {"p_coinc"})
    );

    new TCanvas; 
    gStyle->SetOptStat(0); 
    
    hist_coinc->SetLineColor(kRed);
    hist_coinc->SetFillColor(kRed);
    hist_coinc->SetFillStyle(3004);
    hist_coinc->Draw();
    std::printf("stats in coinc histogram: %.4e\n", hist_coinc->Integral()); 

    hist_accidental->SetLineColor(kBlack);
    hist_accidental->SetFillStyle(0);
    hist_accidental->Draw("SAME");
    std::printf("stats in accidental histogram: %.4e\n", hist_accidental->Integral()); 

    auto legend = new TLegend; 
    legend->SetHeader("frac. of coinc events (p)");
    legend->AddEntry(hist_coinc,      Form("p = %.2f",coinc_frac_coinc));
    legend->AddEntry(hist_accidental, Form("p = %.2f",coinc_frac_accidental));
    legend->Draw(); 
        

    return 0; 
};