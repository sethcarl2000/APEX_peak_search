
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TH1D.h> 
#include <TFile.h>
#include <TVector3.h> 
#include <TAxis.h> 
#include <TCanvas.h>
#include <TStyle.h>  
#include <TLegend.h> 
// stdlib
#include <cmath> 
#include <string> 
#include <vector> 

namespace {
    constexpr double me = 0.501; //MeV 

    constexpr double mass_plot_range[] = { 140., 280. }; 
    constexpr double mass_bin_size = 0.5; 

    constexpr int n_bins = (mass_plot_range[1] - mass_plot_range[0])/mass_bin_size; 

    constexpr double p_promt_cut = 0.4;
    constexpr double p_promt_min = 0.001; 
}

template<typename T> using rptr = ROOT::RDF::RResultPtr<T>; 

void create_prompt_and_accidental_spectra(const std::string& path_infile, const std::string& path_outfile="")
{

    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("track_data", path_infile); 

    ROOT::RDF::Experimental::AddProgressBar(df); 

    auto df_mass = df
        
        .Filter([](double p){ return p > p_promt_min; }, {"p_coinc"})

        .Define("mass", [](const TVector3& Pe, const TVector3& Pp){
    
            //energy of both leptons (lab-frame)
            double E = std::sqrt( Pe.Mag2() + me*me ) + std::sqrt( Pp.Mag2() + me*me ); 

            //spacelike momentum magnitude (again, lab-frame) 
            double P2 = ( Pe + Pp ).Mag2(); 

            return std::sqrt(E*E - P2); 

        }, {"P_e", "P_p"}); 


    


    auto h_all = df_mass.Histo1D({"h_all", ";inv. mass (MeV);", n_bins, mass_plot_range[0], mass_plot_range[1]}, "mass");

    auto h_p   = df_mass.Histo1D({"h_p",   ";P_{prompt};", 150, 0., 1.}, "p_coinc"); 

    auto df_prompt_hi = df_mass
        .Filter([](double p){ return p > p_promt_cut; }, {"p_coinc"});

    auto h_all_hi    = df_prompt_hi.Histo1D({"h_hi", ";inv. mass (MeV);", n_bins, mass_plot_range[0], mass_plot_range[1]}, "mass");
    auto h_prompt_hi = df_prompt_hi.Histo1D({"h_hi", ";inv. mass (MeV);", n_bins, mass_plot_range[0], mass_plot_range[1]}, "mass", "p_coinc");
    auto n_prompt_hi = df_prompt_hi.Sum("p_coinc"); 
    auto n_total_hi  = df_prompt_hi.Count();     


    auto df_prompt_lo = df_mass
        .Filter([](double p){ return p <= p_promt_cut; }, {"p_coinc"});

    auto h_all_lo    = df_prompt_lo.Histo1D({"h_lo", ";inv. mass (MeV);", n_bins, mass_plot_range[0], mass_plot_range[1]}, "mass");
    auto h_prompt_lo = df_prompt_lo.Histo1D({"h_lo", ";inv. mass (MeV);", n_bins, mass_plot_range[0], mass_plot_range[1]}, "mass", "p_coinc");
    auto n_prompt_lo = df_prompt_lo.Sum("p_coinc"); 
    auto n_total_lo  = df_prompt_lo.Count();     

    // this works in the following way: 

    // we divide all events kept (above p_prompt_min) into two categories: those above and below our threshold 'p_prompt_cut'. 
    // 
    // For any given mass bin m \in [m, m+dm), there are a certain number of 'prompt' (np) and 'accidental' (na) events. 
    // The number of events in our mass window above(below) the 'p_prompt_cut' threshold is labeled by N+(N-):  
    // N+ = (np+/(np+ + np-))*np   +   (na+/(na+ + na-))*na 
    // N- = (np-/(np+ + np-))*np   +   (na-/(na+ + na-))*na 
    //
    // We can solve this lin. equation to get a sense of the number of 'prompt' and 'accidental' events in our spectrum. 

    //so, let's count 'em. 
    auto xax = h_all->GetXaxis(); 

    auto h_prompt     = new TH1D("h_prompt",     "Reconstructed Prompt-only spectra",     n_bins, mass_plot_range[0], mass_plot_range[1]); 
    auto h_accidental = new TH1D("h_accidental", "Reconstructed Accidental-only spectra", n_bins, mass_plot_range[0], mass_plot_range[1]); 

    for (int bin=1; bin<=n_bins; bin++) {

        double np_hi = h_prompt_hi->GetBinContent(bin); 
        double na_hi = h_all_hi->GetBinContent(bin) - np_hi; 

        double np_lo = h_prompt_lo->GetBinContent(bin); 
        double na_lo = h_all_lo->GetBinContent(bin) - np_lo; 

        //invert the linear system described above 
        double 
            a{np_hi/(np_hi + np_lo)}, b{na_hi/(na_hi + na_lo)}, 
            c{np_lo/(np_hi + np_lo)}, d{na_lo/(na_hi + na_lo)}; 

        double det = a*d - (b*c); 
        
        double mass = xax->GetBinCenter(bin); 

        double N_hi = h_all_hi->GetBinContent(bin); 
        double N_lo = h_all_lo->GetBinContent(bin); 

        double Np = (+d*N_hi - c*N_lo)/det; 
        double Na = (-b*N_hi + a*N_lo)/det; 

        std::printf(
            "~mass: %.1f ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
            "Np: %+.3e = +(%+.3e)%+.3e - (%+.3e)%+.3e\n"
            "Na: %+.3e = -(%+.3e)%+.3e + (%+.3e)%+.3e\n",
            Np, (d/det),N_hi, (c/det),N_lo, 
            Na, (b/det),N_hi, (a/det),N_lo 
        ); 

        h_prompt    ->Fill( mass, Np ); 
        h_accidental->Fill( mass, Na ); 
    }

    new TCanvas; 
    h_prompt->DrawCopy("HIST"); 

    new TCanvas;
    h_accidental->DrawCopy("HIST"); 

    new TCanvas; 
    auto legend = new TLegend; 
    h_accidental->SetTitle("Reconstructed spectra;inv. mass (MeV);"); 
    h_accidental->SetFillColor(kBlue); 
    h_accidental->SetFillStyle(3004); 
    legend->AddEntry(h_prompt->DrawCopy("HIST"), "Prompt");  
    legend->AddEntry(h_accidental->DrawCopy("SAME HIST"), "Accidental");  
    legend->Draw(); 

    new TCanvas; 
    h_p->DrawCopy(); 

    new TCanvas; 
    h_prompt_hi->SetFillColor(kBlue); 
    h_prompt_hi->SetFillStyle(3004); 
    h_prompt_hi->DrawCopy(); 
    h_prompt_lo->DrawCopy("SAME"); 

    return; 
}