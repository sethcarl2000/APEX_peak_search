#ifndef get_Pprompt_function_h
#define get_Pprompt_function_h

// ROOT headers
#include <TF1.h> 
#include <TH2D.h> 
#include <TH1D.h> 
#include <TAxis.h> 
#include <TFitResult.h>
#include <TFitResultPtr.h> 
#include <TPad.h> 
// stdlib headers
#include <functional> 
#include <vector> 
#include <cmath> 
#include <memory> 
#include <cstdio>


struct PpromptCalculator {
private: 
    friend PpromptCalculator get_PpromptCalculator(TH2D*, double, double, double, double);
    struct Point { double signal, background; }; 
    std::vector<Point> points;
    double m_min{0.}, dm{0.}; 
    static constexpr double norm_const = 0.398942280401; // 1 / sqrt(2*pi)
    
    double t_mean{0.}, t_sigma{1.};

    PpromptCalculator() = default; 

public: 

    //default copy ctor 
    PpromptCalculator(const PpromptCalculator&) = default; 

    /// @brief Compute 'P_prompt' for given mass and t
    /// @param mass mass value 
    /// @param t t (T_R - T_L value)
    /// @return estimated value of 'P_prompt' 
    double operator()(double mass, double t) const {
        const int max_index = points.size()-1; 
        if (max_index < 0) return 0.; 
        int i0 = std::min( (int)std::max( 0., (mass - m_min)/dm ), max_index );
        int i1 = std::min( max_index, i0+1 );
        const auto& pt0 = points.at(i0);
        const auto& pt1 = points.at(i1); 
        double modulo = mass/dm - std::floor( mass/dm ); 
        
        double signal     = pt0.signal     + (pt1.signal     - pt0.signal)*modulo; 
        double background = pt0.background + (pt1.background - pt0.background)*modulo; 
        
        
        double arg = (t - t_mean)/t_sigma; 

        double s = signal * norm_const * std::exp( -arg*arg/2. ) / t_sigma; //both of these have units of MeV^-1 
        double b = background;

        //std::printf(" mass: %5.1f   signal: %+.4e   background: %+.4e   ", mass, s,b);   
        //std::printf(" arg: %.4e = (%.4e - %.4e)/%.4e\n", arg, t, t_mean, t_sigma);   
        
        return s / (s + b); 
    }; 

}; 

PpromptCalculator get_PpromptCalculator(TH2D* mass_vs_t, double t_sigma, double t_mean, double mass_bin_size, double fit_range)
{
    //now, fit each points. 
    PpromptCalculator calc; 

    auto max = mass_vs_t->GetXaxis(); 
    auto tax = mass_vs_t->GetYaxis(); 

    //center of first bin
    calc.m_min = max->GetBinCenter(1); 
    calc.dm = max->GetBinWidth(1); 

    calc.t_mean = t_mean;
    calc.t_sigma = t_sigma; 

    //get an estimate for the number of bins
    int bins_per_fit = mass_bin_size / max->GetBinWidth(1); 

    if (bins_per_fit < 1) {
        Error(__func__, "number of fit slices < 1"); 
        return calc; 
    }

    calc.dm = max->GetBinWidth(1) * ((double)bins_per_fit); 

    const double dt = tax->GetBinWidth(1); 

    new TCanvas; 
    int first_bin =1; 
    while (first_bin < max->GetNbins()) {

        int last_bin = first_bin + bins_per_fit-1; 
        if (last_bin > max->GetNbins()) break;  

        /*std::printf("Fitting bins: [%i, %i]     x-range: [%.1f, %.1f]\n", first_bin, last_bin,
            max->GetBinCenter( first_bin ), 
            max->GetBinCenter( last_bin )
        );*/  

        auto profile = std::unique_ptr<TH1D>(mass_vs_t->ProjectionY("_py", first_bin, last_bin)); 
        profile->SetBit(kMustCleanup);
        profile->ResetBit(kCanDelete); 
        profile->SetDirectory(nullptr); 
        
        auto tf1 = std::make_unique<TF1>("slice_fit", [dt, t_sigma, t_mean](double *x, double *par){
                double signal = par[0]; 
                double background = par[1]; 
                double t = x[0]; 

                double arg = (t - t_mean)/t_sigma; 

                return dt*( background   +   signal*std::exp( -arg*arg/2. )*0.398942280401/t_sigma ); 
            }, 
            t_mean - fit_range*t_sigma, 
            t_mean + fit_range*t_sigma, 
            2
        );
        tf1->SetBit(kMustCleanup); 
        tf1->ResetBit(kCanDelete); 

        //guess the sigma value 
        double background 
            = profile->GetBinContent( tax->FindBin(t_mean - fit_range*t_sigma) )
            + profile->GetBinContent( tax->FindBin(t_mean + fit_range*t_sigma) );
        background *= 0.5/dt; 

        double signal = ( profile->GetMaximum() - background )*t_sigma / (0.398942280401 * dt); 

        tf1->SetParameter(0, background);
        tf1->SetParameter(1, signal); 

        auto fitptr = profile->Fit("slice_fit", "L R Q S"); 

        if (fitptr->IsValid()==false) {
            Error(__func__, "Something went wrong fitting slice: [%.1f, %.1f]", 
                max->GetBinCenter( first_bin ), 
                max->GetBinCenter( last_bin ) 
            );
            return calc; 
        }

        profile->SetMaximum(mass_vs_t->GetMaximum()*1.2); 
        profile->SetMinimum(0.);
        //profile->DrawCopy(); 
        //tf1->Draw("SAME"); 

        //gPad->Modified(); 
        //gPad->Update(); 
        //gPad->SaveAs("plots/slice_fit.gif+20"); 

        first_bin += bins_per_fit; 

        signal = std::max( 0., fitptr->Parameter(0) ); 
        background = std::max( 0., fitptr->Parameter(1) );

        calc.points.emplace_back( signal, background ); 
        
    }
    //std::printf("n. points: %zi\n", calc.points.size()); 
    return calc; 
}


#endif