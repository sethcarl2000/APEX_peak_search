#ifndef bootstrap_sample_hist_C
#define bootstrap_sample_hist_C

#include <TH2D.h>
#include <TAxis.h>
#include <TRandom3.h>
#include <TString.h> 

/// @brief Given 'hist', make another histogram whose bins are randomly chosen from a poisson dist, with the expectation value taken from the 'source' histogram 
/// @param hist source histogram 
/// @param frac fraction of source histogram's total stats to aim fo
/// @return new histogram randomly sampled to mimic given one 
TH2D* bootstrap_sample_hist(TH2D* hist, double frac, TRandom3& rand)
{
    auto xax = hist->GetXaxis();
    auto yax = hist->GetYaxis(); 

    auto hist_cpy = new TH2D(Form("%s_cpy",hist->GetName()), hist->GetTitle(), 
        xax->GetNbins(), xax->GetXmin(), xax->GetXmax(), 
        yax->GetNbins(), yax->GetXmin(), yax->GetXmax()
    );

    for (int ix=1; ix<=xax->GetNbins(); ix++) {
        double x = xax->GetBinCenter(ix); 
        for (int iy=1; iy<=yax->GetNbins(); iy++) {
            double y = yax->GetBinCenter(iy); 

            double n = hist->GetBinContent(ix,iy) * frac; 

            hist_cpy->Fill( x, y, rand.Poisson(n) ); 
        }
    }
    return hist_cpy; 
}

#endif