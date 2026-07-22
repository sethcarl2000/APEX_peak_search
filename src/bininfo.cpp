#include "bininfo.hpp"

#include <stdexcept>

#include <TString.h>
#include <TAxis.h> 

namespace peak_search
{
//___________________________________________________________________________________________
histo_1D_t copy_1D_hist(TH1D* hist)
{   
    if (!hist) {
        throw std::logic_error(Form("in <%s>: hist passed is null",__func__));
        return {}; 
    }

    //make a copy of all bin contents
    auto xax = hist->GetXaxis(); 
    std::vector<bininfo_1D_t> bins; bins.reserve(xax->GetNbins());

    for (int bx=1; bx<=xax->GetNbins(); bx++) {
        bins.emplace_back( xax->GetBinCenter(bx), hist->GetBinContent(bx) );
    }
    return histo_1D_t{ 
        bins, 
        xax->GetXmin(), xax->GetXmax() 
    }; 
};

//___________________________________________________________________________________________
histo_2D_t copy_2D_hist(TH2D* hist)
{
    if (!hist) {
        throw std::logic_error(Form("in <%s>: hist passed is null",__func__));
        return {}; 
    }

    //make a copy of all bin contents
    auto xax = hist->GetXaxis();
    auto yax = hist->GetYaxis();  
    std::vector<bininfo_2D_t> bins; bins.reserve(xax->GetNbins() * yax->GetNbins());

    for (int bx=1; bx<=xax->GetNbins(); bx++) {
        for (int by=1; by<=yax->GetNbins(); by++) {
            bins.emplace_back( 
                xax->GetBinCenter(bx), 
                yax->GetBinCenter(by),
                hist->GetBinContent(bx,by) 
            );
        }
    }
    return histo_2D_t{
        bins, 
        xax->GetXmin(), xax->GetXmax(), 
        yax->GetXmin(), yax->GetXmax()
    }; 
};

};