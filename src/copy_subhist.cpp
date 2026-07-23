#include <copy_subhist.hpp>

#include <TAxis.h>
#include <TString.h> 

#include <stdexcept> 
#include <cmath> 

namespace peak_search
{
   
histo_1D_t copy_subhist(TH1D* hist, double xmin, double xmax)
{
    if (!hist) {
        throw std::logic_error(Form("in <%s>: hist passed is null",__func__));
        return {}; 
    }

    if (xmin >= xmax) {
        throw std::logic_error(Form("in <%s>: end of sub-range requested (%f) is less than start of sub-range (%f).",__func__, xmax, xmin));
        return {}; 
    }

    auto xax = hist->GetXaxis(); 

    xmin = std::max( xmin, xax->GetXmin() );
    xmax = std::min( xmax, xax->GetXmax() );
    

    histo_1D_t ret; 

    //get starting & ending bins
    int bin0 = xax->FindBin(xmin);
    int bin1 = xax->FindBin(xmax);
    
    //set the starting point of the output hist to the leading edge of this bin
    double dx = xax->GetBinWidth(bin0);

    ret.xmin = xax->GetBinCenter(bin0) - dx/2.; 
    ret.xmax = xax->GetBinCenter(bin1) + dx/2.; 

    //now, fill in the bins 
    ret.bins.reserve(bin1 - bin0 + 1); 

    for (int bin=bin0; bin<=bin1; bin++) {

        ret.bins.emplace_back( xax->GetBinCenter(bin), hist->GetBinContent(bin) );
    }

    return ret; 
}   

}; 
