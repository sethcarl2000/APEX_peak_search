
#include <make_histogram_copy.hpp>
//stdlib headers
#include <stdexcept>
//ROOT headers
#include <TError.h>
#include <TAxis.h>

namespace peak_search
{

Histo1D make_histogram_copy(TH1D* hist)
{
    //make a copy of the histogram passed to us 
    Histo1D copy; 

    //check to make sure histogram passed is not null
    if (!hist) {
        throw std::invalid_argument(Form("in <%s>: histogram passed is null!", __func__)); 
        return copy; 
    }
    
    auto xax = hist->GetXaxis();

    const int nbins = xax->GetNbins(); 

    copy.bins.reserve(nbins);

    double x = xax->GetXmin(); 
    double dx = (xax->GetXmax() - xax->GetXmin())/((double)nbins); 

    for (int i=1; i<=nbins; i++) { 
        copy.bins.emplace_back( 
            x,          //xmin
            x += dx,    //xmax
            hist->GetBinContent(i) //bin content 
        );
    }

    return copy; 
}

};