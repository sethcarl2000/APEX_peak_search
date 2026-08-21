
#include <make_histogram_copy.hpp>
//ROOT headers
#include <TError.h>
#include <TAxis.h>
//stdlib headers
#include <stdexcept>
#include <cmath> 

namespace peak_search
{

Histo1D make_histogram_copy(TH1D* hist, double xmin, double xmax)
{
    //make a copy of the histogram passed to us 
    Histo1D copy; 

    //check to make sure histogram passed is not null
    if (!hist) {
        throw std::invalid_argument(Form("in <%s>: histogram passed is null!", __func__)); 
        return copy; 
    }

    
    auto xax = hist->GetXaxis();

    int min_bin, max_bin; 

    //find minimum bin
    if (numbers::is_nan(xmin)) {
        min_bin = 1;  
    } else {
        xmin = std::max(xmin, xax->GetXmin());
        min_bin = xax->FindBin(xmax);
    }

    //find maximum bin
    if (numbers::is_nan(xmax)) {
        max_bin = xax->GetNbins();  
    } else {
        xmax = std::min(xmax, xax->GetXmax());
        max_bin = xax->FindBin(xmax);
    }

    copy.bins.reserve(max_bin - min_bin + 1);

    double x = xax->GetXmin(); 
    double dx = (xax->GetXmax() - xax->GetXmin())/((double)xax->GetNbins()); 

    for (int i=min_bin; i<=max_bin; i++) { 
        copy.bins.emplace_back( 
            x,          //xmin
            x += dx,    //xmax
            hist->GetBinContent(i) //bin content 
        );
    }

    return copy; 
}

};