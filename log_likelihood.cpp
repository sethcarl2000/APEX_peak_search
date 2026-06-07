#include <log_likelihood.hpp>

namespace peak_search
{

//_______________________________________________________________________________________________
double log_likelihood(TH1D* hist, std::function<double(double)> fcn)
{   
    if (!hist) {
        throw std::logic_error("in <log_likelihood>: hist ptr is null"); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }

    auto xax = hist->GetXaxis(); 

    double log_likelihood = 0.; 

    for (int bx=1; bx<xax->GetNbins(); bx++) {

        double expect = fcn( xax->GetBinCenter(bx) );

        //number of events in this bin 
        double ni = hist->GetBinContent(bx);

        log_likelihood += -expect + ni*std::log(expect) - log_factorial(ni);  
    }   
    return log_likelihood; 
}

//_______________________________________________________________________________________________
double log_likelihood(TH2D* hist, std::function<double(double,double)> fcn)
{   
    if (!hist) {
        throw std::logic_error("in <log_likelihood>: hist ptr is null"); 
        return std::numeric_limits<double>::quiet_NaN(); 
    }

    auto xax = hist->GetXaxis(); 
    auto yax = hist->GetYaxis(); 

    double log_likelihood = 0.; 

    for (int bx=1; bx<=xax->GetNbins(); bx++) {
        double x = xax->GetBinCenter(bx);

        for (int by=1; by<=yax->GetNbins(); by++) {
            double y = yax->GetBinCenter(by);        
        
            double expect = fcn( x, y );

            //number of events in this bin 
            double ni = hist->GetBinContent(bx,by);

            log_likelihood += -expect + ni*std::log(expect) - log_factorial(ni);  
        }
    }   
    return log_likelihood; 
}

//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
};