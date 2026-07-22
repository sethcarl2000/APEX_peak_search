#ifndef ToyEventGenerator_h
#define ToyEventGenerator_h

#include <bininfo.hpp>
#include <fit_exponential_poly.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <newton_optimizer.hpp> 
#include <fit_parameter.hpp> 
//ROOT headers
#include <TRandom3.h> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TAxis.h> 
#include <functional> 
//stdlib headers
#include <cstdio> 
#include <iostream> 
#include <cmath> 


class ToyEventGenerator 
{
private: 

    double xmin{-7}, xmax{+7}; 

    double bg_sigma = 2.5;
    double bg_skew = 0.5; 
    double bg_center = 0.; 

    double signal_fraction=1e-4; 

    double signal_sigma = 0.05; 
    double signal_x0 = 2.; 

    TRandom3 my_rand; 

public: 

    void Range(double x0, double x1) { xmin=x0; xmax=x1; }
    void BG_sigma(double x)     { bg_sigma=x; }
    void BG_skew(double x)      { bg_skew=x; }
    void BG_center(double x)    { bg_center=x; }

    //fraction of total data sample that is signal 
    void Signal_fraction(double x)  { signal_fraction=x; }

    void Signal_center(double x)    { signal_x0=x; }
    void Signal_sigma(double x)     { signal_sigma=x; }
    
    double operator()(void) {
        if (my_rand.Rndm() >= signal_fraction) { 
            return Gen_bg(); 
        } 
        return Gen_signal(); 
    }

    double Gen_bg() {
        double ret; 
        do { 
            ret = my_rand.Gaus()*bg_sigma; 
            ret += bg_skew * std::sin( ret*3.1415925636/(bg_sigma*3.) ); 
        } while (ret > xmax || ret < xmin); 
        return ret; 
    }

    double Gen_signal() {
        double ret; 
        do { ret=my_rand.Gaus()*signal_sigma + signal_x0; } while (ret > xmax || ret < xmin); 
        return ret; 
    }

};

#endif