
#include <gauss_integrate.hpp>

namespace {
    //order of quadrature
    constexpr int n_order = 6; 

    constexpr double gauss_weight[] = {
        0.3607615730481386,
        0.3607615730481386, 	
        0.4679139345726910, 	
        0.4679139345726910, 	
        0.1713244923791704, 	
        0.1713244923791704 	
    };

    constexpr double gauss_abscissa[] = {
        +0.6612093864662645,
        -0.6612093864662645,
        -0.2386191860831969,
        +0.2386191860831969,
        -0.9324695142031521,
        +0.9324695142031521
    };
};

namespace peak_search
{

//___________________________________________________________________________________________________________________
double gauss_integrate(const PdfFcn_1D& fcn, const double* params, double xmin,double xmax)
{
    double scale  = (xmax-xmin)/2.;
    double center = (xmax+xmin)/2.; 

    double sum =0.; 
    for (int ix=0; ix<n_order; ix++) {
        double x = scale*gauss_abscissa[ix] + center; 
        sum += fcn(x, params) * gauss_weight[ix];
    }
    return sum; 
}
//___________________________________________________________________________________________________________________
double gauss_integrate(const PdfFcn_2D& fcn, const double* params, double xmin,double xmax, double ymin,double ymax)
{
    double scale_x  = (xmax-xmin)/2.;
    double center_x = (xmax+xmin)/2.; 

    double scale_y  = (ymax-ymin)/2.;
    double center_y = (ymax+ymin)/2.; 

    double sum =0.; 
    for (int ix=0; ix<n_order; ix++) {
        double x = scale_x*gauss_abscissa[ix] + center_x;

        for (int iy=0; iy<n_order; iy++) {
            double y = scale_y*gauss_abscissa[iy] + center_y;
            
            sum += fcn(x,y, params) * gauss_weight[ix]*gauss_weight[iy];
        }
    }
    return sum;
}
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________

};