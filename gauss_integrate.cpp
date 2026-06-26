
#include <gauss_integrate.hpp>

#include <cstdio> 
#include <vector> 

namespace {
    //order of quadrature
    constexpr int n_order = 4; 

    const std::vector<double> gauss_weight{
        0.3478548451374538,
        0.6521451548625461, 	
        0.6521451548625461, 	
        0.3478548451374538
    };

    const std::vector<double> gauss_abscissa{
        -0.8611363115940526,
        -0.3399810435848563,
        +0.3399810435848563,
        +0.8611363115940526
    };
};

namespace peak_search
{

//___________________________________________________________________________________________________________________
double gauss_integrate(const std::function<double(double)>& fcn, double xmin,double xmax)
{
#ifdef DEBUG
    std::printf("in <gauss_integrate>\n");
#endif
    double scale  = (xmax-xmin)/2.;
    double center = (xmax+xmin)/2.; 

    double sum =0.; 
    for (int ix=0; ix<n_order; ix++) {
        double x = scale*gauss_abscissa[ix] + center; 

#ifdef DEBUG
        std::printf("in <gauss_integrate>: calling fcn\n");
#endif
        sum += fcn(x) * gauss_weight[ix];
    }
    return sum * scale; 
}
//___________________________________________________________________________________________________________________
double gauss_integrate(const Fcn1D& fcn, double xmin,double xmax)
{
#ifdef DEBUG
    std::printf("in <gauss_integrate>\n");
#endif
    double scale  = (xmax-xmin)/2.;
    double center = (xmax+xmin)/2.; 

    double sum =0.; 
    for (int ix=0; ix<n_order; ix++) {
        double x = scale*gauss_abscissa[ix] + center; 

#ifdef DEBUG
        std::printf("in <gauss_integrate>: calling fcn\n");
#endif
        sum += fcn(x) * gauss_weight[ix];
    }
    return sum * scale; 
}
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
double gauss_integrate(double (*fcn)(double,double), double xmin,double xmax, double ymin,double ymax)
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
            
            sum += fcn(x,y) * gauss_weight[ix]*gauss_weight[iy];
        }
    }
    return sum * scale_x * scale_y;
}
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________
//___________________________________________________________________________________________________________________

};