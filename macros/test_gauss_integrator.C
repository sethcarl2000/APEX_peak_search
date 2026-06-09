#include <gauss_integrate.hpp>

#include <functional> 
#include <cmath>
#include <cstdio> 

double normal_dist_1d(double x) {
    return std::exp( -0.5*x*x ) / 2.50662827463; 
}

double normal_dist_2d(double x, double y) {
    return normal_dist_1d(x) * normal_dist_1d(y); 
}

int test_gauss_integrator()
{   
    //integrate a normal distribution
    double xmin{-5.}, xmax{+5.}; 

    std::printf("one dimension: \n");
    for (int n_bins=1; n_bins<=10; n_bins++) {

        double sum=0.; 
        double dx = (xmax-xmin)/((double)n_bins); 
        double x=xmin; 
        for (int i=0; i<n_bins; i++) { 
            sum += peak_search::gauss_integrate(normal_dist_1d, x,x+dx); 
            x += dx; 
        }

        std::printf("n. bins: %3i, domain: [%.1f, %.1f], integral = %.8f\n", 
            n_bins, xmin, xmax, 
            sum
        );
    } 

    std::printf("two dimensions: \n");
    for (int n_bins=1; n_bins<=16; n_bins++) {

        double sum=0.; 
        double dx = (xmax-xmin)/((double)n_bins); 
        double x=xmin; 
        for (int i=0; i<n_bins; i++) { 
            double y=xmin; 
            for (int j=0; j<n_bins; j++) {
                sum += peak_search::gauss_integrate(normal_dist_2d, x,x+dx, y,y+dx); 
                y += dx; 
            }
            x += dx; 
        }

        std::printf("n. bins: %3i, domain: x=[%.1f, %.1f], y=[%.1f, %.1f] integral = %.8f\n", 
            n_bins, xmin, xmax, xmin, xmax, 
            sum
        );
    } 


    return 1; 
}
