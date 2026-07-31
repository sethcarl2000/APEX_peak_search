#include "Fcn1D/ExponentialLegendre.hpp"
#include <legendre_polynomial.hpp> 

namespace peak_search
{

//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
ExponentialLegendre::ExponentialLegendre(const std::vector<double>& par, double x_min, double x_max)
    : Fcn1D(par), xmin{x_min}, xmax{x_max}
{
    x_scale  = (x_max - x_min)/2.; 
    x_center = (x_max + x_min)/2.; 
}
//_______________________________________________________________________________________________
double ExponentialLegendre::operator()(double x) const 
{
    //if the only order is '1', or '0', then it's constant
    if (par.empty()) return 0.; 
    
    //first, scale x, so that x_min = -1, and x_max = +1; 
    x = (x - x_center)/x_scale; 

    //now, compute the legendre polynomial, using bonett's recursion formula 
    double Pn_m1  = 0.; //P_{n-1}(x) 
    double Pn     = 1.; //P_{n}(x)
    
    double result = 0.; 

    int n=0;
    int n_max = par.size()-1; 
    
    do {

        //add the current poly to the list of all orders
        result += par[n] * Pn; 

        //compute the next poly, using bonett's recursion formula
        // P_{n+1}(x)
        double N = (double)n; 
        double Pn_p1 = ( (2.*N + 1.)*x*Pn - N*Pn_m1 )/(N + 1.); 

        //now, shift all the definitions up by n => n+1
        ++n; 
        Pn_m1 = Pn; 
        Pn = Pn_p1; 

    } while (n < n_max);

    return std::exp( result ); 
}
//_______________________________________________________________________________________________
double ExponentialLegendre::Di(double x, int i) const 
{
    return (*this)(x) * legendre_polynomial((x - x_center)/x_scale, i); 
}
//_______________________________________________________________________________________________
double ExponentialLegendre::Di_Dj(double x, int i, int j) const 
{
    double ret = (*this)(x);
    
    x = (x - x_center)/x_scale; 
    return ret * legendre_polynomial(x, i) * legendre_polynomial(x, j); 
}
//_______________________________________________________________________________________________

};