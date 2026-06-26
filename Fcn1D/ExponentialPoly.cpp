#include "ExponentialPoly.hpp"

namespace peak_search
{

//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
ExponentialPoly::ExponentialPoly(const std::vector<double>& par, double x_min, double x_max)
    : Fcn1D(par), xmin{x_min}, xmax{x_max}
{
    x_scale  = (x_max - x_min)/2.; 
    x_center = (x_max + x_min)/2.; 
}
//_______________________________________________________________________________________________
double ExponentialPoly::operator()(double x) const 
{
    if (par.empty()) return 1.; 
    
    //scale 'x' 
    x = (x - x_center)/x_scale; 

    //use homer's rule to compute the polynomial efficiently 
    double poly = par.back(); 
    for (int i=((int)par.size())-2; i>=0; i--) { poly = par[i] + poly*x; }
    return std::exp( poly );
}
//_______________________________________________________________________________________________
double ExponentialPoly::Di(double x, int i) const 
{
    x = (x - x_center)/x_scale; 
    return (*this)(x) * numbers::int_pow(x, i); 
}
//_______________________________________________________________________________________________
double ExponentialPoly::Di_Dj(double x, int i, int j) const 
{
    x = (x - x_center)/x_scale; 
    return (*this)(x) * numbers::int_pow(x, i) * numbers::int_pow(x, j); 
}
//_______________________________________________________________________________________________

};