#include "legendre_polynomial.hpp"

#include <cmath> 

namespace peak_search
{

//________________________________________________________________________________________________
double legendre_polynomial(double x, int n)
{
    //this is implemented via bonnet's recursion formula: 
    // 
    // (n+1)*P_{n+1} = (2n+1)x*P_{n} - n*P_{n-1} 
    //
    if (n <= 0) { return 1.; } else { if (n <= 1) return x; }

    double N_end = (double)std::round(n) - 0.1;
    double N = 1.; 
    double Pn{x}, Pn_m1{1.}; 
    do {
        double Pn_p1 = ( (2*N + 1)*x*Pn - N*Pn_m1 )/(N + 1); 

        N += 1.; 
        Pn_m1 = Pn; 
        Pn = Pn_p1; 
    
    } while (N < N_end);
    
    return Pn;
}
//________________________________________________________________________________________________
double legendre_polynomial_dx(double x, int n)
{
    //this is implemented via bonnet's recursion formula (the differential version) 
    // 
    // (n+1)*P_{n+1} = (2n+1)x*P_{n} - n*P_{n-1} 
    //
    double N = (double)n; 

    double Pn{legendre_polynomial(x,n)}, Pn_m1{legendre_polynomial(x,n-1)}; 

    double ret = x*Pn - Pn_m1; 

    return ret * N / ( x*x - 1. ); 
}
//________________________________________________________________________________________________
double legendre_polynomial_ddx(double x, int n)
{
    //using a finite-difference method
    static constexpr double h = 0.01; 

    double Pn_ph = legendre_polynomial_dx(x + 0.5*h, n); 
    double Pn_mh = legendre_polynomial_dx(x - 0.5*h, n); 

    return ( Pn_ph - Pn_mh )/h;    
}
//________________________________________________________________________________________________
//________________________________________________________________________________________________
//________________________________________________________________________________________________
//________________________________________________________________________________________________

}; 