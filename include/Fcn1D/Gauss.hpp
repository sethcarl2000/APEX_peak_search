#ifndef peak_search_Gauss_hpp
#define peak_search_Gauss_hpp

#include <Fcn1D/Fcn1D.hpp>

namespace peak_search
{

/// @brief 1D fcn representing normalized gaussisan, with integral = mu (par[0]). NOTE: only 'mu' (0-th parameter) has derivatives implemented
class Gauss : public Fcn1D {
private: 
    double x0, sigma; 

public: 

    Gauss(double mu, double x0, double sigma);

    Gauss(const Gauss&) = default; 
    
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //mandatory interface: 

    /// @brief evaluate fcn
    /// @param x arg 'x'
    /// @return f(x)
    double operator() (double) const override;

    /// @brief derivative w/r/t parameter i
    /// @param x arg 'x'
    /// @param i index of parameter to take derivative w/r/t 
    /// @return d/d\theta-*_i * f(x)
    double Di (double, int) const override; 

    //double-derivative w/r/t parameters i, j 
    /// @param x arg 'x'
    /// @param i index of 1st-parameter to take derivative w/r/t 
    /// @param j index of 2nd-parameter to take derivative w/r/t 
    /// @return d/d\theta_i * d/d\theta_j * f(x)
    double Di_Dj (double, int,int) const override; 

}; 

}; 

#endif
