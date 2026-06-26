#ifndef peak_search_LegendrePoly_hpp
#define peak_search_LegendrePoly_hpp
//*************************
// -bear's comment. 25 Jun 26

#include <Fcn1D/Fcn1D.hpp>

namespace peak_search
{

class LegendrePoly : public Fcn1D {
private: 
    double xmin, xmax; 
    double x_scale, x_center;

public: 

    LegendrePoly(const std::vector<double>& par={}, double x_min=-1., double x_max=+1.);
    
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
