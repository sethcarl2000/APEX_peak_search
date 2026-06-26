#ifndef peak_search_Gen1DFcn_hpp
#define peak_search_Gen1DFcn_hpp

#include <Fcn1D/Fcn1D.hpp>

#include <numbers.hpp>

#include <vector> 
#include <cmath> 
#include <limits> 
#include <functional> 

namespace peak_search 
{

/// @brief Represents and arbitrary 1D FCN 
class GenFcn1D : public Fcn1D {
    
    using ParFcn        = std::function<double(double,const double*)>; 
    using ParFcnDi      = std::function<double(double,int,const double*)>; 
    using ParFcnDiDj    = std::function<double(double,int,int,const double*)>; 

    ParFcn fFcn; 
    ParFcnDi fFcn_di; 
    ParFcnDiDj fFcn_di_dj; 

public: 
    
    inline GenFcn1D(const std::vector<double>& _par={}, 
        ParFcn _fcn,
        ParFcnDi _fcn_di,
        ParFcnDiDj _fcn_di_dj
    ) : Fcn1D(par), fFcn{_fcn}, fFcn_di{_fcn_di}, fFcn_di_dj{_fcn_di_dj} {};

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //mandatory interface: 

    /// @brief evaluate fcn
    /// @param x arg 'x'
    /// @return f(x)
    inline double operator() (double x) const override {
        return fFcn(x, par.data()); 
    };

    /// @brief derivative w/r/t parameter i
    /// @param x arg 'x'
    /// @param i index of parameter to take derivative w/r/t 
    /// @return d/d\theta-*_i * f(x)
    inline double Di (double x, int i) const override {
        return fFcn_di(x,i, par.data()); 
    }; 

    //double-derivative w/r/t parameters i, j 
    /// @param x arg 'x'
    /// @param i index of 1st-parameter to take derivative w/r/t 
    /// @param j index of 2nd-parameter to take derivative w/r/t 
    /// @return d/d\theta_i * d/d\theta_j * f(x)
    inline double Di_Dj (double x, int i, int j) const override {
        return fFcn_di_dj(x,i,j, par.data()); 
    }; 
};

}; 

#endif