#ifndef peak_search_Fcn1D_hpp
#define peak_search_Fcn1D_hpp

#include <cmath> 
#include <functional> 
#include <vector>
#include <eigen3/Eigen/Core>

namespace peak_search 
{

/// @brief Abstract base-clas for the 
class Fcn1D {
protected: 

    //number of function parameters
    int DoF; 

    //parameters of function
    std::vector<double> par; 

public: 
    
    Fcn1D(int _dof, const std::vector<double>& _par={}) : DoF{_dof}, par{_par} {}; 

    virtual ~Fcn1D() = 0; 

    void SetParams(const double* par); 


    // ~~~~~~~~~ 
    //mandatory interface: 

    /// @brief evaluate fcn
    /// @param x arg 'x'
    /// @return f(x)
    virtual double operator() (double) const = 0;

    /// @brief derivative w/r/t parameter i
    /// @param x arg 'x'
    /// @param i index of parameter to take derivative w/r/t 
    /// @return d/d\theta_i * f(x)
    virtual double Di (double, int) const = 0; 

    //double-derivative w/r/t parameters i, j 
    /// @param x arg 'x'
    /// @param i index of 1st-parameter to take derivative w/r/t 
    /// @param j index of 2nd-parameter to take derivative w/r/t 
    /// @return d/d\theta_i * d/d\theta_j * f(x)
    virtual double Di_Dj (double, int,int) const = 0; 
    
    inline int GetDoF() const { return DoF; }

}; 

}; 

#endif 
