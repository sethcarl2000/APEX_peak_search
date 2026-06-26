#ifndef peak_search_Fcn1D_hpp
#define peak_search_Fcn1D_hpp

#include <fit_parameter.hpp>

#include <eigen3/Eigen/Core>

#include <cmath> 
#include <functional> 
#include <vector>


namespace peak_search 
{

/// @brief Abstract base-clas for the 
class Fcn1D {
protected: 

    //parameters of function
    std::vector<double> par; 

public: 
    
    Fcn1D(const std::vector<double>& _par={}) : par{_par} {}; 

    virtual ~Fcn1D() = 0; 

    //set parameters. 
    // this may be overridden, for example if you want to require that a certain number of parameters are passed.
    virtual void SetParams(const std::vector<double>& par); 
    //set parameters. 
    // this may be overridden, for example if you want to require that a certain number of parameters are passed.
    virtual void SetParams(const Eigen::VectorXd& v); 
    //calls above fcn, but uses
    virtual void SetParams(const std::vector<fit_parameter_t>& par); 
    
    //get copy of parameters
    std::vector<double> GetParamsCpy() const { return par; }

    //get reference to parameters
    std::vector<double>& GetParams() { return par; }

    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
    
    inline int GetDoF() const { return par.size(); }

    // 'nudge' parameters by given value 
    Fcn1D& operator+=(const Eigen::VectorXd& dX); 
}; 

}; 

#endif 
