#ifndef peak_search_FcnSum_hpp
#define peak_search_FcnSum_hpp

#include <Fcn1D/Fcn1D.hpp>
#include <memory>

namespace peak_search
{

/// @brief Sum of two 1D fcn's
class FcnSum : public Fcn1D {
private: 

    Fcn1D *fcnA, *fcnB; 

    size_t n_pars_A{0}, n_pars_B{0}; 

    void UpdateFcn();

public: 

    FcnSum(Fcn1D* _fcnA=nullptr, Fcn1D* _fcnB=nullptr);

    FcnSum& SetFcnA(Fcn1D* ptr); 
    FcnSum& SetFcnB(Fcn1D* ptr); 

    //set parameters. 
    // this may be overridden, for example if you want to require that a certain number of parameters are passed.
    void SetParams(const std::vector<double>& par) override; 
    //set parameters. 
    // this may be overridden, for example if you want to require that a certain number of parameters are passed.
    void SetParams(const Eigen::VectorXd& par) override; 
    //calls above fcn, but uses fit_parameter_t
    void SetParams(const std::vector<fit_parameter_t>& par) override; 

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