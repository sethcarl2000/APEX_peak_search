#ifndef pdf_fcn_hpp
#define pdf_fcn_hpp

#include <functional> 

namespace peak_search
{

class PdfFcn_1D {
private: 
    int n_params; 
    std::function<double(double,const double*)> fcn;

public: 

    PdfFcn_1D(int _n_params, std::function<double(double,const double*)> _fcn)
        : n_params{_n_params}, fcn{_fcn} {}; 

    inline double operator()(double x, const double* params) const { return fcn(x,params); } 

    inline int GetNParams() const { return n_params; }
}; 

class PdfFcn_2D {
private: 
    int n_params; 
    std::function<double(double,double,const double*)> fcn;

public: 

    PdfFcn_2D(int _n_params, std::function<double(double,double,const double*)> _fcn)
        : n_params{_n_params}, fcn{_fcn} {}; 

    inline double operator()(double x, double y, const double* params) const { return fcn(x,y,params); } 

    inline int GetNParams() const { return n_params; }
}; 

};

#endif