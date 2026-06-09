#include "ExponentialPoly.hpp"

namespace peak_search
{

//_______________________________________________________________________________________________
double ExponentialPoly::operator()(double x) const 
{
    if (coeffs.empty()) return numbers::nan; 
    
    double poly = coeffs.back(); 
    for (int i=((int)coeffs.size())-2; i>=0; i--) { poly = coeffs[i] + poly*x; }
    return std::exp( poly );
}
//_______________________________________________________________________________________________
double ExponentialPoly::Eval(double x, const double *coeffs, int n) 
{
    if (!coeffs) return numbers::nan;

    double poly = coeffs[n-1];
    for (int i=n-2; i>=0; i--) { poly = coeffs[i] + poly*x; }
    return std::exp( poly );
}
//_______________________________________________________________________________________________
//_______________________________________________________________________________________________
ExponentialPoly::operator PdfFcn_1D() const 
{
    auto coeffs_cpy = coeffs; 
    int n = coeffs.size();
    return PdfFcn_1D(
        n,
        std::function<double(double,const double*)>{[coeffs_cpy](double x, const double *coeffs) {
            return Eval(x, const_cast<double*>(coeffs_cpy.data()), coeffs_cpy.size());
        }}
    );

} 
//_______________________________________________________________________________________________

};