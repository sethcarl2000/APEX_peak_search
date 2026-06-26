#include <Fcn1D/Fcn1D.hpp>

#include <algorithm> 
#include <sstream>
#include <string> 
#include <stdexcept> 


namespace peak_search
{

//_______________________________________________________________________________________________________________________________
void Fcn1D::SetParams(const std::vector<double>& new_par) { 
    std::copy( new_par.begin(), new_par.end(), par.begin() ); 
}
//_______________________________________________________________________________________________________________________________
void Fcn1D::SetParams(const Eigen::VectorXd& v) { 
    std::copy( v.begin(), v.end(), par.begin() ); 
}
//_______________________________________________________________________________________________________________________________
void Fcn1D::SetParams(const std::vector<fit_parameter_t>& new_par) {
    par.resize(new_par.size()); 
    for (size_t i=0; i<new_par.size(); i++) { par[i] = new_par[i].val; }
}
//_______________________________________________________________________________________________________________________________
Fcn1D& Fcn1D::operator+=(const Eigen::VectorXd& dX) 
{
    //check if the number of parameters for the input vector matches the number of parameters this fcn has
    if (dX.size() != par.size()) {
        std::ostringstream oss; 
        oss << "in <Fcn1D::"<< __func__ << ">: size of operands do not match. Fcn1D n. params: " << par.size() << " vs rhs: " << dX.size();
        throw std::logic_error(oss.str()); 
        return *this; 
    }

    for (size_t i=0; i<par.size(); i++) par[i] += dX(i); 
    
    // qwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwu A
    // -muon's comment 25 jun 26
    return *this; 
}
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________

}; 