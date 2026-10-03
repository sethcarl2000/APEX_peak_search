#include <Fcn1D/Fcn1D.hpp>

#include <algorithm> 
#include <sstream>
#include <string> 
#include <stdexcept> 


namespace peak_search
{

//_______________________________________________________________________________________________________________________________
void Fcn1D::SetParams(std::vector<double>::const_iterator _begin, size_t _length) { 
    par.assign( _begin, _begin + _length ); 
}
//_______________________________________________________________________________________________________________________________
void Fcn1D::SetParams(const Eigen::VectorXd& new_par) { 
    par.assign( new_par.begin(), new_par.end() ); 
}
//_______________________________________________________________________________________________________________________________
void Fcn1D::SetParams(const std::vector<double>& new_par) { 
    par = new_par; 
}
//_______________________________________________________________________________________________________________________________
void Fcn1D::SetParams(const std::vector<fit_parameter_t>& new_par) {
    par.resize(new_par.size()); 
    for (size_t i=0; i<new_par.size(); i++) { par[i] = new_par[i].val; }
}
//_______________________________________________________________________________________________________________________________
Fcn1D& Fcn1D::operator+=(const Eigen::VectorXd& dX) 
{
    std::vector<double> new_params; new_params.reserve(dX.size()); 

    //check if the number of parameters for the input vector matches the number of parameters this fcn has
    if (dX.size() != par.size()) {
        std::ostringstream oss; 
        oss << "in <Fcn1D::"<< __func__ << ">: size of operands do not match. Fcn1D n. params: " << par.size() << " vs rhs: " << dX.size();
        throw std::logic_error(oss.str()); 
        return *this; 
    }

    for (size_t i=0; i<par.size(); i++) new_params.push_back(par[i] + dX(i)); 

    SetParams(new_params); 
    
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