#include <Fcn1D/FcnSum.hpp> 
#include <numbers.hpp>

//ROOT
#include <TString.h> 
//stdlib
#include <stdexcept> 
#include <iostream> 

namespace peak_search
{

//_______________________________________________________________________________________________________________________________
FcnSum& FcnSum::SetFcnA(Fcn1D* ptr) { fcnA=ptr; UpdateFcn(); return *this; }; 
//_______________________________________________________________________________________________________________________________
FcnSum& FcnSum::SetFcnB(Fcn1D* ptr) { fcnA=ptr; UpdateFcn(); return *this; }; 
//_______________________________________________________________________________________________________________________________
FcnSum::FcnSum(Fcn1D* _fcnA, Fcn1D* _fcnB)
    : Fcn1D(), n_pars_A{0}, n_pars_B{0}, fcnA{_fcnA}, fcnB{_fcnB}
{   
    //copy information from our list of parameters to each function 
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return; 
    }

    //make an empty list of parameters 
    auto pA = fcnA->GetParams(); 
    auto pB = fcnB->GetParams();

    par.reserve(pA.size() + pB.size()); 

    size_t i;
    i=0; for (auto val : pA) par.emplace_back(pA[i++]);
    i=0; for (auto val : pB) par.emplace_back(pB[i++]); 

    UpdateFcn(); 
}
//_______________________________________________________________________________________________________________________________
void FcnSum::UpdateFcn()
{
    //copy information from our list of parameters to each function 
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return; 
    }

    n_pars_A = fcnA->GetDoF(); 
    n_pars_B = fcnB->GetDoF();

    if (par.size() != n_pars_A + n_pars_B) {
        throw std::logic_error(Form("in <FcnSum::%s>: size of parameter list in 'sum' fcn (%zi) does not macth size of list in 2 sub-functions (%zi+%zi)", 
            __func__, par.size(), n_pars_A, n_pars_B
        ));
        return; 
    }

    //copy the parameters over to the main list
    size_t i=0; 
    for (auto& val : fcnA->GetParams()) val = par[i++]; 
    for (auto& val : fcnB->GetParams()) val = par[i++];     
}
//_______________________________________________________________________________________________________________________________
void FcnSum::SetParams(const std::vector<double>& new_pars)
{
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return;  
    }

    if (new_pars.size() != n_pars_A + n_pars_B) {
        throw std::logic_error(Form("in <FcnSum::%s>: size of new parameters (%zi) does not match size of both functions together (%i+%i)", 
            __func__, new_pars.size(), fcnA->GetDoF(), fcnB->GetDoF()
        ));
        return; 
    }

    par = new_pars; 

    UpdateFcn(); 
}
//_______________________________________________________________________________________________________________________________
void FcnSum::SetParams(const Eigen::VectorXd& new_pars)
{
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return; 
    }

    if (new_pars.size() != n_pars_A + n_pars_B) {
        throw std::logic_error(Form("in <FcnSum::%s>: size of new parameters (%zi) does not match size of both functions together (%i+%i)", 
            __func__, new_pars.size(), fcnA->GetDoF(), fcnB->GetDoF()
        ));
        return; 
    }
    par.resize(new_pars.size());
    std::copy( new_pars.begin(), new_pars.end(), par.begin() ); 
    
    UpdateFcn(); 
}
//_______________________________________________________________________________________________________________________________
void FcnSum::SetParams(const std::vector<fit_parameter_t>& new_pars)
{
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return; 
    }

    if (new_pars.size() != n_pars_A + n_pars_B) {
        throw std::logic_error(Form("in <FcnSum::%s>: size of new parameters (%zi) does not match size of both functions together (%i+%i)", 
            __func__, new_pars.size(), fcnA->GetDoF(), fcnB->GetDoF()
        ));
        return; 
    }

    //copy list parameters down to each function
    par.resize(new_pars.size()); 
    for (size_t i=0; i<new_pars.size(); i++) { par[i] = new_pars[i].val; }

    UpdateFcn(); 
}
//_______________________________________________________________________________________________________________________________
double FcnSum::operator()(double x) const
{
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return numbers::nan; 
    }

    return (*fcnA)(x) + (*fcnB)(x); 
}
//_______________________________________________________________________________________________________________________________
double FcnSum::Di(double x, int i) const 
{
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return numbers::nan; 
    }

    if (i < (int)n_pars_A) { return fcnA->Di(x, i); }
    return fcnB->Di(x, i-n_pars_A); 
}
//_______________________________________________________________________________________________________________________________
double FcnSum::Di_Dj(double x, int i, int j) const 
{
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA, fcnB
        ));
        return numbers::nan; 
    }
    
   if (i < (int)n_pars_A) {
        if (j <  (int)n_pars_A) { return fcnA->Di_Dj(x, i,j); } else { return 0.; }
   } else {
        if (j >= (int)n_pars_A) { return fcnB->Di_Dj(x, i,j); } else { return 0.; }
   }
}
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________

};