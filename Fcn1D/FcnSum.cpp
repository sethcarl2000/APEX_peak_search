#include <Fcn1D/FcnSum.hpp> 


//ROOT
#include <TString.h> 
//stdlib
#include <stdexcept> 

namespace peak_search
{

//_______________________________________________________________________________________________________________________________
FcnSum& FcnSum::SetFcnA(Fcn1D* ptr) { fcnA=std::unique_ptr<Fcn1D>(ptr); UpdateFcn(); return *this; }; 
//_______________________________________________________________________________________________________________________________
FcnSum& FcnSum::SetFcnB(Fcn1D* ptr) { fcnA=std::unique_ptr<Fcn1D>(ptr); UpdateFcn(); return *this; }; 
//_______________________________________________________________________________________________________________________________
FcnSum::FcnSum(Fcn1D* _fcnA, Fcn1D* _fcnB)
    : Fcn1D(), n_pars_A{0}, n_pars_B{0}, fcnA{std::unique_ptr<Fcn1D>(_fcnA)}, fcnB{std::unique_ptr<Fcn1D>(_fcnB)}
{
    UpdateFcn(); 
}
//_______________________________________________________________________________________________________________________________
void FcnSum::UpdateFcn()
{
    if (fcnA == nullptr || fcnB == nullptr) { 
        throw std::logic_error(Form("in <FcnSum::%s>: one or both functions is nullptr. fcnA=%p, fcnB=%p", 
            __func__, fcnA.get(), fcnB.get()
        ));
        return; 
    }

    n_pars_A = fcnA->GetDoF(); 
    n_pars_B = fcnB->GetDoF();

    par.resize(n_pars_A + n_pars_B); 

    //copy the parameters over to the main list
    std::copy( fcnA->GetParams().begin(), fcnA->GetParams().end(), par.begin() ); 
    std::copy( fcnB->GetParams().begin(), fcnB->GetParams().end(), par.begin()+n_pars_A );
}
//_______________________________________________________________________________________________________________________________
void FcnSum::SetParams(const std::vector<double>& new_pars)
{
    if (new_pars.size() != n_pars_A + n_pars_B) {
        throw std::logic_error(Form("in <FcnSum::%s>: size of new parameters (%zi) does not match size of both functions together (%i+%i)", 
            __func__, new_pars.size(), fcnA->GetDoF(), fcnB->GetDoF()
        ));
        return; 
    }
    par = new_pars; 
    //copy list parameters down to each function
    fcnA->SetParams(new_pars.begin(), n_pars_A); 
    fcnB->SetParams(new_pars.begin()+n_pars_A, n_pars_B); 
}
//_______________________________________________________________________________________________________________________________
void FcnSum::SetParams(const Eigen::VectorXd& new_pars)
{
    if (new_pars.size() != n_pars_A + n_pars_B) {
        throw std::logic_error(Form("in <FcnSum::%s>: size of new parameters (%zi) does not match size of both functions together (%i+%i)", 
            __func__, new_pars.size(), fcnA->GetDoF(), fcnB->GetDoF()
        ));
        return; 
    }
    par.resize(new_pars.size());
    std::copy( new_pars.begin(), new_pars.end(), par.begin() ); 
    
    //copy list parameters down to each function
    fcnA->SetParams(par.begin(), n_pars_A); 
    fcnB->SetParams(par.begin()+n_pars_A, n_pars_B); 
}
//_______________________________________________________________________________________________________________________________
void FcnSum::SetParams(const std::vector<fit_parameter_t>& new_pars)
{
    if (new_pars.size() != n_pars_A + n_pars_B) {
        throw std::logic_error(Form("in <FcnSum::%s>: size of new parameters (%zi) does not match size of both functions together (%i+%i)", 
            __func__, new_pars.size(), fcnA->GetDoF(), fcnB->GetDoF()
        ));
        return; 
    }

    //copy list parameters down to each function
    size_t ii=0; 
    auto& parsA = fcnA->GetParams(); 
    for (size_t i=0; i<n_pars_A; i++) {
        parsA[i] = new_pars[ii++].val; 
    }

    auto& parsB = fcnA->GetParams(); 
    for (size_t i=0; i<n_pars_B; i++) {
        parsB[i] = new_pars[ii++].val; 
    }
}
//_______________________________________________________________________________________________________________________________
double FcnSum::operator()(double x) const
{
    return (*fcnA.get())(x) + (*fcnB.get())(x); 
}
//_______________________________________________________________________________________________________________________________
double FcnSum::Di(double x, int i) const 
{
    if (i < (int)n_pars_A) { return fcnA->Di(x, i); }
    return fcnB->Di(x, i-n_pars_A); 
}
//_______________________________________________________________________________________________________________________________
double FcnSum::Di_Dj(double x, int i, int j) const 
{
   if (i < (int)n_pars_A) {
        if (j <  (int)n_pars_A) { return fcnA->Di_Dj(x, i,j); } else { return 0.; }
   } else {
        if (j >= (int)n_pars_A) { return fcnB->Di_Dj(x, i,j); } else { return 0.; }
   }
}
//_______________________________________________________________________________________________________________________________
//_______________________________________________________________________________________________________________________________

};