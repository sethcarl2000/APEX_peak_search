#ifndef peak_search_FitTestThreadManager_hpp
#define peak_search_FitTestThreadManager_hpp

#include <FitTestFunction.hpp>
#include <Histo1D.hpp>
// ROOT 
#include <TH1D.h>  
#include <TH2D.h> 
#include <TRandom3.h> 
#include <TObject.h> 
// stdlib
#include <thread> 
#include <memory> 
#include <vector>
#include <mutex> 

namespace peak_search
{

class FitTestKernel; 

class FitTestThreadManager {
private: 
    size_t fThreadId; 

    FitTestKernel* fParent; 

    FitTestFunction fTestFcn; 

    //histograms to fill
    std::map<const TH1D*, std::unique_ptr<TH1D>> fTH1D; 
    std::map<const TH2D*, std::unique_ptr<TH2D>> fTH2D; 

    enum class htype { kTH1D, kTH2D }; 

    void make_threadlocal_hist_copy(TObject* source, htype type); 

    //current test mass
    double fMass; 

    //random number generator
    std::unique_ptr<TRandom3> fMyRand; 



public: 
    //we've private-ed the constructor, so only the 'FitTestKernel' can make copies of this object. 
    FitTestThreadManager(
        size_t thread_id, 
        const FitTestFunction& fcn, 
        FitTestKernel* parent, 
        const std::vector<TH1D*>& f_TH1D, 
        const std::vector<TH2D*>& f_TH2D
    );  

    //provide a randomly sampled histogram
    Histo1D get_spectrum(double xmin, double xmax); 

    // get thread-local copy of user TH1D
    TH1D* GetUserTH1D(const TH1D*);

    // get thread-local copy of user TH1D
    TH2D* GetUserTH2D(const TH2D*);

    //run a scan of the mass spectrum at the given mass
    void run_test(double mass); 

    //get the current test mass
    double get_mass() const { return fMass; } 
};

};


#endif