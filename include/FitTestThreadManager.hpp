#ifndef peak_search_FitTestThreadManager_hpp
#define peak_search_FitTestThreadManager_hpp

#include <FitTestFunction.hpp>
#include <Histo1D.hpp>
// ROOT 
#include <TH1D.h>  
#include <TH2D.h> 
#include <TRandom3.h> 
// stdlib
#include <thread> 
#include <memory> 
#include <mutex> 

namespace peak_search
{

class FitTestKernel; 

class FitTestThreadManager {
private: 

    friend class FitTestKernel; 

    size_t fThreadId; 

    FitTestKernel* fParent; 

    FitTestFunction fTestFcn; 

    //histograms to fill
    std::unique_ptr<TH2D> fHist_m_vs_mu, fHist_m_vs_Z, fHist_m_vs_pQ0; 
    std::unique_ptr<TH1D> fHist_pQ0; 

    //current test mass
    double fMass; 

    //random number generator
    std::unique_ptr<TRandom3> fMyRand; 

    //run a scan of the mass spectrum at the given mass
    void run_test(double mass); 

public: 
    //we've private-ed the constructor, so only the 'FitTestKernel' can make copies of this object. 
    FitTestThreadManager(size_t thread_id=0, FitTestKernel* parent=nullptr);  

    ~FitTestThreadManager() = default; 

    //provide a randomly sampled histogram
    Histo1D get_spectrum(double xmin, double xmax); 

    //get the current test mass
    double get_mass() const { return fMass; } 
};

};


#endif