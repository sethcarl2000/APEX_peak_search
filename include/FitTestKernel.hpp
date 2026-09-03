#ifndef peak_search_FitTestKernel_hpp
#define peak_search_FitTestKernel_hpp


#include <Fcn1D/ExponentialPoly.hpp>
#include <SignalFit.hpp>
#include <FitTestThreadManager.hpp>
#include <Histo1D.hpp>
// ROOT
#include <TObject.h> 
// stdlib headers
#include <functional> 
#include <memory>
#include <vector>
#include <map> 

namespace peak_search 
{

class FitTestKernel {
private: 

    static constexpr double fMinMass{140}, fMaxMass{280};
    double fStats{250e6};

    std::unique_ptr<ExponentialPoly> fBackgroundModel{nullptr};     

    //size of a bin 
    double fBinSize{0.5}; 

    double fMinFitMass{150}, fMaxFitMass{270}; 
    size_t fN_steps{100}; 
    
    // append bin contents of 'source' histogram into 'target' histogram
    void copy_histogram(TH1D* source, TH1D* target);

    // append bin contents of 'source' histogram into 'target' histogram
    void copy_histogram(TH2D* source, TH2D* target);

    // make a copy of a given TObjects, and give ourselves control of them. 
    
    // list of user-provided objects 
    std::vector<TH1D*> fUserTH1D{}; 
    std::vector<TH2D*> fUserTH2D{}; 
    
    FitTestKernel(); 
    ~FitTestKernel() = default; 

public: 

    // delete copy ctor & copy assignment operator
    FitTestKernel(const FitTestKernel&) = delete; 
    FitTestKernel& operator=(const FitTestKernel&) = delete; 
    
    static FitTestKernel& Instance(); 


    void RunTest(size_t n_scans, FitTestFunction test_function, size_t n_threads=0);


    void SetMassRange(double min, double max) { fMinFitMass=min; fMaxFitMass=max; }
    void SetNSteps(size_t n_steps) { fN_steps=n_steps; }
    void SetTotalStats(double stats) { fStats =stats; }
    void SetMassBinSize(double binsize) { fBinSize =binsize; }

    Histo1D GetSpectrum(double m_min, double m_max, TRandom3* generator);  

    /// @brief Add a user-created histogram to the list of avialable histograms. 
    /// @param obj User-defined histogram
    /// @return unique ID with which to access histogram 
    size_t AddTH1D(TH1D& obj) { fUserTH1D.emplace_back(&obj); return fUserTH1D.size()-1; } 
    size_t AddTH2D(TH2D& obj) { fUserTH2D.emplace_back(&obj); return fUserTH2D.size()-1; } 

};

};


#endif