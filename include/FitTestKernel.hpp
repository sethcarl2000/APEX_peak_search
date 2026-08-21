#ifndef peak_search_FitTestKernel_hpp
#define peak_search_FitTestKernel_hpp


#include <Fcn1D/ExponentialPoly.hpp>
#include <SignalFit.hpp>
#include <FitTestThreadManager.hpp>
#include <Histo1D.hpp>
// ROOT
#include <ROOT/RDF/HistoModels.hxx>
// stdlib headers
#include <functional> 
#include <memory>
#include <vector>

namespace peak_search 
{

class FitTestKernel {
private: 

    static constexpr double fMinMass{140}, fMaxMass{280};
    static constexpr double fStats{250e6};

    std::unique_ptr<ExponentialPoly> fBackgroundModel{nullptr}; 

    TH2D *fHist_m_vs_mu, *fHist_m_vs_Z, *fHist_m_vs_pQ0;
    TH1D* fHist_pQ0; 

    ROOT::RDF::TH2DModel fModel_m_vs_mu, fModel_m_vs_Z, fModel_m_vs_pQ0;
    ROOT::RDF::TH1DModel fModel_pQ0; 

    //size of a bin 
    double fBinSize{0.5}; 

    double fMinFitMass, fMaxFitMass; 
    size_t fN_steps; 

    TH1D* construct_TH1D(const ROOT::RDF::TH1DModel&);
    TH2D* construct_TH2D(const ROOT::RDF::TH2DModel&);  

    // append bin contents of 'source' histogram into 'target' histogram
    void copy_histogram(TH1D* source, TH1D* target);

    // append bin contents of 'source' histogram into 'target' histogram
    void copy_histogram(TH2D* source, TH2D* target);

public: 

    FitTestKernel(double min_fit_mass, double max_fit_mass, size_t n_steps); 

    void RunTest(size_t n_scans, FitTestFunction test_function, size_t n_threads=0); 

    void DrawResults();

    Histo1D GetSpectrum(double m_min, double m_max, TRandom3* generator);  
};

};


#endif