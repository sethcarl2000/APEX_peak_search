
#include <FitTestThreadManager.hpp>
#include <FitTestKernel.hpp>
#include <compute_Q0.hpp>
#include <generate_toy_events.hpp>
// ROOT 
#include <TError.h>
#include <TString.h> 
#include <TH1.h> 
// stdlib
#include <cmath> 
#include <stdexcept> 
#include <cstdlib> 
#include <cstdio> 
#include <algorithm> 
#include <iostream> 

namespace peak_search 
{

//_________________________________________________________________________________________________________________
FitTestThreadManager::FitTestThreadManager(size_t thread_id, const FitTestFunction& fcn, FitTestKernel* parent, const std::vector<TH1D*>& f_TH1D, const std::vector<TH2D*>& f_TH2D)
    : fThreadId{thread_id}, fParent{parent}, fTestFcn{fcn}
{   
    //seed the random number generator with the current thread-id
    fMyRand = std::make_unique<TRandom3>(thread_id+1);
    
    //initialize thread-local copies of histograms
    
    // TH1D
    for (const auto& h : f_TH1D) {
        if (!h) {
            Fatal(__func__, "Histogram in list passed is null."); 
            return; 
        }
        make_threadlocal_hist_copy(h, htype::kTH1D); 
    }

    // TH2D
    for (const auto& h : f_TH2D) {
        if (!h) {
            Fatal(__func__, "Histogram in list passed is null."); 
            return; 
        }
        make_threadlocal_hist_copy(h, htype::kTH2D); 
    }
}
//_________________________________________________________________________________________________________________
void FitTestThreadManager::run_test(double mass)
{
    fMass = mass; 
    //run the test
    SignalFit result; 
    try {
        fTestFcn(this);
    } catch (const std::exception& e) {
        Error("FitTestThreadManager::run_test", "<thread: %zi> Erorr caught running test. what(): %s", fThreadId, e.what()); 
        std::exit(1); 
    }

/*    //compute Z
    double Q0 = result.Q0;  
    double Z = (Q0 > 0. ? +1 : -1) * std::sqrt( std::fabs(Q0) );

    double pQ0 = compute_Q0_p(Q0);
*/ 
}
//_________________________________________________________________________________________________________________
Histo1D FitTestThreadManager::get_spectrum(double xmin, double xmax)
{
    auto spectrum = fParent->GetSpectrum(xmin, xmax, fMyRand.get()); 

    return spectrum; 
}
//_________________________________________________________________________________________________________________
void FitTestThreadManager::make_threadlocal_hist_copy(TObject* source, htype type)
{
    if (!source) {
        Error(__func__, "source hist is null.\n"); 
        return; 
    }

    const auto copy_name = Form("%s_t%zi", source->GetName(), fThreadId); 

    TH1* copy = dynamic_cast<TH1*>(source->Clone(copy_name)); 

    // tell ROOT that we're going to won this object, and we will manage its memory allocation 
    copy->SetBit(kMustCleanup); 
    copy->ResetBit(kCanDelete); 
    copy->SetDirectory(nullptr); 

    // the maps below index our own thread-local copies of the histograms with the ptrs to the original (which hopefully will be stable over the life of the app...)
    switch (type) {

        case htype::kTH1D : { 
            fTH1D[dynamic_cast<const TH1D*>(source)] = std::unique_ptr<TH1D>( dynamic_cast<TH1D*>(copy) ); 
            break; 
        }

        case htype::kTH2D : { 
            fTH2D[dynamic_cast<const TH2D*>(source)] = std::unique_ptr<TH2D>( dynamic_cast<TH2D*>(copy) ); 
            break;
        }

        default : Break(__func__, "Unknown histogram type"); 
    }
}
//_________________________________________________________________________________________________________________
TH1D* FitTestThreadManager::GetUserTH1D(const TH1D* source)
{
    if (!source) {
        Break(__func__, "Hist passed is null."); 
        return nullptr; 
    }

    auto find_it = fTH1D.find(source); 

    if (find_it == fTH1D.end()) {
        Break(__func__, "Requested local cpy of histogram '%s', but it doesn not exist in our list!", source->GetName()); 
        return nullptr; 
    }
    return find_it->second.get(); 
}
//_________________________________________________________________________________________________________________
TH2D* FitTestThreadManager::GetUserTH2D(const TH2D* source)
{
    if (!source) {
        Break(__func__, "Hist passed is null."); 
        return nullptr; 
    }

    auto find_it = fTH2D.find(source); 

    if (find_it == fTH2D.end()) {
        Break(__func__, "Requested local cpy of histogram '%s', but it doesn not exist in our list!", source->GetName()); 
        return nullptr; 
    }
    return find_it->second.get(); 
}
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________

};