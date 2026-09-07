
#include <FitTest/ThreadManager.hpp>
#include <compute_Q0.hpp>
#include <generate_toy_events.hpp>
// ROOT 
#include <TError.h>
#include <TString.h> 
#include <TH1.h> 
#include <TClass.h> 
// stdlib
#include <cmath> 
#include <stdexcept> 
#include <cstdlib> 
#include <cstdio> 
#include <algorithm> 
#include <iostream> 

namespace peak_search 
{
namespace FitTest
{

//_________________________________________________________________________________________________________________
ThreadManager::ThreadManager(
        size_t thread_id, 
        const Configuration& config,
        const Function& fcn, 
        Fcn1D* background_model, 
        const std::vector<TObject*>& outputs)
    : fThreadId{thread_id}, fTestFcn{fcn}, fParamList{config.params}, fBackgroundModel{background_model}, fStats{config.total_stats}, fConfig{config}
{   
    //seed the random number generator with the current thread-id
    fMyRand = std::make_unique<TRandom3>(thread_id+1);
    
    //initialize thread-local copies of histograms
    for (auto& output : outputs) AddOutput(output); 
}
//_________________________________________________________________________________________________________________
void ThreadManager::AddOutput(TObject* source)
{
    if (!source) {
        throw std::invalid_argument("in <ThreadManager::AddOutput>: source TObject is null."); 
        return; 
    }

    //make a new copy 
    const auto src_name = source->GetName(); 
    const auto cpy_name = Form("%s_t%zi", src_name, fThreadId); 
    fOutputs.emplace_back( source->Clone(cpy_name) ); 

    //set some options with our new copy, to make sure we're the only ones who can delete it 
    auto& copy = fOutputs.back(); 

    copy->SetBit(kMustCleanup); 
    copy->ResetBit(kCanDelete); 

    //if this is a histogram, we need to set it's directory to null, so that it doesn't live and die with any particular TFile 
    if ( copy->IsA()->InheritsFrom( TClass::GetClass<TH1>() ) )
    {
        dynamic_cast<TH1*>(copy.get())->SetDirectory(nullptr); 
    }
}
//_________________________________________________________________________________________________________________
Histo1D ThreadManager::GetSpectrum(size_t n_bins, double m_min, double m_max)
{
    m_min = std::max(m_min, fMinMass);
    m_max = std::min(m_max, fMaxMass);

    double bin_size = (m_max - m_min)/((double)n_bins); 

    double m_center = (m_min + m_max)/2.; 
    double m_span   = ((double)n_bins)*bin_size; 

    Histo1D hist; 
    hist.bins.reserve(n_bins);

    double m = m_center - m_span/2.; 

    for (int i=0; i<n_bins; i++) { 
        hist.bins.emplace_back( m, m+bin_size, 0. ); 
        m += bin_size; 
    }

    //now, generate the toy events 
    generate_toy_events(hist, fBackgroundModel, fStats, *fMyRand.get()); 

    return hist; 
}
//_________________________________________________________________________________________________________________
template<typename T> T* ThreadManager::GetOutput(size_t id)
{
    if (id >= fOutputs.size()) {
        throw std::invalid_argument(Form("in <ThreadManager::GetOutput>: output id %zi requested is invalid; valid range is [0,%zi]", id, fOutputs.size()-1)); 
        return nullptr; 
    }

    return dynamic_cast<T*>(fOutputs[id].get()); 
}
// ----------------------------------------
// explicit template instantiations
template TH1D* ThreadManager::GetOutput(size_t); 
template TH2D* ThreadManager::GetOutput(size_t); 
template TObject* ThreadManager::GetOutput(size_t); 
//_________________________________________________________________________________________________________________
void ThreadManager::ExecuteStepRange(size_t step0, size_t step1)
{
    //first, check to make sure none of the tasks are out-of-bounds
    if (step0 >= fParamList.GetNSteps()) {
        throw std::invalid_argument(Form("in <ThreadManager::ExecuteSteps>: index of first step %lu requested is invalid; valid range is [0,%lu]", step0, fParamList.GetNSteps()-1)); 
        return; 
    }
    if (step1 >= fParamList.GetNSteps()) {
        throw std::invalid_argument(Form("in <ThreadManager::ExecuteSteps>: index of last step %lu requested is invalid; valid range is [0,%lu]", step1, fParamList.GetNSteps()-1)); 
        return; 
    }

    //now, loop over all steps. 
    fStep=step0; 
    for (; fStep<=step1; fStep++) { fTestFcn(this); }
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

}
}