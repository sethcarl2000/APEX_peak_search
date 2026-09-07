#ifndef peak_search_FitTest_ThreadManager_hpp
#define peak_search_FitTest_ThreadManager_hpp

#include <FitTest/Function.hpp>
#include <FitTest/Configuration.hpp>
#include <FitTest/ParameterList.hpp>
#include <FitTest/Outputs.hpp>
#include <FitTest/Run.hpp>
#include <Histo1D.hpp>
#include <Fcn1D/Fcn1D.hpp>
// ROOT 
#include <TRandom3.h> 
#include <TObject.h> 
// stdlib
#include <thread> 
#include <memory> 
#include <vector>
#include <mutex> 

namespace peak_search
{
namespace FitTest
{

class ThreadManager {
private: 

    friend void Run(size_t, Configuration, Outputs, Function, int);  

    size_t fThreadId; 

    Function fTestFcn; 

    ParameterList fParamList; 

    // current step 
    unsigned long fStep{0}; 

    //minimum / maximum mass values in the BG model 
    static constexpr double fMinMass{140}, fMaxMass{280};

    // list of our thread-local copies of user-provided outputs 
    std::vector<std::unique_ptr<TObject>> fOutputs;  

    Configuration fConfig; 

    double fStats; 

    Fcn1D* fBackgroundModel;

    //thread-local random number generator
    std::unique_ptr<TRandom3> fMyRand; 

    /// make thread-local copy of TObject. place it in the list of outputs.  
    void AddOutput(TObject* source); 

    // execute a series of steps, in the inclusive range: [step0, step1-1]
    void ExecuteStepRange(unsigned long step0, unsigned long step1); 

public: 

    //we've private-ed the constructor, so only the 'FitTestManager' can make copies of this object. 
    ThreadManager(
        size_t thread_id, 
        const Configuration& config, 
        const Function& fcn, 
        Fcn1D* background_model, 
        const std::vector<TObject*>& outputs
    );  

    /// @brief Provide a randomly-sampled histogram
    /// @param n_bins number of bins
    /// @param xmin low-edge of generated spectrum
    /// @param xmax hi-edge of generated spectrum 
    /// @return a Histo1D object with randomly-sampled bin values 
    Histo1D GetSpectrum(size_t n_bins, double xmin, double xmax); 

    /// @brief Get access to thread-local copy to user-provided output. 
    /// @tparam T Output type 
    /// @param id 'id' of output in list (See FitTest::Output class)
    /// @return Ptr to thread-local copy of output 
    template<typename T> T* GetOutput(size_t id); 

    /// @brief Get access to a user-defined parameter list for this step. 
    /// @return current value of parameter 
    std::vector<double> GetParamList() const { return fParamList.GetParamList(fStep); }
};

}
}


#endif