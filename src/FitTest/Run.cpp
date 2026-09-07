
#include <FitTest/Run.hpp>
#include <FitTest/ThreadManager.hpp>
#include <Fcn1D/ExponentialPoly.hpp>
#include <read_model_from_file.hpp>
// ROOT
#include <TStopwatch.h> 
#include <TObject.h> 
#include <TH1D.h> 
#include <TAxis.h> 
#include <TH1.h> 
#include <TH2.h> 
#include <TClass.h> 
// stdlib
#include <thread> 
#include <memory>
#include <cstdio> 
#include <iostream> 
#include <cstdlib> 
#include <stdexcept> 
#include <chrono> 
#define DEBUG

namespace peak_search
{
namespace FitTest
{

namespace {
    constexpr char model_path[] = "data/models/exp_poly_19.dat";
    constexpr double fMinMass{140}, fMaxMass{280};

    std::unique_ptr<std::chrono::time_point<std::chrono::system_clock>> start_time; 
    // start system clock timer 
    void start_timer();
    // get elapsed time in seconds since start of timer 
    double get_elapsed_seconds(); 
}; 

/// @brief Draw progress bar 
/// @param progress progress amount (in range [0-1])
/// @param n_ticks number of 'ticks' to drawin the progress bar
/// @return formatted progress bar string
std::string progress_bar(double progress, int n_steps=100); 


/// @brief Copy result from one (thread-local) TObject into another
/// @param source copying data from here
/// @param dest appending data here 
void copy_result(TObject* source, TObject* dest); 


void Run(size_t n_trials, Configuration cfg, Outputs outputs, Function fcn, int run_verbosity)
{
    auto background_model = std::make_unique<ExponentialPoly>(std::vector<double>{}, fMinMass, fMaxMass); 

    try {

        peak_search::read_model_from_file(model_path, background_model.get());

    } catch (const std::exception& e) {

        Error(__func__, "Something went wrong trying to load model from file\n what(): %s", e.what()); 
        return; 
    }

    if (cfg.n_threads<1) { 
        cfg.n_threads = std::thread::hardware_concurrency(); 
    } else { 
        cfg.n_threads = std::min<size_t>(
            cfg.n_threads, 
            (size_t)std::thread::hardware_concurrency()
        ); 
    }

    std::vector<std::unique_ptr<ThreadManager>> thread_managers;

    std::vector<std::thread> threads; 

    thread_managers.reserve(cfg.n_threads); 
    threads        .reserve(cfg.n_threads);
    
    std::mutex scheduler_mutex; 

    size_t scans_done=0; 

    if (run_verbosity>0) std::cout << "\n staring " << n_trials << " trials...\n"; 
    TStopwatch stopwatch; 

    size_t trials_scheduled{0}; 

    
    unsigned long steps_scheduled{0}; 
    const unsigned long steps_per_task = cfg.params.GetNSteps(); 
    
    start_timer(); 
    for (size_t t=0; t<cfg.n_threads; t++) {

        //create the thread manager 
        thread_managers.emplace_back(std::make_unique<ThreadManager>(t, cfg, fcn, background_model.get(), outputs.GetPtrs())); 
        auto& manager = thread_managers.back(); 

        threads.emplace_back([&manager, &scheduler_mutex, &trials_scheduled,n_trials,  &steps_scheduled,steps_per_task,&cfg, t, run_verbosity]{

            //keep running scans until all the scans are done. 
            while (1) {

                //see if there's any new trials left to do.  
                if (trials_scheduled >= n_trials) break; 

                scheduler_mutex.lock();

                unsigned long step_0 = steps_scheduled; 
                unsigned long step_1 = std::min( step_0 + cfg.n_steps_per_task, steps_per_task ); 

                steps_scheduled = step_1; 

                if (run_verbosity>=1) {

                    double fraction_done = ((double)(steps_per_task*trials_scheduled + steps_scheduled))/((double)steps_per_task*n_trials); 
                    if (run_verbosity==1)
                        std::cout << "\r" << progress_bar(fraction_done, 100) << std::flush; 

                    if (run_verbosity>=2) {
                        std::printf("\n"
                            "thread: %2zi Scheduled trial: %3zi, steps [%5lu - %5lu]. %5.1f%%\n",
                            t, 
                            trials_scheduled,
                            step_0, step_1-1, fraction_done*100.
                        );
                    }
                }

                // if all the steps in this trial are done, move on to the next trial. 
                if (steps_scheduled >= steps_per_task) {
                    ++trials_scheduled;  
                    steps_scheduled=0; 
                }
                

                scheduler_mutex.unlock(); 

                try {
                    manager->ExecuteStepRange(step_0, step_1); 
                } catch (const std::exception& e) {
                    Error(__func__, "Exception caugt on trial %zi, step range [%lu-%lu]: %s", trials_scheduled, step_0,step_1, e.what()); 
                    std::exit(1); 
                }
            }
        }); 
    }
    //now, we will wait for all threads to be done. 
    for (auto& thread : threads) thread.join(); 

    if (run_verbosity==1) std::cout << "\r" << progress_bar(1., 100) << "\n"; 

    double cputime = stopwatch.CpuTime(); 
    double realtime = stopwatch.RealTime(); 

    if (run_verbosity>0)
        std::printf("done.\nReal time elapsed: %.3f s, %.3f s cpu time (%.4f ms / step)\n",
            realtime, cputime, 1e3*cputime/((double)n_trials*steps_per_task)
        );


    //now, we can add up sub-results for each histogram. 
    for (auto& manager : thread_managers) {
        
        // for each thread-manager, copy the thread-local results to the global result. 
        //TH1D
        auto& output_list = outputs.GetPtrs(); 
        for (size_t id=0; id<output_list.size(); id++) {
            copy_result(manager->GetOutput<TObject>(id), output_list[id]); 
        }
    }

    //all done! 
    // because the collection of user thread-managers is a vector constructed in the scope of this function, 
    // they will all automatically be deleted as this function exits (right now). 
} 
namespace 
{
//__________________________________________________________________________________________________________________________
void start_timer()
{
    using time_point = std::chrono::time_point<std::chrono::system_clock>; 
    using duration = std::chrono::duration<double, std::milli>; 

    if (!start_time) {
        start_time = std::make_unique<time_point>(std::chrono::system_clock::now()); 
    }
}
//__________________________________________________________________________________________________________________________
double get_elapsed_seconds()
{
    using time_point = std::chrono::time_point<std::chrono::system_clock>; 
    using duration = std::chrono::duration<double, std::milli>; 

    auto now = std::chrono::system_clock::now(); 

    if (!start_time) return numbers::nan; 

    return duration{ now - *start_time.get() }.count() / 1.e3; 
}
}
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
//__________________________________________________________________________________________________________________________
std::string progress_bar(double progress, int n_ticks) 
{

    std::ostringstream oss; 
    oss  << "["; 
    double n_steps_d = (double)n_ticks;
    int n_ticks_full = ((double)n_ticks)*progress;

    for (int i=0; i<n_ticks; i++) { oss << (i<=n_ticks_full ? "=" : " "); }

    oss << Form("]    %4.1f%%", progress*100.); 

    //measure the time, and report it. 
    double seconds = get_elapsed_seconds(); 

    int minutes = std::floor( seconds / 60 ); 
    seconds = seconds - minutes*60; 
    oss << Form("   elapsed: (%02i:%02.0f)", minutes, seconds); 

    return oss.str(); 
}
//__________________________________________________________________________________________________________________________
void copy_result(TObject* source, TObject* dest)
{
    const auto src_class = source->IsA(); 
    const auto dst_class = source->IsA(); 

    if ( src_class != dst_class ) {
        Error(__func__, "Source and dest. classes do not match! source class: %s, dest. class: %s", src_class->GetName(), dst_class->GetName());
        return; 
    }

    if (src_class->InheritsFrom( TClass::GetClass<TH2>() )) {
        auto s = dynamic_cast<TH2*>(source); 
        auto d = dynamic_cast<TH2*>(dest); 
        
        auto xax = s->GetXaxis(); 
        auto yax = s->GetYaxis(); 
        for (int bx=1; bx<=xax->GetNbins(); bx++) {
            for (int by=1; by<=yax->GetNbins(); by++) {
                d->Fill( xax->GetBinCenter(bx), yax->GetBinCenter(by), s->GetBinContent(bx, by) );
            }
        }    
        return; 
    }; 

    if (src_class->InheritsFrom( TClass::GetClass<TH1>() )) {
        auto s = dynamic_cast<TH1*>(source); 
        auto d = dynamic_cast<TH1*>(dest); 
        
        auto xax = s->GetXaxis(); 
        for (int bx=1; bx<=xax->GetNbins(); bx++) {
            d->Fill( xax->GetBinCenter(bx), s->GetBinContent(bx) );
        }    
        return; 
    }; 

    Error(__func__, "Type of output TObject passed is not supported. Type: %s", src_class->GetName()); 
    return; 
}
//__________________________________________________________________________________________________________________________

}
}