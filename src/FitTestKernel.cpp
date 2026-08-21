
#include <FitTestKernel.hpp> 
#include <read_model_from_file.hpp>
#include <generate_toy_events.hpp>
// ROOT
#include <TString.h> 
#include <TAxis.h> 
#include <TCanvas.h> 
#include <TStyle.h>
#include <TStopwatch.h> 
// stdlib
#include <thread> 
#include <cmath> 
#include <mutex> 
#include <string>
#include <sstream> 
#include <iostream> 

namespace {
    constexpr char model_path[] = "data/models/exp_poly_19.dat";
}

namespace peak_search 
{

/// @brief Draw progress bar 
/// @param progress progress amount (in range [0-1])
/// @param n_ticks number of 'ticks' to drawin the progress bar
/// @return formatted progress bar string
std::string progress_bar(double progress, int n_steps=100); 

//_________________________________________________________________________________________________________________
FitTestKernel::FitTestKernel(double min_fit_mass, double max_fit_mass, size_t n_steps)
    : fMinFitMass{min_fit_mass}, fMaxFitMass{max_fit_mass}, fN_steps{n_steps}
{
    const double max_signal_events = 40e3; 
    
    fBackgroundModel = std::make_unique<ExponentialPoly>(std::vector<double>{}, fMinMass, fMaxMass); 

    try {

        peak_search::read_model_from_file(model_path, fBackgroundModel.get());

    } catch (const std::exception& e) {

        Error(__func__, "Something went wrong trying to load model from file\n what(): %s", e.what()); 
        return; 
    }

    fModel_m_vs_mu = ROOT::RDF::TH2DModel{"h_signal", "Best-fit signal parameter '#mu' vs m;signal mass hypothesis (MeV);best-fit #mu", 
        (int)fN_steps/4, fMinFitMass, fMaxFitMass, 
        150, -max_signal_events, max_signal_events
    }; 

    const double max_signal_significance = 6.; 
    fModel_m_vs_Z = ROOT::RDF::TH2DModel{"h_Z", "Significance Z ~ #sqrt{Q0} vs m;signal mass hypothesis (MeV);Significance Z (n. #sigma)", 
        (int)fN_steps/4, fMinFitMass, fMaxFitMass, 
        150, -max_signal_significance, +max_signal_significance
    }; 

    fModel_m_vs_pQ0 = ROOT::RDF::TH2DModel{"h_pZ", "p(Q0) vs m;signal mass hypothesis (MeV);p(Q0)",  
        (int)fN_steps/4, fMinFitMass, fMaxFitMass, 
        50, 0., 1.
    }; 

    fModel_pQ0 = ROOT::RDF::TH1DModel{"h_pQ0", "p(Q0);p(Q0);", 50, 0., 1.}; 

    //construct histograms
    fHist_m_vs_mu   = construct_TH2D(fModel_m_vs_mu);
    fHist_m_vs_Z    = construct_TH2D(fModel_m_vs_Z);
    fHist_m_vs_pQ0  = construct_TH2D(fModel_m_vs_pQ0);

    fHist_pQ0       = construct_TH1D(fModel_pQ0);
}
//_________________________________________________________________________________________________________________
TH1D* FitTestKernel::construct_TH1D(const ROOT::RDF::TH1DModel& model)
{
    return new TH1D(model.fName, model.fTitle, 
        model.fNbinsX, model.fXLow, model.fXUp
    ); 
}
//_________________________________________________________________________________________________________________
TH2D* FitTestKernel::construct_TH2D(const ROOT::RDF::TH2DModel& model)
{
    return new TH2D(model.fName, model.fTitle, 
        model.fNbinsX, model.fXLow, model.fXUp, 
        model.fNbinsY, model.fYLow, model.fYUp
    ); 
}
//_________________________________________________________________________________________________________________
void FitTestKernel::RunTest(size_t n_scans, FitTestFunction test_function, size_t n_threads)
{ 
    if (n_threads<1) { n_threads = std::thread::hardware_concurrency(); }
    else { n_threads = std::min<size_t>(n_threads, (size_t)std::thread::hardware_concurrency); }

    std::vector<std::unique_ptr<FitTestThreadManager>> thread_managers;
    
    std::vector<std::thread> threads; 

    thread_managers.reserve(n_threads); 
    threads.reserve(n_threads);
    
    std::mutex scheduler_mutex; 

    size_t scans_done=0; 

    std::cout << "\n staring " << n_scans << " tests..." << std::flush; 
    TStopwatch stopwatch; 

    for (size_t t=0; t<n_threads; t++) {

        //create the thread manager 
        thread_managers.emplace_back(std::make_unique<FitTestThreadManager>(t, this)); 
        auto& manager = thread_managers.back(); 

        manager->fTestFcn = test_function; 

        //initialize 2D histograms
        auto init_histo_2d = [t,this](const ROOT::RDF::TH2DModel& model){
            auto new_model = model;  
            new_model.fName = Form("%s_%zi",model.fName.Data(),t); 
            auto ptr = std::unique_ptr<TH2D>(construct_TH2D(new_model)); 
            ptr->ResetBit(kCanDelete);
            ptr->SetBit(kMustCleanup);
            ptr->SetDirectory(nullptr);
            return ptr;  
        }; 
        
        manager->fHist_m_vs_mu   = init_histo_2d(fModel_m_vs_mu); 

        manager->fHist_m_vs_Z    = init_histo_2d(fModel_m_vs_Z); 
        
        manager->fHist_m_vs_pQ0  = init_histo_2d(fModel_m_vs_pQ0); 
        

        //initialize 1D histograms
        auto init_histo_1d = [&manager,t,this](const ROOT::RDF::TH1DModel& model){
            auto new_model = model;  
            new_model.fName = Form("%s_%zi",model.fName.Data(),t); 
            auto ptr = std::unique_ptr<TH1D>(construct_TH1D(new_model));
            ptr->ResetBit(kCanDelete);
            ptr->SetBit(kMustCleanup);
            ptr->SetDirectory(nullptr);
            return ptr;  
        }; 

        manager->fHist_pQ0 = init_histo_1d(fModel_pQ0);

        threads.emplace_back([this, &manager, &scheduler_mutex, &scans_done, n_scans, t]{

            //keep running scans until all the scans are done. 
            while (1) {

                scheduler_mutex.lock();
                std::cout << "\r" << progress_bar(((double)scans_done)/((double)n_scans), 100) << std::flush; 
                if (scans_done >= n_scans) { 
                    //all requested scans are already done. 
                    scheduler_mutex.unlock(); 
                    break; 
                } else {
                    //there's at least 1 scan left to do 
                    ++scans_done; 
                    scheduler_mutex.unlock(); 
                }


                double mass_step_size = (fMaxFitMass - fMinFitMass)/((double)fN_steps); 

                double mass = fMinFitMass; 

                for (size_t i=0; i<fN_steps; i++) {

                    manager->run_test(mass); 

                    mass += mass_step_size; 
                }
            }
        }); 
    }
    //now, we will wait for all threads to be done. 
    for (auto& thread : threads) thread.join(); 

    double cputime = stopwatch.CpuTime(); 
    double realtime = stopwatch.RealTime(); 

    std::printf("done.\nReal time elapsed: %.3f s, %.3f s cpu time (%.3f ms / scan)\n",
        realtime, cputime, 1000.*cputime/((double)n_scans)
    );


    //now, we can add up sub-results for each histogram. 
    for (auto& manager : thread_managers) {
        
        //copy 2d histograms
        copy_histogram(manager->fHist_m_vs_mu.get(),     fHist_m_vs_mu);
        copy_histogram(manager->fHist_m_vs_Z.get(),      fHist_m_vs_Z);
        copy_histogram(manager->fHist_m_vs_pQ0.get(),    fHist_m_vs_pQ0);
        
        //copy 1d histograms
        copy_histogram(manager->fHist_pQ0.get(), fHist_pQ0);
    }

    //all done! 
}
//_________________________________________________________________________________________________________________
void FitTestKernel::copy_histogram(TH1D* source, TH1D* target)
{
    auto xax = source->GetXaxis(); 
    for (int bin=1; bin<=xax->GetNbins(); bin++) {
        target->Fill( xax->GetBinCenter(bin), source->GetBinContent(bin) );
    }
}
//_________________________________________________________________________________________________________________
void FitTestKernel::copy_histogram(TH2D* source, TH2D* target)
{
    auto xax = source->GetXaxis(); 
    auto yax = source->GetYaxis(); 
    for (int bx=1; bx<=xax->GetNbins(); bx++) {
        for (int by=1; by<=yax->GetNbins(); by++) {
            target->Fill( xax->GetBinCenter(bx), yax->GetBinCenter(by), source->GetBinContent(bx, by) );
        }
    }    
}
//_________________________________________________________________________________________________________________
void FitTestKernel::DrawResults()
{
    new TCanvas; 
    gStyle->SetOptStat(0);
    fHist_m_vs_mu->Draw("col");

    new TCanvas; 
    fHist_m_vs_Z->Draw("col"); 

    new TCanvas;
    fHist_m_vs_pQ0->Draw("col"); 
    
    new TCanvas; 
    fHist_pQ0->Draw("HIST"); 
}
//_________________________________________________________________________________________________________________
Histo1D FitTestKernel::GetSpectrum(double m_min, double m_max, TRandom3* generator)
{
    m_min = std::max(m_min, fMinMass);
    m_max = std::min(m_max, fMaxMass);

    int n_bins = std::floor((m_max - m_min)/fBinSize);

    double m_center = (m_min + m_max)/2.; 
    double m_span   = ((double)n_bins)*fBinSize; 

    Histo1D hist; 
    hist.bins.reserve(n_bins);

    double m = m_center - m_span/2.; 

    for (int i=0; i<n_bins; i++) { 
        hist.bins.emplace_back( m, m+fBinSize, 0. ); 
        m += fBinSize; 
    }


    //now, generate the toy events 
    generate_toy_events(hist, fBackgroundModel.get(), fStats, *generator); 

    return hist; 
}
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________
//_________________________________________________________________________________________________________________

std::string progress_bar(double progress, int n_ticks) {
    std::ostringstream oss; 
    oss  << "["; 
    double n_steps_d = (double)n_ticks;
    int n_ticks_full = ((double)n_ticks)*progress;

    for (int i=0; i<n_ticks; i++) { oss << (i<=n_ticks_full ? "=" : " "); }

    oss << Form("]    %4.1f%%", progress*100.); 
    return oss.str(); 
}

}; 