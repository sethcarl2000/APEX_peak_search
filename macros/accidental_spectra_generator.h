#ifndef accidental_spectra_generator_h
#define accidental_spectra_generator_h

#include <Fcn1D/ExponentialPoly.hpp>
#include <gauss_integrate.hpp> 
#include <read_model_from_file.hpp>
#include <Histo1D.hpp>
#include <generate_toy_events.hpp> 
// ROOT 
#include <TRandom3.h> 
#include <TString.h> 
// stdlib
#include <cstddef> 
#include <cmath> 
#include <stdexcept> 

namespace peak_search
{

class AccidentalSpectraGenerator {
private: 

    static constexpr char kModelPath[] = "/home/seth-hall/j_research_desktop/APEX_peak_search/data/models/exp_poly_19.dat"; 
    static constexpr double kMass_min{140}, kMass_max{280}; 

    std::unique_ptr<ExponentialPoly> fBackgroundFcn{nullptr}; 

    double fStats{100e6}; 

    TRandom3 fRand; 

public: 

    AccidentalSpectraGenerator(double stats=100e6, int rgen_seed=0) : fStats{stats} {

        fRand = TRandom3(rgen_seed); 
        //try to construct background model 
        try {
            fBackgroundFcn = std::make_unique<ExponentialPoly>(std::vector<double>{}, kMass_min, kMass_max); 

            read_model_from_file(kModelPath, fBackgroundFcn.get());

        } catch (const std::exception& e) {
            throw std::runtime_error(
                "in <"+std::string{__func__}+">: Something went wrong trying to load model from file\n"
                "what(): " + std::string{e.what()}); 
            return; 
        }
    }   

    ~AccidentalSpectraGenerator() = default; 

    /// @brief Generate a brand-new spectra. 
    /// @param n_bins number of (even-sized) bins
    /// @param min low edge of spectrum
    /// @param max hi edge of spectrum
    /// @return Histo1D with randomized data 
    Histo1D MakeNewSpectra(std::size_t n_bins, double min, double max) {

        if (min >= max) {
            throw std::invalid_argument(
                Form("in <%s>: x-range requested is invalid [%.2f,%.2f]",__func__, min, max) 
            );
            return Histo1D{}; 
        }
        // bound the mass
        min = std::max(min, kMass_min);
        max = std::min(max, kMass_max);

        double bin_size = (max-min)/((double)n_bins); 

        double m_center = (max+min)/2.; 
        double m_span   = ((double)n_bins)*bin_size; 

        Histo1D hist; 
        hist.bins.reserve(n_bins);

        double m = m_center - m_span/2.; 

        for (int i=0; i<n_bins; i++) { 
            hist.bins.emplace_back( m, m+bin_size, 0. ); 
            m += bin_size; 
        }

        //now, generate the toy events 
        generate_toy_events(hist, fBackgroundFcn.get(), fStats, fRand); 

        return hist; 
    }

    static constexpr double min_mass() { return kMass_min; }
    static constexpr double max_mass() { return kMass_max; }

    /// @brief Fill an already-created spectrum with events. 
    /// @param hist hist to fill 
    void FillSpectra(Histo1D& hist) {
        generate_toy_events(hist, fBackgroundFcn.get(), fStats, fRand); 
    }
}; 

}

#endif