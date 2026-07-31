
#include "ToyEventGenerator.h"

#include <bininfo.hpp>
#include <copy_subhist.hpp>
#include <fit_exponential_poly.hpp>
#include <fit_exponential_legendre.hpp>
#include <Fcn1D/Gauss.hpp>
#include <Fcn1D/FcnSum.hpp>
#include <newton_optimizer.hpp> 
#include <fit_parameter.hpp> 
#include <compute_Q0.hpp>
#include <numbers.hpp>
//ROOT headers
#include <TRandom3.h> 
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TAxis.h> 
#include <TGraph.h> 
#include <TLine.h> 
#include <TFile.h>
#include <TList.h> 
#include <TObject.h> 
#include <TPad.h> 
#include <TStyle.h> 
//stdlib headers
#include <cstdio> 
#include <functional> 
#include <iostream> 
#include <cmath> 
#include <thread> 
#include <stdexcept> 

#define VERBOSE 0

//this manages our temporary, drawn objects so that each new frame in a gif can have fresh objects. 
class DrawnObjectMgr {
private: 
    TList fDrawnObjects; 
public: 
    
    inline void Draw(TObject* obj, const char* option="") {
        obj->SetBit(kMustCleanup); 
        obj->ResetBit(kCanDelete);
        fDrawnObjects.Add(obj);
        obj->Draw(option);
    }   

    inline void Reset() { fDrawnObjects.Delete(); }
}; 

//get the mass resolution as a fcn of sigma
double GetSigma(double m); 

void scan_spectrum(std::string file_path, std::string output_gif="", std::string histogram_name="h_m")
{   
    //if the user provided a path, then make a gif. 
    const bool make_gif=(output_gif != ""); 

    double m_min{150}, m_max{270};
    double sigma=1.; 
    int n_steps=400; 

    double xmin{140.}, xmax{280};
    double window_size = sigma*7.; 

    //cx 03.22222222222222222222 
    // - Muon's comment (23 Jul 26)
    std::cout << 
        "getting events..." << std::flush;
    
    TH1D* hist;

    try {

        auto file = new TFile(file_path.c_str(), "READ");

        if (!file || file->IsZombie()) {
            Error(__func__, "Unable to open file: %s", file_path.c_str()); 
            return; 
        }

        hist = file->Get<TH1D>(histogram_name.c_str()); 

        if (!hist) {
            Error(__func__, "Could not find inv. mass histogram with name '%s' in file: %s", histogram_name.c_str(), file_path.c_str()); 
            return; 
        }

    } catch (const std::exception& e) {

        Error(__func__, "Something went wrong trying to load data.\n what(): %s", e.what()); 
        return; 
    }

    DrawnObjectMgr draw_mgr; 

    auto c_fit = new TCanvas; 
    gStyle->SetOptStat(0); 
    c_fit->DivideRatios(1,2, {1.}, {0.25, 0.75}, 0.00,0.00); 
    
    auto c_fit_hist = c_fit->cd(2);
    c_fit_hist->SetTopMargin(0.);  
    hist->GetYaxis()->SetRangeUser(0., hist->GetMaximum()*1.1); 
    hist->Draw("E"); 

    //this histogram will track the dist. of p(Q0). 
    auto hist_pQ0 = new TH1D("h_pQ0", "Dist. of p(Q0);p(Q0);", 50, 0., 1.); 


    std::cout << "done.\n" << std::flush; 

    const double dx = (xmax - xmin)/((double)hist->GetXaxis()->GetNbins()); 
    
    std::vector<double> pts_x0, pts_S, pts_sqrtQ0, pts_pQ0;

    pts_x0.reserve(n_steps);
    pts_S.reserve(n_steps);
    pts_sqrtQ0.reserve(n_steps);
    pts_pQ0.reserve(n_steps);

    double x0 = m_min - (m_max - m_min)/((double)n_steps-1);

    auto gaussian_fcn   = peak_search::Gauss(0, x0, sigma); 
    gaussian_fcn.Set_mu(0.); 
    gaussian_fcn.Set_x0(0.);

    for (int i_step=0; i_step < n_steps; i_step++) {

        x0 += (m_max - m_min)/((double)n_steps-1); 

        auto sub_data       = peak_search::copy_subhist(hist, x0-window_size, x0+window_size); 

        if (VERBOSE >= 2) {
            std::cout << "sub histogram xmax=" << x0-window_size << ", " << x0+window_size << "\n"; 
            std::cout << "bins:\n"; 
            for (auto& bin : sub_data.bins) {
                std::cout << "  " << bin.x << " | " << bin.N << "\n"; 
            }
        }

        auto exp_poly_result = peak_search::fit_exponential_legendre(sub_data, 4);

        if (!exp_poly_result) {
            Warning(__func__, "Poly. background fit for x0=%.3f failed.", x0); 
            continue; 
        }

        auto exp_poly = exp_poly_result.data; 

        auto fcn_s_plus_b   = peak_search::FcnSum(&gaussian_fcn, &exp_poly); 

        double Q0 = peak_search::compute_Q0(sub_data, fcn_s_plus_b); 

        if (peak_search::numbers::is_nan(Q0)) {
            Warning(__func__, "computed Q0 is nan"); 
            continue; 
        }

        double pQ0 = peak_search::compute_Q0_p(Q0);

        if (make_gif) {
            auto fcn_b = (std::function<double(double*,double*)>)[exp_poly,dx](double *x, double *par){ 
                //std::cout << "<tf1>: dx = " << dx << "\n"; 
                return dx * exp_poly(x[0]);
            };  
            auto tf1_b = new TF1("tf1_poly", fcn_b, x0-window_size, x0+window_size, 0); 
            tf1_b->SetLineColor(kBlue); 
            tf1_b->SetLineStyle(kDashed); 
            draw_mgr.Draw(tf1_b, "SAME");

            auto fcn_sb = (std::function<double(double*,double*)>)[fcn_s_plus_b,dx](double *x, double *par){ 
                //std::cout << "<tf1>: dx = " << dx << "\n"; 
                return dx * fcn_s_plus_b(x[0]);
            }; 
            auto tf1_sb = new TF1("tf1_poly", fcn_sb, x0-window_size, x0+window_size, 0); 
            tf1_sb->SetLineColor(kRed); 
            tf1_sb->SetLineStyle(kSolid); 
            draw_mgr.Draw(tf1_sb, "SAME");
            
            //draw a vertical line to show where the current best-fit mu is 
            auto vline = new TLine(x0,0., x0,hist->GetMaximum()); 
            vline->SetBit(kMustCleanup); 
            vline->SetLineColor(kBlue);
            vline->SetLineWidth(1);
            vline->SetLineStyle(kBlue);
            draw_mgr.Draw(vline); 

            c_fit->Modified();
            c_fit->Update(); 
            c_fit->SaveAs(std::string{output_gif + "+20"}.c_str()); 

            draw_mgr.Reset(); 
        }

        pts_x0      .emplace_back( x0 ); 
        pts_S       .emplace_back( fcn_s_plus_b.GetParams()[0] * GetSigma(x0) * 2.50662827463 ); 
        pts_sqrtQ0  .emplace_back( (Q0>0.?+1.:-1.) * std::sqrt( std::fabs(Q0) ) ); 
        pts_pQ0     .emplace_back( pQ0 ); 

        hist_pQ0->Fill( pQ0 ); 

        if (VERBOSE >= 1)
            std::printf(" x0=%.2f,  S(mu)=%.1f, sqrt(Q0)=%.2f\n", x0, pts_S[i_step], pts_sqrtQ0[i_step]); 

        
        gaussian_fcn.Set_x0(x0);
        gaussian_fcn.Set_mu(0.);
    }

    TGraph *graph; 
    TLine* line; 

    // 'mu' graph
    auto c_fit_graph = c_fit->cd(1); 
    c_fit_graph->SetBottomMargin(0.); 
    graph = new TGraph(n_steps, pts_x0.data(), pts_S.data());  
    graph->SetTitle("Best-fit S(mu) vs. m;m;best-fit S(mu)"); 
    graph->Draw(); 

    line = new TLine(pts_x0.front(),0., pts_x0.back(),0.); 
    line->SetLineColor(kRed); 
    line->SetLineStyle(kDashed); 
    line->Draw(); 

    // 'Q0' graph
    new TCanvas; 
    graph = new TGraph(n_steps, pts_x0.data(), pts_sqrtQ0.data()); 
    graph->SetTitle("Signficance (Z = #sqrt{|Q0|}) vs. x0;x0;#sqrt{|Q0|}"); 
    graph->Draw(); 

    line = new TLine(pts_x0.front(),0., pts_x0.back(),0.); 
    line->SetLineColor(kRed); 
    line->SetLineStyle(kDashed); 
    line->Draw(); 

    // 'pQ0' graph
    new TCanvas; 
    graph = new TGraph(n_steps, pts_x0.data(), pts_pQ0.data()); 
    graph->SetTitle(" p(Q0) vs. x0;x0;log_{10} p(Q0)"); 
    graph->Draw(); 
    gPad->SetLogy(1); 

    new TCanvas;
    hist_pQ0->SetMinimum(0.); hist_pQ0->SetMaximum( hist_pQ0->GetMaximum()*1.2 ); 
    hist_pQ0->Draw(); 

    return; 
}

double GetSigma(double m) {
    return 1.; 
}
