

//Fcn1D headers
#include <Fcn1D/Fcn1D.hpp>
#include <Fcn1D/ExponentialLegendre.hpp>
#include <Fcn1D/ExponentialPoly.hpp>
//other project headers
#include <bininfo.hpp>
#include <fit_exponential_poly.hpp>
#include <fit_exponential_legendre.hpp>
#include <newton_optimizer.hpp> 
#include <gauss_integrate.hpp>
#include <fit_parameter.hpp> 
#include <compute_Q0.hpp>
#include <numbers.hpp>
#include <read_model_from_file.hpp>
//ROOT headers
#include <TH1D.h> 
#include <TCanvas.h> 
#include <TF1.h> 
#include <TAxis.h> 
#include <TGraph.h> 
#include <TLine.h> 
#include <TFile.h>
#include <TObject.h> 
#include <TPad.h> 
#include <TStyle.h> 
#include <TLegend.h> 
#include <Math/ProbFuncMathCore.h>
//stdlib headers
#include <cstdio> 
#include <functional> 
#include <iostream> 
#include <cmath> 
#include <thread> 
#include <stdexcept> 
#include <vector> 
#include <thread> 
#include <sstream> 
#include <map> 
#include <fstream> 

peak_search::Fcn1D *get_fcn(TH1D* hist, int order, std::string type); 

void fit_full_spectrum(std::string file_path, int model_order=6, std::string type="exp_legendre", std::string path_graphic="", std::string path_polynomial="", std::string histogram_name="h_m")
{
    std::cout << "getting events..." << std::flush;
    
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
    const double xmin{hist->GetXaxis()->GetXmin()}, xmax{hist->GetXaxis()->GetXmax()};


    double *bin_X, *bin_err, *bin_ones; 

    auto poly_model = get_fcn(hist, model_order, type);
    if (!poly_model) {
        Error(__func__, "Fit failed."); 
        return; 
    }

    //save the polynomial 
    if (!path_polynomial.empty()) {

        std::fstream outfile(path_polynomial, std::ios::out | std::ios::trunc);

        if (!outfile.is_open()) {
            Error(__func__, "Something went wrong trying to create polynomial output file in path '%s'. unable to open / create file.", path_polynomial.c_str());
            return; 
        } 

        outfile << "# file data: '"<< path_polynomial <<"'\n";
        outfile << "# model type: " << type << "\n"; 
        outfile << "# model order: " << model_order << "\n"; 

        int i_coeff=0; 
        for (double coefficient : poly_model->GetParams()) {
            outfile << Form("%-3i    %+18.8e\n", i_coeff, coefficient);
            ++i_coeff; 
        }
        outfile.close(); 

    }

    const size_t n_bins = hist->GetNbinsX(); 

    auto hist_p_chi2 = new TH1D("h_p_chi2",  Form("distribution of p(#chi^{2}) for bin residuals (model %s, %i-ord);bin residuals;counts",type.c_str(), model_order), n_bins/10, 0., 1.);
    


    bin_X = new double[n_bins]; 
    bin_err = new double[n_bins]; 
    bin_ones = new double[n_bins]; 

    auto xax = hist->GetXaxis(); 
    const double dx = xax->GetBinWidth(1); 

    for (size_t i=1; i<=n_bins; i++) {

        double x = xax->GetBinCenter(i);
        double expect = peak_search::gauss_integrate(*poly_model, x-dx/2., x+dx/2.);
        double actual = hist->GetBinContent(i); 

        double residual = (actual - expect)/std::sqrt( expect );   

        bin_X[i-1] = x; 
        bin_err[i-1] = residual; 
        //std::printf(" x: %6.1e  actual %6.1e    expect %6.1e    residual: %6.1e\n", x, actual, expect, (actual-expect)/expect); 

        //compute the p-value for the chi2 for each bin
        double p_chi2 = ROOT::Math::chisquared_cdf_c( residual*residual, 1 ); 
        hist_p_chi2->Fill( p_chi2 ); 
    }

    auto c_p = new TCanvas; 
    gStyle->SetOptStat(0); 
    hist_p_chi2->Draw("HIST,E"); 
    hist_p_chi2->GetYaxis()->SetRangeUser( 0., hist_p_chi2->GetMaximum()*2. );


    if (!path_graphic.empty()) {
        c_p->SaveAs(Form("%s_%s_%iord_p_chi2.png", path_graphic.c_str(), type.c_str(), model_order)); 
    }

    auto c_fit = new TCanvas("c_fit", "Fit residuals", 1000, 1200); 
    gStyle->SetOptStat(0); 
    c_fit->DivideRatios(1,2, {1.}, {0.50, 0.50}, 0.00,0.00); 
    auto c_fit_hist = c_fit->cd(1);
    c_fit_hist->SetTopMargin(0.);  
    c_fit_hist->SetLogy(1); 
    c_fit_hist->SetBottomMargin(0); 
    c_fit_hist->SetTopMargin(0.075); 
    

    hist->Scale(2.); 
    hist->GetYaxis()->SetRangeUser(1., hist->GetMaximum()*30); 
    hist->SetBinErrorOption(TH1::kPoisson); 
    hist->SetTitle("Invariant mass (accidental spectrum);;counts / MeV"); 
    hist->Draw("E"); 

    auto fit_fcn = new TF1("fit_fcn", [poly_model](double *X, double *par){ return (*poly_model)(X[0]); }, xmin, xmax, 0); 
    fit_fcn->SetLineColor(kBlue);
    fit_fcn->SetLineStyle(kDashed); 
    fit_fcn->Draw("SAME");

    auto legend = new TLegend(0.4,0.0, 0.8,0.15); 
    legend->AddEntry(fit_fcn, Form("%s, ord: %i",type.c_str(),model_order)); 
    legend->Draw(); 

    //std::vector<double> zeros(n_bins, 0.), ones(n_bins, 1.), bin_X(n_bins,0.), bin_err(n_bins,0.);  


    auto graph = new TGraph(n_bins, bin_X,bin_err); 
    
    auto c_fit_residual = c_fit->cd(2); 
    c_fit_residual->SetTopMargin(0.);
    
    graph->SetTitle(";mass (MeV);Fit Residuals / #sigma"); 
    graph->GetXaxis()->SetRangeUser(xmin, xmax); 
    graph->SetMarkerStyle(kOpenCircle); 
    graph->SetMarkerSize(0.75); 
    graph->Draw(); 

    auto line = new TLine(xmin,0., xmax,0.); 
    line->SetLineColor(kBlue); 
    line->SetLineStyle(kDashed); 
    line->Draw(); 

    if (!path_graphic.empty()) {
        c_fit->SaveAs(Form("%s_%s_%iord.png", path_graphic.c_str(), type.c_str(), model_order)); 
    }

    std::cout << "done.\n" << std::flush; 

    return; 
}

peak_search::Fcn1D* get_fcn(TH1D* hist, int order, std::string type)
{
    //fit a function to the given data. 

    // _____________________________________________________________________________________________________________________
    // Legendre polynomial 
    if (type == "exp_legendre") {

        auto result = peak_search::fit_exponential_legendre(hist, order); 
        if (result.status != peak_search::Status::kSuccess) {
            return nullptr; 
        }
        
        return static_cast<peak_search::Fcn1D*>(new peak_search::ExponentialLegendre(result.data));  
    }
    // _____________________________________________________________________________________________________________________

    // _____________________________________________________________________________________________________________________
    // Legendre polynomial 
    if (type == "exp_poly") {

        auto result = peak_search::fit_exponential_poly(hist, order); 
        if (result.status != peak_search::Status::kSuccess) {
            return nullptr; 
        }
        
        return static_cast<peak_search::Fcn1D*>(new peak_search::ExponentialPoly(result.data));  
    }
    // _____________________________________________________________________________________________________________________

    Error(__func__, "Invalid funciton type: %s", type.c_str()); 
    return nullptr; 
}