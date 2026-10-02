#ifndef make_brazil_flag_plot_h
#define make_brazil_flag_plot_h

// ROOT 
#include <TH2D.h> 
#include <TAxis.h> 
#include <TGraphErrors.h> 
#include <TLegend.h> 
#include <TGraph.h> 
// stdlib 
#include <vector> 
#include <string> 

/// @brief Make classic HEP-style brazilian flag plot, with median and quantiles specified below. 
/// @param hist histogram to convert to brazil-flag plot (quantiles scanned vertically)
/// @param legend_title legend of the title to put on the plot
/// @param cl_1 first central quantile (green)
/// @param cl_2 second central quantile (yellow)
void make_brazil_flag_plot(TH2D* hist,  std::string legend_title="", double cl_1=0.341344746069, double cl_2=0.477249868052)
{
    //go through each bin, and find the cumulative stats corresponding to each cl given. 
    const std::vector<double> levels{ 0.5-cl_2, 0.5-cl_1, 0.5, 0.5+cl_1, 0.5+cl_2 };
    
    auto x_ax = hist->GetXaxis(); const int n_bins_x = x_ax->GetNbins(); 
    auto y_ax = hist->GetYaxis(); const int n_bins_y = y_ax->GetNbins(); 

    const int n_bins = x_ax->GetNbins(); 

    std::vector<double> x; x.reserve(n_bins_x); 
    std::vector<double> y_cl1, y_err_cl1; y_cl1.reserve(n_bins_x); y_err_cl1.reserve(n_bins_x); 
    std::vector<double> y_cl2, y_err_cl2; y_cl2.reserve(n_bins_x); y_err_cl2.reserve(n_bins_x); 
    std::vector<double> y_med; y_med.reserve(n_bins_x); 

    double minval{+1e30}, maxval{-1e30}; 

    for (int bx=1; bx<=n_bins_x; bx++) {

        double integral =0.; 
        for (int by=1; by<=n_bins_y; by++) integral += hist->GetBinContent(bx,by); 
        
        auto find_cumulant = [hist,bx,y_ax,n_bins_y,integral](double p) {

            double cum=0.; 


            double bin_val; 
            for (int by=1; by<=n_bins_y; by++) {

                bin_val = hist->GetBinContent(bx,by) / integral; 
                
                if (cum + bin_val >= p) {
                    
                    double val = y_ax->GetBinCenter(by-1);
                    
                    double remainder = p - cum; 
                    val += ((remainder/bin_val) * y_ax->GetBinWidth(by)); 
                    
                    return val; 

                } else {
                    cum += bin_val; 
                } 
            }
            
            return -1.; 
        };

        double y1_lo = find_cumulant(0.5-cl_1);
        double y1_hi = find_cumulant(0.5+cl_1);
        
        double y2_lo = find_cumulant(0.5-cl_2);
        double y2_hi = find_cumulant(0.5+cl_2);

        minval = std::min(y2_lo, minval); 
        maxval = std::max(y2_hi, maxval); 

        double y_median = find_cumulant(0.5); 

        x.push_back(x_ax->GetBinCenter(bx)); 

        y_cl1    .emplace_back((y1_hi + y1_lo)/2.);
        y_err_cl1.emplace_back((y1_hi - y1_lo)/2.);

        y_cl2    .emplace_back((y2_hi + y2_lo)/2.);
        y_err_cl2.emplace_back((y2_hi - y2_lo)/2.);

        y_med.emplace_back(y_median); 
    }


    auto g2 = new TGraphErrors(n_bins_x, x.data(), y_cl2.data(), nullptr, y_err_cl2.data()); 
    g2->SetTitle(hist->GetTitle()); 
    
    g2->SetFillColor(kYellow);
    g2->Draw("A3"); 

    auto g1 = new TGraphErrors(n_bins_x, x.data(), y_cl1.data(), nullptr, y_err_cl1.data()); 
    g1->SetFillColor(kGreen);
    g1->Draw("3"); 

    auto gmed = new TGraph(n_bins_x, x.data(), y_med.data()); 
    gmed->SetLineStyle(kDashed); 
    gmed->SetLineWidth(2); 
    gmed->Draw("SAME"); 

    auto legend = new TLegend;
    if (!legend_title.empty()) legend->SetHeader(legend_title.c_str());  
    legend->AddEntry(g1, "#pm 1 #sigma");
    legend->AddEntry(g2, "#pm 2 #sigma");
    legend->AddEntry(gmed, "median");
    legend->Draw(); 

    return; 
}

#endif