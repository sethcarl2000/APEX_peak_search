//ROOT
#include <TH1D.h>
#include <ROOT/RDataFrame.hxx>
#include <TGaxis.h> 
#include <TStyle.h>
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TAxis.h> 
#include <TGraph.h> 
//stdlib
#include <string>
#include <vector> 
#include <cstdio> 


/// @brief Draw cumulative stats for full replay, starting at p_coinc = 100%, and moving down towards p_coinc = 0%.
int draw_cumulative_stats(std::string path_infile, const double p_coinc_eval=0.05)
{

    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("track_data", path_infile); 

    const double x_min = 0.;
    const double x_max = 1.; 
    const int n_bins = 200; 

    auto hist = (TH1D*)df.Histo1D({"hh", "", n_bins, x_min, x_max}, "p_coinc")->Clone("h"); 

    auto xax = hist->GetXaxis(); 

    std::vector<double> pts_sum, pts_sum_coinc, pts_p; 
    pts_sum.reserve(xax->GetNbins());  
    pts_sum_coinc.reserve(xax->GetNbins()); 
    pts_p.reserve(xax->GetNbins()); 

    //start at the end
    double cum =0.; 
    double cum_coinc =0.; 
    for (int bx=xax->GetNbins(); bx>=1; bx--) {

        double stats   = hist->GetBinContent(bx); 
        double p_coinc = xax->GetBinCenter(bx); 
        
        cum += stats; 
        cum_coinc += stats * p_coinc; 

        pts_p  .push_back( p_coinc );
        pts_sum.push_back( cum );
        pts_sum_coinc.push_back( cum_coinc );
    }

    auto pts_p_coinc = pts_sum_coinc; 

    //scale the sum of all coinc events 
    for (size_t i=0; i<pts_sum.size(); i++) { 
        pts_p_coinc[i] *= 1./pts_sum[i]; 
    }

    auto graph_sum = new TGraph(
        pts_p.size(), 
        pts_p.data(), 
        pts_sum.data()
    ); 

    auto graph_sum_coinc = new TGraph(
        pts_p.size(), 
        pts_p.data(), 
        pts_sum_coinc.data()
    ); 
    graph_sum->SetTitle("N. stats >= p_coinc;p_coinc;N. stats >= p_coinc"); 
    
    graph_sum->SetLineColor(kBlack);
    graph_sum->SetLineWidth(2);

    graph_sum_coinc->SetLineColor(kRed);
    graph_sum_coinc->SetLineWidth(2);
    
    new TCanvas; 
    graph_sum->Draw(); 
    graph_sum_coinc->Draw("SAME");

    auto legend = new TLegend;
    legend->AddEntry(graph_sum, "All stats");
    legend->AddEntry(graph_sum_coinc, "Coinc");

    legend->Draw(); 


    auto graph_p_coinc = new TGraph(
        pts_p.size(), 
        pts_p.data(), 
        pts_p_coinc.data()
    ); 
    graph_p_coinc->SetTitle("fraction of coinc. events;p_coinc;fraction of stats which are coinc"); 
    
    new TCanvas; 
    graph_p_coinc->SetMinimum(0.);
    graph_p_coinc->SetMaximum(1.1); 
    graph_p_coinc->Draw(); 


    graph_p_coinc->Eval(0.05); 
    
    std::printf(
        "evaluated cut: (p_coinc > %.3f)\n"
        " total stats   %.4e\n"
        " coinc events  %.4e (%.2f%% of total stats)\n",
        p_coinc_eval, 
        graph_sum->Eval(p_coinc_eval),
        graph_sum_coinc->Eval(p_coinc_eval),
        graph_p_coinc->Eval(p_coinc_eval)*100.
    );
        
    return 0; 
}
