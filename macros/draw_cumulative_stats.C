//ROOT
#include <TH1D.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TGaxis.h> 
#include <TStyle.h>
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TAxis.h> 
#include <TGraph.h> 
#include <TString.h> 
//stdlib
#include <string>
#include <vector> 
#include <cstdio> 
#include <fstream> 

/// @brief Draw cumulative stats for full replay, starting at p_coinc = 100%, and moving down towards p_coinc = 0%.
int draw_cumulative_stats(std::string path_infile, std::string path_output="")
{

    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("track_data", path_infile); 

    ROOT::RDF::Experimental::AddProgressBar(df); 

    const double x_min = 0.;
    const double x_max = 1.; 
    const int n_bins = 400; 

    auto hist = (TH1D*)df.Histo1D({"hh", "", n_bins, x_min, x_max}, "p_prompt")->Clone("h"); 

    auto xax = hist->GetXaxis(); 

    std::vector<double> pts_sum, pts_sum_coinc, pts_p, pts_frac_coinc; 
    pts_sum.reserve(xax->GetNbins());  
    pts_sum_coinc.reserve(xax->GetNbins()); 
    pts_p.reserve(xax->GetNbins()); 
    pts_frac_coinc.reserve(xax->GetNbins()); 

    const bool make_csv = !path_output.empty(); 

    struct cum_stats_t { double p_prompt, n_prompt, n_all; }; 

    std::vector<cum_stats_t> pts_cum; pts_cum.reserve(xax->GetNbins()); 

    //start at the end
    double cum =0.; 
    double cum_coinc =0.; 
    for (int bx=xax->GetNbins(); bx>=1; bx--) {

        double stats   = hist->GetBinContent(bx); 
        double p_coinc = xax->GetBinCenter(bx); 
        
        cum += stats; 
        
        if (cum < 1.) continue; 

        cum_coinc += stats * p_coinc; 

        pts_cum.emplace_back( p_coinc, cum_coinc, cum ); 

        pts_p  .push_back( p_coinc - xax->GetBinWidth(bx)/2. );
        pts_sum.push_back( cum );
        pts_sum_coinc.push_back( cum_coinc );
        pts_frac_coinc.push_back( cum_coinc/cum ); 
    }

    if (make_csv) {

        std::fstream outfile; 
        outfile.open(path_output, std::ios::out | std::ios::trunc); 
        
        outfile <<
            "\n" 
            "// Points measured from: '" << path_output << "'\n" 
            "// Each point represents the number of stats above a certain p_prompt threshold.\n"
            "\n"; 
        
        outfile << Form("static constexpr double p_min{%.8f}, p_max{%.8f};\n", pts_cum.back().p_prompt, pts_cum.front().p_prompt); 
        outfile << Form("static constexpr int n_pts{%zu}\n", pts_cum.size());

        outfile <<
            "\n"
            "struct cum_stats_t { double p_prompt, n_prompt, n_all; };\n";  

        outfile << 
            "\n"
            "std::vector<cum_stats_t> data{\n"; 
        
        for (int i=0; i<pts_cum.size(); i++) {
            
            outfile << Form("  {%.8e, %.8e, %.8e}%s\n", 
                pts_cum[i].p_prompt, pts_cum[i].n_prompt, pts_cum[i].n_all,
                (i<pts_cum.size()-1 ? "," : "") 
            );
        }

        outfile << 
            "};\n"; 

        outfile.close(); 
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


    auto graph_frac_coinc = new TGraph(
        pts_p.size(), 
        pts_p.data(), 
        pts_frac_coinc.data()
    ); 
    graph_frac_coinc->SetTitle("fraction of coinc. events;p_coinc;fraction of stats which are coinc"); 
    
    new TCanvas; 
    graph_frac_coinc->SetMinimum(0.);
    graph_frac_coinc->SetMaximum(1.1); 
    graph_frac_coinc->Draw(); 


    graph_frac_coinc->Eval(0.05); 
    
    /*
    std::printf(
        "evaluated cut: (p_coinc > %.3f)\n"
        " total stats   %.4e\n"
        " coinc events  %.4e (%.2f%% of total stats)\n",
        p_coinc_eval, 
        graph_sum->Eval(p_coinc_eval),
        graph_sum_coinc->Eval(p_coinc_eval),
        graph_frac_coinc->Eval(p_coinc_eval)*100.
    );*/ 
     
    

    std::vector<double> rev_sum, rev_sum_coinc, rev_p, rev_frac_coinc; 
    rev_sum.reserve(xax->GetNbins());  
    rev_sum_coinc.reserve(xax->GetNbins()); 
    rev_p.reserve(xax->GetNbins()); 
    rev_frac_coinc.reserve(xax->GetNbins()); 


    //start at the end
    cum =0.; 
    cum_coinc =0.; 
    for (int bx=1; bx<=xax->GetNbins(); bx++) {

        double stats   = hist->GetBinContent(bx); 
        
        if (stats <= 1.) continue; 
        
        double p_coinc = xax->GetBinCenter(bx); 
        
        cum += stats; 
        
        cum_coinc += stats * p_coinc; 

        rev_p  .push_back( p_coinc );
        rev_sum.push_back( cum );
        rev_sum_coinc.push_back( cum_coinc );
        rev_frac_coinc.push_back( cum_coinc/cum ); 
    }


    graph_sum = new TGraph(
        rev_p.size(), 
        rev_p.data(), 
        rev_sum.data()
    ); 
    graph_sum->SetMaximum(cum*1.1);
    graph_sum->SetMinimum(0.);

    graph_sum_coinc = new TGraph(
        rev_p.size(), 
        rev_p.data(), 
        rev_sum_coinc.data()
    ); 
    graph_sum->SetTitle("N. stats < p_coinc;p_coinc;N. stats >= p_coinc"); 
    
    graph_sum->SetLineColor(kBlack);
    graph_sum->SetLineWidth(2);

    graph_sum_coinc->SetLineColor(kRed);
    graph_sum_coinc->SetLineWidth(2);
    
    new TCanvas; 
    graph_sum->Draw(); 
    graph_sum_coinc->Draw("SAME");

    legend = new TLegend;
    legend->AddEntry(graph_sum, "All stats");
    legend->AddEntry(graph_sum_coinc, "Coinc");

    legend->Draw(); 




    return 0; 
}
