//ROOT
#include <TH1D.h>
#include <TH2D.h> 
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <TRandom3.h> 
#include <TGaxis.h> 
#include <TStyle.h>
#include <TCanvas.h> 
#include <TLegend.h> 
#include <TAxis.h> 
#include <TGraph.h> 
#include <TVector3.h>
#include <Math/Vector4D.h> 
#include <TString.h> 
#include <TGraph.h> 
#include <TColor.h> 
#include <TStyle.h> 
#include <TGraphErrors.h> 
//stdlib
#include <string>
#include <vector> 
#include <cstdio> 
#include <cmath> 

namespace {
    constexpr double me2 = 0.511*0.511; 

    constexpr double m_min = 140.; //MeV
    constexpr double m_max = 280.; //MeV 

    constexpr double dm = 1.00; //MeV 
}

/// @brief Draw cumulative stats for full replay, starting at p_coinc = 100%, and moving down towards p_coinc = 0%.
int draw_coinc_and_accidental(std::string path_infile, double threshold_p_accidental=0.25, double threshold_p_coinc=0.95)
{
    ROOT::EnableImplicitMT(); 

    ROOT::RDataFrame df("track_data", path_infile);

    ROOT::RDF::Experimental::AddProgressBar(df); 

    using FourVec = ROOT::Math::XYZTVector; 

    auto df_invmass = df    

        .Define("invariant_mass", [](const TVector3& Pp, const TVector3& Pe)
        { 
            
            double energy = std::sqrt( Pp.Mag2() + me2 ) + std::sqrt( Pe.Mag2() + me2 ); 

            TVector3 P = Pp + Pe; 

            //return (time-like) mass of Pe + PP, using (+---) metric 
            return std::sqrt( energy*energy - P.Mag2() ); 
        }, {"P_p", "P_e"}); 

    auto df_invmass_coinc = df_invmass.Filter([threshold_p_coinc](double p){ return p > threshold_p_coinc; }, {"p_coinc"}); 

    auto df_invmass_accidental = df_invmass.Filter([threshold_p_accidental](double p){ return p < threshold_p_accidental; }, {"p_coinc"}); 

    const int n_bins = std::floor((m_max - m_min)/dm); 

    auto h_coinc = df_invmass_coinc
        .Histo1D<double>({"hh_coinc", "Invariant mass;m_{#pm} (MeV);dP/dm (MeV^{-1})", n_bins, m_min, m_max}, "invariant_mass"); 

    
    auto h_accidental = df_invmass_accidental
        .Histo1D<double>({"hh_accidental", "Invariant mass;m_{#pm} (MeV);dP/dm (MeV^{-1})", n_bins, m_min, m_max}, "invariant_mass"); 

    
    //get *actual* purity of sample (for both)
    auto h_purity_coinc =  
        df.Filter([threshold_p_coinc](double p){ return p > threshold_p_coinc; }, {"p_coinc"})
        .Histo1D({"hhh_c", "", 200, 0., 1.}, "p_coinc");


    auto hist_frac_coinc = df_invmass_coinc
        .Histo1D<double>({"hhh", "", 200, 0., 1.}, "p_coinc"); 

    auto hist_frac_accidental = df_invmass_accidental
        .Histo1D<double>({"hhh", "", 200, 0., 1.}, "p_coinc"); 

    auto get_coinc_fraction = [](ROOT::RDF::RResultPtr<TH1D>& hist)
    {
        double sum = 0.; 
        double coinc = 0.; 
        for (int ib=1; ib<=200; ib++) {
            sum += hist->GetBinContent(ib);
            coinc += hist->GetBinContent(ib) * hist->GetXaxis()->GetBinCenter(ib);
        }
        
        return coinc / sum;
    };

    double coinc_frac_coinc = get_coinc_fraction(hist_frac_coinc); 

    double coinc_frac_accidental = get_coinc_fraction(hist_frac_accidental); 

    auto hist_coinc = dynamic_cast<TH1D*>(h_coinc->Clone("h_coinc")); 
    //normalize the histograms
    double N_coinc = h_coinc->Integral(); 

    auto hist_accidental = dynamic_cast<TH1D*>(h_accidental->Clone("h_accidental")); 
    double N_accidental = h_accidental->Integral(); 

    //now, invert the histograms to find the 'true' coinc and accidental 
    // this is based on the idea that both 'coinc' and 'accidental' histograms are linear combinations of the 'true' coinc and accidental PDFs. 
    // we can combine them this way. 
    // Suppose that histogram 'A' is 75% coinc, and 25% accidental, and hist 'B' is 1% coinc, and 99% accidental. 
    // we can find the 'true' number of coinc and accidental per bin via: 
    // 
    //      nA = 0.75*nc + 0.25*na
    //      nB = 0.01*nc + 0.99*na
    //
    // this is just a linear equation with 2 DoF, we just need to invert the matrix on the rhs to solve for 'nc' and 'na'.  
    std::vector<double> pts_m(n_bins,0.), pts_coinc(n_bins, 0.), pts_coinc_err(n_bins,0.), pts_accidental(n_bins,0.), pts_accidental_err(n_bins,0.); 

    hist_coinc->Scale( 1./(N_coinc*dm) );
    hist_accidental->Scale( 1./(N_accidental*dm) );

    auto xax = hist_coinc->GetXaxis(); 

    double maxval=0.; 

    double norm_coinc{0.}, norm_accidental{0.}; 

    {
        double a = coinc_frac_coinc; 
        double b = 1. - coinc_frac_coinc; 
        double c = coinc_frac_accidental; 
        double d = 1. - coinc_frac_accidental; 

        double denom = a*d - b*c; 

        std::printf(
            "original matrix:\n"
            "   %5.3f   %5.3f\n"
            "   %5.3f   %5.3f\n",
            a, b, 
            c, d
        ); 

        std::printf(
            "inverted matrix:\n"
            "   %5.3f   %5.3f\n"
            "   %5.3f   %5.3f\n",
            d/denom, -c/denom, 
            -b/denom, a/denom
        ); 
    
        for (int i=1; i<=n_bins; i++) {
    
            double m = xax->GetBinCenter(i); 
            pts_m[i-1] = m; 

            double fA = hist_coinc->GetBinContent(i);
            double fB = hist_accidental->GetBinContent(i); 
    
            double fc = (  d*fA / denom ) + ( -c*fB / denom ); 
            double fa = ( -b*fA / denom ) + (  a*fB / denom ); 

            //std::printf(" mass %4.1f    fc = %.3e,   fa = %.3e\n", m, fa, fc); 

            pts_coinc[i-1]      = fc; 
            pts_accidental[i-1] = fa; 

            double variance_A = fA/(N_coinc*dm);  
            double variance_B = fB/(N_accidental*dm);

            pts_coinc_err[i-1]      = 2.*std::sqrt( std::pow(d/denom, 2)*variance_A + std::pow(c/denom,2)*variance_B );
            pts_accidental_err[i-1] = 2.*std::sqrt( std::pow(b/denom, 2)*variance_A + std::pow(a/denom,2)*variance_B );
            
            norm_coinc += fc*dm; 
            norm_accidental += fa*dm; 
        }

        for (auto& x : pts_coinc) { x *= 1./norm_coinc; if (x > maxval) maxval=x; }
        for (auto& x : pts_accidental) { x *= 1./norm_accidental; if (x > maxval) maxval=x; }

        for (auto& x : pts_coinc_err) { x *= 1./norm_coinc; if (x > maxval) maxval=x; }
        for (auto& x : pts_accidental_err) { x *= 1./norm_accidental; if (x > maxval) maxval=x; }
    }

    TLegend* legend; 
    legend = new TLegend; 
    legend->SetHeader("PDFs");

    new TCanvas; 
    auto graph_coinc_err = new TGraphErrors(n_bins, pts_m.data(),pts_coinc.data(), nullptr, pts_coinc_err.data());
    graph_coinc_err->SetFillColor(kGray); 
    graph_coinc_err->GetXaxis()->SetRangeUser(m_min, m_max);
    //graph_coinc_err->Draw("a3"); 
    
    auto graph_coinc = new TGraph(n_bins, pts_m.data(),pts_coinc.data()); 
    graph_coinc->SetLineColor(kBlack); 
    graph_coinc->SetLineStyle(kDashed); 
    graph_coinc->SetTitle("Accidental and Coinc PDF;mass (MeV);dP/dm (MeV^{-1})"); 
    graph_coinc->GetXaxis()->SetRangeUser(m_min, m_max);
    graph_coinc->GetYaxis()->SetRangeUser(0., maxval*1.2);  
    graph_coinc->Draw(); 
    
    legend->AddEntry(graph_coinc, "coinc");

    auto graph_accidental_err = new TGraphErrors(n_bins, pts_m.data(),pts_accidental.data(), nullptr, pts_accidental_err.data());
    graph_accidental_err->SetFillColor(kGray); 
    graph_accidental_err->GetXaxis()->SetRangeUser(m_min, m_max); 
    //graph_accidental_err->Draw("Same3a"); 

    auto graph_accidental = new TGraph(n_bins, pts_m.data(),pts_accidental.data()); 
    graph_accidental->SetLineColor(kBlack); 
    graph_accidental->SetLineStyle(kSolid); 
    graph_accidental->Draw("SAME"); 

    legend->AddEntry(graph_accidental, "accidental");
    
    //legend->AddEntry(graph_coinc_err, "#pm 2#sigma");
    legend->Draw(); 

    //make a graph of the ratio of bin-counts
    std::vector<double> pts_m_ratio, pts_ratio, pts_stddev; 

    for (auto v : {&pts_m_ratio, &pts_ratio, &pts_stddev}) v->reserve(n_bins); 

    auto hist_ratio = new TH2D("h_ratio", "Ratio of coinc / accidental PDF;mass (MeV); dN/dm (coinc) / dN/dm (accidental)", n_bins/2, m_min, m_max, 300, 0., 5.5); 

    maxval=0.;
    double minval=1e9; 
    double min_ratio=1e-9; 
    for (int i=1; i<=n_bins; i++) {

        double m = xax->GetBinCenter(i); 
    
        double coinc{pts_coinc[i-1]}, accident{pts_accidental[i-1]}; 
    
        double ratio = coinc / accident; 
        if (ratio < min_ratio) continue; 

        if (ratio > maxval) maxval=ratio;
        if (ratio < minval) minval=ratio; 

        pts_m_ratio .push_back( m ); 
        pts_ratio   .push_back( ratio ); 
    
        double coinc_err{pts_coinc_err[i-1]/2.}, accident_err{pts_accidental_err[i-1]/2.}; 

        double stddev = std::sqrt( std::pow(coinc_err * 1./accident,2) + std::pow(accident_err * coinc/(accident*accident),2) ); 

        pts_stddev.push_back( stddev ); 
    }
    
    new TCanvas; 
    gStyle->SetPalette(kGreyScale);
    TColor::InvertPalette(); 

    legend = new TLegend; 

    //hist_ratio->Draw("COL2"); 
    auto graph_err = new TGraphErrors(pts_m_ratio.size(), pts_m_ratio.data(), pts_ratio.data(), nullptr, pts_stddev.data()); 
    graph_err->SetFillColor(kGray); 
    graph_err->SetTitle("Ratio of coinc / accidental;mass (MeV);dN/dm (coinc) / dN/dm (accidental)"); 
    graph_err->GetXaxis()->SetRangeUser(m_min, m_max); 
    graph_err->GetYaxis()->SetRangeUser(minval/10., maxval*10.); 
    graph_err->Draw("a3"); 

    auto line = new TLine(m_min,1., m_max,1.); 
    line->SetLineColor(kBlack); 
    line->SetLineStyle(kDashed);
    line->Draw(); 
    
    auto graph_center = new TGraph(pts_m_ratio.size(), pts_m_ratio.data(), pts_ratio.data()); 
    graph_center->Draw("SAME"); 
    
    legend->AddEntry(graph_center, "mean"); 
    legend->AddEntry(graph_err, "#pm 2#sigma"); 
    legend->Draw(); 
    gPad->SetLogy(1); 

    new TCanvas; 
    gStyle->SetOptStat(0); 
    
    hist_coinc->SetLineColor(kRed);
    hist_coinc->SetFillColor(kRed);
    hist_coinc->SetFillStyle(3004);
    hist_coinc->SetMaximum( std::max<double>(hist_coinc->GetMaximum(), hist_accidental->GetMaximum())*1.1 ); 
    hist_coinc->Draw("HIST");
    std::printf("stats in coinc histogram: %.4e\n", N_coinc); 

    hist_accidental->SetLineColor(kBlack);
    hist_accidental->SetFillStyle(0);
    hist_accidental->Draw("HIST, SAME");
    std::printf("stats in accidental histogram: %.4e\n", N_accidental); 

    legend = new TLegend; 
    legend->SetHeader("frac. of coinc events (p)");
    legend->AddEntry(hist_coinc,      Form("p = %.2f",coinc_frac_coinc));
    legend->AddEntry(hist_accidental, Form("p = %.2f",coinc_frac_accidental));
    legend->Draw(); 

    return 0; 
};