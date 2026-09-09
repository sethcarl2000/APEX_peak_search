
// ROOT headers
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h> 
#include <TRandom3.h> 
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <TGraph.h>  
#include <Math/ProbFuncMathCore.h> 
#include <TAxis.h> 
#include <TLine.h> 
// stdlib 
#include <cmath> 
#include <vector> 

namespace {

    //range of drawing
    constexpr double xrange[] = { -5., +5. }; 
    constexpr double trange[] = { -50., +50. }; 
    

    constexpr int n_bins = 100; 

    constexpr double b_rate = 176e4; 
    constexpr double s_rate = 76e4; 

    // background parameters 
    
    // signal paramters (t_spectrum)
    // time 
    constexpr double st_sigma = 4.; 
    constexpr double st_mean  = 15.; 
    // pos
    constexpr double sx_sigma = 1.; 
    constexpr double sx_mean  = -3.; 
    
}

void demonstrate_prompt_extraction()
{

    TRandom3 rand(0);

    struct Event { double x, t; }; 

    auto gen_background = [&rand]() {
        
        Event evt; 
        
        //generate background time 
        evt.t = trange[0] + (trange[1]-trange[0])*rand.Rndm(); 
        
        //generate backgroud x-pos 
        do {
            evt.x = +1.5 + rand.Gaus()*2.;
            evt.x += +0.4*std::sin(evt.x / 5.); 
        } while (evt.x < xrange[0] || evt.x > xrange[1]); 
        
        return evt; 
    };

    auto gen_signal = [&rand]() {
        
        Event evt; 
        
        //generate background time 
        evt.t = st_mean + st_sigma*rand.Gaus(); 
        
        //generate backgroud x-pos 
        do {
            evt.x = +1.5 + rand.Gaus()*2.;
            evt.x += -0.4*std::sin(evt.x / 5.); 
        } while (evt.x < xrange[0] || evt.x > xrange[1]);

        return evt; 
    };

    auto hist_p_s = new TH1D("h_ps", "P_{prompt} - Signal;p_{prompt}", n_bins, 0., 1.); 
    auto hist_p_b = new TH1D("h_pb", "P_{prompt} - Background;p_{prompt}", n_bins, 0., 1.); 

    auto h_m_p = new TH2D("h_m_p", ";x;p_{prompt}", n_bins, xrange[0], xrange[1], n_bins, 0, 1); 

    using ull = unsigned long long; 

    auto h_sx = new TH1D("h_sx", "Prompt;x;", n_bins, xrange[0], xrange[1]); 
    auto h_st = new TH1D("h_st", "Prompt;t;", n_bins, trange[0], trange[1]); 

    auto h_bx = new TH1D("h_bx", "Accidental;x;", n_bins, xrange[0], xrange[1]); 
    auto h_bt = new TH1D("h_bt", "Accidental;t;", n_bins, trange[0], trange[1]); 
    
    auto get_p_prompt = [](double t) {
        using ROOT::Math::normal_cdf; 
        static constexpr double dx = (trange[1] - trange[0])/((double)n_bins); 
        
        double b = b_rate * (dx / (trange[1]-trange[0])); 

        double arg0 = (t-dx/2. - st_mean)/st_sigma; 
        double arg1 = (t+dx/2. - st_mean)/st_sigma; 
        double s = s_rate * ( normal_cdf(arg1, 1.) - normal_cdf(arg0, 1.) );  
        
        return (s / (s + b)); 
    };  

    auto h_lo = new TH1D("hc_lo", "", n_bins, xrange[0], xrange[1]); 
    auto h_hi = new TH1D("hc_hi", "", n_bins, xrange[0], xrange[1]); 

 
    auto h_coinc_x = new TH1D("h_coinc_x", "Expected coinc", n_bins, xrange[0], xrange[1]); 
    auto h_err_x = new TH1D("h_err_x",     "Expected coinc", n_bins, xrange[0], xrange[1]); 

    const double p_cut = 0.40;

    double ns_hi{0.}, ns_lo{0.};

    for (ull i=0; i<b_rate; i++) {

        auto evt = gen_background(); 
        double pi = get_p_prompt(evt.t); 
        
        if (pi > p_cut) { h_hi->Fill(evt.x); ns_hi += pi; }
        else            { h_lo->Fill(evt.x); ns_lo += pi; }

        h_bx->Fill( evt.x ); 
        h_bt->Fill( evt.t ); 
        hist_p_b->Fill( pi );  
        
        h_coinc_x->Fill( evt.x, pi ); 
        h_err_x  ->Fill( evt.x, pi*(2. - pi) ); 
    }
    
    for(ull i=0; i<s_rate; i++) {
        
        auto evt = gen_signal(); 
        double pi = get_p_prompt(evt.t); 
        
        if (pi > p_cut) { h_hi->Fill(evt.x); ns_hi += pi; } 
        else            { h_lo->Fill(evt.x); ns_lo += pi; }
        
        h_sx->Fill( evt.x ); 
        h_st->Fill( evt.t ); 
        hist_p_s->Fill( pi );

        h_coinc_x->Fill( evt.x, pi ); 
        h_err_x  ->Fill( evt.x, pi*(2. - pi) ); 
    }

    double ns_total = ns_lo + ns_hi; 
    double 
        nb_lo{h_lo->Integral() - ns_lo},
        nb_hi{h_hi->Integral() - ns_hi}; 

    double nb_total = nb_lo + nb_hi; 

    //find the p_coinc of both histograms
    double 
        p_lo{ns_lo/h_lo->Integral()}, 
        p_hi{ns_hi/h_hi->Integral()};
    
    auto xax = hist_p_s->GetXaxis(); 

    std::printf("Purity of 'lo' hist: %.4f\n", p_lo);
    std::printf("Purity of 'hi' hist: %.4f\n", p_hi);

    std::vector<double> bin_x; bin_x.reserve(n_bins); 
    xax = h_sx->GetXaxis();
    for (int i=1; i<=n_bins; i++) bin_x.push_back( xax->GetBinCenter(i) );
    
    std::vector<double> bin_avg; bin_avg.reserve(n_bins); 
    std::vector<double> bin_rms; bin_rms.reserve(n_bins);  
    for (int i=2; i<=n_bins; i++) {
        
        double 
            a{ns_hi/ns_total}, b{nb_hi/nb_total},
            c{ns_lo/ns_total}, d{nb_lo/nb_total}; 

        double det = a*d - b*c; 

        double 
            n_hi{h_hi->GetBinContent(i)}, 
            n_lo{h_lo->GetBinContent(i)}; 

        double n_s = ( d*n_hi - c*n_lo)/det; 
        double n_b = (-b*n_hi + a*n_lo)/det; 

        double N = n_s;

        bin_avg.push_back( N );
        bin_rms.push_back( std::sqrt(N) );
    }
    // make signal 
    auto stack_t = new THStack("stack_t", "Time"); 
    auto stack_x = new THStack("stack_x", "Position"); 

    new TCanvas; 
    stack_t->Add(h_bt); 
    h_st->SetFillColor(kBlue);
    h_st->SetFillStyle(3004); 
    stack_t->Add(h_st);  
    stack_t->Draw(); 

    new TCanvas; 
    h_sx->SetFillColor(kBlue);
    h_sx->SetFillStyle(3004); 
    stack_x->Add(h_sx); 
    stack_x->Add(h_bx); 
    stack_x->Draw(); 

    auto g = new TGraph(n_bins, bin_x.data(), bin_avg.data()); 
    g->SetFillColor(kGray); 
    g->Draw("SAME"); 

    new TCanvas; 
    auto stack_p = new THStack("stack_p", "Coinc values"); 
    hist_p_s->SetFillColor(kBlue);
    hist_p_s->SetFillStyle(3004); 
    stack_p->Add(hist_p_s); 
    stack_p->Add(hist_p_b); 
    stack_p->Draw(); 

    std::vector<double> bins_p, bins_frac; 
    for (int i=1; i<=hist_p_s->GetXaxis()->GetNbins(); i++) {

        double s = hist_p_s->GetBinContent(i); 
        double b = hist_p_b->GetBinContent(i); 
        
        if (s + b < 1.) continue; 

        bins_p.push_back(hist_p_s->GetXaxis()->GetBinCenter(i)); 
        bins_frac.push_back(s / (s + b)); 
    }
    new TCanvas; 
    auto gp = new TGraph(bins_p.size(), bins_p.data(), bins_frac.data()); 
    gp->SetTitle("esimated vs observed p_prompt;esimated p_prompt;measured p_prompt"); 
    gp->Draw(); 

    auto line = new TLine(0,0, 1,1); 
    line->SetLineStyle(kDashed); 
    line->SetLineColor(kRed); 
    line->Draw(); 

}
