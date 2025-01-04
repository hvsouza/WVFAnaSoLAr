#include "WVFAnaSoLAr.h"

void fit_sphe_wave0(Int_t ch = 4, Int_t rebin = 1, Double_t deconv = 4){
    
    Calibration Cal;
    Cal.rebin = rebin;
    Cal.channel = ch;
    string histogram = "analyzed";
    Cal.rootFile = "sphe_histograms_Ch"+to_string(ch)+".root";

    Cal.make_free_stddevs = true; // starts with false, if good fitting, change to true
    Cal.searchParameters(histogram.c_str(), deconv, true); // give a first search in the parameters.

    Cal.deltaplus = 1;
    Cal.deltaminus = 0;

    // Cal.drawDebugLines = true;
    Cal.fit_sphe_wave(histogram.c_str(),false); // set true to make if you want to execute "searchParameters" inside here instead
    return;
}
void try_hard(Int_t ch = 0, Int_t minrebin = 1, Int_t maxrebin = 2, Int_t min_sigma = 4 Int_t max_sigma = 5){

    Calibration Cal;
    
    string histogram = "analyzed";
    Cal.rootFile = "sphe_histograms_Ch"+to_string(ch)+".root";

    Cal.make_free_stddevs = false; // starts with false, if good fitting, change to true

    Cal.sphe_fit_try_hard(nullptr, ch, minrebin, maxrebin,  min_sigma, max_sigma);


}
