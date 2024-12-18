// ________________________________________ //
// Author: Henrique Souza
// Filename: giveMeSphe.C
// Created: 2024
// ________________________________________ //
#include "WVFAnaSoLAr.h"



void giveMeSphe(){

    SPHE dark("spe");


    dark.led_calibration = true; // if external trigger + led was used, set true
                                 // start and finish will be the time of integration
    dark.just_a_test     = false; // well.. it is just a test, so `just_this` is the total waveforms analysed
    dark.just_this       = 200;
    dark.channel         = 44;
    dark.rootfile        = "analyzed.root";

    dark.start  = 6400;            // start the search for peaks or start the integration (led)
    dark.finish = 7200;        // fisish the search or finish the integration (led)

    dark.filter         = 8;   // one dimentional denoise filter (0 equal no filder)

    dark.get_wave_form = true; // for getting spe waveforms
    dark.mean_before   = 6000;   // time recorded before and after the peak found
    dark.mean_after    = 10000;
    dark.sphe_charge   = 3062.51; // charge of 1 and 2 p.e. (use fit_sphe.C)
    dark.sphe_charge2  = 7913.06;
    dark.sphe_std      = 1982.13;
    dark.spe_max_val_at_time_cut = 1e12; // after `time_cut`, the signal cannot be higher than this
                                       // this allows to remove after pulses
    dark.time_cut = 2000; // in ns seconds

    // coeficients to surround gaussian of 1 spe.
    // Gain             = (sphe_charge2 - sphe_charge)
    // spe's get events where charge < Gain*deltaplus  and charge < Gain/deltaminus
    // If deltaminus is set to zero, sphe_std*deltaplus will be used instead
    // This value can be checked with fit_sphe.C
    dark.deltaplus  = 1.;
    dark.deltaminus = 0;


    // Not so common to change
    dark.check_selection     = true; // uses(or not) variable `selection` to discard
    dark.withfilter          = true; // Integrate in the filtered waveform
    dark.hnbins = 50000;       // Output histogram bins and limits. Does not change this unless you will analyze alway the same device
    dark.hxmin  = 0;           // The fit function has the `rebin` argument to avoid changing this values
    dark.hxmax  = 0;
    dark.normalize_histogram = false; // will normalize histogram by the average value


    dark.giveMeSphe();



    gROOT->SetBatch(kFALSE);

   
}
