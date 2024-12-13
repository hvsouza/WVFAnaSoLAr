// ________________________________________ //
// Author: Henrique Souza
// Filename: giveMeSphe_darkCount.C
// Created: 2021
// ________________________________________ //
#include "__USER_PATH__/MYCODES.h"



void giveMeSphe(){

    SPHE2 dark("spe");


    dark.led_calibration = true; // if external trigger + led was used, set true
                                 // start and finish will be the time of integration
    dark.just_a_test     = false; // well.. it is just a test, so `just_this` is the total waveforms analysed
    dark.just_this       = 200;
    dark.channel         = 1;
    dark.rootfile        = "analyzed.root";

    dark.nshow_range = {0,100}; // will save some debugging waveforms inside the range.
                                // Example, nshow_range = {0,2} will save waveforms 0, 1 and 2;

    dark.tolerance    = 36;      // n sigmas (smoothed) (not used for led)
                                // if `method` is set to `fix`, the threshold will be the absolute value of tolerance, no baseline is calculated
    dark.baselineTime = 5000;   // is the time to compute the baseline (not used for led)
                                  // If the `method` is set to `dynamic` the baseline is calculated over the range of baselineTime
                                  // and it is updated depending on the next point with a weigth of `influece`
                                  // If `method` is set to `static`, baseline is calculated once using baseLimit as cut
    dark.baseLimit    = 3;      // higher then this wont contribute to the baseline abs(baseLimit) (not used for led)

    dark.start  = 6400;            // start the search for peaks or start the integration (led)
    dark.finish = 7200;        // fisish the search or finish the integration (led)

    dark.timeLow        = 180;   // integration time before peak (not used for led)
    dark.timeHigh       = 640;  // integration time after peak (not used for led)
    dark.lowerThreshold = -999;  // threshold to detect noise (normal waveform) (not used for led)
    dark.maxHits        = 5;    // maximum hit before discarding   (not used for led)
    dark.too_big        = 2000;  // if there is a peak > "too_big" .. wait "waiting" ns for next peak
    dark.waiting        = 1000;
    dark.filter         = 8;   // one dimentional denoise filter (0 equal no filder)
    dark.interactions   = 60;   // for moving avarage filter (not used on led)
    dark.ninter         = 2;    // N times that moving average is applied
    dark.derivate_raw   = true;

    dark.dtime = 4.;            // time step in ns


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
    dark.social_distance = 5.;   // demands that there is a minimum distance of social_distance * timeHigh between 2 consecutive peaks found
    dark.method          = "derivative"; // `dynamic` or `derivative` evaluation of the baseline
                               // `fix` will not evaluate baseline and use raw threshold
                               // See tolerance, baselineTime and baselineLimit above

    dark.check_selection = true; // uses(or not) variable `selection` to discard wvfs
    dark.withfilter      = true; // Integrate in the filtered waveform

    dark.check_selection     = true; // uses(or not) variable `selection` to discard
    dark.withfilter          = true; // Integrate in the filtered waveform
    dark.integrate_from_peak = true; // Set true and the maximum will be searched inside the derivate region.
                                     // Otherwise the integral is done starting from the crossing negative zero of derivative

    dark.hnbins = 50000;       // Output histogram bins and limits. Does not change this unless you will analyze alway the same device
    dark.hxmin  = 0;           // The fit function has the `rebin` argument to avoid changing this values
    dark.hxmax  = 0;
    dark.normalize_histogram = false; // will normalize histogram by the average value


    dark.giveMeSphe();



    gROOT->SetBatch(kFALSE);

   
}
