// ________________________________________ //
// Author: Henrique Souza
// Filename: convert_root_file.C
// Created: 2024
// ________________________________________ //
#include <filesystem>  // C++17
#include "WVFAnaSoLAr.h"


void convert_root_file(string datadir = "./"){

    Read r;

    r.dtime = 16; // step time in ns
    r.nbits = 14; // this is used only for baseline histogram, changing it to 12 might help

    r.baselineTime = 16000; // time limit for baseline in ns
    r.startCharge = 2368;  // start integration time in ns
    r.chargeTime = 2600; // last time to integrate ns
    r.maxRange = 2480; // max range to search for amplitude peak
    r.fast = 200; // fprompt fast integration time
    r.slow = 1700; //fprompt slow integration time
    r.exclusion_baseline = 500; // filtered waveform (ADCs), anything above here will do +exclusion window while evaluating the baseline
    // r.exclusion_baselines = {35,15};
    r.exclusion_window = 1000; // time in ns that it will jump for baseline
    r.filter = 100; // denoise filter.. if filter = 0, no denoise is done.
    r.OnlySomeEvents = false; // Do you want only a few events? Choose it wisely (you can set stopEvent)
    r.stopEvent = 1000; //
    r.noBaseline = false; // if you dont want to calculate baseline
    // r.saveFilter = true;

    // r.channels = {1};
    r.channels = {44, 144};
    // r.channels = {};
    const char *logfilename = r.standard_log_file.c_str();

    system(Form("touch %s && rm %s", logfilename, logfilename));
    std::filesystem::path filePath = datadir;
    if (filePath.extension() == ".root") // Heed the dot.
    {
        string packdata = Form("echo %s > %s", datadir.c_str(), logfilename);
        system(packdata.c_str());
    }
    else
    {
        cout << "HERE..." << endl;
        string packdata = Form("/bin/ls -1 %s/*.root | grep -v analyzed.root > %s", datadir.c_str(), logfilename);
        system(packdata.c_str());
    }
    r.convert_from_root();


}

