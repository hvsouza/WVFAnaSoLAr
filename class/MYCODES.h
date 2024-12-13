// ________________________________________ //
// Author: Henrique Souza
// Filename: MYCODES.h
// Created: 2021
// ________________________________________ //

#ifndef MYCODES
#define MYCODES

#include <fstream>
#include <iostream>
#include "Riostream.h"
#include <time.h>       // time_t, struct tm, difftime, time, mktime
#include <cmath>
#include <numeric>
#include <unistd.h>
#include <chrono>
#include <thread>
#include <cstring> // include the <cstring> header file for memset()
#include <map>
#include <filesystem>  // C++17
#include <stdexcept>


#include "TROOT.h"
#include "TLatex.h"
#include <TMinuit.h>
#include <TStyle.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2D.h>
#include <THStack.h>
#include <TMath.h>
#include <TTree.h>
#include <TRandom3.h>
#include <vector>
#include <TGraph2D.h>
#include <TPolyLine3D.h>
#include <TLine.h>
#include <TTimeStamp.h>
#include <TComplex.h>
#include <TVirtualFFT.h>
#include <TMatrixD.h>
#include "TSystem.h"
#include <TEventList.h>
#include <TPolyMarker.h>
#include <TSpectrum.h>
#include <TGaxis.h>
#include <TLegend.h>
#include <TPaveStats.h>
#include <TMarker.h>
#include <TKey.h>


// #include "DUNEStyle.h"
// // use with:
// //
// // #define DUNESTYLE_ENABLE_AUTOMATICALLY 0
// // dunestyle::SetDuneStyle();

using namespace std;

#include "denoiser.C"
#include "ADC_DATA.h"
#include "readingCodes.C"
#include "wiener_filter.C"
#include "analyzer.C"
#include "analyzer_plots.C"
#include "calibrationCodes.C"
#include "scripts.C"


#endif


