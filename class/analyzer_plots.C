// ________________________________________ //
// Author: Henrique Souza
// Filename: analyzer_plots.C
// Created: 2023-02-08
// ________________________________________ //
//

#include "MYCODES.h"

/**
 *
 * This class holds the methods for ploting for analyzer */
void ANALYZER::sample_plot(Int_t myevent, string opt, Double_t filter, Double_t factor, Int_t mafilter){
  if (opt == "") opt = plot_opt;
  bool state = getWaveformHard(myevent,factor);
  if (!state) return;
  applyMovingAverage(mafilter);
  applyDenoise(filter);

  drawGraph(opt,n_points,&time[0],&ch[kch]->wvf[0]);
}

void ANALYZER::showWaveform(Int_t maxevent, Double_t filter, Int_t dt){

  if (maxevent==0) {
    maxevent = nentries;
  }
  TCanvas *c1 = new TCanvas("c1", "c1",1920,0,1920/2,1080/2);


  TLatex *tx = new TLatex();
  for(Int_t i = 0; i < maxevent; i++){
    if (gSystem->ProcessEvents()) // avoid Canvas freezing and exit with Ctrl+c
      break;
    else{
      if (gROOT->IsBatch() && i!=0) {
        c1->Print("show_wvfs.gif+");
      }
    }
    sample_plot(i,"AL",filter);
    printf("\rEvent %d", i);
    fflush(stdout);
    c1->Modified();
    c1->Update();
    tx->DrawLatex(0.8*c1->GetPad(0)->GetUxmax(), 0.8*c1->GetPad(0)->GetUymax(),Form("Event: %i", i));
    c1->Modified();
    c1->Update();
    if(dt!=0) this_thread::sleep_for(chrono::milliseconds(dt));
  }
  if (gROOT->IsBatch()) {
    c1->Print("show_wvfs.gif++"); //produces infinite loop
  }
}

void ANALYZER::persistence_plot(Int_t nbins, Double_t ymin, Double_t ymax, Double_t filter, string cut, Double_t factor){

  create_persistence(nbins, ymin, ymax, filter, cut, factor);

  hpers->SetStats(kFALSE);
  hpers->Draw("colz");
  hpers->GetXaxis()->SetTitle("Time (ns)");
  hpers->GetYaxis()->SetTitle("Amplitude (ADC Channels)");

  gPad->Update();
  //TPaveStats *ps = (TPaveStats*)hpers->FindObject("stats");
  //ps->SetOptStat(10);

}

void ANALYZER::add_persistence_plot(TH2D *_htemp, Double_t filter, string cut, Double_t factor){
  if(!_htemp) _htemp = hpers;
  getSelection(cut);
  Int_t nev = lev->GetN();
  Int_t iev = 0;
  for(Int_t i = 0; i < nev; i++){
    iev = lev->GetEntry(i);
    getWaveform(iev,kch);
    applyDenoise(filter);
    for (int j = 0; j < n_points; j++) {
      _htemp->Fill(j*dtime,ch[kch]->wvf[j]*factor);
    }
  }
}

TGraph ANALYZER::drawGraph(string opt, Int_t n, Double_t* x, Double_t* y){
  if (opt == "") opt = plot_opt;
  if (x == nullptr) x = time;
  if (y == nullptr) y = ch[kch]->wvf;
  if (n == 0) n = n_points;
  gwvf = new TGraph(n,x,y);
  gwvf->Draw(opt.c_str());
  gwvf->GetXaxis()->SetTitle(xlabel.c_str());
  gwvf->GetYaxis()->SetTitle(ylabel.c_str());

  if(ymax!=0 && ymin!=0){
    gwvf->GetYaxis()->SetRangeUser(ymin,ymax);
  }
  gwvf->GetXaxis()->SetRangeUser(xmin,xmax);
  gwvf->SetLineColor(mylinecolor);
  gwvf->SetMarkerColor(mymarkercolor);
  gwvf->SetEditable(kFALSE);

  TLatex *tx = new TLatex();
  TPad *currentPad = (TPad*)gROOT->GetSelectedPad();
  if(plotShowEvent){
    tx->DrawLatex(0.8*currentPad->GetUxmax(), 0.8*currentPad->GetUymax(),Form("Event: %d", currentEvent));
  }
  return *gwvf;
}

void ANALYZER::showFFT(Int_t naverage, Int_t maxevent, Int_t dt, bool inDecibel){

  if (maxevent==0) {
    maxevent = nentries;
  }
  TH1D *hsf = (TH1D*)w->hfft->Clone("hsf");
  vector<TH1D *> hfft(naverage);
  Int_t k = 0;
  TCanvas *c1 = new TCanvas("c1");
  Double_t maxVal = pow(2,14);
  TH1D *hshowf = nullptr;
  TLatex tx;
  for(Int_t i = 0; i < maxevent; i++){
    getWaveform(i,kch);
    getFFT();
    if (i < naverage) {
      hfft[k] = (TH1D*)w->hfft->Clone(Form("hsf%d",k));
      for (Int_t j = 0; j < n_points/2; j++) hsf->AddBinContent(j+1,hfft[k]->GetBinContent(j+1));
      if (naverage == 1) hsf->Draw();
    }
    else{
      if (k < naverage) {
        for (Int_t j = 0; j < n_points/2; j++) {
          hsf->AddBinContent(j+1,-hfft[k]->GetBinContent(j+1));
          hfft[k]->SetBinContent(j+1,w->hfft->GetBinContent(j+1));
          hsf->AddBinContent(j+1,hfft[k]->GetBinContent(j+1));
        }
        printf("\rEvent %d", i);
        fflush(stdout);
        c1->Modified();
        c1->Update();
        if (gSystem->ProcessEvents())
          break;
        else{
          if (gROOT->IsBatch()) {
            c1->Print("show_ffts.gif+");
          }
        }
        if(dt!=0) this_thread::sleep_for(chrono::milliseconds(dt));
      }
      else{
        if (k == naverage){
          hshowf = (TH1D*)hsf->Clone("hshowf");
          if (inDecibel) w->convertDecibel(hshowf);
          hshowf->Draw();
          tx.DrawLatex(0.8*n_points, 0.8*hshowf->GetMaximum(),Form("Event: %i", i));
        }
        k = 0;
      }
    }

    k++;
    if(naverage == 1) k = 0;
  }

}

void ANALYZER::debugSPE(Int_t event, Int_t moving_average, Int_t n_moving, Double_t xmin, Double_t xmax, vector<Double_t> signal_range, vector<Double_t> not_used,  Double_t filter, Double_t factor, Double_t *SNRs){

  if(!SNRs){
    SNRs = new Double_t[2];
  }
  getWaveform(event,kch);
  applyDenoise(filter);
  TGraph *g = new TGraph(drawGraph());

  g->Draw("AL");
  g->SetLineColor(kGreen);
  g->GetXaxis()->SetRangeUser(xmin,xmax);

  for(Int_t i = 0; i < n_moving; i++){
    applyMovingAverage(moving_average,xmin,xmax);
  }
  // compute signal to noise ratio for amplitude
  // i am avoiding the fist and last 3 points due to residual noise there
  SNRs[0] = computeSNR_simple(xmin+3*dtime, xmax-3*dtime, signal_range);
  drawGraph("SAME");
  Double_t maxav = getMaximum(xmin+8, xmax-8);

  getWaveform(event,kch);
  applyDenoise(filter);
  differenciate(factor);
  for(Int_t i = 0; i < n_moving; i++){
    applyMovingAverage(moving_average,xmin,xmax);
  }
  Double_t maxdif = getMaximum(xmin+8, xmax-8);
  // scaleWvf(maxav/maxdif); // not a good option, it will chage for each peak
  // scaleWvf(1e3);
  // compute signal to noise ratio for amplitude
  // i am avoiding the fist and last 3 points due to residual noise there
  SNRs[1] = computeSNR_simple(xmin+3*dtime, xmax-3*dtime, signal_range);
  drawGraph("SAME");
  gwvf->SetLineColor(kRed);
  // Double_t ymin = -maxav*3;
  // Double_t ymax = maxav*3;
  // g->GetYaxis()->SetRangeUser(ymin,ymax);
}



/**
 * This function will find the best SNR changing the number and range of moving averages applied
 * You need to find an `event` by eye first. Lets say, `event` 3 has 1 p.e. around 6200 ns.
 * Set the xmin and xmax for the computation, ex: 5000 to 10000 ns
 * Set the range of the signal in ns, ex: 6200 to 6700 ns
 * call the function as minimizeParamsSPE(3, 5000, 10000, {6200, 6700})
 **/
void ANALYZER::minimizeParamsSPE(Int_t event, Double_t xmin, Double_t xmax, vector<Double_t> signal_range, vector<Double_t> rangeInter, Double_t filter, Double_t factor){



  // Creating parameters to be tested
  // The moving average window will between 50% to 100% of the signal window
  Double_t signal_min = signal_range[0];
  Double_t signal_max = signal_range[1];
  if(rangeInter[0]==0){
    rangeInter[1] = (signal_range[1] - signal_range[0])/dtime;
    rangeInter[0] = 0.2*rangeInter[1];
  }
  else
    cout << rangeInter[1] << endl;
  Int_t ma_max = rangeInter[1];
  Int_t ma_min = rangeInter[0];
  Int_t nwindows = (ma_max-ma_min)+1;
  vector<Double_t> ma_interactions(nwindows); // compute the number of different m.a. window


  Int_t ma_repeat = 3; // times to repeat the moving average
  Int_t total_tries = ma_repeat * nwindows; // total combinations possible
  vector<vector<Double_t>> parameters(total_tries, std::vector<Double_t>());
  Int_t total_counter = 0;
  // create vector that stores `n. repetitions`, `window`, `SNR` found
  for (Int_t i = 1; i <= ma_repeat; i++) {
    for(Int_t j = 0; j < nwindows; j++){
      ma_interactions[j] = ma_min + j;
      parameters[total_counter] = {(Double_t)i, ma_interactions[j], 0};
      total_counter += 1;
    }
  }

  TCanvas *c1 = new TCanvas("c1", "c1",1920,0,1920,1080);
  vector<Double_t> SNRs = {0,0,0};
  vector<vector<Double_t>> resSNR0(ma_repeat); // stores the values of snr for different iterations
  vector<vector<Double_t>> resSNR1(ma_repeat); // stores the values of snr for different iterations
                                              // We know the length of internal vectors is nwindows, but we push_back

  Double_t max0 = 0;
  Double_t max1 = 0;
  Int_t maxidx0 = 0;
  Int_t maxidx1 = 0;

  getWaveform(event,kch);
  applyDenoise(16);
  SNRs[2] = computeSNR_simple(xmin+3*dtime, xmax-3*dtime, signal_range);

  cout << "idx\tninter\tav. window\tSNR (differential)\tSNR (M.A)\t SNR (normal)" << endl;
  for(Int_t i = 0; i < total_tries; i++){
    debugSPE(event, parameters[i][1], parameters[i][0], xmin, xmax, signal_range, rangeInter, filter, factor, &SNRs[0]);

    c1->Modified();
    c1->Update();
    if (gSystem->ProcessEvents()) // avoid Canvas freezing and exit with Ctrl+c
      break;
    // this_thread::sleep_for(chrono::milliseconds(dt));
    parameters[i][2] = SNRs[1];
    resSNR0[parameters[i][0]-1].push_back(SNRs[0]);
    resSNR1[parameters[i][0]-1].push_back(SNRs[1]);
    cout << i << "\t";
    for (auto j: parameters[i]) {cout << j << "\t";} cout << SNRs[0] << "\t" << SNRs[2]<< endl;

    if (SNRs[0] >= max0) {max0 = SNRs[0]; maxidx0 = i;}
    if (SNRs[1] >= max1) {max1 = SNRs[1]; maxidx1 = i;}

  }

  cout << "\n\n";
  cout << "Maximum SNR (ma)   = " << max0 << " (" << parameters[maxidx0][2] << ") window = " << parameters[maxidx0][1] << " " << " ninter = " << parameters[maxidx0][0] << endl;
  cout << "Maximum SNR (diff) = " << max1 << " (" << parameters[maxidx1][2] << ") window = " << parameters[maxidx1][1] << " " << " ninter = " << parameters[maxidx1][0] << endl;

  vector<TGraph *> gsnr0(3);
  vector<TGraph *> gsnr1(3);
  TMultiGraph *gm0 = new TMultiGraph();
  TMultiGraph *gm1 = new TMultiGraph();

  for (Int_t i = 0; i < 3; i++) {
    gsnr0[i] = new TGraph(nwindows, &ma_interactions[0], &resSNR0[i][0]);
    gsnr0[i]->SetTitle(Form("%d interactions",i+1));
    gsnr0[i]->SetMarkerStyle(21);
    gsnr0[i]->SetMarkerSize(0.7);
    gm0->Add(gsnr0[i]);

    gsnr1[i] = new TGraph(nwindows, &ma_interactions[0], &resSNR1[i][0]);
    gsnr1[i]->SetTitle(Form("%d interactions",i+1));
    gsnr1[i]->GetXaxis()->SetTitle("Moving Average Window (Ticks)");
    gsnr1[i]->GetYaxis()->SetTitle("SNR");
    gsnr1[i]->SetMarkerStyle(21);
    gsnr1[i]->SetMarkerSize(0.7);
    gm1->Add(gsnr1[i]);
  }

  TCanvas *c2 = new TCanvas("c2", "Diff",1920,0,1920,1080);
  gm1->Draw("ap pmc");
  gm1->GetXaxis()->SetTitle("Moving Average Window (Ticks)");
  gm1->GetYaxis()->SetTitle("SNR");
  c2->BuildLegend();


  TCanvas *c3 = new TCanvas("c3", "MA",1920,0,1920,1080);

  gm0->Draw("ap pmc");
  gm0->GetXaxis()->SetTitle("Moving Average Window (Ticks)");
  gm0->GetYaxis()->SetTitle("SNR");
  c3->BuildLegend();

}


void ANALYZER::drawZeroCrossingLines(vector<Int_t> &peaksCross, vector<Int_t> &peaksRise, TCanvas *c, Double_t ymin, Double_t ymax){
  if(!c){
    cout << "We need a Canvas :)" << endl;
    return;
  }
  TPad *p = (TPad*)c->GetPad(0);

  if(ymin == 0 && ymax == 0)
  {
    ymin = p->GetUymin();
    ymax = p->GetUymax();
  }

  Int_t nlines = peaksCross.size();
  vector<TLine *> lns(nlines);
  for(Int_t i = 0; i < nlines; i++){
    Int_t lnx = peaksCross[i]*dtime;
    lns[i] = new TLine(lnx, ymin, lnx, ymax);
    lns[i]->Draw("SAME");
  }
  nlines = peaksRise.size();
  lns.resize(nlines);
  for(Int_t i = 0; i < nlines; i++){
    Int_t lnx = peaksRise[i]*dtime;
    lns[i] = new TLine(lnx, ymin, lnx, ymax);
    lns[i]->Draw("SAME");
  }
}



void ANALYZER::histoTimeTrigger(Int_t nstart, Int_t nfinish, TH1D *_htemp)
{
  if(!_htemp){
    _htemp = new TH1D("","",500,0,0);
  }
  if(nfinish==0) nfinish = nentries;
  Double_t ref = 0;
  for(Int_t i = nstart; i < nfinish; i++){
    getWaveform(i,kch);
    if(i == nstart){
      ref = ch[kch]->time;
    }
    else{
      _htemp->Fill(ch[kch]->time - ref);
      ref = ch[kch]->time;
    }
  }
  _htemp->GetYaxis()->SetTitle("# of events");
  _htemp->GetXaxis()->SetTitle("Time between trigger records (s)");
  _htemp->Draw();
}

void ANALYZER::graphTimeTrigger(Int_t nstart, Int_t nfinish, TGraph *_gtemp)
{
  if(nfinish==0) nfinish = nentries;
  Double_t ref = 0;

  if(!_gtemp){
    _gtemp = new TGraph();
  }
  for(Int_t i = nstart; i < nfinish; i++){
    getWaveform(i,kch);
    if(i == nstart){
      ref = ch[kch]->time;
    }
    else{
      _gtemp->AddPoint(ch[kch]->event, ch[kch]->time - ref);
      ref = ch[kch]->time;
    }
  }
  _gtemp->Draw("ALP");
}


void ANALYZER::check_filtering(vector<Double_t> filter_max_and_step, Int_t event, Int_t rebine, Double_t refFreq){

  Double_t max_filter = filter_max_and_step[0];
  Double_t step_filter = filter_max_and_step[1];

  if(max_filter == 0 || step_filter == 0){
    max_filter = 32;
    step_filter = 8;
  }
  Int_t nf = static_cast<Int_t>(max_filter/step_filter) + 1;

  string reference_filter = filter_type;

  vector<TGraph*> gnf(nf);
  vector<TH1D*> hnf_fft(nf);
  TH1D *href = nullptr;
  TMultiGraph *gm = new TMultiGraph("gm", "gm");
  THStack *hs = new THStack("hs","hs");
  for(Int_t i = 0; i < nf; i++){
    getWaveform(event);
    Double_t filter = 0 + step_filter*i;
    filter_type = reference_filter;
    if(filter_type != "default") setFreqFilter(filter, filter_type);
    if(i!=0) applyDenoise(filter);
    gnf[i] = new TGraph(drawGraph());
    gnf[i]->SetTitle(Form("Filter = %.2f", filter));
    gm->Add(gnf[i]);
    getFFT();
    hnf_fft[i] = (TH1D*)h->Clone(Form("h_nf_%.2f", filter));
    if(i == 0){
      href = (TH1D*)h->Clone("href");
    }
    hnf_fft[i]->Divide(hnf_fft[i], href);
    for(Int_t j = 0; j < hnf_fft[i]->GetNbinsX(); j++){
      hnf_fft[i]->SetBinContent(j+1, 20*log10(hnf_fft[i]->GetBinContent(j+1)));
    }
    hnf_fft[i]->SetTitle(Form("Filter = %.2f", filter));
    hnf_fft[i]->Rebin(rebine);
    hs->Add(hnf_fft[i],"hist");
  }

  getWaveform(event);
  setFreqFilter(refFreq,"low");
  applyDenoise(refFreq);
  TGraph *glow = new TGraph(drawGraph());
  glow->SetTitle(Form("Low pass %0.f MHz", refFreq));
  glow->SetLineColor(kRed);
  getFFT();
  TH1D *hlow_fft = (TH1D*)h->Clone(Form("Low pass %0.f MHz", refFreq));
  hlow_fft->Divide(hlow_fft, href);
  for(Int_t j = 0; j < hlow_fft->GetNbinsX(); j++){
    hlow_fft->SetBinContent(j+1, 20*log10(hlow_fft->GetBinContent(j+1)));
  }
  hlow_fft->SetTitle(Form("Low pass %0.f MHz", refFreq));
  hlow_fft->Rebin(rebine);
  hlow_fft->SetLineColor(kRed);

  TCanvas *c1 = new TCanvas("c1", "c1",1920,0,1920,1080);
  gm->Draw("AL plc");
  glow->Draw("SAME L");

  c1->BuildLegend();

  TCanvas *c2 = new TCanvas("c2", "c2",1920,0,1920,1080);
  c2->SetLogx(1);
  hs->Draw("nostack plc");
  hs->GetXaxis()->SetTitle("Frequency (MHz)");
  hs->GetYaxis()->SetTitle("Gain (dB)");
  hlow_fft->Draw("SAME hist");
  c2->BuildLegend();


  filter_type = reference_filter;

}
