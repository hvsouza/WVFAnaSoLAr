// ________________________________________ //
// Author: Henrique Souza
// Filename: wiener_filter.C
// Created: 2021
// ________________________________________ //
#include "MYCODES.h"


class WIENER{
  public:

    string obj_name;
    Int_t npts = 0;
    Double_t frequency = 250; // in units
    Double_t step = 4; // in ns
    Double_t units_freq = 1e6; // 1e6 = MHz
    Double_t units_step = 1e-9; // 1e9 = ns
    Double_t factor = 1;

    Double_t baseline = 0;

    TH1D *hfft = nullptr;
    TH1D *hPSD = nullptr;
    TH1D *hwvf = nullptr;
    Double_t powerSpectrum = 0;
    TComplex *spec = nullptr;
    Double_t *spec_re = nullptr;
    Double_t *spec_im = nullptr;
    Double_t *re_comp = nullptr;
    Double_t *im_comp = nullptr;
  
    TF1 *f_filter = nullptr;
    TF1 *flar = nullptr;
    TF1 *flar_gaus = nullptr;

    Double_t *res = nullptr;
    Int_t maxBin=0;

    string unit_time;
    string unit_freq;


    void add_fft(WIENER _htemp){
      for(Int_t k=0; k<npts/2+1; k++){
        spec[k] += _htemp.spec[k];
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
      }
    }

    void scale(Double_t fscale){
      for(Int_t k=0; k<npts/2+1; k++){
        spec[k] *= fscale;
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
      }

    }

    // h is the s.p.e. response
    // n is the noise
    // s is the expetect signal (delta-like)
    void wienerGenFilter(WIENER h , WIENER n, WIENER s, Int_t method = 1){
      for(Int_t k=0; k<npts/2+1; k++){
        Double_t S2 = 1;                  // expected signal power spectra (for delta function S2 = 1)
        Double_t H2 = h.spec[k].Rho2();   // single photo-electron power spectra (|H(f)|^2)
        Double_t N2 = n.spec[k].Rho2();   // noise power spectra
        TComplex xH = h.spec[k];          // single photo-electron fourier transformation
        if(method!=1){
          S2 = s.spec[k].Rho2();
        }
        
        spec[k] = TComplex::Conjugate(xH)*S2/(H2*S2+N2);
        // cout << spec[k] << " " << TComplex::Conjugate(xH) << " " << S2 << " " << H2 << " " << N2 << endl;
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
      }

      //Now let's make a backward transform:
    
      TVirtualFFT *fft_final = TVirtualFFT::FFT(1, &npts, "C2R M K");
      fft_final->SetPointsComplex(spec_re,spec_im);
      fft_final->Transform();
      TH1 *hfinal = 0;
      //Let's look at the output
      hfinal = TH1::TransformHisto(fft_final,hfinal,"Ref");
      // hfinal->Scale(factor);

      for(Int_t i = 0; i<npts; i++){
        res[i] = hfinal->GetBinContent(i+1);
        hwvf->SetBinContent(i+1,res[i]);
      }

      fft(hwvf);
      delete hfinal;
   
    }

    void applyWienerFilter(WIENER wg){
      convolve(&wg);
      recompute_hist();
      backfft(*this);
      shift_waveform(hwvf, maxBin);
      hwvf->Sumw2(kFALSE);
      set_flar();
    }

    void build_noise(Double_t sigma = 3.5, Int_t filter = 25,Int_t myseed = 1){
      gRandom->SetSeed(myseed);
      vector<Double_t> mynoise(npts);
      vector<Double_t> mynoise_filtered(npts);
      for(Int_t i = 0; i<npts; i++){
        mynoise[i] = gRandom->Gaus(0,sigma);
      }
      DENOISE dn;
      if(filter > 0) dn.TV1D_denoise<Double_t>(&mynoise[0],&mynoise_filtered[0],npts,filter);
      else mynoise_filtered = mynoise;
      for(Int_t i = 0; i<npts; i++){
        hwvf->SetBinContent(i+1,mynoise_filtered[i]);
      }
      fft(hwvf);
    }

    void apply_filter(){

      // the for is performed in k, need to convert back to frequency
      Double_t convert_freq = frequency/npts;

      Double_t filter = 1;
      for(Int_t k=0; k<npts/2+1; k++){
        filter = f_filter->Eval(convert_freq*k);
        spec[k] = spec[k]*filter;
        // cout << spec[k] << " " << filter << endl;
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
        idx_recompt(k);
      }
    }

    void setBandCut(Double_t lowFreq, Double_t highFreq, WIENER *_temp){
      for(Int_t k=0; k<npts/2+1; k++){
        _temp->spec[k] = 1.;
        // cout << spec[k] << endl;
        _temp->spec_re[k] = _temp->spec[k].Re();
        _temp->spec_im[k] = _temp->spec[k].Im();
      }
      _temp->setFilter(lowFreq, "high");
      _temp->apply_filter();
      _temp->setFilter(highFreq, "low");
      _temp->apply_filter();
      for(Int_t k=0; k<npts/2+1; k++){
        _temp->spec[k] = 1. - _temp->spec[k];
        // _temp->hfft->SetBinContent(k+1,_temp->spec[k].Rho());
        // _temp->hPSD->SetBinContent(k+1,_temp->spec[k].Rho2());
      }
      // _temp->hfft->Draw();
    }
    
    void applyBandCut(WIENER *_temp){
      convolve(_temp);
    }
    void setFilter(Double_t cutoff_frequency, string filter_type){

      // cutoff_frequency is the cutoff frequency in MHz (or your unit set)
      if (filter_type == "gaus")
      {
        cutoff_frequency = sqrt(sqrt(2))*cutoff_frequency;
        f_filter = new TF1("filter","TMath::Gaus(x,[0],[1])",0,frequency);	// A gaussian filter
        // the standard deviation must be sqrt(sqrt(2))* the cutoff frequen, so when x = cutoff frequency Vo/Vi = 0.7
        f_filter->SetParameters(0,cutoff_frequency);
      }
      else if (filter_type == "low")
      {
        f_filter = new TF1("filter","1/sqrt(1+x*x/([0]*[0]))",0,frequency);	// A gaussian filter
        f_filter->SetParameter(0,cutoff_frequency);
      }
      else if (filter_type == "high")
      {
        f_filter = new TF1("filter","1/sqrt(1+[0]*[0]/(x*x))",0,frequency);	// A gaussian filter
        f_filter->SetParameter(0,cutoff_frequency);
      }
      

    }

    void convertDecibel(){
      convertDecibel(hfft);
      convertDecibel(hPSD);
    }
    void convertDecibel(TH1D *htemp, Double_t max = 0){
      if(max==0){
        max = htemp->GetMaximum();
      }
      Int_t ntemp = htemp->GetEntries();
      for(Int_t k=0; k<npts/2+1; k++){
        // cout << htemp->GetBinContent(k+1) << endl;
        htemp->SetBinContent(k+1,20*TMath::Log10(htemp->GetBinContent(k+1)/max));
      }
      htemp->SetEntries(ntemp);
    }
    void fillPlane(){
      for(Int_t k=0; k<npts/2+1; k++){
        spec[k] = 1;
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
        idx_recompt(k);
      }
    }
    void convolve(WIENER *_temp){
      for(Int_t k=0; k<npts/2+1; k++){
        spec[k] = spec[k]*_temp->spec[k];
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
        idx_recompt(k);
      }

    }

    void set_flar(){
      flar = new TF1("flar",Form("[0]*exp(-(x-%f)/[1])+[2]*exp(-(x-%f)/[3])",maxBin*step,maxBin*step),0,npts*step);
      flar->SetParameters(100,10,3,1400);
      flar->SetParName(0,"A_{F}");
      flar->SetParName(1,"#tau_{F}");
      flar->SetParName(2,"A_{S}");
      flar->SetParName(3,"#tau_{S}");

      flar_gaus = new TF1("flar_gaus","([0]*exp(-(x-[5])/[1])*exp([4]*[4]/(2*[1]*[1])))*TMath::Erfc((([5]-x)/[4]+[4]/[1])/TMath::Power(2,0.5))/2. + ([2]*exp(-(x-[5])/[3])*exp([4]*[4]/(2*[3]*[3])))*TMath::Erfc((([5]-x)/[4]+[4]/[3])/TMath::Power(2,0.5))/2.",0,npts*step);
      flar_gaus->SetParameters(100,10,3,1400,20,maxBin*step);
      flar_gaus->SetParName(0,"A_{F}");
      flar_gaus->SetParName(1,"#tau_{F}");
      flar_gaus->SetParName(2,"A_{S}");
      flar_gaus->SetParName(3,"#tau_{S}");
      flar_gaus->SetParName(4,"#sigma");
      flar_gaus->SetParName(5,"t_{0}");
    }
    void frequency_deconv(WIENER y, WIENER G, Double_t cutoff_frequency=0, string filter_type = "gaus"){

      setFilter(cutoff_frequency,filter_type);

      // the for is performed in k, need to convert back to frequency
      Double_t convert_freq = frequency/npts;


      Double_t gaus_blur = 1;
      for(Int_t k=0; k<npts/2+1; k++){
        if(cutoff_frequency!=0) gaus_blur = f_filter->Eval(convert_freq*k);
        spec[k] = y.spec[k]*G.spec[k]*gaus_blur;
        // cout << spec[k] << endl;
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
      }

      //Now let's make a backward transform:
      TVirtualFFT *fft_final = TVirtualFFT::FFT(1, &npts, "C2R M K");
      fft_final->SetPointsComplex(spec_re,spec_im);
      fft_final->Transform();
      TH1 *hfinal = 0;
      //Let's look at the output
      hfinal = TH1::TransformHisto(fft_final,hfinal,"Ref");
      // hfinal->Scale(factor);



      shift_waveform(hfinal,y.maxBin);
      Double_t bl = 0;
      Double_t auxbaseline = 0;
      for(Int_t i=0; i<baseline/step; i++){
        bl+=hfinal->GetBinContent(i+1);
        auxbaseline+=1;
      }
      if(baseline == 0) auxbaseline = 1;
      bl=bl/auxbaseline;

   
      for(Int_t i = 0; i<npts; i++){
        res[i] = hfinal->GetBinContent(i+1);
        hwvf->SetBinContent(i+1,res[i]-bl);
      }

      fft(hwvf);

      maxBin = y.maxBin;
      set_flar();

      delete hfinal;
   
    
    }
    void deconvolve(WIENER y, WIENER h, Double_t cutoff_frequency = 50, string filter_type = "gaus"){
      // y is the signal, h is the device response (a.k.a single photo-electron)

      setFilter(cutoff_frequency,filter_type);
      // the for is performed in k, need to convert back to frequency
      Double_t convert_freq = frequency/npts;


      Double_t gaus_blur = 1;
      for(Int_t k=0; k<npts/2+1; k++){
        if(cutoff_frequency!=0) gaus_blur = f_filter->Eval(convert_freq*k);
        if(!(h.spec[k].Re()==0 && h.spec[k].Im()== 0)){
          spec[k] = y.spec[k]*gaus_blur/h.spec[k];
        }
        else{
          spec[k] = 0;
        }
      
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();

        // or you can do something like this...
        // spec_re[k] = f_filter->Eval(k)*(y.re_comp[k]*h.re_comp[k] + y.im_comp[k]*h.im_comp[k])/(pow(h.re_comp[k],2)+pow(h.im_comp[k],2));
        // spec_im[k] = f_filter->Eval(k)*(h.re_comp[k]*y.im_comp[k] - y.re_comp[k]*h.im_comp[k])/(pow(h.re_comp[k],2)+pow(h.im_comp[k],2));

      }
      //Now let's make a backward transform:
      TVirtualFFT *fft_final = TVirtualFFT::FFT(1, &npts, "C2R M K");
      fft_final->SetPointsComplex(spec_re,spec_im);
      fft_final->Transform();
      TH1 *hfinal = 0;
      //Let's look at the output
      hfinal = TH1::TransformHisto(fft_final,hfinal,"Ref");
      // hfinal->Scale(factor);


      for(Int_t i = 0; i<npts; i++){
        res[i] = hfinal->GetBinContent(i+1);
        hwvf->SetBinContent(i+1,res[i]);
      }

      shift_waveform(hwvf,y.maxBin);
      Double_t bl = 0;
      Double_t auxbaseline = 0;
      for(Int_t i=0; i<baseline/step; i++){
        bl+=hwvf->GetBinContent(i+1);
        auxbaseline+=1;
      }
      if(baseline == 0) auxbaseline = 1;
      bl=bl/auxbaseline;


      for(Int_t i = 0; i<npts; i++){
        res[i] = hwvf->GetBinContent(i+1);
        hwvf->SetBinContent(i+1,res[i]-bl);
      }

      fft(hwvf);

      maxBin = y.maxBin;
      set_flar();
      delete hfinal;
   
    }







    template <class T>
    void shift_waveform(T *h, Int_t new_max, Bool_t rawShift = false){
      Int_t old_max = h->GetMaximumBin();
      if(rawShift) old_max = 0;
      Int_t old_ref = old_max - new_max;
      TH1D *htemp = (TH1D*)h->Clone("htemp");
      Double_t temp;
      if(old_ref<0){
        // cout << " case lower" << endl;
        old_ref = npts-(new_max-old_max);
      }
      for(Int_t i = 1; i<npts-(old_ref); i++){
        temp = htemp->GetBinContent(old_ref+i);
        h->SetBinContent(i,temp);
      }
      Int_t aux = 1;
      for(Int_t i = npts-(old_ref); i<=npts; i++){
        temp = htemp->GetBinContent(aux);
        h->SetBinContent(i,temp);
        aux++;
      }
      delete htemp;
    }



    void idx_recompt(Int_t k){
        hfft->SetBinContent(k+1,spec[k].Rho()); // correction for plot
        hPSD->SetBinContent(k+1,spec[k].Rho2()/2.);
    }
    void recompute_hist(){
      for(Int_t k = 0; k<npts/2+1; k++){
        idx_recompt(k);
      }
    }



  
    void fft(TH1D *hsignal){

      if(maxBin==0) maxBin = hsignal->GetMaximumBin(); //get maximum to realign waveforms later, only one time helps with several waveforms (not all max are the same)
 
      TH1 *hm = 0;
      TVirtualFFT::SetTransform(0);
      hm = hsignal->FFT(hm, "MAG R2C measure");
      //NOTE: for "real" frequencies you have to divide the x-axes range with the range of your function
      //(in this case 4*Pi); y-axes has to be rescaled by a factor of 1/SQRT(n) to be right: this is not done automatically!

      //Look at the DC component and the Nyquist harmonic:
      //That's the way to get the current transform object:
      TVirtualFFT *fft = TVirtualFFT::GetCurrentTransform();
    
      //Use the following method to get the full output:
      fft->GetPointsComplex(re_comp,im_comp);
   
      for(Int_t k = 0; k<npts/2+1; k++){
        spec[k] = TComplex(re_comp[k],im_comp[k])*factor;
        idx_recompt(k);
        // the same as:
        // hfft->SetBinContent(k+1,hm->GetBinContent(k+1)*2*factor);


        // cout << k << " " << re_comp[k] << " " << im_comp[k] << " " << spec[k].Rho2() << endl;
      }
      powerSpectrum = hPSD->Integral();
     
      delete hm;
      delete fft;

    }



    void backfft(WIENER y){

      for(Int_t k=0; k<npts/2+1; k++){
        spec[k] = y.spec[k];
      
        spec_re[k] = spec[k].Re();
        spec_im[k] = spec[k].Im();
        // cout << spec_re[k] << " " << spec_im[k] << endl;

      }
      //Now let's make a backward transform:
      TVirtualFFT *fft_final = TVirtualFFT::FFT(1, &npts, "C2R M K");
      fft_final->SetPointsComplex(spec_re,spec_im);
      fft_final->Transform();
      TH1 *hfinal = 0;
      //Let's look at the output
      hfinal = TH1::TransformHisto(fft_final,hfinal,"Ref");
      // hfinal->Scale(factor);


      hwvf->Reset();
      for(Int_t i = 0; i<npts; i++){
        res[i] = hfinal->GetBinContent(i+1);
        hwvf->SetBinContent(i+1,res[i]/2); // need to correct a sqrt(2) factor that comes around :/
      }

      fft(hwvf);
   
      // flar = new TF1("flar",Form("[0]*exp(-(x-%f)/[1])+[2]*exp(-(x-%f)/[3])",y.maxBin*step,y.maxBin*step),0,npts*step);
      // flar->SetParameters(0.3,10,0.3,1000);
   

      delete hfinal;
      delete fft_final;
   
    }


  
    // the histogram stuff works... deconvolution tested only on LED pulses. Need to check for LAr signal.
    // I need to with LAr signals...
    void rescale_histogram(Int_t newpts){
      npts = newpts;
      TComplex *tmpspec = new TComplex[npts];
      tmpspec = spec;

      Double_t oldfactor = factor;
      startup();
      TH1D *htemp = (TH1D*)hfft->Clone("htemp");
      hfft->SetBins(newpts/2,0,frequency/2);

      Double_t centerHisto = htemp->GetBinCenter(1);
      Double_t widthHisto = htemp->GetBinWidth(1);
      Double_t zeroHisto = centerHisto-widthHisto/2;
      Int_t auxHisto = 1;
      for(Int_t i = 1; i<=npts; i++){
        Double_t refNewHisto = hfft->GetBinCenter(i);
        // cout << refNewHisto << " " << zeroHisto << " " << widthHisto << endl;
        if(refNewHisto <= (zeroHisto+widthHisto)){
          spec[i-1]=tmpspec[auxHisto-1]*factor/oldfactor;
          hfft->SetBinContent(i,htemp->GetBinContent(auxHisto)*factor/oldfactor);
        }
        else{
          auxHisto++;
          zeroHisto = zeroHisto+widthHisto;
          spec[i-1]=tmpspec[auxHisto-1]*factor/oldfactor;
          hfft->SetBinContent(i,htemp->GetBinContent(auxHisto)*factor/oldfactor);
        }
      }
      delete htemp;
    }


    void setFromHist(TH1D *hin){
      npts = hin->GetNbinsX();
      step = 1.;
      Double_t min = hin->GetBinCenter(1) - hin->GetBinWidth(1)/2.;
      Double_t max = hin->GetBinCenter(npts) - hin->GetBinWidth(1)/2.;
      frequency = 1./step;

      units_step = 1;
      units_freq = 1;
      delete hfft;
      delete hPSD;
      delete hwvf;
      hfft = new TH1D(Form("fft_%s",obj_name.c_str()),Form("FFT %s",obj_name.c_str()),npts/2,0,frequency/2);
      hPSD = new TH1D(Form("PSD_%s",obj_name.c_str()),Form("Power Spectral Density %s",obj_name.c_str()),npts/2,0,frequency/2);
      hwvf = new TH1D(Form("wvf_%s",obj_name.c_str()),Form("wvf_%s",obj_name.c_str()),npts,min,max);

      for (Int_t i = 0; i < npts; i++) {
        hwvf->SetBinContent(i+1,hin->GetBinContent(i+1));
      }

    }

  

    WIENER(string myname) : obj_name{myname} {startup();};
    WIENER(string myname, Int_t mynpts) : obj_name{myname},
                                          npts{mynpts}  {startup();}
    WIENER(string myname, Double_t mystep, Double_t myfreq, Double_t myunits_step, Double_t myunits_freq, Int_t mynpts ) : obj_name{myname},
                                                                                                                           frequency{myfreq},
                                                                                                                           step{mystep},
                                                                                                                           units_step{myunits_step},
                                                                                                                           units_freq{myunits_freq},
                                                                                                                           npts{mynpts}
                                                                                                                         
    {startup();}
  

    void startup(){

      hfft = new TH1D(Form("fft_%s",obj_name.c_str()),Form("FFT %s",obj_name.c_str()),npts/2,0,frequency/2);
      hPSD = new TH1D(Form("PSD_%s",obj_name.c_str()),Form("Power Spectral Density %s",obj_name.c_str()),npts/2,0,frequency/2);
      hwvf = new TH1D(Form("wvf_%s",obj_name.c_str()),Form("wvf_%s",obj_name.c_str()),npts,0,npts*step);



      spec = new TComplex[npts];
      spec_re = new Double_t[npts];
      spec_im = new Double_t[npts];
      re_comp = new Double_t[npts];
      im_comp = new Double_t[npts];
      res = new Double_t[npts];
      factor = 1./(npts/2.); // tested and working fine to compare between different time windows
      // on the histograms, a factor 2 is added to match the usual normalization (root does not like this, so I keep complex numbers without)

    
      if(units_step == 1e-9)
        unit_time = "ns";
      else if(units_step == 1e-6)
        unit_time = "#mus";
      else if(units_step == 1e-3)
        unit_time = "ms";
      else if(units_step == 1)
        unit_time = "s";
      else{
        unit_time = "A.U.";
      }

      if(units_freq == 1e9)
        unit_freq = "GHz";
      else if(units_freq == 1e6)
        unit_freq = "MHz";
      else if(units_freq == 1e3)
        unit_freq = "kHz";
      else if(units_freq == 1)
        unit_freq = "Hz";
      else{
        unit_freq = "A.U.";
      }


      hfft->GetXaxis()->SetTitle(Form("Frequency (%s)",unit_freq.c_str()));
      hPSD->GetXaxis()->SetTitle(Form("Frequency (%s)",unit_freq.c_str()));
      hfft->GetYaxis()->SetTitle("ADC");
      hPSD->GetYaxis()->SetTitle("PSD (ADC^{2} Hz^{-1})");

      hwvf->GetXaxis()->SetTitle(Form("Time (%s)",unit_time.c_str()));
      hwvf->GetYaxis()->SetTitle("Amplitude (A.U.)");
      hPSD->SetLineWidth(2);
      hfft->SetLineWidth(2);
    }
  

    // based https://dspcookbook.github.io/optimal-filtering/wiener-filter-2/#3-solution
    // creates filter from source s and reference d (such as noise)
    vector<Double_t> create_filter_from_reference(vector<Double_t> signal, vector<Double_t> reference,Int_t filter_size){
      Int_t n = reference.size();
      vector<Double_t> w(filter_size);
      vector<Double_t> rxx(filter_size);
      vector<Double_t> rxd(filter_size);
      TH1D *hnoise = new TH1D("hnoise","hnoise",500,0,0);
      for(Int_t i = 0; i<n; i++){
        hnoise->Fill(signal[i]);
      }
      Double_t rvv = hnoise->GetStdDev();
      rvv = rvv*rvv;
      // evaluate rxx (and rxd ?)
      for(Int_t i = 0; i<filter_size; i++){

        for(Int_t j = 0; j<n-i; j++){
          rxx[i] += signal[j]*signal[j+i];
          rxd[i] += signal[j]*reference[j+i];
        }
        rxx[i]=rxx[i]/(n-i);
        rxd[i]=rxd[i]/(n-i);
        // if(i==0)rxd[i]=rxx[i]-rvv;
        // else rxd[i] = rxx[i];
      }

      TMatrixD Rxx(filter_size,filter_size);
      for(Int_t i=0; i<filter_size; i++){
        for(Int_t j = 0; j<filter_size; j++){
          Rxx[i][j] = rxx[abs(i-j)];
        }
      }
      cout << "Inverting Rxx matrix... " << endl; Rxx.Invert();
      vector<Double_t> x(filter_size);
      for(Int_t i = 0; i<filter_size; i++){
        x[i] = i;
        for(Int_t j = 0; j<filter_size; j++){
          w[i] += Rxx[i][j]*rxd[j];
        }
        // cout << w[i] << endl;
      }
    
      // TCanvas *c2 = new TCanvas();
      // TGraph *gtest = new TGraph(filter_size,&x[0],&rxd[0]);
      // gtest->Draw("ALP");

      // hnoise->Draw();
    
      return w;
    }

    vector<Double_t> filter(vector<Double_t> signal, vector<Double_t> w){
      Int_t n = signal.size();
      Int_t m = w.size();
      vector<Double_t> res(n);
      for(Int_t i = m-1; i<n; i++){
        for(Int_t j = 0; j<m; j++){
          res[i] += signal[i-j]*w[j];
        }
        // cout << res[i] << endl;
      }
      return res;
    }



  
};
