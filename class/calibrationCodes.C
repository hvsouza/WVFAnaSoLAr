// ________________________________________ //
// Author: Henrique Souza
// Filename: calibrationCodes.C
// Created: 2021
// ________________________________________ //
#include "WVFAnaSoLAr.h"


class  MyFunctionObject{
  public:

    Int_t n_peaks;
    // use constructor to customize your function object
    Double_t operator()(Double_t *x, Double_t *par) {
      Double_t f;
      Double_t xx = x[0];
      f  = abs(par[0])*exp(-0.5*TMath::Power((xx-par[1])/par[2],2)); // first argument
      f = f+abs(par[3])*exp(-0.5*TMath::Power((xx-par[4])/par[5],2));
      f = f+abs(par[6])*exp(-0.5*TMath::Power((xx-par[7])/(TMath::Power((2),0.5)*par[5]),2));
      for(Int_t i = 1; i<n_peaks; i++){
        f = f+ abs(par[i+7])*exp(-0.5*TMath::Power((xx-(par[4]+(i+1)*(par[7]-par[4])))/(TMath::Power((i+2),0.5)*par[5]),2));
      }
      return f;
    }
};


class  MyFunctionFree{
  public:

    Int_t n_peaks;
    // use constructor to customize your function object
    Double_t operator()(Double_t *x, Double_t *par) {
      Double_t f;
      Double_t xx = x[0];
      f  = abs(par[0])*exp(-0.5*TMath::Power((xx-par[1])/par[2],2)); // first argument
      f = f+abs(par[3])*exp(-0.5*TMath::Power((xx-par[4])/par[5],2));
      f = f+abs(par[6])*exp(-0.5*TMath::Power((xx-par[7])/par[8],2));
      for(Int_t i = 1, j = 1; i<n_peaks; i++){
        f = f+ abs(par[j+8])*exp(-0.5*TMath::Power((xx-(par[4]+(i+1)*(par[7]-par[4])))/par[j+8+1],2));
        j+=2;
      }
      return f;
    }
};



































class Calibration

{
    
    
    
  public:
    


    string myname = "c";
    Double_t dtime = 4; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    Int_t nbits = 14; // DIGITIZER bits
    Int_t linhasEvento = 9000;
  
    Int_t channel = 1;
    // _______________ Parameters for giveMeSphe_MeanPulseVersion _______________/
    Double_t time_sample = 6000; // in ns
  
    Double_t timeLow = 3650; // integration limit
    Double_t timeHigh = 4200;
  
    Int_t fn = 0;
  
  
    // _______________ Parameters for fit_sphe_wave _______________/
  
    string rootFile = "";
  
    Int_t n_peaks = 7;
    Double_t peak0 = 0.1;
    Double_t mean0 = -900;
    Double_t sigma0 = 500;
  
    Double_t peak1 = 0.004;
    Double_t mean1 = 5000;
    Double_t sigma1 = 400;
  
    Double_t startpoint = 0.001;
    Double_t poisson_ratio = 3.;
    TF1 *faux = nullptr;
  
    Double_t xmin = -10000;
    Double_t xmax = 40000;

    Double_t deltaplus=1;
    Double_t deltaminus=0;

    Double_t stdVar = 0.5;
  
    Int_t rebin = 4;
  
    Bool_t fixZero = false;
  
    Bool_t make_free_stddevs = true;

    Double_t sphe_charge = 0; // wave0
    Double_t sphe_charge2 = 0; // wave0

  
  
    // This was add here to try fitting dark noise data with already found values
    Bool_t darknoise = false;
  
    Bool_t is_poisson_test = false; // if running tests of poisson statistics

    TH1D *htemp = nullptr;
    TH1D *hcharge = nullptr;


    Bool_t drawDebugLines = false;
    Bool_t show_all_parameters = false;
    Bool_t quitemode = false;

    TTree *thead = nullptr;
    Double_t normfactor = 1;

    Int_t fit_status = -1;
    Double_t snr = 0;
    Double_t chi2 = 0;
    Double_t ndf = 0;
    Int_t lastOneAttempts = 2;

    bool save_plot = false;
    string nameplotpng = "";
    Double_t delta1; // dont change... automatic
    Double_t delta2;

  
    // ____________________________________________________________________________________________________ //
    void fit_sphe_wave(string name, bool optimize = true){

      if (optimize == true) searchParameters(name,2,false);
      makeSphe(name);
    }
    Double_t fact(Double_t n){
      if( n == 1 || n == 0 ){
        return 1;
      }
      else{
        return n*fact(n-1);
      }
    }
    Double_t poisson(Double_t lambda, Double_t k){
      Double_t factk = fact(k);
      Double_t nume = exp(-lambda)*TMath::Power(lambda,k);
      return nume/factk;
    }

    void searchParameters(string histogram, Double_t sigmaSearch = 2, bool debug = false){

      TFile *f1 = nullptr;
      TH1D *h = nullptr;
      if (rootFile != "") {
        f1 = new TFile(rootFile.c_str(),"READ");
        h = (TH1D*)f1->Get(histogram.c_str());
        if(!h) h = (TH1D*)f1->Get("analyzed");
        for (auto fkey: *f1->GetListOfKeys()){
          if (string(fkey->GetName()) == "head")
          {
            thead = (TTree*)f1->Get("head");
            thead->SetBranchAddress("normfactor", &normfactor);
            thead->GetEntry(0);
            thead->ResetBranchAddresses();
            break;
          }
        }
      }
      else{
        h = (TH1D*)htemp->Clone("");
      }


      h->Rebin(rebin);
      Double_t scale = 1/(h->Integral());
      // h->Scale(scale);
      Int_t nbins = h->GetNbinsX();
      Double_t *source = new Double_t[nbins];
      Double_t *destVector = new Double_t[nbins];

      Double_t sigma;
      Int_t fNPeaks = 0;
      Int_t nfound = 0;
      Int_t bin;
      // const Int_t nbins = 1024;
      xmin = h->GetXaxis()->GetXmin();
      xmax = h->GetXaxis()->GetXmax();
      Double_t a;
      h->SetTitle("High resolution peak searching, number of iterations = 3");
      h->GetXaxis()->SetRange(1,nbins);
      TH1D *d = new TH1D("d","",nbins,xmin,xmax);

      for (Int_t i = 0; i < nbins; i++) source[i]=h->GetBinContent(i + 1);


      TSpectrum *s = new TSpectrum();

      Double_t sigmastep = 0.2;
      while(nfound < 3 && sigmaSearch>=1){
        nfound = s->SearchHighRes(source, destVector, nbins, sigmaSearch, 2, kFALSE, 5, kTRUE, 3);
        sigmaSearch-=sigmastep;
      }
      if(nfound < 3){
        cout << "Could not optimize parameters, only " << nfound << " peaks found" << endl;
        debug = true;
      }
      Double_t *xpeaks_t = s->GetPositionX();
      vector<Double_t> xpeaks(nfound);
      Int_t pos0 = 0; // this in case there is some peak before zero
      vector<Double_t> fPositionX(nfound-pos0);
      vector<Double_t> fPositionY(nfound-pos0);
      for (Int_t i = 0; i < nfound; i++) xpeaks[i] = xpeaks_t[i];
      std::sort(xpeaks.begin(),xpeaks.end());
      for (Int_t i = 0; i < nfound-pos0; i++) {
        a=xpeaks[i+pos0];
        bin = 1 + Int_t(a + 0.5);
        fPositionX[i] = h->GetBinCenter(bin);
        fPositionY[i] = h->GetBinContent(bin);
      }


      Double_t upplim_gaus_base = fPositionX[0]+(3./4)*(fPositionX[1]-fPositionX[0])/2;
      Double_t lowlim_gaus_base = fPositionX[0]-(fPositionX[1]-fPositionX[0])*2;
      faux = new TF1("faux","gaus(0)",lowlim_gaus_base,upplim_gaus_base);
      faux->SetParameters(fPositionY[0],fPositionX[0],(fPositionX[1]-fPositionX[0])/sqrt(2));
      h->Fit("faux","R0QW");
      h->Fit("faux","R0QW");

      peak0 = fPositionY[0];
      mean0 = fPositionX[0];
      sigma0 = faux->GetParameter(2);

      peak1 = fPositionY[1];
      mean1 = fPositionX[1];
      sigma1 = 0.8*faux->GetParameter(2);

      Int_t bin_second = 0; // corresponding bin of the second peak
      if (nfound-pos0 >= 3) // in case it was found
      {
        startpoint =  fPositionY[2];
        bin_second = 1 + xpeaks[pos0 + 2] + 0.5;
      }
      else{
        Int_t bin_first_peak = 1 + Int_t(xpeaks[1+pos0]);
        Int_t bin_baseline = 1 + Int_t(xpeaks[0+pos0]);
        bin_second = 2*bin_first_peak-bin_baseline; // same as b + (b-a)
        startpoint = h->GetBinContent(bin_second);
        // cout << startpoint << " " << bin_first_peak << " " << bin_baseline << " " << bin_second << " " << h->GetBinCenter(bin_second) << endl;
      }
      Double_t lowestpt = 0;
      int countlow = 0;
      for (Int_t i = bin_second; i < nbins; i++){
        if(h->GetBinContent(i+1)<= 10*h->GetMinimum(0)*scale){
          lowestpt = h->GetBinCenter(i);
          countlow+=1;
        }
        // else{
        //   countlow=0;
        // }
        if (i == nbins-1){
          lowestpt = h->GetBinCenter(i);
        }
        if (countlow==1) break;
      }
      // for (Int_t i = bin_second; i < nbins; i++){
      //   if(h->GetBinContent(i+1)<= 10*h->GetMinimum(0)*scale){
      //     lowestpt = h->GetBinCenter(i);
      //     break;
      //   }
      //   if (i == nbins-1){
      //     lowestpt = h->GetBinCenter(i);
      //   }
      // }
      n_peaks = Int_t(lowestpt/(fPositionX[1]-fPositionX[0])); // this should probably be -1 !!
      if(n_peaks <= 0 || nfound-pos0 < 3){
        cout << "Something wrong: npeaks = " << n_peaks << endl;
        cout << "Changing to 5" << endl;
        n_peaks = 5;
      }
      else{
        if(n_peaks > 10){ // this value was 18 at some point
          n_peaks = 10;
        }
        xmax = (n_peaks+1)*(fPositionX[1]-fPositionX[0]) + fPositionX[0];
      }


      Double_t total_events = h->Integral("width");
      Double_t zero_events = faux->Integral(xmin, xmax);
      Double_t prob_zeros = zero_events/total_events;

      Double_t width = h->GetBinWidth(1);
      Double_t startbin = h->GetBinCenter(1);

      // bin_center = (nbins-1)*width + startbin;
      Int_t npois = 3;
      Int_t bin_n_peak = round(npois*mean1 - startbin)/width + 1;

      Double_t lambda1 = -TMath::Log(prob_zeros);
      Double_t lambda2 = h->GetMean()/mean1;
      Double_t lambda = lambda1;

      Double_t peak_by_formula1 = startpoint/(pow(poisson(lambda1,npois-1)/poisson(lambda1,npois),npois-2));
      Double_t peak_by_formula2 = startpoint/(pow(poisson(lambda2,npois-1)/poisson(lambda2,npois),npois-2));
      Double_t diff1 = abs(h->GetBinContent(bin_n_peak) - peak_by_formula1);
      Double_t diff2 = abs(h->GetBinContent(bin_n_peak) - peak_by_formula2);


      lambda = diff1 <= diff2 ? lambda1 : lambda2;

      Double_t poisson2 = poisson(lambda,2);
      Double_t poisson3 = poisson(lambda,3);
      Double_t poisson_ratio = poisson2/poisson3;


      if(debug)
      {
        TCanvas *cdb = new TCanvas("cdb", "cdb",1920,0,700,500);

        cdb->cd();
        h->GetYaxis()->SetTitle("# of entries");
        h->GetXaxis()->SetTitle("Charge (ADC*ns)");
        h->Draw("hist");
        TPolyMarker * pm = new TPolyMarker(nfound-pos0, &fPositionX[0], &fPositionY[0]);
        pm->SetMarkerStyle(23);
        pm->SetMarkerColor(kRed);
        pm->SetMarkerSize(1.3);
        pm->Draw("SAME");

        for (Int_t i = 0; i < nbins; i++) d->SetBinContent(i + 1,destVector[i]);
        d->SetLineColor(kRed);
        d->Draw("SAME");

        if (!quitemode)
        {
          printf("Found %d candidate peaks\n",nfound);
          for(Int_t i=0;i<nfound;i++) printf("posx= %f, posy= %f\n",fPositionX[i], fPositionY[i]);
        }
        faux->Draw("SAME");
        if (!quitemode) cout << "npeaks = " << n_peaks << " lowest = " << lowestpt << " spe = " << (fPositionX[1]-fPositionX[0]) << endl;
        if(quitemode){
          delete cdb;
          delete d;
        }
      }



      // delete f1, h, source, destVector;
    }



    // ____________________________________________________________________________________________________ //
    string startingPump(){
      string f = "gaus(0) + gaus(3)";
      Int_t aux = 0;
      for(Int_t i = 0; i<n_peaks; i++){
        f = f + " + gaus(" + to_string(i+6+aux) + ")";
        aux = aux+2;
      }
      return f;
      //         TF1 *func = new TF1("func","gaus(0)+gaus(3)+gaus(6)+gaus(9)+gaus(12)+gaus(15)",-2000,5000);

    }

    // ____________________________________________________________________________________________________ //
    void getMyParameters(Double_t peaks[],Double_t stdpeaks[],TF1 *func){
      for(Int_t i=0; i<n_peaks; i++){
        peaks[i] = (i+2)*(func->GetParameter(4));
        stdpeaks[i] = TMath::Power((i+2),0.5)*func->GetParameter(5);
        //         cout << peaks[i] << endl;
      }
    }

    // ____________________________________________________________________________________________________ //
    void makeSphe(string histogram){

      string histogram_tempo = histogram+"_stat";
      TFile *f1 = nullptr;
      if (rootFile != "") {
        TFile *f1 = new TFile(rootFile.c_str(),"READ");
        hcharge = (TH1D*)f1->Get(histogram.c_str());
        if(!hcharge) hcharge = (TH1D*)f1->Get("analyzed");
      }
      else{
        hcharge = (TH1D*)htemp->Clone("");
      }
      hcharge->Rebin(rebin);
    
      ofstream out;
      out.open(Form("sphe%i.txt", channel), ios::app);

      // ____________________________ Start of sphe fit ____________________________ //
      // hcharge->GetYaxis()->SetTitle("Normalized count");
      hcharge->GetYaxis()->SetTitle("# of events");
      hcharge->GetYaxis()->SetTitleOffset(1.0);
      hcharge->GetXaxis()->SetTitle("Normalized charge [A.U]");

    
      TCanvas *c = new TCanvas("c", "c",1920,0,1920,1080);

      // c1->SetLogy();
      gPad->SetGrid(1,1);
      gPad->SetTicks(1,1);
      if (show_all_parameters){
        gStyle->SetOptFit(1111);
      }
      else{
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
      }

      //First function, will almost fit freely
      TF1 *func = new TF1("func",startingPump().c_str(),xmin,xmax);
      TF1 *fu[2+n_peaks];
      string funame;
      for(Int_t i = 0; i<(2+n_peaks); i++){
        funame = "fu_"+to_string(i);
        fu[i] = new TF1(funame.c_str(),"gaus(0)",xmin,xmax);
      }
    
      Double_t peaks[n_peaks];
      Double_t stdpeaks[n_peaks];
    
      if(peak0==0){
        fixZero = true;
      }
    
      func->SetParameters(peak0,mean0,sigma0,peak1,mean1,sigma1); // this values can change
      // func->SetParLimits(2,0.9*sigma0,2*sigma0);
    
      Int_t aux = 0;
      // TO BE CHECKED ! Maybe it is better to not control the startpoing so hardly
      Double_t temp_startpoint = startpoint;
      for(Int_t i = 0; i<n_peaks; i++){
        func->SetParameter((i+6+aux),temp_startpoint);
        aux++;
        func->SetParameter((i+6+aux),(i+1)*(mean1 - mean0) + mean1);
        aux++;
        func->SetParameter((i+6+aux),sqrt(i+2)*sigma1);
        temp_startpoint = temp_startpoint/poisson_ratio;
      }
      func->SetParName(4,"#mu");
      func->SetParName(5,"#sigma");
    
      func->SetNpx(1000);
    
      if(darknoise){
        func->FixParameter(4,sphe_charge);
        func->FixParameter(7,sphe_charge2);
      }
    
      if(fixZero){
        func->FixParameter(0,0);
        func->FixParameter(1,0);
        func->FixParameter(2,1);
      }
    
      getMyParameters(peaks,stdpeaks,func);
    
    
      Double_t scale = 1/(hcharge->Integral());
      //hcharge->Scale(scale);
      hcharge->Draw("hist");
      // hcharge->Fit("func","R0Q");
      // Debug level:
      // 0 none
      // 1 first general fit
      // 2 second general fit
      // 3 first lastOne fit wiht fix parameters

      Int_t debug_level = 0;
      if (debug_level == 1){
        func->Draw("SAME");
        return;
      }

      // recovery parameters
      getMyParameters(peaks,stdpeaks,func);
      aux=0;
      // use as fixed
      for(Int_t i = 0; i<n_peaks; i++){

        func->FixParameter((i+7+aux),peaks[i]-(i+1)*func->GetParameter(1));
        aux++;
        func->FixParameter((i+7+aux),stdpeaks[i]);
        aux++;
      }

      if(darknoise){
        func->FixParameter(4,sphe_charge);
        func->FixParameter(7,sphe_charge2);
      }

      if(fixZero){
        func->FixParameter(0,0);
        func->FixParameter(1,0);
        func->FixParameter(2,1);
      }
      hcharge->Fit("func","R0Q");
      if (debug_level == 2){
        func->Draw("SAME");
        return;
      }

      // set a new function, now fixed for real

      MyFunctionObject MyFunc;
      MyFunc.n_peaks = n_peaks;

      TF1 *lastOne = new TF1("lastOne",MyFunc,xmin,xmax,7+n_peaks);


      lastOne->SetParameter(0,func->GetParameter(0));
      lastOne->SetParameter(1,func->GetParameter(1));
      lastOne->SetParameter(2,func->GetParameter(2));
      // lastOne->SetParLimits(2,0*func->GetParameter(2),0.5*func->GetParameter(2));

      lastOne->SetParameter(3,func->GetParameter(3));
      lastOne->SetParameter(4,func->GetParameter(4));
      lastOne->SetParameter(5,func->GetParameter(5));

      lastOne->SetParameter(6,func->GetParameter(6));
      lastOne->SetParameter(7,func->GetParameter(7));
      aux = 0;
      for(Int_t i = 1; i<n_peaks; i++){
        lastOne->SetParameter((i+8-1),func->GetParameter(i+8+aux));
        lastOne->SetParName((i+8-1),Form("A_{%d}",i+2));
        aux=aux+2;
      }

      lastOne->SetParName(0,"A_{baseline}");
      lastOne->SetParName(1,"#mu_{baseline}");
      lastOne->SetParName(2,"#sigma_{baseline}");
      lastOne->SetParName(3,"A_{1}");
      lastOne->SetParName(4,"#mu_{1}");
      lastOne->SetParName(5,"#sigma_{1}");
      lastOne->SetParName(6,"A_{2}");
      lastOne->SetParName(7,"#mu_{2}");
      // lastOne->SetParName(8,"#sigma_{2}");
    
      if(darknoise){
        lastOne->FixParameter(4,sphe_charge);
        lastOne->FixParameter(7,sphe_charge2);
      }
      if(fixZero){
        lastOne->FixParameter(0,0);
        lastOne->FixParameter(1,0);
        lastOne->FixParameter(2,1);
      }
      fit_status = -1;

      string lastOneFitOpt = "RQ0";
      if (!quitemode) lastOneFitOpt = "R";
      for(Int_t ktemp = 0; ktemp < lastOneAttempts; ktemp++)
      {
        hcharge->Fit("lastOne",lastOneFitOpt.c_str());
        fit_status = TestFitSuccess();
        if (fit_status == 1) break;
      }
      if(!quitemode) cout << "Fit status before free std dev: " << fit_status << endl;

      if(make_free_stddevs == false){
        fu[0]->SetParameter(0,abs(lastOne->GetParameter(0)));
        fu[0]->SetParameter(1,lastOne->GetParameter(1));
        fu[0]->SetParameter(2,lastOne->GetParameter(2));
    
        fu[1]->SetParameter(0,abs(lastOne->GetParameter(3)));
        fu[1]->SetParameter(1,lastOne->GetParameter(4));
        fu[1]->SetParameter(2,lastOne->GetParameter(5));
    
        fu[2]->SetParameter(0,abs(lastOne->GetParameter(6)));
        fu[2]->SetParameter(1,lastOne->GetParameter(7));
        fu[2]->SetParameter(2,(TMath::Power((2),0.5)*lastOne->GetParameter(5)));
    
        for(Int_t i = 1; i<n_peaks; i++){
          fu[i+2]->SetParameter(0,abs(lastOne->GetParameter(i+7)));
          fu[i+2]->SetParameter(1,(lastOne->GetParameter(4) + (i+1)*(lastOne->GetParameter(7)-lastOne->GetParameter(4))));
          fu[i+2]->SetParameter(2,(TMath::Power((i+2),0.5)*lastOne->GetParameter(5)));
        }
      }
      if (debug_level == 3){
        lastOne->Draw("SAME");
        for(Int_t i = 0; i<(2+n_peaks); i++){
          fu[i]->SetLineColor(kGray+1);
          fu[i]->SetNpx(1000);
          fu[i]->Draw("SAME");
        }
        return;
      }

      MyFunctionFree MyFuncFree;
      MyFuncFree.n_peaks = n_peaks;

      TF1 *lastOneFree = new TF1("lastOneFree",MyFuncFree,xmin,xmax,7+n_peaks*2);

      if(make_free_stddevs == true){
        setParametersFree(lastOneFree,lastOne);

        // Set limit for the gaussians
        Double_t refstd = 0;
        for(Int_t i = 0; i<7+n_peaks*2; i++){
          if (i == 2) refstd = lastOneFree->GetParameter(i); // set reference as from the baseline
          if ((i == 5 || i == 8) | (i > 8 && (i-8)%2 == 0)){ // all other gaussians
            if (lastOneFree->GetParameter(i) < 0.5*refstd | lastOneFree->GetParameter(i) > 5*refstd){
              lastOneFree->SetParameter(i, refstd);
            }
            lastOneFree->SetParLimits(i, 0.5*refstd, 5*refstd);
          }
        }
        // Int_t lastpar = 7+2*n_peaks-1;
        // Double_t lastGausUpperLim = 1.2*lastOneFree->GetParameter(lastpar);
        // Double_t lastGausLowerLim = 0.8*lastOneFree->GetParameter(lastpar);
        // lastOneFree->SetParLimits(lastpar,lastGausLowerLim,lastGausUpperLim);
        // cout << lastGausLowerLim << endl;
        // cout << lastGausUpperLim << endl;

        string lastOneFitFreeOpt = "RQ0S";
        if (!quitemode) lastOneFitFreeOpt = "RS";
        hcharge->Fit("lastOneFree",lastOneFitFreeOpt.c_str());
        chi2 = lastOneFree->GetChisquare();
        ndf = lastOneFree->GetNDF();
        fit_status = TestFitSuccess();
        if(!quitemode) cout << "Fit status with free std dev: " << fit_status << endl;


        fu[0]->SetParameter(0,abs(lastOneFree->GetParameter(0)));
        fu[0]->SetParameter(1,lastOneFree->GetParameter(1));
        fu[0]->SetParameter(2,lastOneFree->GetParameter(2));

        fu[1]->SetParameter(0,abs(lastOneFree->GetParameter(3)));
        fu[1]->SetParameter(1,lastOneFree->GetParameter(4));
        fu[1]->SetParameter(2,lastOneFree->GetParameter(5));

        fu[2]->SetParameter(0,abs(lastOneFree->GetParameter(6)));
        fu[2]->SetParameter(1,lastOneFree->GetParameter(7));
        fu[2]->SetParameter(2,lastOneFree->GetParameter(8));

        for(Int_t i = 1, j = 1; i<n_peaks; i++){
          fu[i+2]->SetParameter(0,abs(lastOneFree->GetParameter(j+8)));
          fu[i+2]->SetParameter(1,(lastOneFree->GetParameter(4) + (i+1)*(lastOneFree->GetParameter(7)-lastOneFree->GetParameter(4))));
          fu[i+2]->SetParameter(2,lastOneFree->GetParameter(j+1+8));
          j+=2;
        }

      }
      else{
        chi2 = lastOneFree->GetChisquare();
        ndf = lastOneFree->GetNDF();
        fit_status = TestFitSuccess();
      }


      lastOne->SetNpx(1000);
      lastOneFree->SetNpx(1000);

      hcharge->Draw("hist");
      

      xmin = fu[0]->GetParameter(1)-5*abs(fu[0]->GetParameter(2));

      hcharge->GetXaxis()->SetRangeUser(1.2*xmin,1.1*xmax);
      // hcharge->StatOverflows(kTRUE);
      hcharge->ClearUnderflowAndOverflow();
      lastOne->SetRange(xmin,xmax);
      lastOneFree->SetRange(xmin,xmax);


      for(Int_t i = 0; i<(2+n_peaks); i++){
        fu[i]->SetLineColor(kGray+1);
        fu[i]->SetNpx(1000);
        fu[i]->Draw("SAME");
      }
      if(make_free_stddevs==false) lastOne->Draw("LP SAME");
      else {
        lastOneFree->Draw("LP SAME");
        lastOne = (TF1*)lastOneFree->Clone("lastOneShow");
      }
      snr = abs((lastOne->GetParameter(4))/lastOne->GetParameter(2));

      if (!quitemode){
        cout << "1th peak = " << lastOne->GetParameter(4) << endl;
        cout << "2th peak = " << lastOne->GetParameter(7) << endl;
        cout << "Std dev  = " << lastOne->GetParameter(5) << endl;
        cout << "sphe charge = " << lastOne->GetParameter(7) - lastOne->GetParameter(4) << endl;
        cout << " SNR = " << (lastOne->GetParameter(4))/sqrt(pow(lastOne->GetParameter(2),2)+pow(lastOne->GetParameter(5),2)) << endl;
        cout << " SNR2 = " << snr << endl;

        out.seekp(0, ios::end);
        if (out.tellp() == 0) {
          out << "# mu1 mu2 m2-m1 sigma1 b sigmab" << endl;
        }
        Double_t _mu1 = abs(lastOne->GetParameter(4))*normfactor;
        Double_t _mu2 = abs(lastOne->GetParameter(7))*normfactor;
        Double_t _mu21 = _mu2-_mu1;
        Double_t _sigma1 = abs(lastOne->GetParameter(5))*normfactor;
        Double_t _b = abs(lastOne->GetParameter(1))*normfactor;
        Double_t _sigmab = abs(lastOne->GetParameter(2))*normfactor;
        out <<  _mu1  << " " << _mu2 << " " << _mu21 << " "
            << _sigma1 << " " << _b << " " << _sigmab << " ";

        out << std::fixed << std::setprecision(2) << snr << endl;
      }

      // ____________________________ Finish of sphe fit ____________________________ //


      //_________________ Drawing lines for the cross-talk probability _________________ //

      sphe_charge = lastOne->GetParameter(4);
      sphe_charge2 = lastOne->GetParameter(7);
      Double_t sphe_std = lastOne->GetParameter(5);

      if(deltaminus!=0)
      {
        delta1 = (sphe_charge2 - sphe_charge)/deltaminus;
        delta2 = deltaplus*(sphe_charge2 - sphe_charge);
      }
      else{
        delta1 = (sphe_charge-deltaplus*sphe_std);
        delta2 = (sphe_charge+deltaplus*sphe_std);
      }
      //     Double_t delta2 = sphe_charge+(sphe_charge2 - sphe_charge)/2;

      Double_t ymax = hcharge->GetMaximum();

      TLine *l1 = new TLine(delta1,0,delta1,ymax);
      TLine *l2 = new TLine(delta2,0,delta2,ymax);
      l1->SetLineWidth(2);
      l2->SetLineWidth(2);
      l1->SetLineColor(kRed);
      l2->SetLineColor(kRed);

      if (drawDebugLines){
        l1->Draw("");
        l2->Draw("");
      }

    
      if(darknoise){
        // ____________________________ Start dark noise CT analysis ____________________________ //
        
        if (!quitemode) cout << "\n\n\nMaking ratio between 1 sphe and 2 sphe: " << endl;
        Double_t zeroAmp = lastOne->GetParameter(0);
        Double_t zeroMean = lastOne->GetParameter(1);
        Double_t zeroSigma = lastOne->GetParameter(2);
        
        TF1 *zeroGaus = new TF1("zeroGaus","gaus(0)",xmin,xmax);
        zeroGaus->SetParameters(zeroAmp,zeroMean,zeroSigma);
        
        Double_t zeroIntegral = zeroGaus->Integral(xmin,xmax);
        Double_t totalIntegral = hcharge->Integral("width");
        
        // ------------ NOTE ------------ //
        // The integral of the histogram
        // could be taken, but the effort
        // is bigger.
        // Also, remember to multiply by
        // the histogram bin width
        
        
        Double_t oneAmp = lastOne->GetParameter(3);
        Double_t oneMean = lastOne->GetParameter(4);
        Double_t oneSigma = lastOne->GetParameter(5);
        
        
        TF1 *oneGaus = new TF1("oneGaus","gaus(0)",xmin,xmax);
        oneGaus->SetParameters(oneAmp,oneMean,oneSigma);
        
        Double_t oneIntegral = oneGaus->Integral(xmin,xmax);
        
        
        
        Double_t twoAmp = lastOne->GetParameter(6);
        Double_t twoMean = lastOne->GetParameter(7);
        Double_t twoSigma = TMath::Power(2,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *twoGaus = new TF1("twoGaus","gaus(0)",xmin,xmax);
        twoGaus->SetParameters(twoAmp,twoMean,twoSigma);
        
        Double_t twoIntegral = twoGaus->Integral(xmin,xmax);

        if (!quitemode)
        {
          cout << "Total 1 = " << oneIntegral << " \t total 2 =  " << twoIntegral << endl;
          cout << "Ratio = " << twoIntegral/oneIntegral << endl;
          cout << "Total events normalized = " << lastOne->Integral(xmin,xmax) << endl;
          cout << "CT probability = " << twoIntegral/(oneIntegral+twoIntegral) << endl;
        }

        
        
        
        
        Double_t threeAmp = lastOne->GetParameter(8);
        Double_t threeMean = 2*twoMean-oneMean;
        Double_t threeSigma = TMath::Power(3,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *threeGaus = new TF1("threeGaus","gaus(0)",xmin,xmax);
        threeGaus->SetParameters(threeAmp,threeMean,threeSigma);
        
        Double_t threeIntegral = threeGaus->Integral(xmin,xmax);
        
        if (!quitemode) cout << "\n\n\n Another possibility would be CT = " << (twoIntegral+threeIntegral)/(oneIntegral+twoIntegral+threeIntegral) << endl;
        
        // ____________________________ Finish Poisson analysis ____________________________ //


      }
    
    
    
      if(is_poisson_test){
        
        Double_t zeroAmp = lastOne->GetParameter(0);
        Double_t zeroMean = lastOne->GetParameter(1);
        Double_t zeroSigma = lastOne->GetParameter(2);
        
        TF1 *zeroGaus = new TF1("zeroGaus","gaus(0)",xmin,xmax);
        zeroGaus->SetParameters(zeroAmp,zeroMean,zeroSigma);
        
        Double_t zeroIntegral = zeroGaus->Integral(xmin,xmax);
        Double_t totalIntegral = hcharge->Integral("width");
        
        // ------------ NOTE ------------ //
        // The integral of the histogram
        // could be taken, but the effort
        // is bigger.
        // Also, remember to multiply by
        // the histogram bin width
        
        
        Double_t oneAmp = lastOne->GetParameter(3);
        Double_t oneMean = lastOne->GetParameter(4);
        Double_t oneSigma = lastOne->GetParameter(5);
        
        
        TF1 *oneGaus = new TF1("oneGaus","gaus(0)",xmin,xmax);
        oneGaus->SetParameters(oneAmp,oneMean,oneSigma);
        
        Double_t oneIntegral = oneGaus->Integral(xmin,xmax);
        
        
        
        Double_t twoAmp = lastOne->GetParameter(6);
        Double_t twoMean = lastOne->GetParameter(7);
        Double_t twoSigma = TMath::Power(2,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *twoGaus = new TF1("twoGaus","gaus(0)",xmin,xmax);
        twoGaus->SetParameters(twoAmp,twoMean,twoSigma);
        
        Double_t twoIntegral = twoGaus->Integral(xmin,xmax);
        
        
        Double_t threeAmp = lastOne->GetParameter(8);
        Double_t threeMean = 2*twoMean-oneMean;
        Double_t threeSigma = TMath::Power(3,0.5)*lastOne->GetParameter(5);
        
        
        TF1 *threeGaus = new TF1("threeGaus","gaus(0)",xmin,xmax);
        threeGaus->SetParameters(threeAmp,threeMean,threeSigma);
        
        Double_t threeIntegral = threeGaus->Integral(xmin,xmax);
        

        Double_t lambda = -TMath::Log(zeroIntegral/totalIntegral);
        
        cout << "Lambda = " << lambda << endl;
        
        Double_t ct = hcharge->GetMean()/(lambda*oneMean);
        
        cout << "Cross-talk = " << ct << endl;
        
        
        
        //         cout << "teste -> " << hcharge->GetMean() << endl;
        // ____________________________ Finish Poisson analysis ____________________________ //
      }
      // ____________________________ Draw legend by hand ____________________________ //

      if (show_all_parameters) return;
      TPaveText *pleg = new TPaveText(0.74,0.29,0.89,0.89, "NDC");
      pleg->SetFillColor(0);
      pleg->SetBorderSize(1);
      TText *ltext = nullptr;
      pleg->SetTextFont(42);
      pleg->SetTextSize(0.03);
      pleg->SetTextAlign(11);
      ltext = pleg->AddText("SPE Calibration");
      ltext->SetTextFont(62);
      ltext->SetTextAlign(22);
      ltext = pleg->AddText(" ");
      ltext = pleg->AddText(Form("Entries: %.0f", hcharge->GetEntries()));
      ltext = pleg->AddText(Form("Mean: %.2f", hcharge->GetMean()));
      ltext = pleg->AddText(Form("Std Dev: %.2f", hcharge->GetStdDev()));
      ltext = pleg->AddText(Form("#chi^{2} / ndf: %.2f/%d", lastOne->GetChisquare(), lastOne->GetNDF()));
      ltext = pleg->AddText(Form("Fitted gaussians: %d", n_peaks+2));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(0),abs(lastOne->GetParameter(0)), lastOne->GetParError(0)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(1),lastOne->GetParameter(1), lastOne->GetParError(1)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(2),abs(lastOne->GetParameter(2)), lastOne->GetParError(2)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(3),abs(lastOne->GetParameter(3)), lastOne->GetParError(3)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(4),abs(lastOne->GetParameter(4)), lastOne->GetParError(4)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(5),abs(lastOne->GetParameter(5)), lastOne->GetParError(5)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(6),abs(lastOne->GetParameter(6)), lastOne->GetParError(6)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(7),abs(lastOne->GetParameter(7)), lastOne->GetParError(7)));
      ltext = pleg->AddText(Form("%s: %.2f #pm %.2f",lastOne->GetParName(8),abs(lastOne->GetParameter(8)), lastOne->GetParError(8)));

      // ((TText*)pleg->GetListOfLines()->Last())->SetTextAlign(23);
      pleg->Draw();
      // ____________________________ FinishDraw legend by hand ____________________________ //

      if (nameplotpng == "") nameplotpng = histogram + ".png";
      if (save_plot)
          c->Print(nameplotpng.c_str());
      if(quitemode)
      {
        delete c;
        if(rootFile=="") delete hcharge;
      }

    }

    void setParametersFree(TF1 *lastOneFree, TF1 *lastOne){
      lastOneFree->SetParameter(0,lastOne->GetParameter(0));
      lastOneFree->SetParameter(1,lastOne->GetParameter(1));
      lastOneFree->SetParameter(2,lastOne->GetParameter(2));
      // lastOneFree->SetParLimits(2,0*lastOne->GetParameter(2),0.5*lastOne->GetParameter(2));

      lastOneFree->SetParameter(3,lastOne->GetParameter(3));
      lastOneFree->SetParameter(4,lastOne->GetParameter(4));
      lastOneFree->SetParameter(5,lastOne->GetParameter(5));

      lastOneFree->SetParameter(6,lastOne->GetParameter(6));
      lastOneFree->SetParameter(7,lastOne->GetParameter(7));
      lastOneFree->SetParameter(8,TMath::Power((2),0.5)*lastOne->GetParameter(5));

      for(Int_t i = 1, j= 1; i<n_peaks; i++){
        lastOneFree->SetParameter((j+9-1),lastOne->GetParameter(i+8-1));
        lastOneFree->SetParameter((j+1+9-1),(TMath::Power((i+2),0.5)*lastOne->GetParameter(5)));
        Double_t lastGausUpperLim = (1+stdVar)*lastOneFree->GetParameter(j+1+9-1);
        Double_t lastGausLowerLim = stdVar*lastOneFree->GetParameter(j+1+9-1);
        lastOneFree->SetParLimits(j+1+9-1,lastGausLowerLim,lastGausUpperLim);
        lastOneFree->SetParName((j+9-1),Form("A_{%d}",i+2));
        lastOneFree->SetParName((j+1+9-1),Form("#sigma_{%d}",i+2));
        // cout << "parameter " << j+9-1 << " = " << lastOneFree->GetParameter(j+9-1);
        // cout << " ... parameter " << j+1+9-1 << " = " << lastOneFree->GetParameter(j+1+9-1) << endl;
        j+=2;
      }


      lastOneFree->SetParName(0,"A_{b}");
      lastOneFree->SetParName(1,"#mu_{b}");
      lastOneFree->SetParName(2,"#sigma_{b}");
      lastOneFree->SetParName(3,"A_{1}");
      lastOneFree->SetParName(4,"#mu_{1}");
      lastOneFree->SetParName(5,"#sigma_{1}");
      lastOneFree->SetParName(6,"A_{2}");
      lastOneFree->SetParName(7,"#mu_{2}");
      lastOneFree->SetParName(8,"#sigma_{2}");

      if(darknoise){
        lastOneFree->FixParameter(4,sphe_charge);
        lastOneFree->FixParameter(7,sphe_charge2);
      }
      if(fixZero){
        lastOneFree->FixParameter(0,0);
        lastOneFree->FixParameter(1,0);
        lastOneFree->FixParameter(2,1);
      }
    }

    bool TestFitSuccess(bool verbose = false)
    {
      std::string minuitstatus = std::string(gMinuit->fCstatu);
      if(minuitstatus.compare("CONVERGED ") != 0 && minuitstatus.compare("OK        ") != 0) //the spaces are important
      {
        if(verbose)
          std::cout << "  Minimization did not converge! (status_\"" << minuitstatus << "\")" << std::endl;
        return false;
      }
      else
        return true;
    }

    vector<Double_t> perform_fit(Int_t ch = 0, Int_t rebin = 1, Double_t deconv = 4, Bool_t quite=true, TH1D *hspe = nullptr){

      this->rebin = rebin;
      this->channel = ch;
      this->quitemode = quite;
      string histogram;
      if (hspe){
        this->rootFile = "";
        histogram = hspe->GetName();
        this->htemp = hspe;
      }
      else{
        histogram = "analyzed";
        if(rootFile == "") this->rootFile = "sphe_histograms_Ch"+to_string(ch)+".root";
      }

      bool tmpfree = this->make_free_stddevs;
      this->make_free_stddevs = true; // starts with false, if good fitting, change to true
      this->searchParameters(histogram.c_str(), deconv, true); // give a first search in the parameters.
      this->make_free_stddevs = tmpfree;

      this->deltaplus = 1;
      this->deltaminus = 0;

      // this->drawDebugLines = true;
      this->fit_sphe_wave(histogram.c_str(),false); // set true to make if you want to execute "searchParameters" inside here instead
      vector<Double_t> ret = {static_cast<Double_t>(this->fit_status), this->snr, this->chi2, this->ndf};
      this->snr = ret[1];
      return ret;
    }

    void sphe_fit_try_hard(TH1D *h, Int_t ch = 0, Int_t minrebin = 1, Int_t maxrebin = 2, Double_t min_sigma = 1, Double_t max_sigma = 8, bool save_last = false, Double_t tolerance = 0.1, Double_t sigma_step = 1){

      gROOT->SetBatch(kTRUE);
      save_plot = false;
      vector<Double_t> vals = {0,0,0,0};
      vector<Double_t> values_best_fit = {0,0,0,0};
      vector<Double_t> values_best_snr = {0,0,0,0};
      vector<Double_t> ref = {0,0,0,0};
      Double_t best_fit = 1e12;
      Double_t best_snr = 0;
      for (Int_t rebin = minrebin; rebin <= maxrebin; rebin*=2){
        for(Double_t deconv = min_sigma; deconv <= max_sigma; deconv+=sigma_step){
          cout << "Trying rebin: " << rebin << " deconv: " << deconv << endl;
          vals = perform_fit(ch,rebin, deconv, true, h);
          if (vals[0] > 0){
            Double_t goodf = vals[2]/vals[3];
            if (goodf < best_fit){
              best_fit = goodf;
              values_best_fit = {static_cast<Double_t>(rebin), static_cast<Double_t>(deconv), vals[1], goodf};
            }
            if (vals[1] > best_snr){
              best_snr = vals[1];
              values_best_snr = {static_cast<Double_t>(rebin), static_cast<Double_t>(deconv), vals[1], goodf};
            }
          }
        }

      }
      cout << "...............: rebin deconv snr chi2/ndf" << endl;
      cout << "values_best_snr:" << " ";
      cout << values_best_snr[0] << " ";
      cout << values_best_snr[1] << " ";
      cout << values_best_snr[2] << " ";
      cout << values_best_snr[3] << endl;
      cout << "values_best_fit:" << " ";
      cout << values_best_fit[0] << " ";
      cout << values_best_fit[1] << " ";
      cout << values_best_fit[2] << " ";
      cout << values_best_fit[3] << endl;
      if (values_best_snr[2] != values_best_fit[2]){
        if (abs(values_best_snr[2] - values_best_fit[2])/values_best_fit[2] < tolerance  || values_best_snr[2] < 0 || abs(values_best_snr[2])>25){
          // nothing to do
        }
        else{
          string Userchoise;
          cout << "Not cool!!" << endl;
          cout << "You have to pick one: SNR [s], Fit[f] ";
          cin >> Userchoise;
          if (Userchoise == "s")
          {
            cout << "Using best for SNR" << endl;
            values_best_fit = values_best_snr;
          }
          else if(Userchoise == "f")
          {
            cout << "Using best for Fit" << endl;
          }
          else{
            return;
          }

        }
      }

      if(values_best_snr == ref){
        cout << "No fit works :( " << endl;
        return;
      }
      if (save_last){
        save_plot = true;
      }
      else{
        gROOT->SetBatch(kFALSE);
      }
      perform_fit(ch, values_best_fit[0], values_best_fit[1], false, h);


    }


    Calibration(string mname = "c"){
      myname = mname;
    }
    
};









// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ //





class SPHE{

  public:


    ///////////////////////////////////////////////////////////////////////////////
    //                      Variables to be changed by user                      //
    ///////////////////////////////////////////////////////////////////////////////

    Bool_t led_calibration = false; // if external trigger + led was used, set true
    // start and finish will be the time of integration

    Bool_t just_a_test = false; // well.. it is just a test, so `just_this` is the total waveforms analysed
    Int_t  just_this   = 200;
    Int_t  channel     = 1;
    string filename    = "analyzed";
    

    vector<Int_t> nshow_range = {0,100}; // will save some debugging waveforms inside the range.
    // Example, nshow_range = {0,2} will save waveforms 0, 1 and 2;

    Double_t tolerance    = 5;    // n sigmas (smoothed) (not used for led)
                                  // if `method` is set to `fix`, the threshold will be the absolute value of tolerance, no baseline is calculated
    Double_t baselineTime = 5000; // is the time to compute the baseline (not used for led)
                                  // If the `method` is set to `dynamic` the baseline is calculated over the range of baselineTime
                                  // and it is updated depending on the next point with a weigth of `influece`
                                  // If `method` is set to `static`, baseline is calculated once using baseLimit as cut
    Double_t baseLimit    = 3;  // higher then this wont contribute to the baseline abs(baseLimit) (not used for led)

    Double_t start  = 0;        // start the search for peaks or start the integration (led)
    Double_t finish = 10250;    // fisish the search or finish the integration (led)

    Double_t timeLow         = 60;   // integration time before peak (not used for led)
    Double_t timeHigh        = 400;  // integration time after peak (not used for led)
    Double_t lowerThreshold  = -1e12;   // threshold to detect noise (normal waveform) (not used for led)
    Double_t maxHits         = 1e12;    // maximum hit before discarding   (not used for led)
    Double_t too_big         = 1000; // if there is a peak > "too_big" .. wait "waiting" ns for next peak
    Double_t waiting         = 1000;
    Double_t filter          = 14;   // one dimentional denoise filter (0 equal no filder)
    Double_t interactions    = 45;   // for moving avarage filter (not used on led)
    Int_t    ninter          = 2;    // N times that moving average is applied
    Double_t diff_multiplier = 1e3;  //derivative can be very small. Use this to make it easier to see
    Bool_t derivate_raw      = true; // If you apply a denoise in the data and want to take derivative on that, set to false

    Double_t dtime = 4.;        // time step in ns


    Double_t get_wave_form = false;    // for getting spe waveforms
    Double_t mean_before   = 0;      // time recorded before and after the peak found
    Double_t mean_after    = 0;
    Double_t sphe_charge   = 1809.52;  // charge of 1 and 2 p.e. (use fit_sphe.C)
    Double_t sphe_charge2  = 3425.95;
    Double_t sphe_std      = 500;      // std dev of the first peak (not needed of deltaminus != 0)
    Double_t spe_max_val_at_time_cut = 20; // after `time_cut`, the signal cannot be higher than this
                                           // this allows to remove after pulses
    Double_t time_cut = 2000; // in ns seconds
    Bool_t cut_with_filtered = true; // will cut using filtered waveform

    // coeficients to surround gaussian of 1 spe.
    // Gain             = (sphe_charge2 - sphe_charge)
    // spe's get events where charge < Gain*deltaplus  and charge < Gain/deltaminus
    // If deltaminus is set to zero, sphe_std*deltaplus will be used instead
    // This value can be checked with fit_sphe.C
    Double_t deltaplus  = 1;
    Double_t deltaminus = 0;


    // Not so common to change
    Double_t social_distance = 2; // demands that there is a minimum distance of social_distance * timeHigh between 2 consecutive peaks found
    string   method          = "derivative"; // `derivative` uses the 1th derivative of the waveform and a fixed threshold
                                             // `dynamic` evaluation of the baseline, search over moving average waveform
                                             // `fix` will not evaluate baseline and use raw threshold over moving average
                                             // See tolerance, baselineTime and baselineLimit above

    bool   check_selection     = true; // uses(or not) variable `selection` to discard
    Bool_t withfilter          = true; // Integrate in the filtered waveform
    Bool_t integrate_from_peak = true; // Set true and the maximum will be searched inside the derivate region.
                                       // Otherwise the integral is done starting from the crossing negative zero of derivative

    Int_t hnbins = 50000;       // Output histogram bins and limits. Does not change this unless you will analyze alway the same device
    Double_t hxmin  = 0;           // The fit function has the `rebin` argument to avoid changing this values
    Double_t hxmax  = 0;
    Bool_t normalize_histogram = false; // will normalize histogram by the average value

    Bool_t do_updateSPEvalues = true; // set to true to automatically search for spheX.txt and retrieve proper spe (need to run twice)


    ///////////////////////////////////////////////////////////////////////////////
    //                         Finish variables for user                         //
    ///////////////////////////////////////////////////////////////////////////////


    // Here are more general ones

    ANALYZER *z      = nullptr; // Main analyser
    TFile    *fout   = nullptr; // File to save histogram
    TFile    *fwvf   = nullptr; // File to save waveforms
    TTree    *twvf   = nullptr; // tree of waveforms

    ADC_DATA sample; // will store the waveforms found
    Double_t *sample_wvf_filtered; // for applying cut with filted wvf


    string myname; // name of the declared class
    Bool_t derivate = false; // control if user has set to derivate
    Int_t kch; // current channel used
    Int_t nshow_finish = 100; // will show `nshow` debugging waveforms up
    Int_t nshow_start = 0; // starting from here
    Bool_t getstaticbase = false; // control if `static` method was set
    Double_t threshold = 0;

    Bool_t get_this_wvf = true;
    Bool_t get_this_charge = true;
    Int_t npts_wvf = 1000;
    Double_t *timeg = nullptr;

    string rootfile    = "analyzed.root";
    string outrootname = "";
    Int_t n_points = 0;


    // ____________________ Variables to calculate and reset ____________________ //

    TH1D *hbase        = nullptr; // not used now
    TH1D *hbase_smooth = nullptr; // compute smoothed baseline
    TH1D *hcharge      = nullptr; // final spe histogram

    Double_t mean = 0; //mean and stddev for hbase
    Double_t stddev = 0;
    vector<Double_t> charge_values;

    vector<Double_t> peakPosition; //store position of the peaks
    vector<Double_t> peakMax; // store peak maximum
    vector<Int_t> peaksFound; // peaks found (threshold apply)
    vector<Int_t> peaksRise; // Rising derivative
    vector<Int_t> peaksCross; // Falling derivative
    vector<Double_t> derivateMax; // Stores the max of the derivative
    vector<Double_t> derivateMaxFound; // Stores the max of the derivative before applying threshold

    vector<Double_t> smooted_wvf; // not needed to reset this
    vector<Double_t> denoise_wvf;

    vector<vector<Double_t>> selected_peaks; // these are for plotting
    vector<vector<Double_t>> selected_time;
    vector<Double_t> temp_selected_peaks;
    vector<Double_t> temp_selected_time;
    vector<Double_t> selected_charge;
    vector<Double_t> selected_charge_time;


    vector<Double_t> discardedTime; //store position of thed
    vector<Double_t> discardedPeak; //store position of thed
    vector<Double_t> discardedCharge; //store position of thed
    vector<Int_t> discarded_idx; //store idx (0, 1 or 2) of the peaks
    const char *select_types[3] = {"Social distance", "Negative hits", "Too big"};
    TH1D *hdiscard = nullptr;
    vector<Double_t> mean_waveform;
    Int_t naverages = 0;

    Double_t delta1 = 0;
    Double_t delta2 = 0;
    // ____________________ ________________________________ ____________________ //

    vector<string> discard_reason = {"Social distance", "Negative hits", "Too big"};

    SPHE(string m_name) : myname{m_name}{
      hbase        = new TH1D(Form("hbase_sphe_%s",myname.c_str()),"histogram for baseline",5*800,-400,400);
      hbase_smooth = new TH1D(Form("hbase_smooth_sphe_%s",myname.c_str()),"histogram for baseline smoothed",5*800,-400,400);
      hcharge      = new TH1D(Form("hcharge_sphe_%s",myname.c_str()),Form("hcharge_sphe_%s",myname.c_str()),hnbins,hxmin,hxmax);
      hdiscard     = new TH1D(Form("hdiscard_sphe_%s",myname.c_str()),Form("hdiscard_sphe_%s",myname.c_str()),3,0,3);
      for(Int_t i = 0; i < (int)(sizeof(select_types)/sizeof(select_types[0])); i++){
        hdiscard->Fill(select_types[i],0);
      }

      smooted_wvf.resize(n_points);
      denoise_wvf.resize(n_points);


    }
    ~SPHE(){
      // if(twvf) delete twvf;
      if(fout) fout->Close();
      if(fwvf) fwvf->Close();
      if(hbase) delete hbase;
      if(hbase_smooth) delete hbase_smooth;
      if(hcharge) delete hcharge;
      if(hdiscard) delete hdiscard;

      // if(z) delete z;
      // if(fout) delete fout;
      // if(fwvf) delete fwvf;
      // if(sample_wvf_filtered) delete sample_wvf_filtered;
      // if(timeg) delete timeg;
    }


    void searchForPeaks(){
    }

    /**
     * This function searches and integrate peaks, used for sphe
     * If `led_calibration` is set to true, integration is made between `start` and `finish`
     * Load everything
     * Create TTree/Histograms necessary
     * Reset area, clean vectors, histograms, etc. for next waveform
     * Process data: moving average, denoise filters, etc.
     * Search for peaks
     * Clean the selection
     * Integrate peaks and store the waveform
     * Plot graphs for debugging
     * Save histogram
     */
    void giveMeSphe(){
      gROOT->SetBatch(kTRUE);

      hcharge->Reset();
      if (hxmax != 0 && hxmin != 0)
        hcharge->SetBins(hnbins, hxmin, hxmax);

      if(!z){ // load analyzer in case it is nullptr
        z = new ANALYZER(myname.c_str());
        z->dtime = dtime; // in case there is no dtime written
        z->setAnalyzer(rootfile);
      }
      if(!z->t1){ // otherwise, check if the ttree was loaded
        z->setAnalyzer(rootfile);
      }
      if(z->setChannel(Form("Ch%d.",channel)) == false) return;

      kch = z->kch;
      dtime = z->dtime; // safety checkpoint
      z->getWaveform(kch);
      n_points = z->n_points;
      smooted_wvf.resize(n_points);
      denoise_wvf.resize(n_points);
      if (outrootname == ""){
        outrootname = Form("sphe_histograms_Ch%i.root",channel);
      }
      fout = new TFile(outrootname.c_str(),"RECREATE");

      // ____________________ Setting up what needs to be set ____________________ //

      if (finish == 0) finish = dtime*n_points;
      if(get_wave_form==true){
        fwvf = new TFile(Form("sphe_waveforms_Ch%i.root",channel),"RECREATE");
        twvf = new TTree("t1","mean waveforms");
        if (led_calibration && mean_before > start) mean_before = start;
        if (led_calibration && mean_after < finish) mean_after = finish+dtime;
      }
      else{ // in the case we are not taking waveform, I change this values for the same setup of integration
        mean_before = (led_calibration) ? start : timeLow;
        mean_after  = (led_calibration) ? finish : timeHigh;
        
      }
      // if the user set mean_before = 24 and mean_after = 100
      // we need to get  6 + 25 points + 1 = 32  (+1 because of the peak start)
      // if I found a point in around 100 ns, that is i = 25, I perform a while that goes
      // from i-mean_before = 19 up to i+mean_after = 50 (included) = 32 pts
      npts_wvf = (mean_before/dtime + mean_after/dtime) + 1;
      if (led_calibration){
        if (mean_after > n_points*z->dtime) mean_after = n_points*dtime;
        npts_wvf = (mean_after/dtime - mean_before/dtime);
      }

      mean_waveform.clear();
      mean_waveform.resize(npts_wvf);
      naverages = 0;

      sample.Set_npts(npts_wvf);
      sample_wvf_filtered = new Double_t[npts_wvf];
      if (get_wave_form) twvf->Branch(Form("Ch%i.",channel),&sample);

      if (do_updateSPEvalues && get_wave_form) updateSPEvalues();


      nshow_start = nshow_range[0];
      nshow_finish = nshow_range[1];
      n_points = z->n_points;
      timeg = new Double_t[n_points];
      for(int i = 0; i < n_points; i++){
        timeg[i] = i*dtime;
      }
      if(led_calibration==true){
        method = "led";
      }
      else{
        // this might be too much, time will tell. Keep an eye in the RAM memory
        discardedTime.reserve(n_points);
        discardedPeak.reserve(n_points);
        discardedCharge.reserve(n_points);
        discarded_idx.reserve(n_points);
        peakPosition.reserve(n_points);
        peakMax.reserve(n_points);
        peaksFound.reserve(n_points);
        derivateMax.reserve(n_points);
        derivateMaxFound.reserve(n_points);
        temp_selected_peaks.reserve(n_points);
        temp_selected_time.reserve(n_points);
        selected_charge.reserve(n_points);
        selected_charge_time.reserve(n_points);
      }

      if(method == "static") getstaticbase = true;
      else if(method == "fix"){
      }
      else if(method == "derivative"){
        derivate = true;
        peaksRise.reserve(n_points);
        peaksCross.reserve(n_points);
      }

      Double_t delta = sphe_charge2 - sphe_charge;
      if(deltaminus == 0){
        delta1 = sphe_charge - deltaplus*sphe_std;
        delta2 = sphe_charge + deltaplus*sphe_std;
      }
      else{
        delta1 = delta/deltaminus;
        delta2 = delta*deltaplus;
      }

      // ____________________ ___________________________ ____________________ //



      Int_t nentries = z->nentries;
      if(just_a_test){nentries = just_this;}
      for(Int_t i = 0; i<nentries; i++){
        theGreatReset();
        z->getWaveform(i,kch); // get waveform by memory
        if(check_selection && z->ch[kch]->selection != 0) continue;
        
        if(i>nshow_finish && static_cast<Int_t>(i)%200==0){
          cout << i << " out of " << nentries << "\r" << flush;
        }
        processData();

        if(method == "static"){
          searchForPeaks();
        }
        else if(derivate){
          // Search for crossing points in the derivated waveform
          // positive crossing stored in peaksRise
          // negative crossing stored in peaksCross
          // from start to finish, discounting what was for getting integral with an extra point (not waveform, this filter is done later)
          Double_t start_of_search_offset = (timeLow>interactions*dtime) ? timeLow : interactions*dtime;
          Double_t finish_of_search_offset = (timeHigh>interactions*dtime) ? timeHigh : interactions*dtime;
          z->zeroCrossSearch(&smooted_wvf[0], peaksRise, peaksCross, start+start_of_search_offset+dtime, finish-finish_of_search_offset-dtime);
          derivateApplyThreshold();
        }
        cleanPeaks();
        integrateSignals(sample);
        if(snap()){
          drawMySamples();
          if(!led_calibration){
            cout << "Event " << z->ch[kch]->event << " total of peaks: " << peaksFound.size() << ", Valid = " << selected_charge.size() << endl;
            for(unsigned int j = 0; j < selected_charge.size(); j++){
              cout << "\t\t Charge = " << selected_charge[j] << " at " << selected_charge_time[j] << " ns " << endl;
            }
          }
        } // draw sample graphs
        
      }
      cout << "\n\n" << endl;


      if(get_wave_form==false){
        // fwvf->Close();
        // system(Form("rm sphe_waveforms_Ch%i.root",kch));
      }
      else{
        fwvf->WriteObject(twvf,"t1","TObject::kOverwrite");
        Int_t ndiv = (naverages == 0 ? 1 : naverages);
        for(Int_t i = 0; i<npts_wvf; i++){
          mean_waveform[i] = mean_waveform[i]/ndiv;
        }

        TGraph *gmean = new TGraph(npts_wvf,timeg,&mean_waveform[0]);
        fwvf->WriteObject(gmean,"mean","TObject::kOverwrite");
        if (z->thead){
          auto newtree = z->thead->CloneTree();
          fwvf->WriteObject(newtree,"head","TObject::kOverwrite");
        }
        cout << "A total of " << naverages << " waveforms where found "<< endl;
      }

      Double_t normfactor = 1.;
      if(normalize_histogram){
        nth_element(charge_values.begin(), charge_values.begin()+charge_values.size()*0.16, charge_values.end());
        nth_element(charge_values.begin(), charge_values.begin()+charge_values.size()*0.84, charge_values.end());
        Double_t minusstd = charge_values[charge_values.size()*0.16];
        Double_t plusstd = charge_values[charge_values.size()*0.84];
        Double_t _mean = 0;
        Int_t _count = 0;
        for (auto v: charge_values){
          if (v >= minusstd && v <= plusstd){
            _mean += v;
            _count += 1;
          }
        }
        // normfactor = abs(std::reduce(.begin(), charge_values.end()) / charge_values.size());
        normfactor = abs(_mean)/_count;
      }
      for(auto v: charge_values){
        hcharge->Fill(v/normfactor);
      }
      charge_values.clear();

      // fout->WriteObject(hcharge,Form("%s_%i",filename.c_str(),z->getIdx()),"TObject::kOverwrite");
      fout->cd();
      fout->WriteObject(hcharge,Form("%s",filename.c_str()),"TObject::kOverwrite");
      if(!led_calibration) fout->WriteObject(hdiscard,"hdiscard","TObject::kOverwrite");


      TTree *theadspe = new TTree("head","header");
      theadspe->Branch("tolerance", &tolerance);
      theadspe->Branch("baselineTime", &baselineTime);
      theadspe->Branch("baseLimit", &baseLimit);
      theadspe->Branch("start", &start);
      theadspe->Branch("finish", &finish);
      theadspe->Branch("timeLow", &timeLow);
      theadspe->Branch("timeHigh", &timeHigh);
      theadspe->Branch("lowerThreshold", &lowerThreshold);
      theadspe->Branch("maxHits", &maxHits);
      theadspe->Branch("too_big", &too_big);
      theadspe->Branch("waiting", &waiting);
      theadspe->Branch("filter", &filter);
      theadspe->Branch("interactions", &interactions);
      theadspe->Branch("ninter", &ninter);
      theadspe->Branch("diff_multiplier", &diff_multiplier);
      theadspe->Branch("derivate_raw", &derivate_raw);
      theadspe->Branch("dtime", &dtime);
      theadspe->Branch("social_distance", &social_distance);
      theadspe->Branch("method", &method);

      theadspe->Branch("check_selection", &check_selection);
      theadspe->Branch("withfilter", &withfilter);
      theadspe->Branch("integrate_from_peak", &integrate_from_peak);
      theadspe->Branch("normfactor", &normfactor);
      theadspe->Fill();
      fout->WriteObject(theadspe, "head", "TObject::kOverwrite");



    }
    void theGreatReset(){
      hbase->Reset();
      if(method != "led"){
        if(derivate){
          peaksRise.clear();
          peaksCross.clear();
        }
        else{
          hbase_smooth->Reset();
        }
        peakPosition.clear();
        peakMax.clear();
        selected_peaks.clear();
        selected_time.clear();
        temp_selected_peaks.clear();
        temp_selected_time.clear();
        selected_charge.clear();
        selected_charge_time.clear();
        discardedTime.clear();
        discardedPeak.clear();
        discardedCharge.clear();
        discarded_idx.clear();
        peaksFound.clear();
        derivateMax.clear();
        derivateMaxFound.clear();
      }
      
      denoise_wvf.clear();
      denoise_wvf.resize(n_points);
      smooted_wvf.clear();
      smooted_wvf.resize(n_points);
    }

    void processData(){
      Double_t *vec = nullptr;

      z->applyDenoise(filter, z->ch[kch]->wvf, &denoise_wvf[0]); // denoise is stored at denoise_wvf. If no filter, they are equal

      /**
       * From https://terpconnect.umd.edu/~toh/spectrum/Differentiation.html
       *It makes no difference whether the smooth operation is applied before or after the differentiation. What is important, however,
       *is the nature of the smooth, its smooth ratio (ratio of the smooth width to the width of the original peak), and the number of
       *times the signal is smoothed. The optimum values of smooth ratio for derivative signals is approximately 0.5 to 1.0. For a first
       *derivative, two applications of a simple rectangular smooth (or one application of a triangular smooth) is adequate. For a second
       *derivative, three applications of a simple rectangular smooth or two applications of a triangular smooth is adequate. The general rule is:
       *for the nth derivative, use at least n+1 applications of a rectangular smooth. (The Matlab signal processing program iSignal automatically
       *provides the desired type of smooth for each derivative order).
       *Usediscarded_idx to check the best option
       **/
      if(derivate){
        if(derivate_raw) z->differenciate(diff_multiplier,z->ch[kch]->wvf,&smooted_wvf[0]); // multiply by 1e3 so we can see something :)
        else z->differenciate(diff_multiplier,&denoise_wvf[0],&smooted_wvf[0]); // multiply by 1e3 so we can see something :)
        vec = &smooted_wvf[0];
      }
      else{
        vec = z->ch[kch]->wvf;
      }
      for(Int_t i = 0; i < ninter; i++){
        if (i == interactions-1 && method == "static") getstaticbase = true;
        smooted_wvf = movingAverage(vec, interactions, getstaticbase);
        vec = &smooted_wvf[0];
      }

    }

    bool goodSocialDistance(Int_t id2, Int_t id1){
      if(abs(id2-id1) <= social_distance*timeHigh/dtime){
        // if(snap()) cout << id2 << " " << id1 << " " << abs(id2-id1) << " " <<  social_distance*timeHigh/dtime << endl;
        return false;
      }
      else{
        return true;
      }
    }
    bool checkTooBig(Bool_t &wait_now, Double_t &refWait, Int_t pos){
      Double_t peakmax;
      if(derivate){
        peakmax = denoise_wvf[pos];
      }
      else{
        cout << "WARNING!!! ................................ DONT!" << endl;
        peakmax = peakMax[pos];
      }
      if(peakmax > too_big){
        wait_now = true;
        refWait = pos*dtime;
        return true;
      }
      else{

        return false;
      }

    }

    void derivateApplyThreshold(){
      Int_t ntotal = (int)peaksCross.size();

      for(Int_t i = 0; i < ntotal; i++){
        Double_t crossPositive = peaksRise[i];
        Double_t candidatePosition = peaksCross[i]; // peaks previously selected

        Double_t dmax = z->getMaximum(crossPositive*dtime, candidatePosition*dtime, &smooted_wvf[0]);
        // if(snap()) cout << crossPositive*dtime << " " << candidatePosition*dtime << " " << dmax  << "\n";
        if(dmax < tolerance) {
          continue;
        }
        if(integrate_from_peak == true){
          z->getMaximum(crossPositive*dtime, candidatePosition*dtime, &denoise_wvf[0]);
          peaksFound.push_back(z->temp_pos/dtime);
        }
        else{
          peaksFound.push_back(peaksCross[i]) ;
        }

        derivateMaxFound.push_back(dmax);
      }

    }
    void cleanPeaks(){
      Int_t ntotal = (int)peaksFound.size();
      Bool_t wait_now = false;
      Double_t refWait = 0;// reference time that started waiting from a big pulse

      Bool_t discard_by_distance = false;
      Bool_t current_good_distance = false;

      for(Int_t i = 0; i < ntotal; i++){
        Double_t candidatePosition = peaksFound[i]; // peaks previously selected

        if(checkTooBig(wait_now, refWait, candidatePosition)){
          discard_by_distance = false;
          // if(snap()){
          //   // cout << wait_now << " " << refWait << " " << candidatePosition << " " << waiting << endl;
          // }
          fill_discarded(candidatePosition, smooted_wvf[candidatePosition], 0, 2);
          continue;
        }

        if(wait_now){ // I am waiting for a big pulse to finish
          if(candidatePosition*dtime >= (refWait+waiting)){
            // if(snap()){
            //   cout << "Got out of big pulse: " << candidatePosition << endl;
            // }
            wait_now = false;
          }
          else{
            fill_discarded(candidatePosition, smooted_wvf[candidatePosition], 0, 2);
            continue;
          }
        }

        //request social distance
        if( (ntotal > 1  && !(current_good_distance = goodSocialDistance(peaksFound[i+1],candidatePosition))) || discard_by_distance){
          // if(snap()) cout << discard_by_distance << endl;
          if (!current_good_distance) discard_by_distance = true;
          else discard_by_distance = false;

          // discarding current peak by social
          fill_discarded(candidatePosition, smooted_wvf[candidatePosition], 0, 0);
          continue;
        }
        peakPosition.push_back(peaksFound[i]);
        derivateMax.push_back(derivateMaxFound[i]);
      }
    }

    void fill_discarded(Int_t pidx, Double_t val, Double_t charge, Int_t discardidx){
      if(snap())
      {
        discardedTime.push_back(pidx*dtime);
        discardedPeak.push_back(val);
        discardedCharge.push_back(charge);
        discarded_idx.push_back(discardidx);
      }
      hdiscard->Fill(select_types[discardidx], 1);
    }

    vector<Double_t> delay_line(vector<Double_t> v, Double_t delay_time){
      if(delay_time==0) return v;
      vector<Double_t> res(v.size());
      for(int i=0; i<(int)v.size(); i++){
        res[i]=v[i] - (i-delay_time>=0 ? v[i-delay_time] : 0);
      }
      return res;
    }

    template<class T>
    vector<Double_t> movingAverage(T* v, Int_t myinte, Bool_t eval_baseline){
      Int_t midpoint = 0;
      Int_t interactions = 0;
      Double_t width = 0;
      Double_t sum = 0;
      Int_t n = n_points;

      vector<Double_t> res(n,0);
      if(myinte==0) { // nothing to do
        for (Int_t i = 0; i < n; i++) {
          res[i] = v[i];
        }
        return res;
      }
      else{
        interactions = myinte;
      }

      if(interactions%2==0){ // if it is even, we insert the middle point, e.g. 8 interactions takes 4 before, mid, 4 later
        midpoint = interactions/2+1;    //midpoint will be 5 here
        width = interactions+1;
      }
      else{
        midpoint = (interactions-1)/2 + 1; // e.g. 9 interactions the midpoint will be 5
        width = interactions;
      }

      for(Int_t i = 0; i < n; i++){
        if(i<midpoint || i>(n-midpoint) || interactions == 0){ // make it to start at i = 5 and finish at i = (3000-5) = 2995
          res[i] = v[i];
        }
        else if(i*dtime > start && i*dtime < finish){
          for(Int_t j = (i-midpoint); j < (i+midpoint); j++) { //first one: from j = (5-5); j<(5+5)
            sum = sum+v[j];
          }
          res[i] = (sum/width);

          if(eval_baseline){ // in case we are computing the baseline
            if(i*dtime<=baselineTime && abs(res[i])<baseLimit){
              hbase_smooth->Fill(res[i]);
            }
          }
        }
        else{

          res[i] = 0;
        }
        sum=0;
      }
      if(eval_baseline){ // in case we are computing the baseline
        mean = hbase_smooth->GetMean();
        stddev = hbase_smooth->GetStdDev();
      }
      return res;
    }

    void integrateSignals(ADC_DATA &sample){
      unsigned int npeaks = peakPosition.size();
      if (led_calibration){
        npeaks = 1;
        peakPosition.push_back(0);
      }

      for(unsigned int i = 0; i < npeaks; i++){
        Int_t peakPosIdx = peakPosition[i]; // position of the peak as idx ( don't worry for led )
        temp_selected_peaks.clear();
        temp_selected_time.clear();
        getIntegral(peakPosIdx, sample); // If get_mean_wvf is set to false, there is no problem!
        Double_t charge = z->temp_charge;
        Double_t peak = z->temp_max;
        if(get_this_charge){
          charge_values.push_back(charge);
          if(snap())
          {
            selected_charge.push_back(charge);
            selected_charge_time.push_back(peakPosIdx*dtime);
            selected_peaks.push_back(temp_selected_peaks);
            selected_time.push_back(temp_selected_time);
          }
        }
        else{
          continue; // If I didn't get the charge, I will not take the waveform
        }
        if(get_this_wvf){
          sample.selection = 0;
          if(charge>=delta1 && charge <= delta2){
            sample.selection = 1;
            for(Int_t j = 0; j < npts_wvf; j++){
              mean_waveform[j] += sample.wvf[j];
              if(j*dtime >= time_cut && ((!cut_with_filtered && sample.wvf[j] >= spe_max_val_at_time_cut) || (cut_with_filtered && sample_wvf_filtered[j] >= spe_max_val_at_time_cut))){
                sample.selection = 2;
                for(Int_t k = j; k >= 0; k--){
                  mean_waveform[k] -= sample.wvf[k];
                }
                break;
              }
            }
            if(sample.selection == 1) naverages+=1;
          }
          else{
          }
          sample.peak = peak;
          sample.peakpos = peakPosIdx*dtime;
          sample.charge = charge;
          sample.event = z->currentEvent;
          if(derivate) sample.fprompt = derivateMax[i];
          if(get_wave_form) twvf->Fill();
        }
      }
    }

    void getIntegral(Int_t peakPosIdx, ADC_DATA &sample){
      get_this_charge = true;
      if(get_wave_form) get_this_wvf = true; // in the case not, I dont take the waveform
      else get_this_wvf = false;
      Double_t from = peakPosIdx - mean_before/dtime;
      Double_t to = peakPosIdx + mean_after/dtime;
      
      if(from < 0 || to >= n_points){ // check if boundaries are being respected
        from = peakPosIdx - timeLow/dtime;
        to = peakPosIdx + timeHigh/dtime;
        if(get_wave_form) get_this_wvf = false; // in the case not, I dont take the waveform
      }
      Int_t iwvf = 0;

      Double_t res = 0;
      Double_t max = -1e12;
      Int_t negativeHits = 0;
      Double_t integralfrom = peakPosIdx - timeLow/dtime;
      Double_t integralto = peakPosIdx + timeHigh/dtime;

      if(led_calibration){
        from = mean_before/dtime;
        to = mean_after/dtime-1;
        integralfrom = start/dtime;
        integralto = finish/dtime;
        if (from > integralfrom || from < 0) from = integralfrom;
        if (to  < integralto || to > n_points) to = n_points-1;
        if(get_wave_form) get_this_wvf = true;
      }
      // cout << "Cycling over " << from*dtime << " to " << to*dtime << " ns " << endl;
      // cout << "Integrating over " << integralfrom*dtime << " to " << integralto*dtime << " ns " << endl;
      for(Int_t i = from; i <= to; i++, iwvf++){
        Double_t val = denoise_wvf[i];
        if(get_this_wvf){
          sample.wvf[iwvf] = z->ch[kch]->wvf[i];
          sample_wvf_filtered[iwvf] = val;
          sample.time = z->ch[kch]->time;
        }
        if(i >= integralfrom && i <= integralto){
          if (val < lowerThreshold){
            negativeHits++;
          }
          
          if(negativeHits >= maxHits || val > too_big){
            get_this_wvf = false;
            get_this_charge = false;
            fill_discarded(peakPosIdx, smooted_wvf[peakPosIdx], res*dtime, 1);
            break;
          }
          if(snap())
          {
            temp_selected_peaks.push_back(val);
            temp_selected_time.push_back(i*dtime);
          }
          if(withfilter == false){ val = z->ch[kch]->wvf[i]; }
          res += val;
          if(val>=max){
            max = val;
          }
        }
      }
      // store values here because of lazyness
      z->temp_charge = res*dtime;
      z->temp_max = max;

    }

    // ____________________________________________________________________________________________________ //
    Bool_t snap(){
      if(z->currentEvent >= nshow_start && z->currentEvent < nshow_finish && !led_calibration){
        return true;
      }
      else{
        return false;
      }
      return false; // just to be sure oO
    }
    void drawMySamples(){
      unsigned int npeaksselected = selected_charge.size();
      unsigned int npeakstotal = peaksFound.size();
      string sampleName = "ev_" + to_string(z->currentEvent)+"_"+to_string(npeaksselected)+"_of_"+to_string(npeakstotal);

      TGraph g_smooth(n_points, timeg, &smooted_wvf[0]);
      TGraph g_normal(n_points, timeg, &denoise_wvf[0]);
      // TCanvas *c1 = new TCanvas(sampleName.c_str(),sampleName.c_str(),1920,0,700,500);
      // this is not working when saving
      TCanvas *c1 = new TCanvas(sampleName.c_str(),sampleName.c_str());


      c1->cd(1);
      g_smooth.SetLineColor(kRed);
      g_smooth.SetLineWidth(3);

      g_normal.SetLineColor(kBlue);
      g_normal.SetTitle(" ");

      g_normal.SetEditable(kFALSE);
      g_smooth.SetEditable(kFALSE);

      g_normal.Draw("AL");
      g_smooth.Draw("L SAME");

      if (derivate || getstaticbase ) {
        threshold = tolerance;
        mean = 0;
      }

      TLine lmean(start,mean,finish,mean);
      TLine ldev(start,mean+threshold,finish,mean+threshold);

      lmean.SetLineColor(kGreen);
      ldev.SetLineColor(kGreen);

      lmean.SetLineWidth(2);
      ldev.SetLineWidth(2);

      lmean.Draw();
      ldev.Draw();
      Int_t n = selected_charge.size();
      vector<TGraph*> g_selected(n);
      for(Int_t i = 0; i < n; i++){
        g_selected[i] = new TGraph(selected_time[i].size(), &selected_time[i][0], &selected_peaks[i][0]);
        g_selected[i]->SetMarkerColor(kBlack);
        stringstream stream;
        stream << std::fixed << std::setprecision(2) << selected_charge[i];
        string gtitle = "charge = " + stream.str();
        g_selected[i]->SetTitle(gtitle.c_str());
        g_selected[i]->SetEditable(kFALSE);
        g_selected[i]->Draw("P* SAME");
      }

      Int_t ndis = discardedPeak.size();
      TGraph *g_discarded = nullptr;
      Int_t marker_style = 21;
      Color_t marker_color = kRed;
      for(Int_t i = 0; i < ndis; i++){
        g_discarded = new TGraph(1, &discardedTime[i], &discardedPeak[i]);
        stringstream stream_discarded;
        stream_discarded << std::fixed << std::setprecision(2) << discardedCharge[i];
        string gtitle_discarded = "charge = " + stream_discarded.str();
        g_discarded->SetTitle(gtitle_discarded.c_str());
        if(discarded_idx[i] == 0){ // social distance
          marker_style = 22; // triangle
          marker_color = kYellow;
        }
        else if(discarded_idx[i] == 1){ // negative hits
          marker_style = 29; // filled star
          marker_color = kGreen;
        }
        else if(discarded_idx[i] == 2){ // Too big
          marker_style = 21; // filled square
          marker_color = kRed;
        }
        g_discarded->SetMarkerSize(2);
        g_discarded->SetMarkerStyle(marker_style);
        g_discarded->SetMarkerColor(marker_color);
        g_discarded->SetEditable(kFALSE);
        g_discarded->Draw("SAME P");
      }

      // if(derivate) z->drawZeroCrossingLines(peaksCross, peaksRise, c1,0,tolerance);
      fout->WriteObject(c1,(sampleName.c_str()),"TObject::kOverwrite");
      for(Int_t i = 0; i < n; i++){
        delete g_selected[i];
      }
      g_selected.clear();

      delete c1;


    }

    void updateSPEvalues(){

      ifstream fsphe;
      string fsphe_name = Form("sphe%d.txt", channel);
      fsphe.open(fsphe_name, ios::in);
      if(fsphe.good() && fsphe.is_open()){ // Ok
        // cout << "Reading file " << files << " ... " << endl;
      }
      else{
        cout << "The file " << fsphe_name << " did not open!!" << endl;
        cout << "Keeping standard values... " << endl;
        return;
      }
      fsphe.seekg(-1,ios_base::end);                // go to one spot before the EOF

      bool keepLooping = true;
      bool first_one = true;
      while(keepLooping) {
        char ch;
        fsphe.get(ch);                            // Get current byte's data

        if((int)fsphe.tellg() <= 1) {             // If the data was at or before the 0th byte
          fsphe.seekg(0);                       // The first line is the last line
          keepLooping = false;                // So stop there
        }
        else if(ch == '\n' && first_one==false) {                   // If the data was a newline
          keepLooping = false;                // Stop at the current position.
        }
        else {                                  // If the data was neither a newline nor at the 0 byte
          fsphe.seekg(-2,ios_base::cur);        // Move to the front of that data, then to the front of the data before it
        }
        first_one = false;
      }

      string lastLine;
      Double_t dummy;
      fsphe >> sphe_charge >> sphe_charge2 >> dummy >> sphe_std;
      // cout << sphe_charge << " " << sphe_charge2 << " " << dummy << " " << sphe_std << endl;

      fsphe.close();

      if(deltaminus!=0)
      {
        delta1 = (sphe_charge2 - sphe_charge)/deltaminus;
        delta2 = deltaplus*(sphe_charge2 - sphe_charge);
      }
      else{
        delta1 = (sphe_charge-deltaplus*sphe_std);
        delta2 = (sphe_charge+deltaplus*sphe_std);
      }

    }
};

              
