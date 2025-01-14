// ________________________________________ //
// Author: Henrique Souza
// Filename: analyzer.C
// Created: 2022
// ________________________________________ //
//
#include "WVFAnaSoLAr.h"



class ANALYZER{

  public:
    vector<vector<Double_t>> _baseIntervals;

    TFile *f = nullptr;
    TTree *t1 = nullptr;
    TTree *thead = nullptr;
    Int_t nentries = 0;
    vector<TBranch*> b;
    vector<ADC_DATA*> ch;
    Int_t nchannels = 1;
    vector<Int_t> channels = {1,2};
    vector<string> schannel;
    Double_t dtime = 2;
    string myname;
    string filename = "analyzed.root";
    TEventList *lev = nullptr;

    Int_t kch = 0;
    Int_t currentEvent = 0;
    DENOISE dn;
    WIENER *w = nullptr;
    string filter_type = "default";

    Int_t n_points;
    Bool_t invert = false;
    Bool_t put_my_offset_back = false;
    Double_t default_scale = 1;
    vector<vector<Double_t>> raw;
    vector<vector<Double_t>> wvf;
    vector<TH1D*> haverage;
    vector<TH1D*> hfft;
    TH1D *h = nullptr;
    Double_t *time = nullptr;
    Double_t *tempraw = nullptr;

    string plot_opt = "AL";
    TGraph *gwvf;
    string xlabel = "Time (ns)";
    string ylabel = "Amplitude (ADC Channels)";
    TH2D *hpers = nullptr;
    bool plotShowEvent = false;

    Double_t ymin = 0;
    Double_t ymax = 0;
    Double_t xmin = 0;
    Double_t xmax = 0;
    Double_t temp_charge = 0;
    Double_t temp_max = 0;
    Double_t temp_min = 0;
    Double_t temp_pos = 0;
    Double_t temp_pos_min = 0;
    Double_t temp_dummy = 0;

    Color_t mylinecolor = kBlack;
    Color_t mymarkercolor = kBlack;


    TCanvas *cpers = nullptr;


    struct Baselineparameters {
      Double_t baselineStart;
      Double_t baselineTime;
      string selection = "";
      Double_t filter = 0;
      Double_t exclusion_baseline = 0;
      Double_t exclusion_window = 1000;
      Double_t validfraction = 3.;
      Int_t selection_base = 0;
      Bool_t use_bits = true;
    };
    Baselineparameters baselineparameters;

    struct Header {
      Double_t dtime;
      Double_t startcharge;
      Double_t endcharge;
      Double_t baselineStart;
      Double_t baselineTime;
      Double_t maxRange;
      Double_t fast_time;
      Double_t slow_time;
      vector<Double_t> *exclusion_baselines = nullptr;
      Double_t exclusion_window;
      Double_t filter;
    };
    Header head;
    

    // _________________________ Methods related to the class maintenance _________________________ //

    // This allows to create a file, a tree and a branch outside the class
    // The reference type will allow us to change the pointer address
    void setAnalyzer(string m_filename = "analyzed.root"){
      filename = m_filename;
      if(f == nullptr) f = new TFile(filename.c_str(),"READ");
      if(t1 == nullptr) t1 = (TTree*)f->Get("t1");
      TList *lb = (TList*)t1->GetListOfBranches();

      if(lev==nullptr) lev = new TEventList(Form("lev_%s",myname.c_str()),Form("lev_%s",myname.c_str()));

      nchannels = lb->GetEntries();

      b.resize(nchannels);
      schannel.resize(nchannels);
      channels.resize(nchannels);
      ch.resize(nchannels);


      if(thead == nullptr)
      {
        bool has_head = false;
        for (auto fkey: *f->GetListOfKeys()){
          if (string(fkey->GetName()) == "head")
          {
            thead = (TTree*)f->Get("head");
            thead->SetBranchAddress("dtime", &head.dtime);
            thead->SetBranchAddress("startcharge", &head.startcharge);
            thead->SetBranchAddress("endcharge", &head.endcharge);
            thead->SetBranchAddress("baselineStart", &head.baselineStart);
            thead->SetBranchAddress("baselineTime", &head.baselineTime);
            thead->SetBranchAddress("maxRange", &head.maxRange);
            thead->SetBranchAddress("fast_time", &head.fast_time);
            thead->SetBranchAddress("slow_time", &head.slow_time);
            thead->SetBranchAddress("exclusion_baselines", &head.exclusion_baselines);
            thead->SetBranchAddress("exclusion_window", &head.exclusion_window);
            thead->SetBranchAddress("filter", &head.filter);
            thead->GetEntry(0);
            dtime = head.dtime;
            thead->ResetBranchAddresses();
            break;
          }
        }
      }

      for (Int_t i = 0; i < nchannels; i++) {
        schannel[i] = lb->At(i)->GetName();
        // Find the position of the first digit
        size_t digitPos = schannel[i].find_first_of("0123456789");

        // If a digit is found, extract the numeric part
        if (digitPos != std::string::npos) {
          // Use substr to get the numeric part from the string
          std::string numberString = schannel[i].substr(digitPos);

          // Convert the substring to an integer
          try {
            channels[i] =  std::stoi(numberString);
          } catch (const std::invalid_argument& e) {
            // Handle the case where conversion to integer fails
            std::cerr << "Error: " << e.what() << std::endl;
          }
        }
        ch[i] = new ADC_DATA();
        string temp_ch_name = schannel[i];
        b[i] = t1->GetBranch(temp_ch_name.c_str());
        b[i]->SetAddress(&ch[i]);
      }
      b[0]->GetEntry(0);
      n_points = ch[0]->npts;
      setEmpty();
      nentries = t1->GetEntries();
    }

    void setEmpty(){
      string wienername = myname + "_w";
      if(!w){
        w = new WIENER(wienername.c_str(),dtime,1./dtime/1e-9/1e6,1e-9,1e6,n_points);
      }
      raw.resize(nchannels);
      wvf.resize(nchannels);
      haverage.resize(nchannels);
      hfft.resize(nchannels);
      for(Int_t i = 0; i < nchannels; i++){
        raw[i].resize(n_points);
        wvf[i].resize(n_points);
        hfft[i] = nullptr;
      }
      time = new Double_t[n_points];
      tempraw = new Double_t[n_points];
      for (int j = 0; j < n_points; j++) {
        time[j] = j*dtime;
      }
      xmin = 0;
      xmax = n_points*dtime;
    }

    Int_t getIdx(){
      return channels[kch];
    }
    void print(){
      t1->Print();
    }
    Bool_t setChannel(string mych = "Ch1."){
      for(Int_t i = 0; i < nchannels; i++){
        if(mych == schannel[i])
        {
          kch = i;
          return true;
        }
      }

      printf("%s not found, run s.print() to check the branches\n",mych.c_str());
      return false;
    }


    void setAnalyzerExt(TFile *&ft, TTree *&tr, vector<TBranch *> &bt){
      f = ft;
      tr = (TTree*)ft->Get("t1");
      t1 = tr;
      bt.resize(nchannels);
      setAnalyzer();
      for (Int_t i = 0; i < nchannels; i++) bt[i] = b[i];
    }


    // _________________________ ________________________________________ _________________________ //


    // _________________________ Methods that get me something _________________________ //
    
    void getFFT(Double_t *_v = nullptr, bool inDecibel = false){
      if(_v == nullptr) _v = ch[kch]->wvf;
      for(Int_t i = 0; i < n_points; i++){
        w->hwvf->SetBinContent(i+1,_v[i]);
      }
      w->fft(w->hwvf);
      if(inDecibel) w->convertDecibel();
      h = w->hfft;
    }

    Double_t compute_true_baseline(Double_t baseline, Double_t offset){
      return (baseline/pow(2,14) - offset)*2000;
    }
    Double_t getMean(Double_t from, Double_t to){
      if (from == to) return 0;
      Double_t res = 0;
      Double_t totalev = 0;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        res += ch[kch]->wvf[i];
        totalev+=1;
      }
      return res/totalev;
    }

    Double_t getMaximum(Double_t from, Double_t to, Double_t *v = nullptr){
      if (!v) v = ch[kch]->wvf;
      Double_t max = -1e12;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        if(v[i]>max){
          max = v[i];
          temp_max = max;
          temp_pos = i*dtime;
        }
      }
      return max;
    }

    Double_t getMinimum(Double_t from, Double_t to, Double_t *v = nullptr){
      if (!v) v = ch[kch]->wvf;
      Double_t min = 1e12;
      if (to > n_points*dtime) to = n_points*dtime;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        if(v[i]<min){
          min = v[i];
          temp_min = min;
          temp_pos_min = i;
        }
      }
      return min;
    }

    void getMaxMin(Double_t from, Double_t to, Double_t *v = nullptr){
      if (!v) v = ch[kch]->wvf;
      Double_t max = -1e12;
      Double_t min = 1e12;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        if(v[i]>max){
          max = v[i];
          temp_max = max;
          temp_pos = i;
        }
        if(v[i]<min){
          min = v[i];
          temp_min = min;
          temp_pos_min = i;
        }
      }

    }

    Double_t getRMS(Double_t from, Double_t to, Double_t *v = nullptr){
      if (!v) v = ch[kch]->wvf;
      Double_t rms = 0;
      Double_t n = 0;
      for (Int_t i = from/dtime; i < to/dtime; i++) {
        rms += v[i]*v[i];
        n+=1;
      }
      return sqrt(rms/n);
    }

    void gen_RMS(Int_t channel = 0, Double_t from = 0, Double_t to = 0, string selection="", Double_t filter=0, TH1D *htemp = nullptr){
      getSelection(selection);
      kch = channel;
      htemp->GetYaxis()->SetTitle("# of events");
      htemp->GetXaxis()->SetTitle("RMS (ADCs)");
      if (to == 0) to = n_points*dtime;

      for(Int_t i = 0; i < lev->GetN(); i++){
        getWaveform(lev->GetEntry(i),kch);
        applyDenoise(filter);
        Double_t val = getRMS(from,to);
        htemp->Fill(val);
      }
    }


    void gen_rise_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}, Double_t filter = 0, TH1D *htemp = nullptr, string selection = ""){
      getSelection(selection);
      kch = channel;
      htemp->GetYaxis()->SetTitle("# of events");
      htemp->GetXaxis()->SetTitle("Rise time_{90} (ns)");

      for(Int_t i = 0; i < lev->GetN(); i++){
        getWaveform(lev->GetEntry(i),kch);
        applyDenoise(filter);
        Double_t val = rise_time(kch, baseline_range, ispulse, peak_range);
        htemp->Fill(val);
      }
    }

    void gen_fall_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}, Double_t filter = 0, TH1D *htemp = nullptr, string selection = ""){
      getSelection(selection);
      kch = channel;
      htemp->GetYaxis()->SetTitle("# of events");
      htemp->GetXaxis()->SetTitle("Fall time_{90} (ns)");

      for(Int_t i = 0; i < lev->GetN(); i++){
        getWaveform(lev->GetEntry(i),kch);
        applyDenoise(filter);
        Double_t val = fall_time(kch, baseline_range, ispulse, peak_range);
        htemp->Fill(val);
      }
    }
    Double_t linear_interpole_tot(Int_t i, Double_t th){
      if (i - 1 < 0) {
        return i;
      }
      Double_t t1 = i-1;
      Double_t t2 = i;
      Double_t A1 = ch[kch]->wvf[i-1];
      Double_t A2 = ch[kch]->wvf[i];

      return ((t2-t1)*th - (A1*t2- A2*t1))/(A2-A1);
    }

    Double_t getTOT(Int_t channel = 0, Double_t from = 0, Double_t to = 0, Double_t time_trigger = 20, Double_t minimum_gap_ns = 12){
      kch = channel;
      Bool_t triggered = false;
      Double_t time_mark = 0;
      for(Int_t i = from/dtime; i<to/dtime; i++){
        if(ch[channel]->wvf[i]>time_trigger && triggered == false){
          triggered = true;
          time_mark = linear_interpole_tot(i, time_trigger)*dtime;
          temp_pos = time_mark;
          // break;
        }
        else if(triggered==true && ch[channel]->wvf[i]<time_trigger && (i*dtime-time_mark)>=minimum_gap_ns){
          time_mark = linear_interpole_tot(i, time_trigger)*dtime - time_mark;
          return time_mark;
        }
      }
      return 0;
    }

    void integrate(Double_t from = 0, Double_t to = 0, Double_t percent = 0, Bool_t percent_withlimit = false, Double_t offset = 0, bool energy = false){
      Double_t res = 0;
      Double_t max = -1e12;
      if (to == 0) to = n_points*dtime;
      if(percent == 0){ // normal integration
        for(Int_t i = from/dtime; i < to/dtime; i++){
          if (!energy)
            res += ch[kch]->wvf[i]+offset;
          else
            res += ch[kch]->wvf[i]*ch[kch]->wvf[i]+offset;
          if(ch[kch]->wvf[i]>=max){
            max = ch[kch]->wvf[i];
            temp_pos = i*dtime;
          }

        }
      }
      else{
        max = getMaximum(from, to);
        Int_t imin = 0;
        Int_t imax = n_points;
        if(percent_withlimit){
          imin = from/dtime;
          imax = to/dtime;
        }
        for(Int_t i = temp_pos/dtime; i >= imin ;i--){
          Double_t val = ch[kch]->wvf[i];
          if(val >= percent*max){
            res += val+offset;
          }
          else{
            break;
          }
        }
        for(Int_t i = temp_pos/dtime+1; i < imax ;i++){
          Double_t val = ch[kch]->wvf[i];
          if(val >= percent*max){
            res += val+offset;
          }
          else{
            break;
          }
        }
      }
      temp_charge = res*dtime;
      temp_max = max;
    }
    void getIntegral(TH1D *_htemp = nullptr, Double_t from = 0, Double_t to = 0, string selection = "", Double_t filter = 0, Double_t sphe = 1., Double_t percent = 0){
      getSelection(selection);
      Int_t nev = lev->GetN();
      Int_t iev = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        integrate(from, to, percent);
        _htemp->Fill(temp_charge/sphe);
      }
    }


    void cumulative_integral(vector<Double_t> &output, Double_t from = 0, Double_t to = 0){
      Double_t res = 0;
      Double_t max = -1e12;
      Double_t maxcharge = -1e12;
      Int_t idxmaxcharge = 0;
      Double_t pos = 0;
      if (to == 0) to = n_points*dtime;
      output.resize(n_points);
      // Assigns value 0 to all the elements in the vector
      std::fill(output.begin(), output.end(), 0);
      for(Int_t i = from/dtime; i < to/dtime; i++){
        res += ch[kch]->wvf[i]*dtime;
        output[i] = res;
        if(ch[kch]->wvf[i]>=max){
          max = ch[kch]->wvf[i];
          pos = i*dtime;
        }
        if(res >= maxcharge){
          maxcharge = res;
          idxmaxcharge = i;
        }
      }
      temp_charge = res;
      temp_max = max;
      temp_pos = pos;

      for (Int_t i = 0; i < (Int_t)output.size(); i++){
        output[i] = output[i]/output[idxmaxcharge];
      }

    }

    void exportAsHistogram(TH1D *hexport = nullptr, Double_t x_start = 0){

      if (hexport == nullptr){
        hexport = new TH1D("","",n_points, x_start-dtime/2, x_start+n_points*dtime - dtime/2);
        h = hexport;
      }
      for(Int_t i = 0; i < n_points; i++){
        hexport->SetBinContent(i+1, ch[kch]->wvf[i]);
      }
    }

    void getWaveFromHistogram(TH1D *htemp){
      if (htemp->GetNbinsX() != n_points){
        cout << "Not same amount of samples! Graph has " << htemp->GetNbinsX() << endl;
        return;
      }
      for(Int_t i = 0; i < htemp->GetNbinsX(); i++){
        ch[kch]->wvf[i] = htemp->GetBinContent(i+1);
      }
    }
    void getWaveFromGraph(TGraph *gtemp){
      Double_t *xtemp = nullptr;
      Int_t ntemp = gtemp->GetN();
      if (ntemp != n_points){
        cout << "Not same amount of samples! Graph has " << ntemp << endl;
        return;
      }
      for(Int_t i = 0; i < ntemp; i++){
        ch[kch]->wvf[i] = *(gtemp->GetY()+i);
      }
    }


    void getWaveform(Int_t myevent = 0, Int_t k = -1){
      if (k == -1) k = kch;
      if (k>=nchannels){
        cout << "There are only " << nchannels << " in the TTree, execute print() to check channels" << endl;
        return;
      }
      kch = k;
      b[k]->GetEvent(myevent);
      n_points = ch[k]->npts;
      currentEvent = ch[k]->event;
      if (default_scale!=1) scaleWvf(default_scale);

      if (!invert && !put_my_offset_back) return;
      else if (invert && put_my_offset_back) addOffsetWithScale(0,0,0,-1);
      else if (invert) scaleWvf(-1);
      else if (put_my_offset_back) addOffet();


    }

    bool getWaveformHard(Int_t myevent = 0, Double_t factor = 1){
      if (kch>=nchannels){
        cout << "There are only " << nchannels << " in the TTree, execute print() to check channels" << endl;
        return false;
      }
      b[kch]->GetEvent(myevent);
      n_points = ch[kch]->npts;
      currentEvent = ch[kch]->event;
      if (!invert && !put_my_offset_back){}
      else if (invert && put_my_offset_back) addOffsetWithScale(0,0,0,-1);
      else if (invert) scaleWvf(-1);
      else if (put_my_offset_back) addOffet();

      for (int j = 0; j < n_points; j++) {
        raw[kch][j] = ch[kch]->wvf[j]*factor;
        ch[kch]->wvf[j] = ch[kch]->wvf[j]*factor;
        wvf[kch][j] = raw[kch][j];
        time[j] = j*dtime;
      }
      return true;

    }

    void draw_rise_lines(Double_t start, Double_t finish, Double_t baseline_level, Double_t peak_level){
      Double_t gxmin = gPad->GetUxmin();
      Double_t gxmax = gPad->GetUxmax();
      TLine *lbase = new TLine(gxmin, baseline_level, gxmax, baseline_level);
      TLine *lpeak = new TLine(gxmin, peak_level, gxmax, peak_level);
      TLine *l90 = new TLine(gxmin, 0.9*(peak_level-baseline_level)+baseline_level, gxmax, 0.9*(peak_level-baseline_level)+baseline_level);
      TLine *l10 = new TLine(gxmin, 0.1*(peak_level-baseline_level)+baseline_level, gxmax, 0.1*(peak_level-baseline_level)+baseline_level);

      TLine *lstart = new TLine(start, baseline_level, start, peak_level);
      TLine *lfinish = new TLine(finish, baseline_level, finish, peak_level);

      set_up_lines(lbase,kGreen+2);
      set_up_lines(lpeak,kGreen+2);
      set_up_lines(l90, kRed);
      set_up_lines(l10, kRed);
      set_up_lines(lstart, kRed);
      set_up_lines(lfinish, kRed);

    }

    Double_t triggerTime(vector<Double_t> baseline_range = {0,0}, vector<Double_t> peak_range = {0,0}, Double_t threshold = 0, bool debug = false){
      Double_t baseline_level = getMean(baseline_range[0],baseline_range[1]);

      Int_t startpoint = peak_range[0]/dtime;
      Int_t endpoint = peak_range[1]/dtime;

      Double_t time_mark;
      for(Int_t i = startpoint; i < endpoint; i++){
        if (ch[kch]->wvf[i] >= threshold) {
          if(i == startpoint) break;
          time_mark = linear_interpole_tot(i, threshold - baseline_level)*dtime;
          temp_pos = time_mark;
          return time_mark;
        }
      }
      temp_pos = 0;
      return 0;
    }


    Double_t rise_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}, bool debug = false){
      kch = channel;
      Double_t baseline_level = getMean(baseline_range[0],baseline_range[1]);
      temp_dummy = baseline_level;

      Double_t peak_level = 0;
      if (ispulse){
        peak_level = getMaximum(peak_range[0],peak_range[1]);
      }
      else{
        peak_level = getMean(peak_range[0],peak_range[1]);
      }
      Int_t startpoint = peak_range[0]/dtime;
      Int_t endpoint = peak_range[1]/dtime;
      Bool_t triggered = false;
      Double_t time_mark = 0;

      for(Int_t i = startpoint; i < endpoint; i++){
        if (ch[kch]->wvf[i] >= 0.1*(peak_level-baseline_level)+baseline_level && triggered==false) {
          triggered = true;
          time_mark = linear_interpole_tot(i, 0.1*(peak_level-baseline_level)+baseline_level)*dtime;
        }
        else if(triggered == true && ch[kch]->wvf[i]>=0.9*(peak_level - baseline_level)+baseline_level){
          if(debug) draw_rise_lines(time_mark, i*dtime, baseline_level, peak_level);
          temp_pos = linear_interpole_tot(i, 0.9*(peak_level - baseline_level)+baseline_level)*dtime;
          // time_mark = i*dtime - time_mark;
          temp_max = peak_level;
          return temp_pos - time_mark;
        }
      }
      return 0;
    }

    Double_t fall_time(Int_t channel = 0, vector<Double_t> baseline_range = {0,0}, Bool_t ispulse = true, vector<Double_t> peak_range = {0,0}, bool debug = false){
      kch = channel;
      Double_t baseline_level = getMean(baseline_range[0],baseline_range[1]);
      temp_dummy = baseline_level;

      Double_t peak_level = 0;
      if (ispulse){
        peak_level = getMaximum(peak_range[0],peak_range[1]);
      }
      else{
        peak_level = getMean(peak_range[0],peak_range[1]);
      }
      Int_t startpoint = temp_pos/dtime;
      Int_t endpoint = peak_range[1]/dtime;
      Bool_t triggered = false;
      Double_t time_mark = 0;

      for(Int_t i = startpoint; i < endpoint; i++){
        if (ch[kch]->wvf[i] <= 0.9*(peak_level-baseline_level)+baseline_level && triggered==false) {
          triggered = true;
          time_mark = linear_interpole_tot(i, 0.9*(peak_level-baseline_level)+baseline_level)*dtime;
          // time_mark = i*dtime;
        }
        else if(triggered == true && ch[kch]->wvf[i]<=0.1*(peak_level - baseline_level)+baseline_level){
          if(debug) draw_rise_lines(time_mark, i*dtime, baseline_level, peak_level);
          temp_pos = linear_interpole_tot(i, 0.1*(peak_level - baseline_level)+baseline_level)*dtime;
          // time_mark = i*dtime - time_mark;
          temp_max = peak_level;
          return temp_pos - time_mark;
        }
      }
      return 0;
    }

    void getBackFFT(Double_t *_filtered = nullptr){
      if (_filtered == nullptr){_filtered = ch[kch]->wvf;}
      w->backfft(*w);
      for(Int_t i = 0; i < n_points; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }

    }
    void averageFFT(Int_t maxevent = 0, Int_t startevent = 0, string selection = "", bool inDecibel = true, Double_t factor = 0, Double_t filter = 0){
      getSelection(selection);
      Int_t nev = lev->GetN();
      if (maxevent==0) {
        maxevent = nev;
      }
      if (startevent+maxevent <= nev) {
        nev = startevent+maxevent;
      }
      else{
        cout << "Maximum number of events is " << nev << endl;
        return;
      }
      Int_t iev = 0;
      hfft[kch] = (TH1D*)w->hfft->Clone(Form("hfft_%s_ch%d",myname.c_str(),kch));
      hfft[kch]->Reset();
      cout << "\n";
      Int_t total = 0;
      for(Int_t i = startevent; i < nev; i++){
        iev = lev->GetEntry(i);
        printev(i,nev);
        getWaveform(iev,kch);
        applyDenoise(filter);
        // applyFreqFilter();
        getFFT();
        for (Int_t j = 0; j < n_points/2; j++) hfft[kch]->AddBinContent(j+1,w->hfft->GetBinContent(j+1));
        total++;
      }
      cout << "\n";
      hfft[kch]->Scale(1./total);
      hfft[kch]->SetEntries(total);

      if (inDecibel){
        w->convertDecibel(hfft[kch], factor);
        hfft[kch]->GetYaxis()->SetTitle("Magnitude (dB)");
      }
      h = hfft[kch];
    }

    void averageWaveform(Int_t maxevent = 0, string selection = "", Double_t filter =  0){
      if (maxevent==0) {
        maxevent = nentries;
      }
      if(!haverage[kch]) haverage[kch] = new TH1D(Form("haverage_%s_Ch%d",myname.c_str(),channels[kch]),"Averaged waveform",n_points,0,n_points*dtime);
      else haverage[kch]->Reset();
      haverage[kch]->GetYaxis()->SetTitle("Amplitude (ADC Channels)");
      haverage[kch]->GetXaxis()->SetTitle("Time (ns)");
      Int_t total = 0;
      getSelection(selection);
      Int_t nev = lev->GetN();
      if (maxevent < nev) {
        nev = maxevent;
      }
      Int_t iev = 0;
      for(Int_t i = 0; i < nev; i++){
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        total += 1;
        for (Int_t j = 0; j < n_points; j++){
          haverage[kch]->AddBinContent(j+1,ch[kch]->wvf[j]);
        }

      }
      haverage[kch]->Scale(1./total);
      haverage[kch]->SetEntries(total);
      haverage[kch]->Sumw2(kFALSE);
      h = haverage[kch];
    }

    void zeroCrossSearch(Double_t *derwvf, vector<Int_t> &peaksRise, vector<Int_t> &peaksCross, Int_t start, Int_t finish){
      if(start == 0) start = dtime;
      bool lastIsPositive = false; // I want to always start with a positive crossing
      for(Int_t i = start/dtime; i < finish/dtime-1; i++){
        if(lastIsPositive==false && derwvf[i] >= 0 && derwvf[i-1] <=0 && derwvf[i+1]>0){
          peaksRise.push_back(i);
          lastIsPositive = true;
          i = i+1;
        }
        else if(lastIsPositive && derwvf[i] <= 0 && derwvf[i-1] >=0 && derwvf[i+1] < 0){
          peaksCross.push_back(i);
          lastIsPositive = false;
        }
      }
    }


    void convert_1D_histogram(TH2D *horiginal, TH1D **hnew, vector<vector<Double_t>> &_hy, vector<vector<Double_t>> &_hx, Int_t rebin = 1){
      Int_t nx = horiginal->GetNbinsX();
      Int_t ny = horiginal->GetNbinsY();
      Double_t hmin = horiginal->GetXaxis()->GetBinLowEdge(1);
      Double_t hmax = horiginal->GetXaxis()->GetBinLowEdge(nx) + horiginal->GetXaxis()->GetBinWidth(nx);

      *hnew = new TH1D("hnew","hnew",nx/rebin,hmin,hmax);
      Double_t sum = 0;
      Double_t div = 0;
      Double_t stddev = 0;
      Double_t n = 0;
      Int_t nrebins = 0;
      int ndiv = 0;
      for(Int_t i = 0; i < nx; i++) { // for loop in the X axis of the histogram;
        for(Int_t j = 0; j < ny; j++) { // for loop in the Y axis of the histogram;
          Double_t weight = horiginal->GetBinContent(i+1,j+1);
          if (weight != 0){
            Double_t bnpes = horiginal->GetYaxis()->GetBinCenter(j);
            sum += weight*bnpes;
            div += weight;
          }
        }
        ndiv++;
        if(ndiv == rebin){
          if(div == 0){sum = 0; div = 1;}
          for(Int_t j = 0; j < ny; j++) { // for loop in the Y axis of the histogram;
            Double_t weight = horiginal->GetBinContent(i+1,j+1);
            if (weight != 0){
              Double_t bnpes = horiginal->GetYaxis()->GetBinCenter(j);
              stddev += weight*(bnpes-sum/div)*(bnpes-sum/div);
              n += 1;
            }
          }
          if (n == 1){ n = 2;}
          stddev = (stddev/div)*n/(n-1);
          (*hnew)->SetBinContent(nrebins+1, sum/div);

          Double_t err = sqrt(abs(sum)/div);
          (*hnew)->SetBinError(nrebins+1,sqrt(stddev + err*err));
          if(sum!=0){
            _hy[0].push_back(sum/div);
            _hy[1].push_back((*hnew)->GetBinError(nrebins+1));
            _hx[0].push_back((*hnew)->GetBinCenter(nrebins+1));
            _hx[1].push_back((*hnew)->GetBinWidth(nrebins+1)/2.);
          }

          nrebins++;

          sum = 0;
          div = 0;
          stddev = 0;
          n = 0;
          ndiv = 0;
        }


      }

    }


    void getSPEIntervals(Double_t clockstart,
                         Double_t clockperiod,
                         Double_t integration_time,
                         string selection = "",
                         Double_t filter = 0,
                         Double_t percent = 0){
      getSelection(selection);
      Int_t nev = lev->GetN();
      Int_t iev = 0;
      Double_t from = 0;
      Double_t to = 0;
      Bool_t correct_baseline = false;

      if (_baseIntervals.size()!=0){
        correct_baseline = true;
        if ((int)_baseIntervals.size() != nev){
          cout << "Baseline correction not possible, different size of events... " << endl;
        }
      }
      for(Int_t i = 0; i < nev; i++){
        printev(i, nev);
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        from = clockstart;
        to = clockstart + integration_time;
        int count = 0;
        Double_t newBase = 0;
        while(from < n_points*dtime && to < n_points*dtime && to < xmax){
          if (correct_baseline){
            newBase = _baseIntervals[i][count];
          }
          integrate(from, to, percent, false, -newBase);
          if (newBase != -999){
            h->Fill(temp_charge - newBase);
          }
          from += clockperiod;
          to = from + integration_time;
          count+=1;
        }
      }
      cout << "\n";
      TFile *fout = new TFile(Form("sphe_histograms_Ch%i.root",channels[kch]),"RECREATE");
      fout->WriteObject(h,"analyzed","TObject::kOverwrite");
    }



    void getBaselineIntervals(Double_t clockstart,
                              Double_t clockperiod,
                              Double_t baseline_time,
                              string selection = "",
                              Double_t filter = 0,
                              Double_t exclusion_baseline = 0,
                              Double_t exclusion_window = 0,
                              Double_t validfraction = 3.,
                              TH1D *hbase = nullptr){
      getSelection(selection);
      Int_t nev = lev->GetN();
      Int_t iev = 0;
      Double_t from = 0;
      Double_t to = 0;
      _baseIntervals.clear();
      _baseIntervals.resize(nev);
      for(Int_t i = 0; i < nev; i++){
        printev(i, nev);
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        from = clockstart;
        to = clockstart + baseline_time;
        while(from < n_points*dtime && to < n_points*dtime && to < xmax){
          baselineparameters.baselineStart = from;
          baselineparameters.baselineTime = to;
          baselineparameters.exclusion_baseline = exclusion_baseline;
          baselineparameters.exclusion_window = exclusion_window;
          baselineparameters.validfraction = validfraction;
          Double_t base = compute_baseline(hbase);
          from += clockperiod;
          to = from + baseline_time;
          _baseIntervals[i].push_back(base);
        }
      }
    }

    void extract_waveforms(Double_t from,
                           Double_t to,
                           string selection = "",
                           Bool_t correct_baseline = true,
                           string filename="sample.root",
                           Double_t clockstart = 0,
                           Double_t clockperiod = 0,
                           TH1D *hbase = nullptr
                           ){

      TFile *fout = new TFile(filename.c_str(), "READ");
      if (!fout->IsZombie()){
        cout << "File " << filename << " already exists!" << endl;
        cout << "Do you want to overwrite it? (yes/no): " << endl;
        fout->Close();

        // Prompt user for overwrite confirmation
        std::string answer;
        std::cin >> answer;

        if (answer != "yes") {
            std::cout << "Exiting without overwriting the file." << std::endl;
            return;
        }
      }
      ADC_DATA newch;

      gROOT->cd();
      fout = new TFile(filename.c_str(), "RECREATE");
      TTree *tout = new TTree("t1", "ADC processed waveform");
      // tout->SetAutoSave(false);
      TBranch *newb = tout->Branch(Form("Ch%d.",this->getIdx()), &newch);

      TTree *theadout = new TTree("head","header");
      theadout->Branch("dtime", &head.dtime);
      theadout->Branch("startcharge", &head.startcharge);
      theadout->Branch("endcharge", &head.endcharge);
      theadout->Branch("baselineStart", &head.baselineStart);
      theadout->Branch("baselineTime", &head.baselineTime);
      theadout->Branch("maxRange", &head.maxRange);
      theadout->Branch("fast_time", &head.fast_time);
      theadout->Branch("slow_time", &head.slow_time);
      theadout->Branch("exclusion_baselines", &head.exclusion_baselines);
      theadout->Branch("exclusion_window", &head.exclusion_window);
      theadout->Branch("filter", &head.filter);
      theadout->Fill();

      if (to/dtime == this->n_points){
        cout << "Limit reached " << "changing to " << this->n_points*dtime << endl;
        to = this->n_points*dtime-1;
      }
      Int_t npts = (to/dtime - from/dtime);
      // cout << npts << " " << to << " " << from << " " << dtime << " " << to/dtime << " " << from/dtime << endl;
      newch.Set_npts(npts);

      getSelection(selection);
      Int_t nev = lev->GetN();
      Int_t iev = 0;
      from += clockstart;
      to += clockstart;
      baselineparameters.baselineStart += clockstart;
      baselineparameters.baselineTime += clockstart;
      Double_t fromstart = from;
      Double_t tostart = to;
      Double_t baseStartstart = baselineparameters.baselineStart;
      Double_t baseTimestart = baselineparameters.baselineTime;

      Double_t filter = (correct_baseline) ? baselineparameters.filter : 0;
      Double_t newbase = 0;
      for(Int_t i = 0; i < nev; i++){
        printev(i, nev);
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        from = fromstart;
        to = tostart;
        baselineparameters.baselineStart = baseStartstart;
        baselineparameters.baselineTime = baseTimestart;
        while(from < n_points*dtime && to < n_points*dtime && to < xmax){
          if (correct_baseline){
            newbase = compute_baseline(hbase);
            getWaveform(iev,kch);
          }
          for (Int_t k = 0; k < npts; k++){
              newch.wvf[k] = this->ch[kch]->wvf[(int)(k + from/dtime)] - newbase;
            // aux++;
          }
          // cout << aux << " " << npts << " " << to << " " << from << " " << dtime << " " << to/dtime << " " << from/dtime << endl;
          newch.peak = this->ch[kch]->peak;
          newch.peakpos = this->ch[kch]->peakpos;
          newch.charge = this->ch[kch]->charge;
          newch.fprompt = this->ch[kch]->fprompt;
          newch.event = this->ch[kch]->event;
          newch.time = this->ch[kch]->time;
          newch.base = this->ch[kch]->base;
          newch.selection = baselineparameters.selection_base;

          tout->Fill();

          if (clockperiod != 0){
            from += clockperiod;
            to += clockperiod;
            baselineparameters.baselineStart += clockperiod;
            baselineparameters.baselineTime += clockperiod;
            if (correct_baseline  && (baselineparameters.baselineStart >= n_points*dtime || baselineparameters.baselineTime >= n_points*dtime)){
              break;
            }
          }
          else{
            break;
          }
        }

      }

      cout << "\n\n";
      fout->WriteTObject(tout, "t1");
      fout->WriteTObject(theadout, "head");
      delete newb;
      delete tout;
      fout->Close();


    }
    void create_persistence(Int_t nbins, Double_t ymin, Double_t ymax, Double_t filter, string cut, Double_t factor){

      Int_t nbinsx = (xmax-xmin)/dtime;
      if(!hpers) hpers = new TH2D("hpers","hpers",nbinsx,xmin,xmax,nbins,ymin,ymax);
      else{
        hpers->Reset();
        hpers->SetBins(nbinsx, xmin, xmax, nbins, ymin, ymax);
      }
      if(!cpers) cpers = new TCanvas(Form("cpers_%s", myname.c_str()), Form("cpers_%s", myname.c_str()),1920,0,1920,1080);
      else{cpers->cd();}

      getSelection(cut);
      Int_t nev = lev->GetN();
      Int_t iev = 0;

      for(Int_t i = 0; i < nev; i++){

        printev(i,nev);
        iev = lev->GetEntry(i);
        getWaveform(iev,kch);
        applyDenoise(filter);
        // applyFreqFilter();

        for (int j = 0; j < n_points; j++) {
          hpers->Fill(j*dtime,ch[kch]->wvf[j]*factor);
        }

      }
      cout << "\n";
    }
    // _________________________ _____________________________ _________________________ //


    // _________________________ Filters and processing _________________________ //
    void scaleWvf(Double_t factor, Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch]->wvf;
      for (Int_t i = 0; i < n_points; i++) {
        _filtered[i] = _filtered[i]*factor;
      }
    }

    void addOffet(Double_t offset = 0, Double_t from = 0, Double_t to = 0){
      if(offset == 0){
        offset = ch[kch]->base;
      }
      if(to == 0) to = n_points*dtime;
      for(Int_t i = from/dtime; i < to/dtime; i++){
        ch[kch]->wvf[i] = ch[kch]->wvf[i] + offset;
      }
    }

    void addOffsetWithScale(Double_t offset = 0, Double_t from = 0, Double_t to = 0, Double_t factor = 1){
      if(offset == 0){
        offset = ch[kch]->base;
      }
      if(to == 0) to = n_points*dtime;
      for(Int_t i = from/dtime; i < to/dtime; i++){
        ch[kch]->wvf[i] = ch[kch]->wvf[i]*factor + offset;
      }
    }
    void checkSignals(Double_t **_raw, Double_t **_filtered){
      if (*_raw == nullptr  && *_filtered == nullptr){
        *_filtered = ch[kch]->wvf;
        *_raw = tempraw;
          for (Int_t i = 0; i < n_points; i++) {
            (*_raw)[i] = ch[kch]->wvf[i];
          }
      }
      else if(*_raw == *_filtered){
        for (Int_t i = 0; i < n_points; i++) {
          tempraw[i] = (*_raw)[i];
        }
        *_raw = tempraw;
      }
    }
    void differenciate(Double_t factor = 1, Double_t *_raw = nullptr, Double_t *_shifted = nullptr){
      checkSignals(&_raw,&_shifted);
      _shifted[0] = 0;
      _shifted[n_points-1] = 0;
      for(int i=1; i<n_points-1; i++){
        _shifted[i] = factor*(_raw[i+1] - _raw[i-1])/(2*dtime);
        // _shifted[i]=_raw[i] - (i-delay_time>=0 ? _raw[i-delay_time] : 0);
      }
    }

    void applyMovingAverage(Int_t mafilter = 0, Double_t start = 0, Double_t finish = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      if (finish == 0) finish = n_points*dtime;
      if(mafilter!=0) {
        checkSignals(&_raw,&_filtered);
        dn.movingAverage(_raw,_filtered,mafilter,start/dtime,finish/dtime);
      }
    }

    void applyCumulativeAverage(Int_t mafilter = 0, Double_t start = 0, Double_t finish = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      if (finish == 0) finish = n_points*dtime;
      if(mafilter!=0) {
        checkSignals(&_raw,&_filtered);
        dn.cumulativeAverage(_raw,_filtered,mafilter,start/dtime,finish/dtime);
      }
    }

    void applyDenoise(Double_t filter = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      checkSignals(&_raw,&_filtered);
      if (filter == 0 && _filtered == ch[kch]->wvf) return;
      if (filter_type == "default"){
        applyTV1D(filter, _raw, _filtered);
      }
      else if (filter_type == "ma"){
        applyMovingAverage(filter);
      }
      else if(filter_type == "band")
        applyBandCut();
      else{
        applyFreqFilter();
      }
    }

    void applyTV1D(Double_t filter = 0, Double_t *_raw = nullptr, Double_t *_filtered = nullptr){
      checkSignals(&_raw,&_filtered);
      if (filter == 0 && _filtered == ch[kch]->wvf) return;
      dn.TV1D_denoise(_raw,_filtered,n_points,filter);
    }

    void setFreqFilter(Double_t frequency_cut, string filter_type = "gaus"){
      this->filter_type = filter_type;
      w->setFilter(frequency_cut,filter_type);
    }
    void applyFreqFilter(Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch]->wvf;
      for(Int_t i = 0; i < n_points; i++){
        w->hwvf->SetBinContent(i+1,_filtered[i]);
      }
      w->fft(w->hwvf);
      w->apply_filter();
      w->backfft(*w);
      for(Int_t i = 0; i < n_points; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }
      // w->hwvf->Draw("");
    }


    void applyBandCut(WIENER *_temp = nullptr, Double_t *_filtered = nullptr){
      if(_filtered == nullptr) _filtered = ch[kch]->wvf;
      if(_temp == nullptr) _temp = w;
      for(Int_t i = 0; i < n_points; i++){
        w->hwvf->SetBinContent(i+1,_filtered[i]);
      }
      w->fft(w->hwvf);
      w->applyBandCut(_temp);
      w->backfft(*w);
      for(Int_t i = 0; i < n_points; i++){
        _filtered[i] = w->hwvf->GetBinContent(i+1);
      }
    }

    void makeCopy(Double_t *cpy, Double_t *original = nullptr){
      if (original == nullptr) original = ch[kch]->wvf;
      for (Int_t i = 0; i < n_points; i++) {
        cpy[i] = original[i];
      }
    }

    void shift_waveform(TH1D *h, Int_t new_max, Bool_t rawShift = false){
      Int_t old_max = h->GetMaximumBin();
      if(rawShift) old_max = 0;
      Int_t old_ref = old_max - new_max;
      TH1D *htemp = (TH1D*)h->Clone("htemp");
      Double_t temp;
      if(old_ref<0){
        // cout << " case lower" << endl;
        old_ref = n_points-(new_max-old_max);
      }
      for(Int_t i = 1; i<n_points-(old_ref); i++){
        temp = htemp->GetBinContent(old_ref+i);
        h->SetBinContent(i,temp);
      }
      Int_t aux = 1;
      for(Int_t i = n_points-(old_ref); i<=n_points; i++){
        temp = htemp->GetBinContent(aux);
        h->SetBinContent(i,temp);
        aux++;
      }
      delete htemp;
    }


    Double_t compute_baseline(TH1D *hbase = nullptr){

      vector<Double_t> range_base = {baselineparameters.baselineStart, baselineparameters.baselineTime};
      Double_t exclusion_baseline = baselineparameters.exclusion_baseline;
      Double_t exclusion_window = baselineparameters.exclusion_window;
      Double_t validfraction = baselineparameters.validfraction;

      
      Double_t result = 0;
      if(!hbase){
        if(!h) {
          getMaxMin(range_base[0], range_base[1]);
          
          Int_t min = (int)temp_min;
          Int_t max = (int)temp_max;
          Int_t nbins = (max-min);
          if (baselineparameters.use_bits)
            h = new TH1D("hbase","finding baseline",TMath::Power(2,14),0,TMath::Power(2,14));
          else
            h = new TH1D("hbase","finding baseline",nbins, min+0.5, max+0.5);
        }
        hbase = h;
      }
      hbase->Reset();
      // cout << "Computing from " << range_base[0] << " to " << range_base[1] << " " << exclusion_baseline << " " << exclusion_window << endl;
      for(Int_t i=range_base[0]/dtime; i<range_base[1]/dtime; i++) hbase->Fill(ch[kch]->wvf[i]);
      Double_t res0 = hbase->GetBinCenter(hbase->GetMaximumBin());
      Double_t hmean = hbase->GetMean();
      Double_t hstd = hbase->GetStdDev();
      Double_t bins=0;
      // cout << "Maximum bin: " << res0 << endl;
      for(Int_t i=range_base[0]/dtime; i<range_base[1]/dtime;){
        if(ch[kch]->wvf[i] > res0 + exclusion_baseline || ch[kch]->wvf[i]<res0 - exclusion_baseline) {
          // cout << "Too big... wvf[" << i << "] = " << ch[kch]->wvf[i] << endl;
          i+=exclusion_window/dtime;
        }
        else{
          result += ch[kch]->wvf[i];
          bins+=1;
          i++;
        }
      }
      if(bins>0)result/=bins;
      if(bins > ((range_base[1]-range_base[0])/dtime)/validfraction){
        baselineparameters.selection_base = 0;
      }
      else{
        // result = res0;
        baselineparameters.selection_base = 1;
        // cout << "Not enough points probably... valid points = " << bins << endl;
      }

      delete h;
      h = nullptr;
      return result;
    }

    Double_t reval_baseline(TH1D *hbase = nullptr){
      Double_t result = compute_baseline(hbase);

      addOffet(-result);
      return result;
    }

    void computeSlope(Int_t gap = 1, Double_t *_raw = nullptr, Double_t *_out = nullptr){
      checkSignals(&_raw,&_out);
      for (Int_t i = gap; i < n_points; i++){
        _out[i] = (_raw[i]-_raw[i-gap])/(gap);
      }
    }

    // _________________________ ______________________ _________________________ //


    // _________________________ Methods for selection _________________________ //
    bool checkHigher(Double_t a, Double_t b){
      if (a > b){return true;}
      else{return false;}
    }
    bool checkLower(Double_t a, Double_t b){
      if (a < b){return true;}
      else{return false;}
    }

    void invertSelection(){
      TEventList *ttemp = new TEventList("ttemp", "ttemp");
      t1->Draw(">>ttemp","");
      ttemp->Subtract(lev);
      lev->Reset();
      for (Int_t i = 0; i < ttemp->GetN(); i++) {
        lev->Enter(ttemp->GetEntry(i));
      }
      delete ttemp;
    }

    void excludeByAmplitude(Double_t filter = 0, Double_t xmin = 0, Double_t xmax = 0, Double_t limit = 100, string type = "higher"){
      f->cd();
      lev = (TEventList*)gDirectory->Get(Form("lev_%s",myname.c_str()));
      if(lev->GetN() == 0){
        getSelection("");
      }
      TEventList *ttemp = new TEventList("ttemp", "ttemp");

      if (xmax == 0 || xmax>n_points*dtime) {
        xmax = n_points*dtime;
      }
      for(Int_t i = 0; i < lev->GetN(); i++){
        Int_t ev = lev->GetEntry(i);
        getWaveform(ev, kch);
        applyDenoise(filter);
        bool contact = false;
        for(Int_t j = xmin/dtime; j < xmax/dtime; j++){
          if(type == "higher") contact = checkHigher(ch[kch]->wvf[j], limit);
          else contact = checkLower(ch[kch]->wvf[j], limit);
          if(contact){
            ttemp->Enter(ev);
            break;
          }
        }
      }
      lev->Subtract(ttemp);
      delete ttemp;
    }

    void genSelectByAmplitude(Double_t filter = 0, Double_t xmin = 0, Double_t xmax = 0, Double_t limit = 100, string type = "higher"){
      f->cd();
      lev = (TEventList*)gDirectory->Get(Form("lev_%s",myname.c_str()));
      lev->Reset();
      TEventList *ttemp = new TEventList("ttemp", "ttemp");

      if (xmax == 0 || xmax>n_points*dtime) {
        xmax = n_points*dtime;
      }
      for(Int_t i = 0; i < nentries; i++){
        // Int_t ev = lev->GetEntry(i);
        Int_t ev = i;
        printev(ev, nentries);
        getWaveform(ev, kch);
        applyDenoise(filter);
        bool contact = false;
        for(Int_t j = xmin/dtime; j < xmax/dtime; j++){
          if(type == "higher") contact = checkHigher(ch[kch]->wvf[j], limit);
          else contact = checkLower(ch[kch]->wvf[j], limit);
          if(contact){
            ttemp->Enter(ev);
            break;
          }
        }
      }
      cout << "\n";
      cout << ttemp->GetN() << " events found. " << endl;
      lev->Add(ttemp);
      delete ttemp;
    }


    void getSelection(string selection)
    {
      f->cd();
      if(selection!="use_mine"){
        lev->Resize(nentries);
        t1->Draw(Form(">>lev_%s",myname.c_str()),selection.c_str());
      }
      else{
        if (lev->GetN() == 0){
          cout << "No event selectected... " << endl;
        }
      }
    }
    // _________________________ _____________________ _________________________ //

    void set_up_lines(TLine *l, Color_t color){
      l->SetLineColorAlpha(color,0.6);
      l->SetLineStyle(9);
      l->Draw("SAME");

    }

    Double_t computeSNR_simple(Double_t xmin, Double_t xmax, vector<Double_t> signal_range){
      Double_t snr = 0;
      Double_t avg = 0;
      Double_t avgsqr = 0;
      Double_t sum = 0;
      Double_t sumsqr = 0;
      Double_t navg = 0;
      Double_t stddev = 0;

      for(Int_t i = xmin/dtime; i < xmax/dtime; i++){
        if(i*dtime > signal_range[1] || i*dtime < signal_range[0]){
          sum += ch[kch]->wvf[i];
          sumsqr += ch[kch]->wvf[i]*ch[kch]->wvf[i];
          navg += 1;
        }
      }
      avg = sum/navg;
      avgsqr = sumsqr/navg;
      stddev = sqrt(avgsqr - avg*avg);
      Double_t max = getMaximum(signal_range[0], signal_range[1]);
      snr = max/stddev;
      return snr;
    }


    void printev(Int_t i, Int_t nev, Int_t reference=200){
      if(nev<2*reference){
        reference = 20;
      }
      if (reference==1){
        cout << i << " out of " << nev << "\r" << flush;
      }
      else if(static_cast<Int_t>(i)%reference==0){
        cout << i << " out of " << nev << "\r" << flush;
      }
    }



    void showFFT(Int_t naverage = 10, Int_t maxevent = 0, Int_t dt = 100, bool inDecibel = true);
    void sample_plot(Int_t myevent = 0, string opt = "", Double_t filter = 0, Double_t factor = 1., Int_t mafilter = 0);
    void showWaveform(Int_t maxevent = 0, Double_t filter = 0, Int_t dt = 0);
    void persistence_plot(Int_t nbins = 500, Double_t ymin = -500, Double_t ymax = 500, Double_t filter = 0, string cut="", Double_t factor = 1);
    void add_persistence_plot(TH2D *_htemp = nullptr, Double_t filter = 0, string cut = "", Double_t factor = 1);
    TGraph drawGraph(string opt = "", Int_t n = 0, Double_t* x = nullptr, Double_t* y = nullptr);
    void debugSPE(Int_t event, Int_t moving_average, Int_t n_moving, Double_t xmin, Double_t xmax, vector<Double_t> signal_range, vector<Double_t> not_used, Double_t filter = 16, Double_t factor = 200, Double_t *SNRs = nullptr);
    void minimizeParamsSPE(Int_t event, Double_t xmin, Double_t xmax, vector<Double_t> signal_range, vector<Double_t> rangeInter = {0,0}, Double_t filter = 16, Double_t factor = 200);
    void drawZeroCrossingLines(vector<Int_t> &peaksCross, vector<Int_t> &peaksRise, TCanvas *c = nullptr, Double_t ymin = 0, Double_t ymax = 0);
    void histoTimeTrigger(Int_t nstart = 0, Int_t nfinish = 0, TH1D *_htemp = nullptr);
    void graphTimeTrigger(Int_t nstart = 0, Int_t nfinish = 0, TGraph *_gtemp = nullptr);
    void check_filtering(vector<Double_t> filter_max_and_step = {0,0}, Int_t event = 0, Int_t rebine = 1, Double_t refFreq = 20);

    ANALYZER(string m_myname = "z") : myname{m_myname}{
    }

};



