#include "WVFAnaSoLAr.h"

class IntegrationScan{
  public:
    bool check_selection = false; // uses(or not) variable `selection` to discard wvfs
    bool withfilter      = true; // Integrate in the filtered waveform
    Double_t filter = 12;
    Double_t hnbins = 2000;
    Double_t hxmin = -2.5;
    Double_t hxmax = 40;
    bool normalize_histogram = true;
    vector<Double_t> start_first_stop_step = {0,0,0,1};

    Int_t minrebin = 1;
    Int_t maxbins = 2;
    Int_t min_sigma = 1;
    Int_t max_sigma = 8;


    Double_t too_big = 15000;

    string rootfile = "analyzed.root";


    void scan_integration_windows(Int_t ch = 0){

      SPHE dark("spe_scan");
      dark.led_calibration = true; // if external trigger + led was used, set true
      dark.channel         = ch;
      dark.filter = filter;   // one dimentional denoise filter (0 equal no filder)
      dark.check_selection = check_selection; // uses(or not) variable `selection` to discard wvfs
      dark.too_big = too_big;
      dark.withfilter      = true; // Integrate in the filtered waveform
      dark.hnbins = hnbins;
      dark.hxmin = hxmin;
      dark.hxmax = hxmax;
      dark.normalize_histogram = normalize_histogram;
      dark.rootfile = rootfile;


      Double_t start = start_first_stop_step[0];
      Double_t first = start_first_stop_step[1];
      Double_t step = start_first_stop_step[3];
      string output = "";
      for (Double_t finish = first; finish <= start_first_stop_step[2]; finish+=step){
        Calibration cal;
        cal.make_free_stddevs = true;
        string hresname = Form("hcharge_%d_%d_Ch%d", static_cast<Int_t>(start), static_cast<Int_t>(finish), ch);
        dark.outrootname = hresname+".root";
        dark.start = start;
        dark.finish = finish;
        dark.giveMeSphe();
        dark.fout->Close();

        cal.rootFile = dark.outrootname;
        cal.nameplotsave = hresname + ".png";
        cal.sphe_fit_try_hard(nullptr, ch, minrebin, maxbins, min_sigma, max_sigma, true);
        string tmp = Form("%.0f %.0f %.2f\n", start, finish, cal.snr);
        output += tmp;
      }
      cout << output << endl;

    }

};



class WCuts{
  public:
    Double_t filter;
    Double_t start;
    Double_t finish;
    Double_t cut;
    string type = "higher";
    string perform = "include";
    WCuts( Double_t _filter, Double_t _start, Double_t _finish, Double_t _cut, string _type, string _perform):
      filter{_filter}, start{_start}, finish{_finish}, cut{_cut}, type{_type}, perform{_perform} {} 
};
class WaveSelector{
  public:
    ANALYZER *zWS = nullptr;
    string selection = "";

    Double_t filter = 100;
    Double_t start = 100;
    vector<WCuts> cuts;// always:
                                 // filter, start, stop, cut,
                                 // "higher/lower", 1, -1
                                 // include/exclude
    bool isSPE = false;
    void waveselector(int ch,
                      string filename,
                      bool debug = true
                      ){

      cout << "Computing waveform..." << endl;

      zWS = new ANALYZER("zWS");
      zWS->setAnalyzer(filename);

      bool chok = false;
      if (zWS->channels.size()>1) chok = zWS->setChannel(Form("Ch%d.", ch));
      else chok = true;
      if (!chok) return;

      TGraph *gsphe = nullptr;
      vector<Double_t> avgspe(zWS->n_points,0);
      vector<Double_t> time(zWS->n_points,0);

      bool stop = false;
      if (isSPE){
        zWS->getSelection(Form("Ch%d.selection == 1", ch));
      }
      else{
        zWS->getSelection(selection);
      }

      for (const WCuts &cut: cuts)
      {
        if (cut.perform == "include")
          zWS->genSelectByAmplitude(cut.filter,cut.start, cut.start, cut.cut, cut.type);
        else
          zWS->excludeByAmplitude(cut.filter,cut.start, cut.finish, cut.cut, cut.type);
      }

      Double_t navg = 0;
      for (int i = 0; i < zWS->lev->GetN(); i++) {
        Int_t iev = zWS->lev->GetEntry(i);
        zWS->getWaveform(iev);
        for(Int_t j = 0; j < zWS->n_points; j++){
          avgspe[j] += zWS->ch[zWS->kch]->wvf[j];
        }
        navg += 1;

      }
      cout << "Valid Waveforms = " << navg << endl;
      for(Int_t i = 0; i < zWS->n_points; i++){
        avgspe[i] = avgspe[i]/navg;
        time[i] = i*zWS->dtime;
      }
      gsphe = new TGraph(zWS->n_points, &time[0], &avgspe[0]);
      zWS->baselineparameters.baselineStart      = 0;
      zWS->baselineparameters.baselineTime       = 2300;
      zWS->baselineparameters.exclusion_baseline = 20;
      zWS->baselineparameters.exclusion_window = 400;
      zWS->reval_baseline();
      if (debug)
      {
        zWS->persistence_plot(2000,-2000,30000,10,"use_mine");
        // zWS->getFFT();
        zWS->getWaveFromGraph(gsphe);
        zWS->drawGraph("SAME");
      }
    }
};


