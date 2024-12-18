
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
    Int_t max_sigma = 8;

    string rootfile = "analyzed.root";


    void scan_integration_windows(Int_t ch = 0){

      SPHE dark("spe_scan");
      dark.led_calibration = true; // if external trigger + led was used, set true
      dark.channel         = ch;
      dark.filter = filter;   // one dimentional denoise filter (0 equal no filder)
      dark.check_selection = check_selection; // uses(or not) variable `selection` to discard wvfs
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
        cal.nameplotpng = hresname + ".png";
        cal.sphe_fit_try_hard(nullptr, ch, minrebin, maxbins, max_sigma,true);
        string tmp = Form("%.0f %.0f %.2f\n", start, finish, cal.snr);
        output += tmp;
      }
      cout << output << endl;

    }

};
