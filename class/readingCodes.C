// ________________________________________ //
// Author: Henrique Souza
// Filename: readingCodes.C
// Created: 2024
// ________________________________________ //
#include "WVFAnaSoLAr.h"
#include <arpa/inet.h>


class Read{

  private:

    // Variables in the TTree
    UInt_t m_event = 0;
    UInt_t m_sn = 0;
    UInt_t m_ch = 0;
    ULong64_t m_ntpm_ms = 0;
    ULong64_t m_taim_ns = 0;
    ULong64_t m_tslppsm_ns = 0;
    TH1S *m_h = nullptr;

    string fmapname;
    UInt_t m_totalevents = 0;
  public:



    Int_t n_points;

    Double_t dtime = 16; // steps (ADC's MS/s, 500 MS/s = 2 ns steps)
    Int_t nbits = 14;
    Int_t basebits = nbits;
    Bool_t saveFilter       = false;

    Double_t startCharge = 3300;
    Double_t maxRange = 5000;
    Double_t fast = 500;
    Double_t slow = 17000;
    Double_t filter = 9;
    Double_t exclusion_baseline = 30;
    map<UInt_t, Double_t> map_exclusion_threshold_baselines = {{0,0}};
    Double_t exclusion_window = 500;
    Double_t currentTime = 0;

    Bool_t OnlySomeEvents = false; // Do you want only one event? Choose it wisely
    Int_t stopEvent = 2000;

    Double_t baselineTime = 10000; // time limit to start looking for baseline
    Double_t baselineStart = 0;
    Double_t chargeTime = 18000; // last time to integrate
    Bool_t noBaseline=false;
    Double_t baselineFraction = 1/3.;

    Double_t nbinsbase = TMath::Power(2,basebits);
    Double_t minbase = -TMath::Power(2, basebits-1);
    Double_t maxbase =  TMath::Power(2,basebits-1);


    string standard_log_file = "fileslog.log";

    Double_t polarity=-1; // set -1 invert to negative pulses
    vector<UInt_t> channels = {};
    vector<double> exclusion_baselines = {};

    TH1D *hbase;
    Double_t max_bits = 0;

    std::map<vector<uint>, uint> chmap; // map board, ch -> chid
    std::map<uint, vector<uint>> channels_to_boardch_map; // map chid -> board, ch
    std::map<vector<uint>, uint> channels_to_get_map; // map of channels to pick, gotta go fast
    std::map<uint, uint> boardorder_map; // map board serial in order...


    DENOISE dn;



    std::map<std::vector<uint>, uint> load_channel_map(){
      std::map<vector<uint>, uint> tmpchmap;
      std::filesystem::path filePath = __FILE__; // this file path
      fmapname = filePath.parent_path();
      fmapname += "/channelmaps/solars_channel_map.csv";
      ifstream fmap;
      fmap.open(fmapname.c_str(), ifstream::in);
      UInt_t board;
      UInt_t ch;
      UInt_t fch;
      while(!fmap.fail()){
        fmap >> board;
        if(fmap.bad() || fmap.fail()){
          break;
        }
        fmap >> ch;
        fmap >> fch;
        vector<uint> v = {board, ch};
        tmpchmap[v] = fch;
        boardorder_map[board] = 1;
      }
      fmap.close();
      return tmpchmap;
    }

    void read_info_from_root(string logfilename, bool load_ch_map_from_root, bool update_map){

      auto tmpmap = this->load_channel_map();
      cout << "Getting channels maps..." << endl;
      if (!load_ch_map_from_root){
        chmap = tmpmap;
      }

      ifstream logfile;
      logfile.open(logfilename.c_str(),ios::in);

      if(logfile.good() && logfile.is_open()){ // Ok
      }
      else{
        cout << "The file " << logfilename << " did not open!!" << endl;
        return;
      }
      string dataname;
      while(!logfile.fail()){
        logfile >> dataname;
        if(logfile.bad() || logfile.fail()){
          break;
        }
        TFile *ftmp = new TFile(dataname.c_str(), "READ");
        TTree *ttmp = (TTree*)ftmp->Get("rwf");
        TBranch *bevent = ttmp->GetBranch("event");
        TBranch *bsn = ttmp->GetBranch("sn");
        TBranch *bch = ttmp->GetBranch("ch");
        bevent->SetAddress(&m_event);
        bsn->SetAddress(&m_sn);
        bch->SetAddress(&m_ch);

        auto nentries = ttmp->GetEntries();
        for (auto i = 0; i < nentries; i++){
          bevent->GetEntry(i);
          bsn->GetEntry(i);
          bch->GetEntry(i);

          vector<uint> board_ch = {m_sn, m_ch};
          if (load_ch_map_from_root)
          {
            if (chmap[board_ch] == 1){
              // All scanned...
              break;
            }
            chmap[board_ch] += 1;
            boardorder_map[m_sn] = 1;
          }

        }

        ftmp->Close();
        delete ftmp;
        auto nchannels = chmap.size();
        if (nentries%nchannels != 0){
          throw std::invalid_argument("WTF!");
        }
        auto events_per_channel = nentries/nchannels;
        m_totalevents += events_per_channel;
      }
      int aux = 0;
      for(const auto& boards : boardorder_map)
      {
        boardorder_map[boards.first] = aux;
        aux+=1;
      }

      logfile.close();

      if (!load_ch_map_from_root){
        return;
      }
      else{
        for(const auto& elem : chmap)
        {
          // std::cout << elem.first[0] << " " << elem.first[1] << " " << elem.second << std::endl;
          UInt_t fboard = elem.first[0];
          UInt_t fBoardF = boardorder_map[fboard];
          UInt_t ch = elem.first[1];
          UInt_t fch = fBoardF*100 + ch;
          chmap[elem.first] = fch;
          if (!update_map && tmpmap.count(elem.first)==0){
            cout << "WARNING!!! .... Channel map in repository needs to be updated..." << endl;
          }
        }
      }
    }
    void get_ch_info(string logfilename, bool load_ch_map_from_root, bool update_map=false){

      read_info_from_root(logfilename, load_ch_map_from_root, update_map);
      bool fill_channels = channels.empty();
      // creating inverted map
      for(const auto& elem : chmap){
          UInt_t fboard = elem.first[0];
          UInt_t fch = elem.first[1];
          UInt_t fBoardF = elem.second;
          channels_to_boardch_map[fBoardF] = {fboard, fch};

          if (fill_channels){
            channels.push_back(fBoardF);
          }
      }
      for (const uint &chnumber: channels ){
        channels_to_get_map[channels_to_boardch_map[chnumber]] = 1;
      }

      if (update_map && load_ch_map_from_root){
        ofstream fout;
        fout.open(fmapname, std::ofstream::out);
        for(const auto& elem : chmap)
        {
          // std::cout << elem.first[0] << " " << elem.first[1] << " " << elem.second << std::endl;
          UInt_t fboard = elem.first[0];
          UInt_t fBoardF = boardorder_map[fboard];
          UInt_t fch = fBoardF*100 + elem.first[1];
          fout << fboard << " " << elem.first[1] << " " << fch << "\n";
        }
      }
      setup_baseline_thresholds();

    }

    void convert_from_root(Bool_t load_ch_map_from_root=true, bool update_map=false){
      max_bits = TMath::Power(2,nbits);
      if (maxRange < startCharge) maxRange = startCharge + startCharge;
      get_ch_info(standard_log_file, load_ch_map_from_root, update_map);
      readFiles(standard_log_file); //use it like this
      return;
    }


    void readFiles(string files){

      ifstream filesdata;
      filesdata.open(files.c_str(),ios::in);
      string rootfile;
      TFile *f1;
      TTree *t1;
      TTree *thead;
      hbase = new TH1D("hbase","finding baseline", nbinsbase, minbase, maxbase);
      Double_t tEvent = 0;
      vector<ADC_DATA*> ch(channels.size());
      vector<TBranch*> bch(channels.size());


      Bool_t first_file = true;

      if(filesdata.good() && filesdata.is_open()){ // Ok
        cout << "Reading file " << files << " ... " << endl;
      }
      else{
        cout << "The file " << files << " did not open!!" << endl;
        return;
      }
      string temp = "";
      // While for read the entire file
      while(!filesdata.fail()){
        filesdata >> temp;
        if(filesdata.bad() || filesdata.fail()){
          break;
        }

        // ______________________________ Create root files at first files only __________________________________
        if(first_file){
          first_file = false;
          rootfile = "analyzed.root";

          // string erase = "rm " + rootfile;
          // system(erase.c_str());

          f1 = new TFile(rootfile.c_str(),"RECREATE");
          t1 = new TTree("t1","ADC processed waveform");
          t1->SetEntries(m_totalevents);

          for (const uint &chnumber: channels ){
            auto i = &chnumber - &channels[0];
            bch[i] = t1->Branch(Form("Ch%d.",chnumber),&ch[i]);
          }

          thead = new TTree("head","header");
          thead->Branch("dtime", &dtime);
          thead->Branch("startcharge", &startCharge);
          thead->Branch("endcharge", &chargeTime);
          thead->Branch("baselineStart", &baselineStart);
          thead->Branch("baselineTime", &baselineTime);
          thead->Branch("maxRange", &maxRange);
          thead->Branch("fast_time", &fast);
          thead->Branch("slow_time", &slow);
          thead->Branch("exclusion_baselines", &exclusion_baselines);
          thead->Branch("exclusion_window", &exclusion_window);
          thead->Branch("filter", &filter);
          thead->Fill();

          f1->Write();

          f1->Close();
        }
        // _______________________________________________________________________________________________________

        readData(temp,rootfile,tEvent);

      }

    }






    // This function will read your data and create a root file with the same name
    void readData(string ifile,string rootfile, Double_t &tEvent){

      TFile *f1 = new TFile(rootfile.c_str(),"UPDATE");
      TTree *t1 = (TTree*)f1->Get("t1");

      map<uint, ADC_DATA*> ch;
      map<uint, TBranch*> bcho;


      for (const uint &i: channels ){
        bcho[i] = t1->GetBranch(Form("Ch%d.",i));
        bcho[i]->SetAddress(&ch[i]);
      }


      TFile *fin = new TFile(ifile.c_str(), "READ");
      TTree *trwf = (TTree*)fin->Get("rwf");
      TBranch *bevent = trwf->GetBranch("event");
      TBranch *bsn = trwf->GetBranch("sn");
      TBranch *bch = trwf->GetBranch("ch");
      // TBranch *bntpm_ms = trwf->GetBranch("ntpm_ms");
      // TBranch *btaim_ns = trwf->GetBranch("taim_ns");
      // TBranch *btslppsm_ns = trwf->GetBranch("tslppsm_ns");
      TBranch *bth1s_ptr = trwf->GetBranch("th1s_ptr");

      bevent->SetAddress(&m_event);
      bsn->SetAddress(&m_sn);
      bch->SetAddress(&m_ch);
      // bntpm_ms->SetAddress(&m_ntpm_ms);
      // btaim_ns->SetAddress(&m_taim_ns);
      // btslppsm_ns->SetAddress(&m_tslppsm_ns);
      bth1s_ptr->SetAddress(&m_h);

      auto nentries = trwf->GetEntries();

      for (auto i = 0; i < nentries; i++){
        bevent->GetEntry(i);
        bsn->GetEntry(i);
        bch->GetEntry(i);
        if(i%200==0){
          cout << i << " out of " << nentries << "\r" << flush;
        }
        vector<uint> board_ch = {m_sn, m_ch};
        if (channels_to_get_map.count(board_ch)>0)
        {
          uint thechannel = chmap[board_ch];
          bth1s_ptr->GetEntry(i);
          ch[thechannel]->event = m_event;
          ch[thechannel]->board = m_sn;
          getvalues(thechannel, *ch[thechannel], m_h);
          bcho[thechannel]->Fill();
        }

      }
      cout << "\n";

      f1->WriteObject(t1,"t1","TObject::kOverwrite");

      f1->Close();

    }



    void getvalues(uint nch, ADC_DATA &ch, TH1S *htmp){
      Int_t n_points = htmp->GetNbinsX();
      ch.Set_npts(n_points);
      vector<Double_t> raw(n_points);
      vector<Double_t> filtered(n_points);
      for (Int_t i = 0; i < n_points; i++){
        ch.wvf[i] = htmp->GetBinContent(i+1)*polarity;
        filtered[i] = ch.wvf[i];
      }
      if(filter>0) dn.TV1D_denoise<Double_t>(&ch.wvf[0],&filtered[0],n_points,filter);
      Double_t bl = baseline(&filtered[0],ch.selection, map_exclusion_threshold_baselines[nch]);

      ch.peak = 0;
      ch.base = bl;
      Double_t fastcomp = 0;
      Double_t slowcomp = 0;
      ch.charge=0;
      for(Int_t i = 0; i<n_points; i++){
        ch.wvf[i] = ch.wvf[i]-bl;
        filtered[i] = filtered[i]-bl;
        if(i>=startCharge/dtime && i<chargeTime/dtime){
          ch.charge+=filtered[i]*dtime;
          if(i <= maxRange/dtime){
            if(ch.peak==0){ ch.peak = filtered[i]; ch.peakpos = i*dtime;}
            else if(ch.peak<filtered[i]){ch.peak=filtered[i];ch.peakpos = i*dtime;}
          }
          if(i<(startCharge+fast)/dtime){
            fastcomp+=ch.wvf[i];
          }
          if(i<(startCharge+slow)/dtime){
            slowcomp+=ch.wvf[i];
          }

        }
      }
//     cout << fastcomp << " " << slowcomp << endl;
      ch.fprompt = fastcomp/slowcomp;
    }


    Double_t baseline(Double_t v[], UInt_t &selection, Double_t exclusion_threshold){
      if(noBaseline) return 0;
      Double_t result = 0;
      hbase->Reset();
      for(Int_t i=baselineStart/dtime; i<baselineTime/dtime; i++) hbase->Fill(v[i]);
      Double_t res0 = hbase->GetBinCenter(hbase->GetMaximumBin());
      Double_t hmean = hbase->GetMean();
      Double_t hstd = hbase->GetStdDev();


      Bool_t changed_mean = false;


      Double_t bins=0;
      for(Int_t i=baselineStart/dtime; i<baselineTime/dtime;){
        if(v[i] > res0 + exclusion_threshold || v[i]<res0 - exclusion_threshold) {
          i+=exclusion_window/dtime;
        }
        else{
          result += v[i];
          bins+=1;
          i++;
        }
      }
      if(bins>0) result/=bins;
      if(bins > ((baselineTime - baselineStart)/dtime)*baselineFraction){
        selection = 0;
        return result;
      }
      else{
        if(changed_mean==false)selection = 1;
        else selection = 2;
        return res0;
      }
    }

    void setup_baseline_thresholds(){
      Int_t nchannels = (int)channels.size();
      auto mapIter = map_exclusion_threshold_baselines.find(0);

      if (mapIter != map_exclusion_threshold_baselines.end()){ // in case the map was not set
        map_exclusion_threshold_baselines.erase(0);
        for (const uint &chnumber: channels ){
          map_exclusion_threshold_baselines[chnumber] = exclusion_baseline;
        }
      }
      else{
        for (const uint &chnumber: channels ){
          if (map_exclusion_threshold_baselines.count(chnumber) == 0)
          {
            map_exclusion_threshold_baselines[chnumber] = exclusion_baseline;
          }
          exclusion_baselines.push_back(map_exclusion_threshold_baselines[chnumber]);
        }

      }
    }


};


