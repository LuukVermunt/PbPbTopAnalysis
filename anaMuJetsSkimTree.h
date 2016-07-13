#ifndef anaMuJetsSkimTree_h
#define anaMuJetsSkimTree_h

#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>



class anaMuJetsSkimTree{

public:

  anaMuJetsSkimTree();
  ~anaMuJetsSkimTree();

  //Interface functions
  void OpenFiles(TString FileName, int i);
  void BuildHistograms();
  void CalculateNormalizationHistograms(TH1F* h[4], int option);
  void NormalizeHistograms(int option, int Drawoption);
  double CalculateDeltaPhi(int ilep, int index_b);
  void SetAddressBranches(int i);
  std::vector<int> BuildCSVVectorLeadingbJets(std::vector<int> indexJets);
  void FillHistogramsBeforeCuts(int option, int i);
  void FillHistogramsAfterCuts(int option, int i);
  void LayoutHistograms(int option);
  TCanvas* PlotHistograms(int option2);

private:
  //Constant variables
  const Float_t muM = .1056583715;

  //Objects
  TFile* f_[21];      //Change, or think of something, when more files are included
  TTree* tr_[21];     //Change, or think of something, when more files are included

  //4 types of files: Data, MC tt, MC W, MC DY
  TH1F* h_events_[4];
  TH1F* h_totevents[4];

  TH1F* h_lepPt_[4];
  TH1F* h_lepPseu_[4];
  TH1F* h_jetHt_[4];
  TH1F* h_jetCSV_[4];
  //TH1F* h_lepbjetMinv_[4];
  TH1F* h_lepbjetMinv_min_[4];
  //TH1F* h_ttMinv_[4];
  TH1F* h_ttMinv_min_[4];
  TH1F* h_Phi_1stb_[4];
  TH1F* h_Phi_allb_[4];
  TH1F* h_Phi_b_minpi_[4];
  TH2F* h_2dCSV_[4];
  TH2F* h_2dCSV2_[4];
  TGraph* CSV_effcut;
  TGraph* CSV_effcut2;

  THStack* lepPt_Stack;
  THStack* lepPseu_Stack;
  THStack* jetHt_Stack;
  THStack* jetCSV_Stack;
  THStack* lepbjetMinv_Stack;
  THStack* ttMinv_Stack;


  //Global variables
  unsigned int run_, lumi_;
  ULong64_t evt_;
  int hiBin;
  float vz;
  
  //Max leptons in skimFile = 2, to be sure I used 4
  int nLep;
  int lepID[4];
  float lepPt[4];
  float lepPhi[4];
  float lepEta[4];
  int lepChg[4];
  float lepIso[4];

  //Max jets = 500
  int nJt;
  double sum_jtPt;
  int numb_bJets;
  float jetPt[500];
  float jtPhi[500];
  float jtEta[500];
  float jtM[500];
  float discr_csvV1[500];


  //Cut values
  double CSVCut1_;
  double CSVCut2_;
};

#endif