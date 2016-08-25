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

  TTree* tr_[5];     //Change, or think of something, when more files are included
  TRandom3* rand = new TRandom3();
  bool JetSmearing = true;

  anaMuJetsSkimTree();
  ~anaMuJetsSkimTree(); 

  //Interface functions
  void SetCSVCut(double CSVcut) {/*cout << "Changing CSVCut from " << CSVCut_ << " to " << CSVcut << endl*/CSVCut_ = CSVcut; };
  void SetnumbB(double numbB) {cout << "Changing number b-jets Cut from " << numbB_ << " to " << numbB << endl; numbB_ = numbB; };
  void SetnumbJ(double numbJets) {cout << "Changing number jets Cut from " << numbJets_ << " to " << numbJets << endl; numbJets_ = numbJets; };
  void SetjetEtaCut(double jtEtaCut) {cout << "Changing jtEta Cut from " << jtEtaCut_ << " to " << jtEtaCut << endl; jtEtaCut_ = jtEtaCut; }
  
  double GetCSVCut() {return CSVCut_; };
  int GethiBin() {return hiBin; };
  double GetjetPt(int index) {return jetPt[index]; };
  double GetlepPt(int index) {return lepPt[index]; };
  double GetlepIso(int index) {return lepIso[index]; };
  double GetjetCSVv1(int index) {return discr_csvV1[index]; };
  double GetjetCSVv2(int index) {return discr_csvV2[index]; };
  int GetnLep() {return nLep; };

  void OpenFiles(TString FileName, int i);
  void BuildHistograms();
  int FindLeadingMuonAfterCuts(int option);
  std::vector<int> FindJetsAfterCuts(int indexMuon);
  std::vector<int> FindbJets(std::vector<int> indexJets);
  std::vector<int> FindbJetsCSVv2(std::vector<int> indexJets);
  void CalculateNormalizationHistograms(TH1F* h[4], int option, int nBjets);
  void NormalizeHistograms(int option, int nBjets);
  void NormalizeHistogramsToOne(int i, int nBjets);
  double CalculateDeltaPhi(int ilep, int index_b);
  void SetAddressBranches(int i);
  std::vector<int> BuildCSVVectorLeadingbJets(std::vector<int> indexJets);
  std::vector<int> BuildCSVVectorLeadingbJetsCSVv2(std::vector<int> indexJets);
  std::vector<int> BuildpTVectorLeadingbJets(std::vector<int> indexbJets);
  void FillHistogramsBeforeCuts(int option);
  void FillHistogramsAfterCuts(int option);

private:
  //Constant variables
  const Float_t muM = .1056583715;

  //Objects
  TFile* f_[5];      //Change, or think of something, when more files are included

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
  float discr_csvV2[500];

  //Cut values
  int numbB_ = 1;
  int numbJets_ = 4;
  double CSVCut_ = 0.75;
  double mupTCut_ = 18;
  double jtpTCut_ = 30.;
  double jtEtaCut_ = 2.;
  double drJetToMuonCut_ = 0.3;
  double WMassHighCut_ = 170.;
};

#endif
