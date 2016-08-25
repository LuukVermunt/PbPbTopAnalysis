#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

void OpenFiles(TString FileName, int i);
void SetAddressBranches(int i);
bool checkSelectedMuon(int indexSelectedMuon, int entry);
bool checkMultipleGenParticles(std::vector<int> index);
std::vector<double> CalculateDeltaEtaAndPhi(bool possibleMultW, bool possibleMultAntiW, int w1, TLorentzVector TL_munu);
std::vector<int> FindbJets(std::vector<int> indexJets);
std::vector<int> BuildCSVVectorLeadingbJets(std::vector<int> indexJets);

const Float_t muM = .1056583715;
const Float_t nuM = .00017;

int nJets = 4;
int nbJets = 1;
const int nFiles = 1;
const bool isDebug = false;
//if(isDebug) 

TTree* tr_[nFiles];
TFile* f_[nFiles];

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
float jtPt[500];
float jtPhi[500];
float jtEta[500];
float jtM[500];
float discr_csvV1[500];
float discr_csvV2[500];
float discr_tcHighEff[500];
float discr_tcHighPur[500];
int refparton_flavorForB[500];

//Max generated particles = 5000;
int nGen;
int genPdg[5000];
float genPt[5000];
float genPhi[5000];
float genEta[5000];
float genChg[5000];


void checkOriginSelectedMuon(){

  //if(!strcmp(outFileName.c_str(), "")){
  //  std::cout << "No output specified. return" << std::endl;
  //  return;
  //}
  OpenFiles("~/Documents/MCtt_tcHigh_All.root", 0);            //NoLepIso cut done
  //TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  int test = 0, yes = 0;
  for(int i = 0; i < nFiles; i++){

    SetAddressBranches(i);

    for(int entry = 0; entry < (int)tr_[i]->GetEntries(); entry++){
      
      tr_[i]->GetEntry(entry);

      if(hiBin < 160) continue;

      int indexMuon = -999;
      double Ptlepfirst = 0.;
      for(int ilep = 0; ilep<nLep; ilep++){
        if(lepPt[ilep] > Ptlepfirst){
          Ptlepfirst = lepPt[ilep];
          indexMuon = ilep;
        }
      }

      if(lepPt[indexMuon] < 18.) continue;

      if(hiBin<20 && lepIso[indexMuon]>0.58) continue;
      else if(hiBin>=20 && hiBin<60 && lepIso[indexMuon]>0.45) continue;
      else if(hiBin>=60 && hiBin<100 && lepIso[indexMuon]>0.3) continue;
      else if(hiBin>=100 && hiBin<140 && lepIso[indexMuon]>0.24) continue;
      else if(hiBin>=140 && lepIso[indexMuon]>0.18) continue;

      std::vector<int> indexJets;
      for(int ij = 0; ij < nJt; ij++){

        if(jtPt[ij] < 30.) continue;
        if(fabs(jtEta[ij]) > 2.) continue;

        double drJetToMuon = sqrt( pow( TMath::Abs(TVector2::Phi_mpi_pi(lepPhi[indexMuon] - jtPhi[ij]) ) ,2) + pow(jtEta[ij]-lepEta[indexMuon],2) ); 
        if(drJetToMuon < 0.3)  continue;

        indexJets.push_back(ij);
      }

      std::vector<int> indexbJets = FindbJets(indexJets);
      std::vector<int> srt_indexJets = BuildCSVVectorLeadingbJets(indexJets);

      if((int)indexJets.size() < nJets && (int)indexbJets.size() < nbJets) continue;
/*
      for(int iMC = 0; iMC < nGen; iMC++){
        if(fabs(genPdg[iMC]) < 10 || fabs(genPdg[iMC]) == 24){
          cout << iMC << ": " << genPdg[iMC] << " " << genEta[iMC] << " " << genPhi[iMC] << endl;
        }
      }
      cout << endl << endl;
*/
      int iWplus(-99), iWmin(-99);
      for(int iMC = 0; iMC < nGen; iMC++){
        if(genPdg[iMC] == -24) iWmin = iMC;
        if(genPdg[iMC] == 24) iWplus = iMC;
      }

      TLorentzVector j1, j2, j12;
      std::vector<int> vector_j1, vector_j2;
      std::vector<TLorentzVector> vector_j12;
      for(int i = 0; i < nJt; i++){
        j1.SetPtEtaPhiM(jtPt[i],jtEta[i],jtPhi[i],jtM[i]);
        for(int j = i+1; j < nJt; j++){
          j2.SetPtEtaPhiM(jtPt[j],jtEta[j],jtPhi[j],jtM[j]);
          j12 = j1 + j2;

          double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(j12.Phi() - genPhi[iWmin]) );
          double deltaeta = TMath::Abs(j12.Eta() - genEta[iWmin] );

          //cout << j12.M() << " " ;
          if(j12.M() < 90. && j12.M() > 70. /*&& deltaeta < 0.2 && deltaphi < 0.2 */){
            vector_j1.push_back(i);
            vector_j2.push_back(j);
            vector_j12.push_back(j12);
          }
        }
      }

      cout << vector_j1.size() << " " << vector_j2.size() << endl;
      if((int)vector_j1.size() == 1){
        //cout << jtPt[vector_j1[0]]<<" " << jtEta[vector_j1[0]]<<" "<<jtPhi[vector_j1[0]]<<" "<<jtM[vector_j1[0]]<<endl;
        cout << vector_j12[0].Pt() <<" "<<vector_j12[0].Eta()<<" "<<vector_j12[0].Phi()<<" "<<vector_j12[0].M()<<endl;
        cout << genPt[iWplus]<<" "<<genEta[iWplus]<<" "<<genPhi[iWplus]<<endl;
        cout << genPt[iWmin]<<" "<<genEta[iWmin]<<" "<<genPhi[iWmin]<<endl;
      }

	    //if(checkSelectedMuon(indexMuon, entry)) yes++;
      test++;
  	}
    cout << yes << " good matches of " << test << " in total of " << tr_[i]->GetEntries() << " events" << endl;
  }
}



void OpenFiles(TString FileName, int i){
  f_[i] = TFile::Open(FileName);
  tr_[i] = dynamic_cast<TTree*>(f_[i]->Get("skimTree"));

  if(f_[i] == 0){
    cout << "ERROR: Failed to open file " << FileName << endl;
  } else {
    cout << "Opened file: " << FileName << endl;
  }
}

void SetAddressBranches(int i){

  tr_[i]->SetBranchAddress("run", &run_);
  tr_[i]->SetBranchAddress("evt", &evt_);
  tr_[i]->SetBranchAddress("lumi", &lumi_);
  tr_[i]->SetBranchAddress("hiBin", &hiBin);
  tr_[i]->SetBranchAddress("vz", &vz);

  tr_[i]->SetBranchAddress("nLep", &nLep);
  tr_[i]->SetBranchAddress("lepID", lepID);
  tr_[i]->SetBranchAddress("lepPt", lepPt);
  tr_[i]->SetBranchAddress("lepPhi", lepPhi);
  tr_[i]->SetBranchAddress("lepEta", lepEta);
  tr_[i]->SetBranchAddress("lepChg", lepChg);
  tr_[i]->SetBranchAddress("lepIso", lepIso);

  tr_[i]->SetBranchAddress("nJt", &nJt);
  tr_[i]->SetBranchAddress("jtPt", jtPt);
  tr_[i]->SetBranchAddress("jtPhi", jtPhi);
  tr_[i]->SetBranchAddress("jtEta", jtEta);
  tr_[i]->SetBranchAddress("jtM", jtM);
  tr_[i]->SetBranchAddress("discr_csvV1", discr_csvV1); 
  tr_[i]->SetBranchAddress("discr_csvV2", discr_csvV2); 
  tr_[i]->SetBranchAddress("discr_tcHighEff", discr_tcHighEff); 
  tr_[i]->SetBranchAddress("discr_tcHighPur", discr_tcHighPur); 
  tr_[i]->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);

  tr_[i]->SetBranchAddress("nGen", &nGen);
  tr_[i]->SetBranchAddress("genPdg", genPdg);
  tr_[i]->SetBranchAddress("genPt", genPt);
  tr_[i]->SetBranchAddress("genPhi", genPhi);
  tr_[i]->SetBranchAddress("genEta", genEta);
  tr_[i]->SetBranchAddress("genChg", genChg);

}


bool checkSelectedMuon(int indexSelectedMuon, int entry){

/////////////////////////////////////  
//Save indexes muons, muon-neutrinos and W's
/////////////////////////////////////
  std::vector<int> indexmuon, indexantimuon, indexmuonnu, indexantimuonnu, indextaunu, indexantitaunu, indexW, indexantiW;
  for(int igen = 0; igen < nGen; igen++){
    if( genPdg[igen] == 13) 		  indexmuon.push_back(igen);
    else if(genPdg[igen] == -13)	indexantimuon.push_back(igen);
    else if(genPdg[igen] == 14)		indexmuonnu.push_back(igen);
    else if(genPdg[igen] == -14)	indexantimuonnu.push_back(igen);
    else if(genPdg[igen] == 16)   indextaunu.push_back(igen);
    else if(genPdg[igen] == -16)  indexantitaunu.push_back(igen);
    else if(genPdg[igen] == 24)		indexW.push_back(igen);
    else if(genPdg[igen] == -24)	indexantiW.push_back(igen);
  }
/////////////////////////////////////

/////////////////////////////////////
//Check if there are multiple of the same particles produced
/////////////////////////////////////
  bool possibleMultMuon = checkMultipleGenParticles(indexmuon);
  bool possibleMultAntiMuon = checkMultipleGenParticles(indexantimuon);
  bool possibleMultMuonNu = checkMultipleGenParticles(indexmuonnu);
  bool possibleMultAntiMuonNu = checkMultipleGenParticles(indexantimuonnu);
  bool possibleMultTauNu = checkMultipleGenParticles(indextaunu);
  bool possibleMultAntiTauNu = checkMultipleGenParticles(indexantitaunu);
  bool possibleMultW = checkMultipleGenParticles(indexW);
  bool possibleMultAntiW = checkMultipleGenParticles(indexantiW);

/*
  cout << std::boolalpha << possibleMultMuon << " " << possibleMultAntiMuon << " ";
  cout << std::boolalpha << possibleMultMuonNu << " " << possibleMultAntiMuonNu << " ";
  cout << std::boolalpha << possibleMultW << " " << possibleMultAntiW << endl;

  if(possibleMultMuon == true){
    for(int igen = 0; igen < nGen; igen++){
      if( genPdg[igen] == 13) cout << " " << igen << " Mu (pT, Eta, Phi) (" << genPt[igen] << ", " << genEta[igen] << ", " << genPhi[igen] << ")" << endl;
    }
  }

  if(possibleMultMuonNu == true){
    for(int igen = 0; igen < nGen; igen++){
      if( genPdg[igen] == 14) cout << " " << igen << " Nu (pT, Eta, Phi) (" << genPt[igen] << ", " << genEta[igen] << ", " << genPhi[igen] << ")" << endl;
    }
  }

  if(possibleMultW == true){
    for(int igen = 0; igen < nGen; igen++){
      if( genPdg[igen] == 24) cout << " " << igen << " W (pT, Eta, Phi) (" << genPt[igen] << ", " << genEta[igen] << ", " << genPhi[igen] << ")" << endl;
    }
  }
*/
/////////////////////////////////////

/////////////////////////////////////
//Select the final particles
/////////////////////////////////////
  int mu1(-999), mu2(-999), antimu1(-999), antimu2(-999), nu1(-999), nu2(-999), antinu1(-999), antinu2(-999), w1(-999), w2(-999), antiw1(-999), antiw2(-999);
  int taunu1(-999), taunu2(-999), antitaunu1(-999), antitaunu2(-999);
  if(indexmuon.size()>0){
    if(possibleMultMuon == true){ mu1 = indexmuon[(int)indexmuon.size() - 1];  mu2 = indexmuon[(int)indexmuon.size() - 2]; }
    else{                          mu1 = indexmuon[(int)indexmuon.size() - 1]; }
  }
  if(indexantimuon.size()>0){
    if(possibleMultAntiMuon == true){  antimu1 = indexantimuon[(int)indexantimuon.size() - 1]; antimu2 = indexantimuon[(int)indexantimuon.size() - 2]; }
    else{                             antimu1 = indexantimuon[(int)indexantimuon.size() - 1]; }
  }
  if(indexmuonnu.size()>0){
    if(possibleMultMuonNu == true){ nu1 = indexmuonnu[(int)indexmuonnu.size() - 1]; nu2 = indexmuonnu[(int)indexmuonnu.size() - 2]; }
    else{                           nu1 = indexmuonnu[(int)indexmuonnu.size() - 1]; }
  }
  if(indexantimuonnu.size()>0){
    if(possibleMultAntiMuonNu == true){ antinu1 = indexantimuonnu[(int)indexantimuonnu.size() - 1];  antinu2 = indexantimuonnu[(int)indexantimuonnu.size() - 2]; }
    else{                               antinu1 = indexantimuonnu[(int)indexantimuonnu.size() - 1]; }
  }
  if(indextaunu.size()>0){
    if(possibleMultTauNu == true){  taunu1 = indextaunu[(int)indextaunu.size() - 1]; taunu2 = indextaunu[(int)indextaunu.size() - 2]; }
    else{                           taunu1 = indextaunu[(int)indextaunu.size() - 1]; }
  }
  if(indexantitaunu.size()>0){
    if(possibleMultAntiTauNu == true){  antitaunu1 = indexantitaunu[(int)indexantitaunu.size() - 1];  antitaunu2 = indexantitaunu[(int)indexantitaunu.size() - 2]; }
    else{                               antitaunu1 = indexantitaunu[(int)indexantitaunu.size() - 1]; }
  }
  if(indexW.size() > 0){
    if(possibleMultW == true){ w1 = indexW[(int)indexW.size() - 1];  w2 = indexW[(int)indexW.size() - 2]; }
    else{                      w1 = indexW[(int)indexW.size() - 1]; }
  }
  if(indexantiW.size() > 0){
    if(possibleMultAntiW == true){ antiw1 = indexantiW[(int)indexantiW.size() - 1];  antiw2 = indexantiW[(int)indexantiW.size() - 2]; }
    else{                          antiw1 = indexantiW[(int)indexantiW.size() - 1]; }
  }      
/*
  cout << "Event " << entry << endl;
  if((int)indexmuon.size()>0){ 
    cout << " Mu " << mu1 << " (pT, Eta, Phi) (" << genPt[mu1] << ", " << genEta[mu1] << ", " << genPhi[mu1] << ")" << endl;
    if(possibleMultMuon == true) cout << " Mu " << mu2 << " (pT, Eta, Phi) (" << genPt[mu2] << ", " << genEta[mu2] << ", " << genPhi[mu2] << ")" << endl;
  }
  if((int)indexantimuon.size()>0){ 
    cout << " antiMu " << antimu1 << " (pT, Eta, Phi) (" << genPt[antimu1] << ", " << genEta[antimu1] << ", " << genPhi[antimu1] << ")" << endl;
    if(possibleMultAntiMuon == true) cout << " Mu " << antimu2 << " (pT, Eta, Phi) (" << genPt[antimu2] << ", " << genEta[antimu2] << ", " << genPhi[antimu2] << ")" << endl;
  }
  if((int)indexmuonnu.size()>0){ 
    cout << " Nu " << nu1 << " (pT, Eta, Phi) (" << genPt[nu1] << ", " << genEta[nu1] << ", " << genPhi[nu1] << ")" << endl;
    if(possibleMultMuonNu == true) cout << " Nu " << nu2 << " (pT, Eta, Phi) (" << genPt[nu2] << ", " << genEta[nu2] << ", " << genPhi[nu2] << ")" << endl;
  }
  if((int)indexantimuonnu.size()>0){ 
    cout << " antiNu " << antinu1 << " (pT, Eta, Phi) (" << genPt[antinu1] << ", " << genEta[antinu1] << ", " << genPhi[antinu1] << ")" << endl;
    if(possibleMultAntiMuonNu == true) cout << " AntiNu " << antinu2 << " (pT, Eta, Phi) (" << genPt[antinu2] << ", " << genEta[antinu2] << ", " << genPhi[antinu2] << ")" << endl;
  }
  if((int)indexW.size()>0){ 
    cout << " W " << w1 << " (pT, Eta, Phi) (" << genPt[w1] << ", " << genEta[w1] << ", " << genPhi[w1] << ")" << endl;
    if(possibleMultW == true) cout << " W " << w2 << " (pT, Eta, Phi) (" << genPt[w2] << ", " << genEta[w2] << ", " << genPhi[w2] << ")" << endl;
  }
  if((int)indexantiW.size()>0){ 
    cout << " antiW " << antiw1 << " (pT, Eta, Phi) (" << genPt[antiw1] << ", " << genEta[antiw1] << ", " << genPhi[antiw1] << ")" << endl;
    if(possibleMultAntiW == true) cout << " antiW " << antiw2 << " (pT, Eta, Phi) (" << genPt[antiw2] << ", " << genEta[antiw2] << ", " << genPhi[antiw2] << ")" << endl;
  }
  cout << endl;
*/
/////////////////////////////////////

/////////////////////////////////////
//Match selected muon with generated muon
/////////////////////////////////////
  int musel(-999), antimusel(-999);
  if(possibleMultAntiMuon == false){
    //cout << " antiMu (pT, Eta, Phi) (" << lepPt[indexSelectedMuon] << ", " << lepEta[indexSelectedMuon] << ", " << lepPhi[indexSelectedMuon] << ")" << endl;
    //cout << " Gen antiMu (pT, Eta, Phi) (" << genPt[antimu1] << ", " << genEta[antimu1] << ", " << genPhi[antimu1] << ")" << endl;
    antimusel = antimu1;
  } else {
    double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[antimu1] - lepPhi[indexSelectedMuon]) );
    double deltaeta = TMath::Abs( genEta[antimu1] - lepEta[indexSelectedMuon] );

    if(deltaeta > 0.3 || deltaphi > 0.3) antimusel = antimu2;
    else antimusel = antimu1;
  }
  if(possibleMultMuon == false){
    //cout << " Mu (pT, Eta, Phi) (" << lepPt[indexSelectedMuon] << ", " << lepEta[indexSelectedMuon] << ", " << lepPhi[indexSelectedMuon] << ")" << endl;
    //cout << " Gen Mu (pT, Eta, Phi) (" << genPt[mu1] << ", " << genEta[mu1] << ", " << genPhi[mu1] << ")" << endl;
    musel = mu1;
  } else {
    double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[mu1] - lepPhi[indexSelectedMuon]) );
    double deltaeta = TMath::Abs( genEta[mu1] - lepEta[indexSelectedMuon] );

    if(deltaeta > 0.3 || deltaphi > 0.3) musel = mu2;
    else musel = mu1;
  }

  /*if(musel != -999){
    cout << " Mu (pT, Eta, Phi) (" << lepPt[indexSelectedMuon] << ", " << lepEta[indexSelectedMuon] << ", " << lepPhi[indexSelectedMuon] << ")" << endl;
    cout << " Gen Mu (pT, Eta, Phi) (" << genPt[musel] << ", " << genEta[musel] << ", " << genPhi[musel] << ")" << endl;
  }
  if(antimusel != -999){
    cout << " Mu (pT, Eta, Phi) (" << lepPt[indexSelectedMuon] << ", " << lepEta[indexSelectedMuon] << ", " << lepPhi[indexSelectedMuon] << ")" << endl;
    cout << " Gen Mu (pT, Eta, Phi) (" << genPt[antimusel] << ", " << genEta[antimusel] << ", " << genPhi[antimusel] << ")" << endl;
  }*/


  if(musel == -999 && antimusel == -999){
    cout << "No match found with muon" << endl;
    return false;
  }
/////////////////////////////////////

/////////////////////////////////////
//Match with generated neutrino
/////////////////////////////////////
  int nusel, antinusel;
  TLorentzVector TL_mu, TL_nu, TL_munu, TL_munu_temp;
  if(lepChg[indexSelectedMuon]>0){
    if(possibleMultMuonNu == false) nusel = nu1;
    else nusel = nu1;

    TL_mu.SetPtEtaPhiM(genPt[antimusel], genEta[antimusel], genPhi[antimusel], muM);
    TL_nu.SetPtEtaPhiM(genPt[nusel], genEta[nusel], genPhi[nusel], nuM);
  } else {
    if(possibleMultAntiMuonNu == false) antinusel = antinu1;
    else antinusel = antinu1;

    TL_mu.SetPtEtaPhiM(genPt[musel], genEta[musel], genPhi[musel], muM);
    TL_nu.SetPtEtaPhiM(genPt[antinusel], genEta[antinusel], genPhi[antinusel], nuM);
  }
  TL_munu = TL_mu + TL_nu;
  TL_munu_temp = TL_munu;

  std::vector<double> deltaEtaPhi;
  if(lepChg[indexSelectedMuon]>0) deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, w1, TL_munu);
  else                            deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, antiw1, TL_munu);

  bool accepted = false;
  if( (TL_munu.M() < 90. && TL_munu.M() > 70.) || (deltaEtaPhi[0] < 0.3 && deltaEtaPhi[1] < 0.3) ){
    accepted = true;
  }
  else{
    if(lepChg[indexSelectedMuon]>0){
      if(possibleMultMuonNu == true)      TL_nu.SetPtEtaPhiM(genPt[nu2], genEta[nu2], genPhi[nu2], nuM);
    } else {
      if(possibleMultAntiMuonNu == true)  TL_nu.SetPtEtaPhiM(genPt[antinu2], genEta[antinu2], genPhi[antinu2], nuM);
    }
  }
  TL_munu = TL_mu + TL_nu;

  deltaEtaPhi.clear();
  if(lepChg[indexSelectedMuon]>0) deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, w1, TL_munu);
  else                            deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, antiw1, TL_munu);

  if(accepted == false){
    if( (TL_munu.M() < 90. && TL_munu.M() > 70.) || (deltaEtaPhi[0] < 0.3 && deltaEtaPhi[1] < 0.3) ){
      accepted = true;
    } else {

      if(taunu1 == -999 && antitaunu1 == -999) return false;

      TLorentzVector TL_taunu, TL_taunu2, TL_mununu, TL_mununu_temp;
      if(lepChg[indexSelectedMuon]>0){   
        TL_taunu.SetPtEtaPhiM(genPt[taunu1], genEta[taunu1], genPhi[taunu1], 100*nuM);
        if(possibleMultTauNu == true) TL_taunu2.SetPtEtaPhiM(genPt[taunu2], genEta[taunu2], genPhi[taunu2], 100*nuM);
      }
      else{ 
        TL_taunu.SetPtEtaPhiM(genPt[antitaunu1], genEta[antitaunu1], genPhi[antitaunu1], 100*nuM);
        if(possibleMultAntiTauNu == true) TL_taunu2.SetPtEtaPhiM(genPt[antitaunu2], genEta[antitaunu2], genPhi[antitaunu2], 100*nuM);
      }
      TL_mununu = TL_munu_temp + TL_taunu;
      TL_mununu_temp = TL_mununu;

      deltaEtaPhi.clear();
      if(lepChg[indexSelectedMuon]>0) deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, w1, TL_mununu);
      else                            deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, antiw1, TL_mununu);

      if( (TL_mununu.M() < 100. && TL_mununu.M() > 60.) || (deltaEtaPhi[0] < 0.4 && deltaEtaPhi[1] < 0.4) ){
        //Muon from W-decay. Larger interval because three particles instead of 2
      }
      else{
        TL_mununu = TL_munu + TL_taunu;
      }

      deltaEtaPhi.clear();
      if(lepChg[indexSelectedMuon]>0) deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, w1, TL_mununu);
      else                            deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, antiw1, TL_mununu);
      
      if( (TL_mununu.M() < 100. && TL_mununu.M() > 60.) || (deltaEtaPhi[0] < 0.4 && deltaEtaPhi[1] < 0.4) ){
        //Muon from W-decay
      } else if(possibleMultTauNu == true || possibleMultAntiTauNu == true) {
        TL_mununu = TL_munu_temp + TL_taunu2;
      } 

      deltaEtaPhi.clear();
      if(lepChg[indexSelectedMuon]>0) deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, w1, TL_mununu);
      else                            deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, antiw1, TL_mununu);


      if( (TL_mununu.M() < 100. && TL_mununu.M() > 60.) || (deltaEtaPhi[0] < 0.4 && deltaEtaPhi[1] < 0.4) ){
        //Muon from W-decay
      } else if(possibleMultTauNu == true || possibleMultAntiTauNu == true) {
        TL_mununu = TL_munu + TL_taunu2;
      } 

      deltaEtaPhi.clear();
      if(lepChg[indexSelectedMuon]>0) deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, w1, TL_mununu);
      else                            deltaEtaPhi = CalculateDeltaEtaAndPhi(possibleMultW, possibleMultAntiW, antiw1, TL_mununu);

      if( (TL_mununu.M() < 100. && TL_mununu.M() > 60.) || (deltaEtaPhi[0] < 0.4 && deltaEtaPhi[1] < 0.4) ){
        //Muon from W-decay
      } else {
        cout << "LV_MuNu (pT, Eta, Phi, M) (" << TL_munu.Pt() << ", " << TL_munu.Eta() << ", " << TL_munu.Phi() << ", " << TL_munu.M() << ")" << endl;
        cout << "LV_MuNu_temp (pT, Eta, Phi, M) (" << TL_munu_temp.Pt() << ", " << TL_munu_temp.Eta() << ", " << TL_munu_temp.Phi() << ", " << TL_munu_temp.M() << ")" << endl;
        cout << "LV_MuNuNu (pT, Eta, Phi, M) (" << TL_mununu.Pt() << ", " << TL_mununu.Eta() << ", " << TL_mununu.Phi() << ", " << TL_mununu.M() << ")" << endl;
        cout << "LV_MuNuNu_temp (pT, Eta, Phi, M) (" << TL_mununu_temp.Pt() << ", " << TL_mununu_temp.Eta() << ", " << TL_mununu_temp.Phi() << ", " << TL_mununu_temp.M() << ")" << endl;
        if(lepChg[indexSelectedMuon]>0) cout << " Gen W (pT, Eta, Phi) (" << genPt[w1] << ", " << genEta[w1] << ", " << genPhi[w1] << ")" << endl;
        else                            cout << " Gen W (pT, Eta, Phi) (" << genPt[antiw1] << ", " << genEta[antiw1] << ", " << genPhi[antiw1] << ")" << endl;
        cout << "W " << std::boolalpha << possibleMultW << " " << possibleMultAntiW << endl;
        cout << "TauNu " << std::boolalpha << possibleMultTauNu << " " << possibleMultAntiTauNu << endl;
        cout << "      Selected muon doesn't correspond to W-decay" << endl;
        return false;
      }
    }
  }
  
  /*if(lepChg[indexSelectedMuon]>0){
    cout << " LV_MuNu (pT, Eta, Phi, M) (" << TL_munu.Pt() << ", " << TL_munu.Eta() << ", " << TL_munu.Phi() << ", " << TL_munu.M() << ")" << endl;
    cout << " Gen W (pT, Eta, Phi) (" << genPt[w1] << ", " << genEta[w1] << ", " << genPhi[w1] << ")" << endl;
    if(possibleMultW == true) cout << " Gen W (pT, Eta, Phi) (" << genPt[w2] << ", " << genEta[w2] << ", " << genPhi[w2] << ")" << endl;
  } else {
    cout << " LV_MuNu (pT, Eta, Phi, M) (" << TL_munu.Pt() << ", " << TL_munu.Eta() << ", " << TL_munu.Phi() << ", " << TL_munu.M() << ")" << endl;
    cout << " Gen W (pT, Eta, Phi) (" << genPt[antiw1] << ", " << genEta[antiw1] << ", " << genPhi[antiw1] << ")" << endl;
    if(possibleMultAntiW == true) cout << " Gen W (pT, Eta, Phi) (" << genPt[antiw2] << ", " << genEta[antiw2] << ", " << genPhi[antiw2] << ")" << endl;
  }
  cout << endl;*/
  cout << "------Selected muon comes from W-decay" << endl;
  return true;
}

bool checkMultipleGenParticles(std::vector<int> index){
	
  bool multiple = false;
  int i = (int)index.size() - 2;
  if((int)index.size() > 1){
    if(genPt[index[i]] - genPt[index[i+1]] < -10.)              multiple = true;
    if(TMath::Abs(genPhi[index[i+1]] - genPhi[index[i]]) > 0.5) multiple = true;
    if(TMath::Abs(genEta[index[i+1]] - genEta[index[i]]) > 0.5) multiple = true;
  }

  return multiple;
}


std::vector<double> CalculateDeltaEtaAndPhi(bool possibleMultW, bool possibleMultAntiW, int w1, TLorentzVector TL_munu){

  //Don't take into account the events with multiple W's. Think of something

  double deltaphi, deltaeta;
  if(possibleMultW == false && possibleMultAntiW == false){
    deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(TL_munu.Phi() - genPhi[w1]) );
     deltaeta = TMath::Abs(TL_munu.Eta() - genEta[w1] );
  } else {
    deltaphi = 99.;
    deltaeta = 99.;
  }

  std::vector<double> EtaAndPhi;

  EtaAndPhi.push_back(deltaeta);
  EtaAndPhi.push_back(deltaphi);

  return EtaAndPhi;

}

std::vector<int> FindbJets(std::vector<int> indexJets){
  
  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(discr_csvV1[indexJets[ij]] > 0.75)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}

std::vector<int> BuildCSVVectorLeadingbJets(std::vector<int> indexJets){
/*Returns a vector with the indices of the jets from highest CSV to lowest CSV*/

  const int array_size = (int) indexJets.size();

  int sortedindex_indexJets[array_size];
  double CSV_array[array_size];

  for(int i = 0; i < (int)indexJets.size(); i++){ 
    sortedindex_indexJets[i] = i;
    CSV_array[i] = discr_csvV1[indexJets[i]];
  }
  TMath::Sort(array_size, CSV_array, sortedindex_indexJets, true);

  std::vector<int> strind_indJets(sortedindex_indexJets, sortedindex_indexJets + array_size);
  return strind_indJets;
}