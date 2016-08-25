#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

void OpenFiles(TString FileName, int i);
void SetAddressBranches(int i);
std::vector<int> FindbJetsCSVv1(std::vector<int> indexJets, double cutValue);
std::vector<int> FindbJetstcHighEff(std::vector<int> indexJets, double cutValue);
std::vector<int> FindbJetstcHighPur(std::vector<int> indexJets, double cutValue);
std::vector<int> FindbJetssvtxm(std::vector<int> indexJets, double cutValue);
std::vector<int> FindbJetssvtxpt(std::vector<int> indexJets, double cutValue);
std::vector<int> FindbJetsprob(std::vector<int> indexJets, double cutValue);
std::vector<int> ownFlavourMatchingAlgo(int entry);

int nJets = 4;
int nbJets = 1;
//const int nDiscrCuts = 180; //tcHigh
//const int nDiscrCuts = 8; //svtxm
//const int nDiscrCuts = 15; //svtxpt
const int nDiscrCuts = 20; int nCut = 15;//discr_prob
const int nFiles = 3;
const bool isDebug = false;
//if(isDebug) std::cout << __LINE__ << std::endl;

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
float svtxm[500];
float svtxpt[500];
float discr_prob[500];
int refparton_flavorForB[500];

//Max generated particles = 5000;
int nGen;
int genPdg[5000];
float genPt[5000];
float genPhi[5000];
float genEta[5000];
float genChg[5000];


void plotDiscr_tcHigh(std::string outFileName = "", bool ownAlgo = false){

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  OpenFiles("~/Documents/MCtt_tcHigh_All2.root", 0);            //NoLepIso cut done
  OpenFiles("~/Documents/MCW_tcHigh_All2.root", 1);            //NoLepIso cut done
  OpenFiles("~/Documents/MCDY_tcHigh_All2.root", 2);          //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  TH1F* h_tcHigh[4];      //1 and 2 are b-jets, 3 and 4 are non-bjets (refParton used)
  for(int i = 0; i < 4; i++){
    if(i == 0 || i == 2) h_tcHigh[i] = new TH1F(Form("h_tcHigh_%d",i),"tcHighEff distribution;tcHighEff;Counts",120,-10,/*30*/50);
    else h_tcHigh[i] = new TH1F(Form("h_tcHigh_%d",i),"tcHighPur distribution;tcHighPur;Counts",120,-10,/*30*/50);
  } 	
  TH1F* h_svtx[4];        //1 and 2 are b-jets, 3 and 4 are non-bjets (refParton used)
  for(int i = 0; i < 4; i++){
    if(i == 0 || i == 2) h_svtx[i] = new TH1F(Form("h_svtx_%d",i),"svtxm distribution;svtxm;Counts",100,0,7);
    else h_svtx[i] = new TH1F(Form("h_svtx_%d",i),"svtxpt distribution;svtxpt;Counts",100,0,200);
  }
  TH1F* h_discr_prob[2];
  for(int i = 0; i < 2; i++) h_discr_prob[i] = new TH1F(Form("h_discr_prob_%d",i),"discr_prob distribution;discr_prob;Counts",100,0,3.5);


  int nPassedCSVCut010[nFiles][nDiscrCuts], nPassedCSVCut1030[nFiles][nDiscrCuts], nPassedCSVCut30100[nFiles][nDiscrCuts];
  int nDidntPassCSVCut010[nFiles][nDiscrCuts], nDidntPassCSVCut1030[nFiles][nDiscrCuts], nDidntPassCSVCut30100[nFiles][nDiscrCuts];
  double CSVcuts[nDiscrCuts];
  for(int i = 0; i < nFiles; i++){
    for(int j = 0; j < nDiscrCuts; j++){
      nPassedCSVCut010[i][j] = nPassedCSVCut1030[i][j] = nPassedCSVCut30100[i][j] = 0;
      nDidntPassCSVCut010[i][j] = nDidntPassCSVCut1030[i][j] = nDidntPassCSVCut30100[i][j] = 0;
    }
  }


  for(int i = 0; i < nFiles; i++){

    SetAddressBranches(i);

    for(int entry = 0; entry < (int)tr_[i]->GetEntries(); entry++){
      
      tr_[i]->GetEntry(entry);

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
    
      if((int)indexJets.size() < nJets) continue;

      if(i == 0 && ownAlgo == true){
        std::vector<int> indexMCbjets = ownFlavourMatchingAlgo(entry);
        if((int)indexMCbjets.size()>0 && indexMCbjets[0] == -999){ 
          cout << "Drop event........................" << endl; continue;
          indexMCbjets.clear(); //empty vector and test for Wjets and DY
        }
        
        for(int ij = 0; ij < (int)indexJets.size(); ij++){
        //for(int ij = 0; ij < nJt; ij++){
          
          bool inMCvector(false);
          for(int ijMC = 0; ijMC < (int)indexMCbjets.size(); ijMC++){
            if( indexJets[ij] == indexMCbjets[ijMC]) inMCvector = true;
          }
          if(inMCvector == true){
            h_tcHigh[0]->Fill(discr_tcHighEff[indexJets[ij]]);
            h_tcHigh[1]->Fill(discr_tcHighPur[indexJets[ij]]);
            h_svtx[0]->Fill(svtxm[indexJets[ij]]);
            h_svtx[1]->Fill(svtxpt[indexJets[ij]]);
            h_discr_prob[0]->Fill(discr_prob[indexJets[ij]]);
          } else {
            h_tcHigh[2]->Fill(discr_tcHighEff[indexJets[ij]]);
            h_tcHigh[3]->Fill(discr_tcHighPur[indexJets[ij]]);
            h_svtx[2]->Fill(svtxm[indexJets[ij]]);
            h_svtx[3]->Fill(svtxpt[indexJets[ij]]);
            h_discr_prob[1]->Fill(discr_prob[indexJets[ij]]);
          }
        }
      } else {
        if(i == 0 && hiBin >= 60){
          for(int ij = 0; ij < (int)indexJets.size()/*nJt*/; ij++){
            if(fabs(refparton_flavorForB[indexJets[ij]]) == 5 /*|| fabs(refparton_flavorForB[indexJets[ij]]) == 4*/){ 
              h_tcHigh[0]->Fill(discr_tcHighEff[indexJets[ij]]);
              h_tcHigh[1]->Fill(discr_tcHighPur[indexJets[ij]]);
              h_svtx[0]->Fill(svtxm[indexJets[ij]]);
              h_svtx[1]->Fill(svtxpt[indexJets[ij]]);
              h_discr_prob[0]->Fill(discr_prob[indexJets[ij]]);
            } else{
              h_tcHigh[2]->Fill(discr_tcHighEff[indexJets[ij]]);
              h_tcHigh[3]->Fill(discr_tcHighPur[indexJets[ij]]);
              h_svtx[2]->Fill(svtxm[indexJets[ij]]);
              h_svtx[3]->Fill(svtxpt[indexJets[ij]]);
              h_discr_prob[1]->Fill(discr_prob[indexJets[ij]]);
            }
          }
        }
      }

      double csvv1Cut = 0.;
      double tcHighCut = 30.; //(till some negative values)
      double svtxmCut = 0.; //(till 3)
      double svtxptCut = 0.; //(till 100)
      double probCut = -0.0001; //(till 2.5)
      for(int l = 0; l < nDiscrCuts; l++){
        
        double CSVv1CUTforROC = csvv1Cut + l * 0.05;
        double tcCUTforROC = tcHighCut - l * 0.25;
        double svtxmCUTforROC = svtxmCut + l * 0.33;
        double svtxptCUTforROC = svtxptCut + l * 3.;
        double probCUTforROC = probCut + l * 0.08;

        std::vector<int> indexbJets = FindbJetsCSVv1(indexJets, CSVv1CUTforROC);
        CSVcuts[l] = CSVv1CUTforROC;
        //std::vector<int> indexbJets = FindbJetstcHighPur(indexJets, tcCUTforROC);
        //CSVcuts[l] = tcCUTforROC;
        //std::vector<int> indexbJets = FindbJetstcHighEff(indexJets, tcCUTforROC);
        //CSVcuts[l] = tcCUTforROC;
        //std::vector<int> indexbJets = FindbJetssvtxm(indexJets, svtxmCUTforROC);
        //CSVcuts[l] = svtxmCUTforROC;
        //std::vector<int> indexbJets = FindbJetssvtxpt(indexJets, svtxptCUTforROC);
        //CSVcuts[l] = svtxptCUTforROC;
        //std::vector<int> indexbJets = FindbJetsprob(indexJets, probCUTforROC);
        //CSVcuts[l] = probCUTforROC;

        if(hiBin < 20){
          if((int)indexbJets.size() < nbJets) nDidntPassCSVCut010[i][l]++;
          else                                nPassedCSVCut010[i][l]++;
        } else if(hiBin >= 20 && hiBin < 60){
          if((int)indexbJets.size() < nbJets) nDidntPassCSVCut1030[i][l]++;
          else                                nPassedCSVCut1030[i][l]++;
        } else {
          if((int)indexbJets.size() < nbJets) nDidntPassCSVCut30100[i][l]++;
          else                                nPassedCSVCut30100[i][l]++;
        }
      }
    }
  }



/////////////////////////////////////  
//Plots 1
/////////////////////////////////////
  Color_t colors[4] = {2,2,4,4};
  for(int i = 0; i < 4; i++){
    //h_tcHigh[i]->Sumw2();
    //h_tcHigh[i]->Scale(1./h_tcHigh[i]->Integral() );
    //h_svtx[i]->Sumw2();
    //h_svtx[i]->Scale(1./h_svtx[i]->Integral() );
    h_tcHigh[i]->SetLineColor(colors[i]);
    h_tcHigh[i]->SetLineWidth(2);
    h_svtx[i]->SetLineColor(colors[i]);
    h_svtx[i]->SetLineWidth(2);
  }
  Color_t colors2[2] = {2,4};
  for(int i = 0; i < 2; i++){
    //h_discr_prob[i]->Sumw2();
    //h_discr_prob[i]->Scale(1./h_discr_prob[i]->Integral() );
    h_discr_prob[i]->SetLineWidth(2);
    h_discr_prob[i]->SetLineColor(colors2[i]);
  }

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(h_tcHigh[0],"MC b-jet","l");
  leg->AddEntry(h_tcHigh[2],"No MC b-jet","l");

  TCanvas* cst = new TCanvas("cst","tcHigh (non)b-jets",10,10,1400,500);
  cst->Divide(2);
  cst->cd(1);
  gStyle->SetOptStat(0); 
  //h_tcHigh[2]->Draw("HIST");
  //h_tcHigh[0]->Draw("same, HIST");
  h_svtx[2]->Draw("HIST");
  h_svtx[0]->Draw("same, HIST");
  leg->Draw();
  cst->cd(2);
  //h_tcHigh[3]->Draw("HIST");
  //h_tcHigh[1]->Draw("same, HIST");
  h_svtx[3]->Draw("HIST");
  h_svtx[1]->Draw("same, HIST");
  leg->Draw();
/*
  TCanvas* cst2 = new TCanvas("cst2","tcHigh (non)b-jets",10,10,1400,500);
  cst2->cd();
  gStyle->SetOptStat(0); 
  h_discr_prob[1]->Draw("HIST");
  h_discr_prob[0]->Draw("same, HIST");
  leg->Draw();
*/  

/////////////////////////////////////  
//Plots 2
/////////////////////////////////////
  double Normaliz[3] = {1, 0.71091, 0.084058};
  for(int i = 0; i < nFiles; i++){
    cout << "Going to normalize File " << i << " using " << Normaliz[i] << endl;
    for(int j = 0; j < nDiscrCuts; j++){
      nPassedCSVCut010[i][j] = Normaliz[i] * nPassedCSVCut010[i][j];
      nDidntPassCSVCut010[i][j] = Normaliz[i] * nDidntPassCSVCut010[i][j];
      nPassedCSVCut1030[i][j] = Normaliz[i] * nPassedCSVCut1030[i][j];
      nDidntPassCSVCut1030[i][j] = Normaliz[i] * nDidntPassCSVCut1030[i][j];
      nPassedCSVCut30100[i][j] = Normaliz[i] * nPassedCSVCut30100[i][j];
      nDidntPassCSVCut30100[i][j] = Normaliz[i] * nDidntPassCSVCut30100[i][j];
    }
  }

  double x010[nDiscrCuts], y010[nDiscrCuts], x1030[nDiscrCuts], y1030[nDiscrCuts], x30100[nDiscrCuts], y30100[nDiscrCuts];
  for(int i = 0; i < nFiles; i++){
    cout << "File: " << i << endl;
    for(int j = 0; j < nDiscrCuts; j++){

      cout << "  CSV cut : " << CSVcuts[j] << endl;
      cout << "     HiBin 0-10 : " << nPassedCSVCut010[i][j] << " +- " << sqrt(nPassedCSVCut010[i][j]) << " , " << nDidntPassCSVCut010[i][j] << " +- " << sqrt(nDidntPassCSVCut010[i][j])<< endl;
      cout << "     HiBin 10-30 : " << nPassedCSVCut1030[i][j] << " +- " << sqrt(nPassedCSVCut1030[i][j]) << " , " << nDidntPassCSVCut1030[i][j]<< " +- " << sqrt(nDidntPassCSVCut1030[i][j]) << endl;
      cout << "     HiBin 30-100 : " << nPassedCSVCut30100[i][j] << " +- " << sqrt(nPassedCSVCut30100[i][j]) << " , " << nDidntPassCSVCut30100[i][j] << " +- " << sqrt(nDidntPassCSVCut30100[i][j])<< endl;

      if(i==0){
        x010[j] = 100.* nPassedCSVCut010[i][j]/(nPassedCSVCut010[i][j] + nDidntPassCSVCut010[i][j]);
        x1030[j] = 100.* nPassedCSVCut1030[i][j]/(nPassedCSVCut1030[i][j] + nDidntPassCSVCut1030[i][j]);
        x30100[j] = 100.* nPassedCSVCut30100[i][j]/(nPassedCSVCut30100[i][j] + nDidntPassCSVCut30100[i][j]);
      }
      else if(i == 1){ 
        y010[j] = 100.* (nPassedCSVCut010[i][j] + nPassedCSVCut010[2][j])/(nPassedCSVCut010[i][j] + nDidntPassCSVCut010[i][j] + nPassedCSVCut010[2][j] + nDidntPassCSVCut010[2][j]);
        y1030[j] = 100.* (nPassedCSVCut1030[i][j] + nPassedCSVCut1030[2][j])/(nPassedCSVCut1030[i][j] + nDidntPassCSVCut1030[i][j] + nPassedCSVCut1030[2][j] + nDidntPassCSVCut1030[2][j]);
        y30100[j] = 100.* (nPassedCSVCut30100[i][j] + nPassedCSVCut30100[2][j])/(nPassedCSVCut30100[i][j] + nDidntPassCSVCut30100[i][j] + nPassedCSVCut30100[2][j] + nDidntPassCSVCut30100[2][j]);
        
        //y010[j] = 100.* (nPassedCSVCut010[i][j] + nPassedCSVCut010[2][j] + nPassedCSVCut010[3][j])/(nPassedCSVCut010[i][j] + nDidntPassCSVCut010[i][j] + nPassedCSVCut010[2][j] + nDidntPassCSVCut010[2][j] + nPassedCSVCut010[3][j] + nDidntPassCSVCut010[3][j]);
        //y1030[j] = 100.* (nPassedCSVCut1030[i][j] + nPassedCSVCut1030[2][j] + nPassedCSVCut1030[3][j])/(nPassedCSVCut1030[i][j] + nDidntPassCSVCut1030[i][j] + nPassedCSVCut1030[2][j] + nDidntPassCSVCut1030[2][j] + nPassedCSVCut1030[3][j] + nDidntPassCSVCut1030[3][j]);
        //y30100[j] = 100.* (nPassedCSVCut30100[i][j] + nPassedCSVCut30100[2][j] + nPassedCSVCut30100[3][j])/(nPassedCSVCut30100[i][j] + nDidntPassCSVCut30100[i][j] + nPassedCSVCut30100[2][j] + nDidntPassCSVCut30100[2][j] + nPassedCSVCut30100[3][j] + nDidntPassCSVCut30100[3][j]);

        y010[j] = 100. - y010[j];
        y1030[j] = 100. - y1030[j];
        y30100[j] = 100. - y30100[j];
      }
    }
  }

  TGraph* CSV_effcut010 = new TGraph(nDiscrCuts,x010,y010);
  TGraph* CSV_workp010 = new TGraph(1, &x010[nCut], &y010[nCut]);
  TGraph* CSV_effcut1030 = new TGraph(nDiscrCuts,x1030,y1030);
  TGraph* CSV_workp1030 = new TGraph(1, &x1030[nCut], &y1030[nCut]);
  TGraph* CSV_effcut30100 = new TGraph(nDiscrCuts,x30100,y30100);
  TGraph* CSV_workp30100 = new TGraph(1, &x30100[nCut], &y30100[nCut]);
  TGraph* CSV_workp30100leg = new TGraph(1, &x30100[nCut], &y30100[nCut]);

  //CSV_effcut010->SetNameTitle("g_CSV_effcut010","Efficiency tcHigh cut, Centrality 0-10;Sig eff; Back rejection eff (W + DY)");
  //CSV_effcut1030->SetNameTitle("g_CSV_effcut1030","Efficiency tcHigh cut, Centrality 10-30;Sig eff; Back rejection eff (W + DY)");
  //CSV_effcut30100->SetNameTitle("g_CSV_effcut30100","Efficiency tcHigh cut, Centrality 30-100;Sig eff; Back rejection eff (W + DY)");
  CSV_effcut010->SetNameTitle("g_CSV_effcut010","Efficiency tcHigh cut, Centrality 0-20;Sig eff; Back rejection eff (W + DY)");
  CSV_effcut1030->SetNameTitle("g_CSV_effcut1030","Efficiency tcHigh cut, Centrality 20-80;Sig eff; Back rejection eff (W + DY)");
  CSV_effcut30100->SetNameTitle("g_CSV_effcut30100","Efficiency tcHigh cut, Centrality 80-100;Sig eff; Back rejection eff (W + DY)");

  CSV_effcut010->SetMarkerStyle(20);
  CSV_workp010->SetMarkerStyle(29);
  CSV_effcut1030->SetMarkerStyle(21);
  CSV_workp1030->SetMarkerStyle(29);
  CSV_effcut30100->SetMarkerStyle(22);
  CSV_workp30100->SetMarkerStyle(29);
  CSV_workp30100leg->SetMarkerStyle(29);

  CSV_effcut010->SetMarkerColor(8);
  CSV_effcut010->SetLineColor(8);
  CSV_effcut010->SetLineWidth(2);
  CSV_effcut010->SetMarkerSize(1.25);
  CSV_workp010->SetMarkerColor(8);
  CSV_workp010->SetMarkerSize(3.5);
  CSV_effcut1030->SetMarkerColor(kBlue);
  CSV_effcut1030->SetLineColor(kBlue);
  CSV_effcut1030->SetLineWidth(2);
  CSV_effcut1030->SetMarkerSize(1.25);
  CSV_workp1030->SetMarkerColor(kBlue);
  CSV_workp1030->SetMarkerSize(3.5);
  CSV_effcut30100->SetMarkerColor(kRed);
  CSV_effcut30100->SetLineColor(kRed);
  CSV_effcut30100->SetLineWidth(2);
  CSV_effcut30100->SetMarkerSize(1.25);
  CSV_workp30100->SetMarkerColor(kRed);
  CSV_workp30100->SetMarkerSize(3.5);
  CSV_workp30100leg->SetMarkerSize(2.5);

  TMultiGraph *mg = new TMultiGraph();
  //mg->SetNameTitle("mg","ROC: tcHighPur Cut;Sig eff; Back rejection eff (W + DY)");
  //mg->SetNameTitle("mg","ROC: tcHighEff Cut;Sig eff; Back rejection eff (W + DY)");
  //mg->SetNameTitle("mg","ROC: svtxm Cut;Sig eff; Back rejection eff (W + DY)");
  //mg->SetNameTitle("mg","ROC: svtxpt Cut;Sig eff; Back rejection eff (W + DY)");
  //mg->SetNameTitle("mg","Sig vs Back event selection: discr_prob;Sig event eff; Back event rejection eff (W + DY)");
  //mg->SetNameTitle("mg","Sig vs Back event selection: svtxpt;Sig event eff; Back event rejection eff (W + DY)");
  mg->SetNameTitle("mg","Sig vs Back event selection: CSVv1;Sig event eff; Back event rejection eff (W + DY)");
  mg->Add(CSV_effcut010, "lp");
  mg->Add(CSV_effcut1030, "lp");
  mg->Add(CSV_effcut30100, "lp");
  mg->Add(CSV_workp010,"p");
  mg->Add(CSV_workp1030,"p");
  mg->Add(CSV_workp30100,"p");

  TLegend* leg2 = new TLegend(0.7,0.7,0.9,0.9);
  leg2->AddEntry(CSV_effcut010,"PbPb 0-10%","lp");
  leg2->AddEntry(CSV_effcut1030,"PbPb 10-30%","lp");
  leg2->AddEntry(CSV_effcut30100,"PbPb 30-100%","lp");
  leg2->AddEntry(CSV_workp30100leg,Form("Cutvalue = %.2f",CSVcuts[nCut]),"p");

  TCanvas *cst2 = new TCanvas("cst2","Histograms",10,10,1800,1000);//1800,1000);
  cst2->cd();
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(0,100);
  mg->GetYaxis()->SetRangeUser(0,100);
  leg2->Draw();

  cout << endl << endl << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  cst->Write();
  cst2->Write();
  outFile->Write();
  delete outFile;
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
  tr_[i]->SetBranchAddress("svtxm", svtxm); 
  tr_[i]->SetBranchAddress("svtxpt", svtxpt);
  tr_[i]->SetBranchAddress("discr_prob", discr_prob); 
  tr_[i]->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);

  tr_[i]->SetBranchAddress("nGen", &nGen);
  tr_[i]->SetBranchAddress("genPdg", genPdg);
  tr_[i]->SetBranchAddress("genPt", genPt);
  tr_[i]->SetBranchAddress("genPhi", genPhi);
  tr_[i]->SetBranchAddress("genEta", genEta);
  tr_[i]->SetBranchAddress("genChg", genChg);

}

std::vector<int> FindbJetsCSVv1(std::vector<int> indexJets, double cutValue){

  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(discr_csvV1[indexJets[ij]] > cutValue)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}

std::vector<int> FindbJetstcHighPur(std::vector<int> indexJets, double cutValue){

  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(discr_tcHighPur[indexJets[ij]] > cutValue)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}

std::vector<int> FindbJetstcHighEff(std::vector<int> indexJets, double cutValue){

  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(discr_tcHighEff[indexJets[ij]] > cutValue)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}

std::vector<int> FindbJetssvtxm(std::vector<int> indexJets, double cutValue){

  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(svtxm[indexJets[ij]] > cutValue)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}

std::vector<int> FindbJetssvtxpt(std::vector<int> indexJets, double cutValue){

  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(svtxpt[indexJets[ij]] > cutValue)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}


std::vector<int> FindbJetsprob(std::vector<int> indexJets, double cutValue){

  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(discr_prob[indexJets[ij]] > cutValue)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}

std::vector<int> ownFlavourMatchingAlgo(int entry){

/////////////////////////////////////  
//Drop strange events without b-bbar or t-ttbar pair
/////////////////////////////////////
  bool bquark(false), tquark(false);
  std::vector<int> indexb, indexbbar;
  std::vector<int> returnvector;
  for(int igen = 0; igen < nGen; igen++){
    if( genPdg[igen] == 5){ bquark = true; indexb.push_back(igen); }
    if( genPdg[igen] == 6){ tquark = true; }
  }
  if(bquark == false || tquark == false){
    cout << "Event " << entry << ". No b-quark = " << std::boolalpha << bquark << " or t-quark " << tquark << endl;
    returnvector.push_back(-999);
    return returnvector;
  }

  bool bbarquark(false), tbarquark(false);
  for(int igen = 0; igen < nGen; igen++){
    if( genPdg[igen] == -5){ bbarquark = true; indexbbar.push_back(igen); }
    if( genPdg[igen] == -6){ tbarquark = true; }
  }
  if(bbarquark == false || tbarquark == false){
    cout << "Event " << entry << ". No bbar-quark = " << std::boolalpha << bquark << " or tbar-quark " << tquark << endl;
    returnvector.push_back(-999);
    return returnvector;
  }
/////////////////////////////////////


/////////////////////////////////////
//Check if there are more b or bbar quarks
/////////////////////////////////////
  bool possible2b = false;
  int i = (int)indexb.size() - 2;
  if((int)indexb.size() > 1){
    if(genPt[indexb[i]] - genPt[indexb[i+1]] < -10.)              possible2b = true;
    if(TMath::Abs(genPhi[indexb[i+1]] - genPhi[indexb[i]]) > 0.5) possible2b = true;
    if(TMath::Abs(genEta[indexb[i+1]] - genEta[indexb[i]]) > 0.5) possible2b = true;
  }

  bool possible2bbar = false;
  i = (int)indexbbar.size() - 2;
  if((int)indexbbar.size() > 1){
    if(genPt[indexbbar[i]] - genPt[indexbbar[i+1]] < -10.)              possible2bbar = true;
    if(TMath::Abs(genPhi[indexbbar[i+1]] - genPhi[indexbbar[i]]) > 0.5) possible2bbar = true;
    if(TMath::Abs(genEta[indexbbar[i+1]] - genEta[indexbbar[i]]) > 0.5) possible2bbar = true;
  }
/////////////////////////////////////


/////////////////////////////////////
//Select the final quarks (after gluon scattering)
/////////////////////////////////////
  int bquark1, bquark2, bbarquark1, bbarquark2;
  if(possible2b == true){  bquark1 = indexb[(int)indexb.size() - 1];              bquark2 = indexb[(int)indexb.size() - 2]; } 
  else {                   bquark1 = indexb[(int)indexb.size() - 1]; }
  if(possible2bbar == true){  bbarquark1 = indexbbar[(int)indexbbar.size() - 1];  bbarquark2 = indexbbar[(int)indexbbar.size() - 2]; }
  else {                      bbarquark1 = indexbbar[(int)indexbbar.size() - 1]; }
/////////////////////////////////////

/////////////////////////////////////
//Matching jet to b-quarks
/////////////////////////////////////
  int bjet1(-999), bjet2(-999), bbarjet1(-999), bbarjet2(-999);
  double temp_dR = 999.;
  for(int ijet = 0; ijet < nJt; ijet++){

    double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bquark1] - jtPhi[ijet]) );
    double deltaeta = TMath::Abs( genEta[bquark1] - jtEta[ijet] );
    double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bquark1]-jtEta[ijet],2) );

    if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;

    if(drJetToGen < temp_dR){ bjet1 = ijet; temp_dR = drJetToGen; }
  }
  temp_dR = 999.;
  if(possible2b == true){
    for(int ijet = 0; ijet < nJt; ijet++){

      if(ijet == bjet1) continue;

      double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bquark2] - jtPhi[ijet]) );
      double deltaeta = TMath::Abs( genEta[bquark2] - jtEta[ijet] );
      double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bquark2]-jtEta[ijet],2) );

      if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;

      if(drJetToGen < temp_dR){ bjet2 = ijet; temp_dR = drJetToGen; }
    }
  }
  temp_dR = 999.;
  for(int ijet = 0; ijet < nJt; ijet++){

    if(ijet == bjet1 || ijet == bjet2) continue;

    double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bbarquark1] - jtPhi[ijet]) );
    double deltaeta = TMath::Abs( genEta[bbarquark1] - jtEta[ijet] );
    double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bbarquark1]-jtEta[ijet],2) );

    if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;
      
    if(drJetToGen < temp_dR){ bbarjet1 = ijet; temp_dR = drJetToGen; }
  }
  temp_dR = 999.;
  if(possible2bbar == true){
    for(int ijet = 0; ijet < nJt; ijet++){

      if(ijet == bjet1 || ijet == bjet2 || ijet == bbarjet1) continue;

      double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bbarquark2] - jtPhi[ijet]) );
      double deltaeta = TMath::Abs( genEta[bbarquark2] - jtEta[ijet] );
      double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bbarquark2]-jtEta[ijet],2) );

      if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;
       
      if(drJetToGen < temp_dR){ bbarjet2 = ijet; temp_dR = drJetToGen; }
    }
  }
/////////////////////////////////////

/////////////////////////////////////
//Building MC b-jet vector
/////////////////////////////////////
  std::vector<int> indexMCbjets;
  if(bjet1 != -999)     indexMCbjets.push_back(bjet1);
  if(bbarjet1 != -999)  indexMCbjets.push_back(bbarjet1);
  if(possible2b == true && bjet2 != -999)       indexMCbjets.push_back(bjet2);
  if(possible2bbar == true && bbarjet2 != -999) indexMCbjets.push_back(bbarjet2);

  return indexMCbjets;
}