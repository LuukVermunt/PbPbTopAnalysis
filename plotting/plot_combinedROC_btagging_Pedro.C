#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>

void OpenFiles(TString FileName, int i);
void SetAddressBranches(int i);

int nJets = 4;
int nbJets = 1;
int hiBinLow = 0;
int hiBinHigh = 200;
const int nDiscrCuts = 15; int nCut = 1;
const int nFiles = 1;
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



void plot_combinedROC_btagging_Pedro(std::string outFileName = ""){

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  OpenFiles("~/Documents/MCtt_tcHigh_All2.root", 0);            //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  TGraph* g_hiBinint[4];
  Color_t colors[4] = {1,2,4,6};
  Style_t styles[4] = {20,21,22,23};

  double CSVcuts[nDiscrCuts];
  int hiBinint[nDiscrCuts], hiBin010[nDiscrCuts], hiBin1030[nDiscrCuts], hiBin30100[nDiscrCuts];
  int hiBinintnotag[nDiscrCuts], hiBin010notag[nDiscrCuts], hiBin1030notag[nDiscrCuts], hiBin30100notag[nDiscrCuts];
  int hiBinintfake[nDiscrCuts], hiBin010fake[nDiscrCuts], hiBin1030fake[nDiscrCuts], hiBin30100fake[nDiscrCuts];
  int hiBinintlightgood[nDiscrCuts], hiBin010lightgood[nDiscrCuts], hiBin1030lightgood[nDiscrCuts], hiBin30100lightgood[nDiscrCuts];

  for(int jj = 0; jj < 4; jj++){

    for(int i = 0; i < nDiscrCuts; i++){
      hiBinint[i] = hiBin010[i] = hiBin1030[i] = hiBin30100[i] = 0;
      hiBinintnotag[i] = hiBin010notag[i] = hiBin1030notag[i] = hiBin30100notag[i] = 0;
      hiBinintfake[i] = hiBin010fake[i] = hiBin1030fake[i] = hiBin30100fake[i] = 0;
      hiBinintlightgood[i] = hiBin010lightgood[i] = hiBin1030lightgood[i] = hiBin30100lightgood[i] = 0;
    }

    SetAddressBranches(0);

    for(int entry = 0; entry < (int)tr_[0]->GetEntries(); entry++){
      
      tr_[0]->GetEntry(entry);

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

      double csvCut = 0.; //(till 1)
      double tcHighCut = -5.; //(till some negative values)
      double svtxmCut = 0.; //(till 3)
      double svtxptCut = 0.; //(till 100)
      double probCut = -0.0001; //(till 2.5)
      for(int l = 0; l < nDiscrCuts; l++){
        
        double CSVcut = csvCut + l * 0.1;
        double probCUTforROC = probCut + l * 0.08;//0.12;
        double svtxptCUTforROC = svtxptCut + l * 3;
        double tcCUTforROC = tcHighCut + l * 3;
        
        if(jj == 0)      CSVcuts[l] = CSVcut;
        else if(jj == 1) CSVcuts[l] = probCUTforROC;
        else if(jj == 2) CSVcuts[l] = svtxptCUTforROC;
        else            CSVcuts[l] = tcCUTforROC;
        
				for(int ij = 0; ij < (int)indexJets.size(); ij++){

          double value;
          if(jj == 0)      value = discr_csvV1[indexJets[ij]];
          else if(jj == 1) value = discr_prob[indexJets[ij]];
          else if(jj == 2) value = svtxpt[indexJets[ij]];
          else            value = discr_tcHighPur[indexJets[ij]];

          if(value > CSVcuts[l]){
          
          	if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBinint[l]++;
            else                             											hiBinintfake[l]++;

            if(hiBin < 40/*20*/){ 
              if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBin010[l]++;
              else                             											hiBin010fake[l]++;
            }
            else if(hiBin >= 40/*20*/ && hiBin < 160/*60*/){
              if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBin1030[l]++;
              else                             											hiBin1030fake[l]++;
            }
            else{
              if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBin30100[l]++;
              else						                            					hiBin30100fake[l]++;
            }

          } else {
			
          	if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBinintnotag[l]++;
            else                             											hiBinintlightgood[l]++;

						if(hiBin < 40/*20*/){ 
              if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBin010notag[l]++;
              else                             											hiBin010lightgood[l]++;
            }
            else if(hiBin >= 40/*20*/ && hiBin < 160/*60*/){
              if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBin1030notag[l]++;
              else                             											hiBin1030lightgood[l]++;
            }
            else{
              if( fabs(refparton_flavorForB[indexJets[ij]]) == 5 )  hiBin30100notag[l]++;
              else						                         						  hiBin30100lightgood[l]++;
            }
          }
        }
      }
    }

    double btageffint[nDiscrCuts], lighteffint[nDiscrCuts];
    double btageff010[nDiscrCuts], lighteff010[nDiscrCuts];
    double btageff1030[nDiscrCuts], lighteff1030[nDiscrCuts];
    double btageff30100[nDiscrCuts], lighteff30100[nDiscrCuts];
    for(int i = 0; i < nDiscrCuts; i++){
      btageffint[i] = (double) hiBinint[i] / (hiBinintnotag[i] + hiBinint[i]);
      btageff010[i] = (double) hiBin010[i] / (hiBin010notag[i] + hiBin010[i]);
      btageff1030[i] = (double) hiBin1030[i] / (hiBin1030notag[i] + hiBin1030[i]);
      btageff30100[i] = (double) hiBin30100[i] / (hiBin30100notag[i] + hiBin30100[i]);

      lighteffint[i] = 1. - (double) hiBinintfake[i] / (hiBinintlightgood[i] + hiBinintfake[i]);
      lighteff010[i] = 1. - (double) hiBin010fake[i] / (hiBin010lightgood[i] + hiBin010fake[i]);
      lighteff1030[i] = 1. - (double) hiBin1030fake[i] / (hiBin1030lightgood[i] + hiBin1030fake[i]);
      lighteff30100[i] = 1. - (double) hiBin30100fake[i] / (hiBin30100lightgood[i] + hiBin30100fake[i]);
    }

    g_hiBinint[jj] = new TGraph(nDiscrCuts,btageffint,lighteffint);
    g_hiBinint[jj]->SetNameTitle("g_hiBinint","Performance of b-tagging (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
    g_hiBinint[jj]->SetMarkerStyle(styles[jj]);
    g_hiBinint[jj]->SetMarkerColor(colors[jj]);
    g_hiBinint[jj]->SetLineColor(colors[jj]);
    g_hiBinint[jj]->SetLineWidth(2);
    g_hiBinint[jj]->SetMarkerSize(1.25);
  }


  TMultiGraph *mg = new TMultiGraph();
  mg->SetNameTitle("mg","Performance of b-tagging (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  mg->Add(g_hiBinint[0], "lp");
  mg->Add(g_hiBinint[1], "lp");
  mg->Add(g_hiBinint[2], "lp");
  mg->Add(g_hiBinint[3], "lp");
  
  TLegend* leg99 = new TLegend(0.7,0.7,0.9,0.9);
  leg99->AddEntry(g_hiBinint[0],"CSVv1","lp");
  leg99->AddEntry(g_hiBinint[1],"discr_prob","lp");
  leg99->AddEntry(g_hiBinint[2],"svtxpt","lp");
  leg99->AddEntry(g_hiBinint[3],"tcHighEff","lp");
  
  TCanvas *cst99 = new TCanvas("cst99","Histograms",10,10,1800,1000);
  cst99->cd();
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(0,1);
  mg->GetYaxis()->SetRangeUser(0,1);
  leg99->Draw();

  cout << endl << endl << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  cst99->Write();
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