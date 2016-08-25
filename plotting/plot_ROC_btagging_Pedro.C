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
//const int nDiscrCuts = 11; int nCut = 1; //CSV
//const int nDiscrCuts = 20; int nCut = 1; //tcHigh
//const int nDiscrCuts = 8; int nCut = 1; //svtxm
//const int nDiscrCuts = 11; int nCut = 1; //svtxpt
const int nDiscrCuts = 20; int nCut = 15; //discr_prob
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



void plot_ROC_btagging_Pedro(std::string outFileName = "", bool ownAlgo = false){

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  OpenFiles("~/Documents/MCtt_tcHigh_All2.root", 0);            //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  TH1F* h_tcHigh[4];      //1 and 2 are b-jets, 3 and 4 are non-bjets (refParton used)
  for(int i = 0; i < 4; i++){
    if(i == 0 || i == 2) h_tcHigh[i] = new TH1F(Form("h_tcHigh_%d",i),Form("tcHighEff distribution (PbPb %d-%d);tcHighEff;Counts",hiBinLow/2,hiBinHigh/2),120,-10,/*30*/50);
    else h_tcHigh[i] = new TH1F(Form("h_tcHigh_%d",i),Form("tcHighPur distribution (PbPb %d-%d);tcHighPur;Counts",hiBinLow/2,hiBinHigh/2),120,-10,/*30*/50);
  } 	
  TH1F* h_svtx[4];        //1 and 2 are b-jets, 3 and 4 are non-bjets (refParton used)
  for(int i = 0; i < 4; i++){
    if(i == 0 || i == 2) h_svtx[i] = new TH1F(Form("h_svtx_%d",i),Form("svtxm distribution (PbPb %d-%d);svtxm;Counts",hiBinLow/2,hiBinHigh/2),100,0,7);
    else h_svtx[i] = new TH1F(Form("h_svtx_%d",i),Form("svtxpt distribution (PbPb %d-%d);svtxpt;Counts",hiBinLow/2,hiBinHigh/2),100,0,200);
  }
  TH1F* h_discr_prob[2];
  TH1F* h_discr_CSVv1[2];
  TH1F* h_discr_CSVv2[2];
  for(int i = 0; i < 2; i++) h_discr_prob[i] = new TH1F(Form("h_discr_prob_%d",i),Form("discr_prob distribution (PbPb %d-%d);discr_prob;Counts",hiBinLow/2,hiBinHigh/2),100,0,3.5);
  for(int i = 0; i < 2; i++) h_discr_CSVv1[i] = new TH1F(Form("h_discr_CSVv1_%d",i),Form("discr_CSVv1 distribution (PbPb %d-%d);discr_CSVv1;Counts",hiBinLow/2,hiBinHigh/2),20,0,1);
  for(int i = 0; i < 2; i++) h_discr_CSVv2[i] = new TH1F(Form("h_discr_CSVv2_%d",i),Form("discr_CSVv2 distribution (PbPb %d-%d);discr_CSVv2;Counts",hiBinLow/2,hiBinHigh/2),20,0,1);


  double CSVcuts[nDiscrCuts];
  int hiBinint[nDiscrCuts], hiBin010[nDiscrCuts], hiBin1030[nDiscrCuts], hiBin30100[nDiscrCuts];
  int hiBinintnotag[nDiscrCuts], hiBin010notag[nDiscrCuts], hiBin1030notag[nDiscrCuts], hiBin30100notag[nDiscrCuts];
  int hiBinintfake[nDiscrCuts], hiBin010fake[nDiscrCuts], hiBin1030fake[nDiscrCuts], hiBin30100fake[nDiscrCuts];
  int hiBinintlightgood[nDiscrCuts], hiBin010lightgood[nDiscrCuts], hiBin1030lightgood[nDiscrCuts], hiBin30100lightgood[nDiscrCuts];

  for(int i = 0; i < nDiscrCuts; i++){
    hiBinint[i] = hiBin010[i] = hiBin1030[i] = hiBin30100[i] = 0;
    hiBinintnotag[i] = hiBin010notag[i] = hiBin1030notag[i] = hiBin30100notag[i] = 0;
    hiBinintfake[i] = hiBin010fake[i] = hiBin1030fake[i] = hiBin30100fake[i] = 0;
    hiBinintlightgood[i] = hiBin010lightgood[i] = hiBin1030lightgood[i] = hiBin30100lightgood[i] = 0;
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

      if(hiBin >= hiBinLow && hiBin <= hiBinHigh){
        for(int ij = 0; ij < (int)indexJets.size(); ij++){
          if(fabs(refparton_flavorForB[indexJets[ij]]) == 5){ 
            h_tcHigh[0]->Fill(discr_tcHighEff[indexJets[ij]]);
            h_tcHigh[1]->Fill(discr_tcHighPur[indexJets[ij]]);
            h_svtx[0]->Fill(svtxm[indexJets[ij]]);
            h_svtx[1]->Fill(svtxpt[indexJets[ij]]);
            h_discr_prob[0]->Fill(discr_prob[indexJets[ij]]);
            h_discr_CSVv1[0]->Fill(discr_csvV1[indexJets[ij]]);
            h_discr_CSVv2[0]->Fill(discr_csvV2[indexJets[ij]]);
          } else{
            h_tcHigh[2]->Fill(discr_tcHighEff[indexJets[ij]]);
            h_tcHigh[3]->Fill(discr_tcHighPur[indexJets[ij]]);
            h_svtx[2]->Fill(svtxm[indexJets[ij]]);
            h_svtx[3]->Fill(svtxpt[indexJets[ij]]);
            h_discr_prob[1]->Fill(discr_prob[indexJets[ij]]);
            h_discr_CSVv1[1]->Fill(discr_csvV1[indexJets[ij]]);
            h_discr_CSVv2[1]->Fill(discr_csvV2[indexJets[ij]]);
          }
        }
      }

      double csvCut = 0.; //(till 1)
      double tcHighCut = 30.; //(till some negative values)
      double svtxmCut = 0.; //(till 3)
      double svtxptCut = 0.; //(till 100)
      double probCut = -0.0001; //(till 2.5)
      for(int l = 0; l < nDiscrCuts; l++){
        
        double CSVcut = csvCut + l * 0.05;
        double tcCUTforROC = tcHighCut - l * 2.5;
        double svtxmCUTforROC = svtxmCut + l * 0.33;
        double svtxptCUTforROC = svtxptCut + l * 3;
        double probCUTforROC = probCut + l * 0.08;//0.12;
        
        CSVcuts[l] = CSVcut;
        //CSVcuts[l] = tcCUTforROC;
        //CSVcuts[l] = svtxmCUTforROC;
        //CSVcuts[l] = svtxptCUTforROC;
        //CSVcuts[l] = probCUTforROC;

				for(int ij = 0; ij < (int)indexJets.size(); ij++){

          if(discr_csvV1[indexJets[ij]] > CSVcut){
          //if(discr_csvV2[indexJets[ij]] > CSVcut){
          //if(discr_tcHighPur[indexJets[ij]] > tcCUTforROC){
          //if(discr_tcHighEff[indexJets[ij]] > tcCUTforROC){
          //if(svtxm[indexJets[ij]] > svtxmCUTforROC){
          //if(svtxpt[indexJets[ij]] > svtxptCUTforROC){
          //if(discr_prob[indexJets[ij]] > probCUTforROC){
 
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
    //h_discr_CSVv1[i]->Sumw2();
    //h_discr_prob[i]->Scale(1./h_discr_prob[i]->Integral() );
    //h_discr_CSVv1[i]->Scale(1./h_discr_CSVv1[i]->Integral() );
    h_discr_prob[i]->SetLineWidth(2);
    h_discr_CSVv1[i]->SetLineWidth(2);
    h_discr_CSVv2[i]->SetLineWidth(2);
    h_discr_prob[i]->SetLineColor(colors2[i]);
    h_discr_CSVv1[i]->SetLineColor(colors2[i]);
    h_discr_CSVv2[i]->SetLineColor(colors2[i]);
  }

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(h_tcHigh[0],"MC b-jet","l");
  leg->AddEntry(h_tcHigh[2],"MC light-jet","l");

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

  TCanvas* cst2 = new TCanvas("cst2","tcHigh (non)b-jets",10,10,700,500);
  cst2->cd();
  gStyle->SetOptStat(0); 
  //h_discr_prob[1]->Draw("HIST");
  //h_discr_prob[0]->Draw("same, HIST");
  h_svtx[2]->Draw("HIST");
  h_svtx[0]->Draw("same, HIST");
  h_tcHigh[2]->Draw("HIST");
  h_tcHigh[0]->Draw("same, HIST");
  //h_svtx[3]->Draw("HIST");
  //h_svtx[1]->Draw("same, HIST");
  //h_discr_CSVv1[1]->Draw("HIST");
  //h_discr_CSVv1[0]->Draw("same, HIST");
  //h_discr_CSVv2[1]->Draw("HIST");
  //h_discr_CSVv2[0]->Draw("same, HIST");
  leg->Draw();


/////////////////////////////////////  
//Plots 2
/////////////////////////////////////
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

  TGraph* g_hiBinint = new TGraph(nDiscrCuts,btageffint,lighteffint);
  TGraph* g_workpint = new TGraph(1, &btageffint[nCut], &lighteffint[nCut]);
  g_hiBinint->SetNameTitle("g_hiBinint","Performance of b-tagging: CSVv1 (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  //g_hiBinint->SetNameTitle("g_hiBinint","Performance of b-tagging: CSVv2 (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  //g_hiBinint->SetNameTitle("g_hiBinint","Performance of b-tagging: discr_prob (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  //g_hiBinint->SetNameTitle("g_hiBinint","Performance of b-tagging: svtxm (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  //g_hiBinint->SetNameTitle("g_hiBinint","Performance of b-tagging: tcHighEff (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  //g_hiBinint->SetNameTitle("g_hiBinint","Performance of b-tagging: svtxpt (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  g_hiBinint->SetMarkerStyle(20);
  g_workpint->SetMarkerStyle(29);
  g_hiBinint->SetMarkerColor(8);
  g_hiBinint->SetLineColor(8);
  g_workpint->SetMarkerColor(8);
  g_workpint->SetMarkerSize(3.5);

  TGraph* g_hiBin010 = new TGraph(nDiscrCuts,btageff010,lighteff010);
  TGraph* g_workp010 = new TGraph(1, &btageff010[nCut], &lighteff010[nCut]);
  g_hiBin010->SetNameTitle("g_hiBin010","Performance of b-tagging (PbPb 0-10%);b-jet tagging efficiency; light-jet rejection efficiency");
  g_hiBin010->SetMarkerStyle(20);
  g_workp010->SetMarkerStyle(29);
  g_hiBin010->SetMarkerColor(8);
  g_hiBin010->SetLineColor(8);
  g_hiBin010->SetLineWidth(2);
  g_hiBin010->SetMarkerSize(1.25);
  g_workp010->SetMarkerColor(8);
  g_workp010->SetMarkerSize(3.5);
  TGraph* g_hiBin1030 = new TGraph(nDiscrCuts,btageff1030,lighteff1030);
  TGraph* g_workp1030 = new TGraph(1, &btageff1030[nCut], &lighteff1030[nCut]);
  g_hiBin1030->SetNameTitle("g_hiBin1030","Performance of b-tagging (PbPb 10-30%);b-jet tagging efficiency; light-jet rejection efficiency");
  g_hiBin1030->SetMarkerStyle(21);
  g_workp1030->SetMarkerStyle(29);
  g_hiBin1030->SetMarkerColor(kBlue);
  g_hiBin1030->SetLineColor(kBlue);
  g_hiBin1030->SetLineWidth(2);
  g_hiBin1030->SetMarkerSize(1.25);
  g_workp1030->SetMarkerColor(kBlue);
  g_workp1030->SetMarkerSize(3.5);
  TGraph* g_hiBin30100 = new TGraph(nDiscrCuts,btageff30100,lighteff30100);
  TGraph* g_workp30100 = new TGraph(1, &btageff30100[nCut], &lighteff30100[nCut]);
  g_hiBin30100->SetNameTitle("g_hiBin30100","Performance of b-tagging (PbPb 30-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  g_hiBin30100->SetMarkerStyle(22);
  g_workp30100->SetMarkerStyle(29);
  g_hiBin30100->SetMarkerColor(kRed);
  g_hiBin30100->SetLineColor(kRed);
  g_hiBin30100->SetLineWidth(2);
  g_hiBin30100->SetMarkerSize(1.25);
  g_workp30100->SetMarkerColor(kRed);
  g_workp30100->SetMarkerSize(3.5);
  TGraph* g_workp30100leg = new TGraph(1, &btageff30100[nCut], &lighteff30100[nCut]);
  g_workp30100leg->SetMarkerStyle(29);
  g_workp30100leg->SetMarkerSize(2.5);

  TMultiGraph *mg = new TMultiGraph();
  //mg->SetNameTitle("mg","Performance of b-tagging: discr_prob;b-jet tagging efficiency; light-jet rejection efficiency");
  //mg->SetNameTitle("mg","Performance of b-tagging: svtxpt;b-jet tagging efficiency; light-jet rejection efficiency");
  mg->SetNameTitle("mg","Performance of b-tagging: CSVv1;b-jet tagging efficiency; light-jet rejection efficiency");
  mg->Add(g_hiBin010, "lp");
  mg->Add(g_hiBin1030, "lp");
  mg->Add(g_hiBin30100, "lp");
  mg->Add(g_workp010,"p");
  mg->Add(g_workp1030,"p");
  mg->Add(g_workp30100,"p");

  TLegend* leg99 = new TLegend(0.7,0.7,0.9,0.9);
  leg99->AddEntry(g_hiBin010,"PbPb 0-20%","lp");//"PbPb 0-10%","p");
  leg99->AddEntry(g_hiBin1030,"PbPb 20-80%","lp");//PbPb 10-30%","p");
  leg99->AddEntry(g_hiBin30100,"PbPb 80-100%","lp");//PbPb 30-100%","p");
  leg99->AddEntry(g_workp30100leg,Form("CSV = %.2f",CSVcuts[nCut]),"p");

  TCanvas *cst99 = new TCanvas("cst99","Histograms",10,10,1800,1000);
  cst99->cd();
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(0,1);
  mg->GetYaxis()->SetRangeUser(0,1);
  leg99->Draw();

	TCanvas *cst98 = new TCanvas("cst98","Histograms",10,10,700,500);
  cst98->cd();
  g_hiBinint->Draw("ALP");
  g_hiBinint->GetXaxis()->SetLimits(0,1);
  g_hiBinint->GetYaxis()->SetRangeUser(0,1);


  cout << endl << endl << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  cst->Write();
  cst2->Write();
  cst98->Write();
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