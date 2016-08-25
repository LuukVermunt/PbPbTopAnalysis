#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "anaMuJetsSkimTree.h"
#include "histControlDistributions.h"

#include <string>
#include <vector>

TCanvas* PlotEfficiencyHighInterval();
TCanvas* PlotHistograms(int opt);

TH2F* ABCD_plot_[4][6];
int hiBinLow[6] = {0,20,60,100,140,0};
int hiBinHigh[6] = {20,60,100,140,200,200};
double lepIsoCut[6] = {0.58, 0.45, 0.3, 0.24, 0.18, 0.58};
double lepIsoCutHigh[6] = {2, 2, 2, 2, 2, 2};
double lepPtCut[6] = {18., 18., 18., 18., 18., 18.};

void plotABCD_HighInterval(std::string outFileName = "", int nJets = 4, int nbJets = 1){

  anaMuJetsSkimTree ana;
  //ana.rand->SetSeed(0);

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }

  ana.OpenFiles("~/Documents/MCtt_tcHigh_All2.root", 1);          //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCW_tcHigh_All2.root", 2);           //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCDY_tcHigh_All2.root", 3);          //NoLepIso cut done
  ana.OpenFiles("~/Documents/DataOld_tcHigh_All2.root", 4);          //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  for(int j = 0; j < 4; j++){
    for(int i = 0; i < 6; i++){
      ABCD_plot_[j][i] = new TH2F(Form("ABCD_plot_%d_%d",j,i),Form("ABCD (hiBin: %d - %d) 4j%db;lepIso;pT",hiBinLow[i], hiBinHigh[i],nbJets),8001,-0.01,80,3940,15,1000);
    }
  }
  ana.BuildHistograms();

  cout << endl << "Making ABCD plot for " << nJets << "j" << nbJets << "b." << endl;
  cout << "Selecting " << nbJets << " b-jet with CSVv1 > " << ana.GetCSVCut() << endl;

  for(int j = 1; j < 5; j++){
    ana.SetAddressBranches(j);

    int eventsPassedCuts = 0;

    for(int entry = 0; entry < (int)ana.tr_[j]->GetEntries(); entry++){
      ana.tr_[j]->GetEntry(entry);
      int indexMuon;
      double Ptlepfirst = 0.;
      if(ana.GetnLep() > 1){
        for(int ilep = 0; ilep<ana.GetnLep(); ilep++){
          if(ana.GetlepPt(ilep) >= Ptlepfirst){
            Ptlepfirst = ana.GetlepPt(ilep);
            indexMuon = ilep;
          }
        }
      } 
      else  indexMuon = 0;

      //Store the jets that pass the basic cuts
      std::vector<int> indexJets; 
      if(j != 4){ ana.JetSmearing = true; indexJets = ana.FindJetsAfterCuts(indexMuon); }
      else{       ana.JetSmearing = false; indexJets = ana.FindJetsAfterCuts(indexMuon); }

      std::vector<int> indexbJets = ana.FindbJets(indexJets);

      if((int)indexJets.size() < nJets || (int)indexbJets.size() < nbJets) continue;
        
      for(int i = 0; i < 6; i++){
        if(ana.GethiBin() >= hiBinLow[i] && ana.GethiBin() < hiBinHigh[i]){
          ABCD_plot_[j-1][i]->Fill( ana.GetlepIso(indexMuon) , ana.GetlepPt(indexMuon));

          if(i != 5 && ana.GetlepIso(indexMuon) <= lepIsoCut[i] && ana.GetlepPt(indexMuon) >= 18.)  eventsPassedCuts++;
        }  
      }
    }

    cout << " File " << j << " , Events in area C = " << eventsPassedCuts << endl;
  }

  TCanvas* cst = PlotEfficiencyHighInterval();

  //TCanvas* cst = PlotHistograms(3);
  //TCanvas* cst2 = PlotHistograms(1);
  //TCanvas* cst99 = PlotHistograms(99);

  cout << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  //gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  gSystem->cd("~/Documents");
  cst->Write();
  //cst2->Write();
  //cst99->Write();
  //outFile->Write();
  delete outFile;
}

TCanvas* PlotEfficiencyHighInterval(){

  double lepIsoHigh1stInt = 0.58;
  int binx1, binx2, binx2high, binx3, biny1, biny2, biny3;
  double higherbound[14], EWoverData[14], errorEWoverData[14], errorLepIso[14];

  for(int i = 1; i < 15; i++){

    double lepIsoHigh = lepIsoHigh1stInt + i * 0.1;
    higherbound[i-1] = lepIsoHigh;

    binx1 = ABCD_plot_[3][0]->GetXaxis()->GetFirst();
    binx2 = ABCD_plot_[3][0]->GetXaxis()->FindBin( lepIsoCut[0] );
    binx2high = ABCD_plot_[3][0]->GetXaxis()->FindBin( lepIsoHigh );
    binx3 = ABCD_plot_[3][0]->GetXaxis()->GetLast();
    biny1 = ABCD_plot_[3][0]->GetYaxis()->GetFirst();
    biny2 = ABCD_plot_[3][0]->GetYaxis()->FindBin( lepPtCut[i] );
    biny3 = ABCD_plot_[3][0]->GetYaxis()->GetLast();

    double nData = ABCD_plot_[3][0]->Integral(binx2, binx2high, biny1, biny3);
    double ntt = /*0.00431057 * */ABCD_plot_[0][0]->Integral(binx2, binx2high, biny1, biny3);
    double nWjets = /*0.75377 * */ABCD_plot_[1][0]->Integral(binx2, binx2high, biny1, biny3);
    double nDY = /*0.0891257 * */ABCD_plot_[2][0]->Integral(binx2, binx2high, biny1, biny3);

    EWoverData[i-1] = (0.00431057 * ntt + 0.75377 * nWjets + 0.0891257 * nDY)/nData;
    //EWoverData[i-1] = (ntt + nWjets + nDY)/nData;
    errorEWoverData[i-1] = EWoverData[i-1] * sqrt(1./nData + (0.00431057 * 0.00431057 * ntt + 0.75377 * 0.75377 * nWjets + 0.0891257 * 0.0891257 * nDY)/((0.00431057 * ntt + 0.75377 * nWjets + 0.0891257 * nDY)*(0.00431057 * ntt + 0.75377 * nWjets + 0.0891257 * nDY)) );
    //errorEWoverData[i-1] = EWoverData[i-1] * sqrt(1./nData + (ntt + nWjets + nDY)/( (ntt + nWjets + nDY)*(ntt + nWjets + nDY) ) );
    errorLepIso[i-1] = 0;
    cout << " " << i << ": Bound = " << higherbound[i-1] << " , EW/Data = " << EWoverData[i-1] << " +- " << errorEWoverData[i-1] << endl;
  }

  TGraphErrors* g_hist1 = new TGraphErrors(14,higherbound,EWoverData,errorLepIso,errorEWoverData);
  g_hist1->SetNameTitle("g_hist1",Form("EW / Data in sideband (PbPb %d-%d%%);lepIso; EW/Data",hiBinLow[0]/2,hiBinHigh[0]/2));
  g_hist1->SetMarkerStyle(20);
  g_hist1->SetMarkerColor(2);
  g_hist1->SetLineColor(2);
  
  lepIsoHigh1stInt = 0.45;
  for(int i = 1; i < 15; i++){

    double lepIsoHigh = lepIsoHigh1stInt + i * 0.14;
    higherbound[i-1] = lepIsoHigh;

    binx1 = ABCD_plot_[3][1]->GetXaxis()->GetFirst();
    binx2 = ABCD_plot_[3][1]->GetXaxis()->FindBin( lepIsoCut[1] );
    binx2high = ABCD_plot_[3][1]->GetXaxis()->FindBin( lepIsoHigh );
    binx3 = ABCD_plot_[3][1]->GetXaxis()->GetLast();
    biny1 = ABCD_plot_[3][1]->GetYaxis()->GetFirst();
    biny2 = ABCD_plot_[3][1]->GetYaxis()->FindBin( lepPtCut[i] );
    biny3 = ABCD_plot_[3][1]->GetYaxis()->GetLast();

    double nData = ABCD_plot_[3][1]->Integral(binx2, binx2high, biny1, biny3);
    double ntt = 0.00431057 * ABCD_plot_[0][1]->Integral(binx2, binx2high, biny1, biny3);
    double nWjets = 0.75377 * ABCD_plot_[1][1]->Integral(binx2, binx2high, biny1, biny3);
    double nDY = 0.0891257 * ABCD_plot_[2][1]->Integral(binx2, binx2high, biny1, biny3);

    EWoverData[i-1] = (ntt + nWjets + nDY)/nData;
    cout << " " << i << ": Bound = " << higherbound[i-1] << " , EW/Data = " << EWoverData[i-1] << endl;
  }

  TGraph* g_hist2 = new TGraph(14,higherbound,EWoverData);
  g_hist2->SetNameTitle("g_hist1",Form("EW / Data in sideband (PbPb %d-%d%%);lepIso; EW/Data",hiBinLow[1]/2,hiBinHigh[1]/2));
  g_hist2->SetMarkerStyle(20);
  g_hist2->SetMarkerColor(4);
  g_hist2->SetLineColor(4);

  //TCanvas *cst = new TCanvas("cst","Histograms",10,10,1400,500);
  /*cst->Divide(2);
  cst->cd(1);
  g_hist1->Draw("ALP");
  cst->cd(2);
  g_hist2->Draw("ALP");*/
  TCanvas *cst = new TCanvas("cst","Histograms",10,10,700,500);
  cst->Divide();
  cst->cd(1);
  g_hist1->Draw("ALP");
  
  return cst;
}




TCanvas* PlotHistograms(int opt){

  int nEntriesMJ[4], nEntries[4], nEntriestt[4], nEntriesW[4], nEntriesDY[4], nSumEvents[5];
  int binx1, binx2, binx3, biny1, biny2, biny3;

  TCanvas *cst = new TCanvas(Form("cst_%d",opt),"Histograms",10,10,1800,1000);
  cst->Divide(3,2);

  bool MinusWjets = false;
  if(opt == 99){
    MinusWjets = true;
    opt = 3;
  }

  TLine *line1[6];
  TLine *line2[6];

  TText *t = new TText();
  t->SetTextSize(0.03);

  for(int i = 0; i < 5/*6*/; i++){
    cst->cd(i+1);
    line1[i] = new TLine(lepIsoCut[i],15,lepIsoCut[i],100);
    line2[i] = new TLine(0,lepPtCut[i],6,lepPtCut[i]);
    line1[i]->SetLineColor(kRed);
    line2[i]->SetLineColor(kRed);
    ABCD_plot_[opt][i]->Draw("COLZ");
    line1[i]->Draw();
    line2[i]->Draw();
    t->DrawText(0.04,15.5,"A");
    t->DrawText(2.5,15.5,"B");
    t->DrawText(0.04,55,"C");
    t->DrawText(2.5,55,"D");

    binx1 = ABCD_plot_[opt][i]->GetXaxis()->GetFirst();
    binx2 = ABCD_plot_[opt][i]->GetXaxis()->FindBin( lepIsoCut[i] );
    binx3 = ABCD_plot_[opt][i]->GetXaxis()->GetLast();
    biny1 = ABCD_plot_[opt][i]->GetYaxis()->GetFirst();
    biny2 = ABCD_plot_[opt][i]->GetYaxis()->FindBin( lepPtCut[i] );
    biny3 = ABCD_plot_[opt][i]->GetYaxis()->GetLast();

    if(MinusWjets == true){
      nEntries[0] = ABCD_plot_[opt][i]->Integral(binx1, binx2, biny1, biny2);
      nEntries[1] = ABCD_plot_[opt][i]->Integral(binx2, binx3, biny1, biny2);
      nEntries[2] = ABCD_plot_[opt][i]->Integral(binx1, binx2, biny2, biny3);
      nEntries[3] = ABCD_plot_[opt][i]->Integral(binx2, binx3, biny2, biny3);
      nEntriestt[0] = 0.005349285 * ABCD_plot_[0][i]->Integral(binx1, binx2, biny1, biny2);
      nEntriestt[1] = 0.005349285 * ABCD_plot_[0][i]->Integral(binx2, binx3, biny1, biny2);
      nEntriestt[2] = 0.005349285 * ABCD_plot_[0][i]->Integral(binx1, binx2, biny2, biny3);
      nEntriestt[3] = 0.005349285 * ABCD_plot_[0][i]->Integral(binx2, binx3, biny2, biny3);
      nEntriesW[0] = 0.71091 * ABCD_plot_[1][i]->Integral(binx1, binx2, biny1, biny2);
      nEntriesW[1] = 0.71091 * ABCD_plot_[1][i]->Integral(binx2, binx3, biny1, biny2);
      nEntriesW[2] = 0.71091 * ABCD_plot_[1][i]->Integral(binx1, binx2, biny2, biny3);
      nEntriesW[3] = 0.71091 * ABCD_plot_[1][i]->Integral(binx2, binx3, biny2, biny3);
      nEntriesDY[0] = 0.084058 * ABCD_plot_[2][i]->Integral(binx1, binx2, biny1, biny2);
      nEntriesDY[1] = 0.084058 * ABCD_plot_[2][i]->Integral(binx2, binx3, biny1, biny2);
      nEntriesDY[2] = 0.084058 * ABCD_plot_[2][i]->Integral(binx1, binx2, biny2, biny3);
      nEntriesDY[3] = 0.084058 * ABCD_plot_[2][i]->Integral(binx2, binx3, biny2, biny3);

      for(int j = 0; j < 4; j++){
        //nEntriesMJ[j] = nEntries[j] - nEntriesW[j] - nEntriesDY[j] - nEntriestt[j];
        nEntriesMJ[j] = nEntries[j];
      }

      if(i == 0){
        cout << "The expected events for Multijets in 0.404 nb^(-1) is: " << endl;
        nSumEvents[0] = (double)nEntriesMJ[0] * nEntriesMJ[3] / (double)nEntriesMJ[1];
        cout << "   hiBin 0-20: " << nSumEvents[0] << endl;
      } else if(i == 1){ 
        nSumEvents[1] = (double)nEntriesMJ[0] * nEntriesMJ[3] / (double)nEntriesMJ[1];
        cout << "   hiBin 20-60: " << nSumEvents[1] << endl;
      } else if(i == 2){ 
        nSumEvents[2] = (double)nEntriesMJ[0] * nEntriesMJ[3] / (double)nEntriesMJ[1];
        cout << "   hiBin 60-100: " << nSumEvents[2] << endl;
      } else if(i == 3){ 
        nSumEvents[3] = (double)nEntriesMJ[0] * nEntriesMJ[3] / (double)nEntriesMJ[1];
        cout << "   hiBin 100-140: " << nSumEvents[3] << endl;
      } else if(i == 4){ 
        nSumEvents[4] = (double)nEntriesMJ[0] * nEntriesMJ[3] / (double)nEntriesMJ[1];
        cout << "   hiBin 140-200: " << nSumEvents[4] << endl;
      }

      t->DrawText(0.75,90,Form("lepPt Cut: >=%.0f, lepIso Cut: <=%.2f", lepPtCut[i], lepIsoCut[i]));
      t->DrawText(0.75,85,Form("Data:     A: %d, B: %d, C: %d, D: %d", nEntries[0], nEntries[1], nEntries[2], nEntries[3]));
      t->DrawText(0.75,82,Form("ttbar:    A: %d, B: %d, C: %d, D: %d", nEntriestt[0], nEntriestt[1], nEntriestt[2], nEntriestt[3]));
      t->DrawText(0.75,79,Form("Wjets:    A: %d, B: %d, C: %d, D: %d", nEntriesW[0], nEntriesW[1], nEntriesW[2], nEntriesW[3]));
      t->DrawText(0.75,76,Form("DY:       A: %d, B: %d, C: %d, D: %d", nEntriesDY[0], nEntriesDY[1], nEntriesDY[2], nEntriesDY[3]));
      //t->DrawText(0.75,70,Form("Data - ttbar - Wjets - DY:   A: %d, B: %d, C: %d, D: %d", nEntriesMJ[0], nEntriesMJ[1], nEntriesMJ[2], nEntriesMJ[3]));
      t->DrawText(0.75,64,Form("MJ: N_c = N_a N_d/N_b = %f", (double)nEntriesMJ[0] * nEntriesMJ[3] / (double)nEntriesMJ[1]));
      t->DrawText(0.75,61,Form("N_c - N'_c = %f",(double)nEntriesMJ[0] * nEntriesMJ[3] /( (double)nEntriesMJ[1]) - nEntriestt[2] - nEntriesW[2] - nEntriesDY[2] ));
    } else {
      nEntries[0] = ABCD_plot_[opt][i]->Integral(binx1, binx2, biny1, biny2);
      nEntries[1] = ABCD_plot_[opt][i]->Integral(binx2, binx3, biny1, biny2);
      nEntries[2] = ABCD_plot_[opt][i]->Integral(binx1, binx2, biny2, biny3);
      nEntries[3] = ABCD_plot_[opt][i]->Integral(binx2, binx3, biny2, biny3);

      t->DrawText(1,70,Form("Entries plot, area A: %d, B: %d, C: %d, D: %d", nEntries[0], nEntries[1], nEntries[2], nEntries[3]));
      t->DrawText(1,67,Form("N_c = N_a N_d/N_b = %f", (double)nEntries[0] * nEntries[3] / (double)nEntries[1]));
    }

    ABCD_plot_[opt][i]->GetXaxis()->SetRangeUser(0,6);
    ABCD_plot_[opt][i]->GetYaxis()->SetRangeUser(0,100);
  }

  int nTotal = 0;
  for(int i = 0; i < 4; i++){
    if(nSumEvents[i] > 0) nTotal += nSumEvents[i];
  }

  if(MinusWjets == true){
    cout << endl << endl << "Total number of expected Multijets events is : " << nTotal << endl << endl;
  }
 
  cst->Update();
  return cst;
}
