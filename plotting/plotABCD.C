#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "anaMuJetsSkimTree.h"
#include "histControlDistributions.h"

#include <string>
#include <vector>

TCanvas* PlotHistograms(int opt);

TH2F* ABCD_plot_[4][6];
int hiBinLow[6] = {0,20,60,100,140,0};
int hiBinHigh[6] = {20,60,100,140,200,200};
double lepIsoCut[6] = {0.58, 0.45, 0.3, 0.24, 0.18, 0.58};
double lepIsoCutHigh[6] = {1.1, 1.1, 1.1, 1.1, 1.1, 1.1};
//double lepIsoCutHigh[6] = {0.58, 0.45, 0.3, 0.24, 0.18, 0.58};
double lepPtCut[6] = {18., 18., 18., 18., 18., 18.};

void plotABCD(std::string outFileName = "", int nJets = 4, int nbJets = 1){

  anaMuJetsSkimTree ana;
  //ana.rand->SetSeed(0);

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }

  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_W/MCW.root", 2);            //NoLepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_DY/MCDY.root", 3);          //NoLepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/DataMultijets.root", 4);	//NoLepIso cut done
  //ana.OpenFiles("~/Documents/MCtt.root", 1);            //NoLepIso cut done
  //ana.OpenFiles("~/Documents/MCW.root", 2);            //NoLepIso cut done
  //ana.OpenFiles("~/Documents/MCDY.root", 3);          //NoLepIso cut done
  //ana.OpenFiles("~/Documents/DataMultijets.root", 4);  //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCtt_tcHigh_All2.root", 1);            //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCW_tcHigh_All2.root", 2);            //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCDY_tcHigh_All2.root", 3);          //NoLepIso cut done
  ana.OpenFiles("~/Documents/DataOld_tcHigh_All2.root", 4);  //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  for(int j = 0; j < 4; j++){
    for(int i = 0; i < 6; i++){
      ABCD_plot_[j][i] = new TH2F(Form("ABCD_plot_%d_%d",j,i),Form("ABCD 4j%db (PbPb %d-%d%%);muon Isolation;muon p_{T}",nbJets,hiBinLow[i]/2, hiBinHigh[i]/2),8001,-0.01,80,3940,15,1000);
    }
  }
  ana.BuildHistograms();

  cout << endl << "Making ABCD plot for " << nJets << "j" << nbJets << "b." << endl;
  cout << "Selecting " << nbJets << " b-jet with CSVv1 > " << ana.GetCSVCut() << endl;

  for(int j = 1; j < 5; j++){
    ana.SetAddressBranches(j);

    int eventsPassedCuts = 0;

    int t1(0), t2(0);

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
      t1++;
      if((int)indexJets.size() < nJets || (int)indexbJets.size() < nbJets) continue;
      t2++;
      for(int i = 0; i < 6; i++){
        if(ana.GethiBin() >= hiBinLow[i] && ana.GethiBin() < hiBinHigh[i]){
          ABCD_plot_[j-1][i]->Fill( ana.GetlepIso(indexMuon) , ana.GetlepPt(indexMuon));

          if(i != 5 && ana.GetlepIso(indexMuon) <= lepIsoCut[i] && ana.GetlepPt(indexMuon) >= 18.)  eventsPassedCuts++;
        }  
      }
    }

    cout << " File " << j << " , Events in area C = " << eventsPassedCuts << endl;
    cout << "   t1 = " << t1 << " , t2 = " << t2 << endl;
  }

  //TCanvas* cst = PlotHistograms(3);
  //TCanvas* cst2 = PlotHistograms(1);
  TCanvas* cst99 = PlotHistograms(99);

  cout << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  //gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  gSystem->cd("~/Documents");
  //cst->Write();
  //cst2->Write();
  cst99->Write();
  outFile->Write();
  delete outFile;
}

TCanvas* PlotHistograms(int opt){

  int nEntriesMJ[4], nEntries[4], nEntriestt[4], nEntriesW[4], nEntriesDY[4], nSumEvents[5];
  int binx1, binx2, binx22, binx3, biny1, biny2, biny3;

  TCanvas *cst = new TCanvas(Form("cst_%d",opt),"Histograms",10,10,1800,1000);
  cst->Divide(3,2);

  bool MinusWjets = false;
  if(opt == 99){
    MinusWjets = true;
    opt = 3;
  }

  TLine *line1[6];
  TLine *line2[6];
  TLine *line22[6];
  TLine *line3[6];

  TText *t = new TText();
  t->SetTextSize(0.03);

  for(int i = 0; i < 5/*6*/; i++){
    cst->cd(i+1);
    line1[i] = new TLine(lepIsoCut[i],15,lepIsoCut[i],100);
    line2[i] = new TLine(0,lepPtCut[i],lepIsoCut[i],lepPtCut[i]);
    line22[i] = new TLine(lepIsoCutHigh[i],lepPtCut[i],6,lepPtCut[i]);
    line3[i] = new TLine(lepIsoCutHigh[i],15,lepIsoCutHigh[i],100);
    line1[i]->SetLineColor(kRed);
    line2[i]->SetLineColor(kRed);
    line22[i]->SetLineColor(kRed);
    line3[i]->SetLineColor(kRed);
    line1[i]->SetLineWidth(2);
    line2[i]->SetLineWidth(2);
    line22[i]->SetLineWidth(2);
    line3[i]->SetLineWidth(2);
    ABCD_plot_[opt][i]->Draw("COLZ");
    line1[i]->Draw();
    line2[i]->Draw();
    line22[i]->Draw();
    line3[i]->Draw();
    t->DrawText(0.04,15.5,"A");
    t->DrawText(2.5,15.5,"B");
    t->DrawText(0.04,55,"C");
    t->DrawText(2.5,55,"D");

    binx1 = ABCD_plot_[opt][i]->GetXaxis()->GetFirst();
    binx2 = ABCD_plot_[opt][i]->GetXaxis()->FindBin( lepIsoCut[i] );
    binx22 = ABCD_plot_[opt][i]->GetXaxis()->FindBin( lepIsoCutHigh[i] );
    binx3 = ABCD_plot_[opt][i]->GetXaxis()->GetLast();
    biny1 = ABCD_plot_[opt][i]->GetYaxis()->GetFirst();
    biny2 = ABCD_plot_[opt][i]->GetYaxis()->FindBin( lepPtCut[i] );
    biny3 = ABCD_plot_[opt][i]->GetYaxis()->GetLast();

    if(MinusWjets == true){
      nEntries[0] = ABCD_plot_[opt][i]->Integral(binx1, binx2, biny1, biny2);
      nEntries[1] = ABCD_plot_[opt][i]->Integral(binx22, binx3, biny1, biny2);
      nEntries[2] = ABCD_plot_[opt][i]->Integral(binx1, binx2, biny2, biny3);
      nEntries[3] = ABCD_plot_[opt][i]->Integral(binx22, binx3, biny2, biny3);
      nEntriestt[0] = 0.00431057 * ABCD_plot_[0][i]->Integral(binx1, binx2, biny1, biny2);
      nEntriestt[1] = 0.00431057 * ABCD_plot_[0][i]->Integral(binx22, binx3, biny1, biny2);
      nEntriestt[2] = /*0.00431057 */ ABCD_plot_[0][i]->Integral(binx1, binx2, biny2, biny3);
      nEntriestt[3] = 0.00431057 * ABCD_plot_[0][i]->Integral(binx22, binx3, biny2, biny3);
      nEntriesW[0] = 0.75377 * ABCD_plot_[1][i]->Integral(binx1, binx2, biny1, biny2);
      nEntriesW[1] = 0.75377 * ABCD_plot_[1][i]->Integral(binx22, binx3, biny1, biny2);
      nEntriesW[2] = /*0.75377 */ ABCD_plot_[1][i]->Integral(binx1, binx2, biny2, biny3);
      nEntriesW[3] = 0.75377 * ABCD_plot_[1][i]->Integral(binx22, binx3, biny2, biny3);
      nEntriesDY[0] = 0.0891257 * ABCD_plot_[2][i]->Integral(binx1, binx2, biny1, biny2);
      nEntriesDY[1] = 0.0891257 * ABCD_plot_[2][i]->Integral(binx22, binx3, biny1, biny2);
      nEntriesDY[2] = /*0.0891257 */ ABCD_plot_[2][i]->Integral(binx1, binx2, biny2, biny3);
      nEntriesDY[3] = 0.0891257 * ABCD_plot_[2][i]->Integral(binx22, binx3, biny2, biny3);

      for(int j = 0; j < 4; j++){
        nEntriesMJ[j] = nEntries[j] - nEntriesW[j] - nEntriesDY[j] - nEntriestt[j];
        //nEntriesMJ[j] = nEntries[j];
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
      //t->DrawText(0.75,61,Form("N_c - N'_c = %f",(double)nEntriesMJ[0] * nEntriesMJ[3] /( (double)nEntriesMJ[1]) - nEntriestt[2] - nEntriesW[2] - nEntriesDY[2] ));
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
  

  line1[0]->SetY2(70);// = new TLine(lepIsoCut[0],15,lepIsoCut[0],70);
  //line2[0] = new TLine(0,lepPtCut[0],lepIsoCut[0],lepPtCut[0]);
  line22[0]->SetX2(4);// = new TLine(lepIsoCutHigh[0],lepPtCut[0],4,lepPtCut[0]);
  line3[0]->SetY2(70);// = new TLine(lepIsoCutHigh[0],15,lepIsoCutHigh[0],70);
  t->SetTextSize(0.05);
  t->SetTextColor(kRed);
  TCanvas *cst2 = new TCanvas("cst_solo","Histograms",10,10,800,500);
  cst2->cd();
  ABCD_plot_[3][0]->Draw("COLZ");
  line1[0]->Draw();
  line2[0]->Draw();
  line22[0]->Draw();
  line3[0]->Draw();
  t->DrawText(0.25,15.49,"A");
  t->DrawText(2.5,15.49,"B");
  t->DrawText(0.25,55,"C");
  t->DrawText(2.5,55,"D");  
  ABCD_plot_[3][0]->GetXaxis()->SetRangeUser(0,4);
  ABCD_plot_[3][0]->GetYaxis()->SetRangeUser(0,70);

  cst2->Update();

  return cst2;
}
