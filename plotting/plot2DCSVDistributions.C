#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "anaMuJetsSkimTree.h"
#include "histControlDistributions.h"

#include <string>
#include <vector>

void LayoutHistograms();
TCanvas* PlotHistograms(int nbJets, double usedCSVcut);


void plot2DCSVDistributions(std::string outFileName = "", int nJets = 4, int nbJets = 0){

  anaMuJetsSkimTree ana;

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/Data.root",0);      //LepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_tt/MCtt.root", 1);      //NoLepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_W/MCW.root", 2);      //NoLepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_DY/MCDY.root", 3);      //NoLepIso cut done
  ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/DataMultijets.root", 4);  //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  double usedCSVcut = ana.GetCSVCut();

  //cout << "Making Control Distribution plots in the following pT interval." << endl;
  //cout << "     ptLow = " << pTLow << " , pTHigh = " << pTHigh << endl;
  cout << endl << "Making 2D CSV-Distribution plots for " << nJets << "j" << nbJets << "b." << endl;
  if(nbJets == 0)  cout << "   So we are not placing a cut on the CSVv1 variable" << endl << endl;
  else if(nbJets == 1) cout << "   We require 1 bjet by placing a cut on CSVv1 > " << usedCSVcut << endl << endl;
  else cout << "   We require 2 bjets by placing a cut on CSVv1 > " << usedCSVcut << endl << endl;
    

  ana.BuildHistograms();
  for(int i = 0; i < 5; i++)  ana.FillHistogramsBeforeCuts(i);
  LayoutHistograms();
  TCanvas* cst = PlotHistograms(nbJets, usedCSVcut);
  
  cout << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  outFile->Write();
  cst->Write();
  delete outFile;
}


void LayoutHistograms(){

  h_2dCSV_[0]->SetTitle("CSV of 1st vs 2nd b-jet (Data)");
  h_2dCSV_[1]->SetTitle("CSV of 1st vs 2nd b-jet (MC tt)");
  h_2dCSV_[2]->SetTitle("CSV of 1st vs 2nd b-jet (MC W)");
  h_2dCSV_[3]->SetTitle("CSV of 1st vs 2nd b-jet (MC DY)");
  h_2dCSV_[4]->SetTitle("CSV of 1st vs 2nd b-jet (Data Multijets)");

  h_2dCSV2_[0]->SetTitle("CSV of 1st vs 3rd b-jet (Data)");
  h_2dCSV2_[1]->SetTitle("CSV of 1st vs 3rd b-jet (MC tt)");
  h_2dCSV2_[2]->SetTitle("CSV of 1st vs 3rd b-jet (MC W)");
  h_2dCSV2_[3]->SetTitle("CSV of 1st vs 3rd b-jet (MC DY)");
  h_2dCSV2_[4]->SetTitle("CSV of 1st vs 3rd b-jet (Data Multijets)");
}

TCanvas* PlotHistograms(int nbJets, double usedCSVcut){

  TCanvas *cst = new TCanvas("cst","Histograms",10,10,1800,1000);
  TLine *line1 = new TLine(usedCSVcut,-0.1,usedCSVcut,1);
  TLine *line2 = new TLine(-0.1,usedCSVcut,1,usedCSVcut);
  line1->SetLineColor(kRed);
  line2->SetLineColor(kRed);

  TText *t = new TText();
  t->SetTextSize(0.03);

  int binx1, binx2, biny1, biny2;
  double perc_cutCSV[5], perc_cutCSV2[5];
  int TotEntries[5], TotEntries2[5];
  if(nbJets != 0){
    if(nbJets == 1){
      binx1 = h_2dCSV_[0]->GetXaxis()->FindBin( usedCSVcut );
      binx2 = h_2dCSV_[0]->GetXaxis()->FindBin(1);
      biny1 = h_2dCSV_[0]->GetYaxis()->GetFirst();
      biny2 = h_2dCSV_[0]->GetYaxis()->FindBin(1);
    } else {
      binx1 = h_2dCSV_[0]->GetXaxis()->FindBin( usedCSVcut );
      binx2 = h_2dCSV_[0]->GetXaxis()->FindBin(1);
      biny1 = h_2dCSV_[0]->GetYaxis()->FindBin( usedCSVcut );
      biny2 = h_2dCSV_[0]->GetYaxis()->FindBin(1);
    }

    for(int i = 0; i < 5; i++){
      TotEntries[i] = (int)h_2dCSV_[i]->GetEntries();
      TotEntries2[i] = (int)h_2dCSV2_[i]->GetEntries();
      perc_cutCSV[i] = 100.*h_2dCSV_[i]->Integral(binx1, binx2, biny1, biny2)/h_2dCSV_[i]->GetEntries();
      perc_cutCSV2[i] = 100.*h_2dCSV2_[i]->Integral(binx1, binx2, biny1, biny2)/h_2dCSV2_[i]->GetEntries();
    }
  }

  cst->Divide(5,2);
  gStyle->SetOptStat(0);
  for(int i = 0; i < 5; i++){
    cst->cd(i+1);
    h_2dCSV_[i]->Draw("COLZ");
    if(nbJets != 0){
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV[i]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries[i]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV_[i]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      if(nbJets != 1) line2->Draw();
    }
  }
  for(int i = 0; i < 5; i++){
    cst->cd(i+6);
    h_2dCSV2_[i]->Draw("COLZ");
    if(nbJets != 0){
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV2[i]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries2[i]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV2_[i]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      if(nbJets != 1) line2->Draw();
    }
  }

  cst->Update();
  return cst;
}
