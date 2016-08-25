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
TCanvas* PlotHistograms(int nbJets, int nJets);

void plotControlDistributions(std::string outFileName = "", int nJets = 4, int nbJets = 1){

  anaMuJetsSkimTree ana;
  ana.rand->SetSeed(0);

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/Data.root",0);			//LepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_tt/MCtt.root", 1);			//NoLepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_W/MCW.root", 2);			//NoLepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_DY/MCDY.root", 3);			//NoLepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/DataMultijets.root", 4);	//NoLepIso cut done
  //ana.OpenFiles("~/Documents/Data.root", 0);            //LepIso cut done
  //ana.OpenFiles("~/Documents/MCtt_tcHigh_All2.root", 1);          //NoLepIso cut done
  //ana.OpenFiles("~/Documents/MCW_tcHigh_All2.root", 2);            //NoLepIso cut done
  //ana.OpenFiles("~/Documents/MCDY_tcHigh_All2.root", 3);          //NoLepIso cut done
  //ana.OpenFiles("~/Documents/DataMultijets.root", 4);  //NoLepIso cut done
  ana.OpenFiles("~/Documents/DataOld_tcHigh_All2.root", 0);            //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCtt_tcHigh_All2.root", 1);          //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCW_tcHigh_All2.root", 2);            //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCDY_tcHigh_All2.root", 3);          //NoLepIso cut done
  ana.OpenFiles("~/Documents/DataOld_tcHigh_All2.root", 4);  //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  //cout << "Making Control Distribution plots in the following pT interval." << endl;
  //cout << "     ptLow = " << pTLow << " , pTHigh = " << pTHigh << endl;
  cout << endl << "Making Control Distribution plots for " << nJets << "j" << nbJets << "b." << endl;
  ana.SetnumbB(nbJets);
  ana.SetnumbJ(nJets);

  if(nbJets == 0){

  	cout << "   So we are not placing a cut on the CSVv1 variable" << endl << endl;
  	ana.BuildHistograms();
  	for(int i = 0; i < 5; i++)  ana.FillHistogramsBeforeCuts(i);

  } else {
  	
  	cout << "   We are placing a cut on the CSVv1 variable to tag b-jets" << endl << endl;
  	ana.BuildHistograms();
  	for(int i = 0; i < 5; i++)  ana.FillHistogramsAfterCuts(i);

  }

  LayoutHistograms();
  for(int i = 0; i < 5; i++)  ana.NormalizeHistograms(i, nbJets);
  //for(int i = 0; i < 4; i++) NormalizeHistogramsToOne(i, nbJets);
  TCanvas* cst = PlotHistograms(nbJets, nJets);

  cout << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  cst->Write();
  outFile->Write();
  delete outFile;
}




void LayoutHistograms(){

  Color_t colors[4] = {920,800,42,32};

  h_lepPt_[0]->SetMarkerStyle(20);
  h_lepPseu_[0]->SetMarkerStyle(20);
  h_jetHt_[0]->SetMarkerStyle(20);
  h_jetCSV_[0]->SetMarkerStyle(20);
  h_lepbjetMinv_min_[0]->SetMarkerStyle(20);
  h_ttMinv_min_[0]->SetMarkerStyle(20);
  h_qqMinv_min_[0]->SetMarkerStyle(20);
  h_pTselectbjet_[0]->SetMarkerStyle(20);
  h_Phi_1stb_[0]->SetMarkerStyle(20);
  h_Phi_allb_[0]->SetMarkerStyle(20);
  h_Phi_b_minpi_[0]->SetMarkerStyle(20);

  for(int i = 1; i < 5; i++){
  	h_lepPt_[i]->SetFillColor(colors[i-1]);
    h_lepPseu_[i]->SetFillColor(colors[i-1]);
    h_jetHt_[i]->SetFillColor(colors[i-1]);
    h_jetCSV_[i]->SetFillColor(colors[i-1]);
    h_lepbjetMinv_min_[i]->SetFillColor(colors[i-1]);
    h_ttMinv_min_[i]->SetFillColor(colors[i-1]);
    h_qqMinv_min_[i]->SetFillColor(colors[i-1]);
    h_pTselectbjet_[i]->SetFillColor(colors[i-1]);
    h_Phi_1stb_[i]->SetFillColor(colors[i-1]);
    h_Phi_allb_[i]->SetFillColor(colors[i-1]);
    h_Phi_b_minpi_[i]->SetFillColor(colors[i-1]);
  }
}


TCanvas* PlotHistograms(int nbJets, int nJets){

  //TCanvas *cst = new TCanvas("cst","Histograms",10,10,1800,1000);
  TCanvas *cst = new TCanvas("cst","Histograms",10,10,1800,500);
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  TLegend* leg2 = new TLegend(0.12,0.68,0.32,0.88);

  lepPt_Stack = new THStack("lepPt_Stack",Form("p_{T}(muon) distributions (%dj%db);p_{T};Counts",nJets, nbJets));
  lepPseu_Stack = new THStack("lepPseu_Stack",Form("|#eta|(muon) distributions (%dj%db);|#eta|;Counts",nJets, nbJets));
  jetHt_Stack = new THStack("jetHt_Stack",Form("H_{T} distributions (%dj%db);H_{T};Counts",nJets, nbJets));
  jetCSV_Stack = new THStack("jetCSV_Stack",Form("CSV(all jets) distributions (%dj%db);CSVv1;Counts",nJets, nbJets));
  pTselectbjet_Stack = new THStack("pTselectbjet_Stack",Form("p_{T}(b-jets) distributions (%dj%db);p_{T};Counts",nJets, nbJets));
  lepbjetMinv_min_Stack = new THStack("lepbjetMinv_min_Stack",Form("M_{lb} distributions (%dj%db);M_{lb};Counts",nJets, nbJets));
  ttMinv_min_Stack = new THStack("ttMinv_min_Stack",Form("M_{ttbar} distributions (%dj%db);M_{ttbar};Counts",nJets, nbJets));
  qqMinv_min_Stack = new THStack("qqMinv_min_Stack",Form("M_{qq} distributions (%dj%db);M_{qq};Counts",nJets, nbJets));
  h_Phi_1stb_Stack = new THStack("h_Phi_1stb_Stack",Form("#Delta#phi>0.5#pi muon vs leading-CSV-bjet (%dj%db);#Delta#phi;Counts",nJets, nbJets));
  h_Phi_allb_Stack = new THStack("h_Phi_allb_Stack",Form("#Delta#phi>0.5#pi muon vs all b-jets (%dj%db);#Delta#phi;Counts",nJets, nbJets));
  h_Phi_b_minpi_Stack = new THStack("h_Phi_b_minpi_Stack",Form("Min(#Delta#phi - #pi) of muon vs b-jet (%dj%db);#Delta#phi;Counts",nJets, nbJets));

  //for(int j = 3; j > 0; j--){
  for(int j = 1; j < 6; j = j+2){

  	if(j == 5) j = 2;

    lepPt_Stack->Add(h_lepPt_[j]);
    lepPseu_Stack->Add(h_lepPseu_[j]);
    jetHt_Stack->Add(h_jetHt_[j]);
    
    if(nbJets == 0) jetCSV_Stack->Add(h_jetCSV_[j]);
    else { 
      pTselectbjet_Stack->Add(h_pTselectbjet_[j]);
      lepbjetMinv_min_Stack->Add(h_lepbjetMinv_min_[j]);
      ttMinv_min_Stack->Add(h_ttMinv_min_[j]);
      qqMinv_min_Stack->Add(h_qqMinv_min_[j]);
      h_Phi_1stb_Stack->Add(h_Phi_1stb_[j]);
      h_Phi_allb_Stack->Add(h_Phi_allb_[j]);
      h_Phi_b_minpi_Stack->Add(h_Phi_b_minpi_[j]);
    }

  }

  leg->AddEntry(h_lepPt_[0],"Data","p");
  leg->AddEntry(h_lepPt_[1],"ttbar (MC)","f");
  leg->AddEntry(h_lepPt_[2],"Wjets (MC)","f");
  leg->AddEntry(h_lepPt_[3],"DY (MC)","f");
  leg->AddEntry(h_lepPt_[4],"Multijets (Data)","f");
  leg2->AddEntry(h_lepPt_[0],"Data","p");
  leg2->AddEntry(h_lepPt_[1],"ttbar (MC)","f");
  leg2->AddEntry(h_lepPt_[2],"Wjets (MC)","f");
  leg2->AddEntry(h_lepPt_[3],"DY (MC)","f");
  leg2->AddEntry(h_lepPt_[4],"Multijets (Data)","f");

  //cst->Divide(4,3);
  //cst->Divide(3,2);
  cst->Divide(3);
  cst->cd(1);
  gPad->SetLogy();
  lepPt_Stack->Draw(); 
  h_lepPt_[0]->Draw("same,ep");
  if(lepPt_Stack->GetMaximum() < h_lepPt_[0]->GetMaximum()){
    lepPt_Stack->SetMaximum((int)( h_lepPt_[0]->GetMaximum() + 0.15 * h_lepPt_[0]->GetMaximum() ));
  }
  lepPt_Stack->SetMinimum(0.1);
  leg->Draw();
 
  cst->cd(2);
/*  gPad->SetLogy();
  lepPseu_Stack->Draw();
  h_lepPseu_[0]->Draw("same,ep");
  if(lepPseu_Stack->GetMaximum() < h_lepPseu_[0]->GetMaximum()){
    lepPseu_Stack->SetMaximum((int)( h_lepPseu_[0]->GetMaximum() + 0.15 * h_lepPseu_[0]->GetMaximum() ));
  }
  lepPseu_Stack->SetMinimum(0.1);
  leg->Draw();
  
  cst->cd(3);
  //gPad->SetLogy();
  jetHt_Stack->Draw();
  h_jetHt_[0]->Draw("same,ep");
  if(jetHt_Stack->GetMaximum() < h_jetHt_[0]->GetMaximum()){
    jetHt_Stack->SetMaximum((int)( h_jetHt_[0]->GetMaximum() + 0.15 * h_jetHt_[0]->GetMaximum() ));
  }
  //jetHt_Stack->SetMinimum(0.1);
  leg->Draw();
  
  cst->cd(4);
  if(nbJets == 0){
    //gPad->SetLogy();
    jetCSV_Stack->Draw();
    h_jetCSV_[0]->Draw("same,ep");
    if(jetCSV_Stack->GetMaximum() < h_jetCSV_[0]->GetMaximum()){
      jetCSV_Stack->SetMaximum((int)( h_jetCSV_[0]->GetMaximum() + 0.15 * h_jetCSV_[0]->GetMaximum() ));
    }
    //jetCSV_Stack->SetMinimum(0.1);
    leg->Draw();
  } 
  else {    
*/    gPad->SetLogy();
    pTselectbjet_Stack->Draw();
    h_pTselectbjet_[0]->Draw("same,ep");
    if(pTselectbjet_Stack->GetMaximum() < h_pTselectbjet_[0]->GetMaximum()){
      pTselectbjet_Stack->SetMaximum((int)( h_pTselectbjet_[0]->GetMaximum() + 0.15 * h_pTselectbjet_[0]->GetMaximum() ));
    }
    pTselectbjet_Stack->SetMinimum(0.1);
    leg->Draw();
/*
    //cst->cd(9);
    cst->cd(1);
    gPad->SetLogy();
	  lepbjetMinv_min_Stack->Draw();
    h_lepbjetMinv_min_[0]->Draw("same,ep");
    if(lepbjetMinv_min_Stack->GetMaximum() < h_lepbjetMinv_min_[0]->GetMaximum()){
      lepbjetMinv_min_Stack->SetMaximum((int)( h_lepbjetMinv_min_[0]->GetMaximum() + 0.15 * h_lepbjetMinv_min_[0]->GetMaximum() ));
    }
    lepbjetMinv_min_Stack->SetMinimum(0.1);
    leg->Draw();
*/    
    //cst->cd(10);
    cst->cd(3);
    gPad->SetLogy();
    qqMinv_min_Stack->Draw();
    h_qqMinv_min_[0]->Draw("same,ep");
    if(qqMinv_min_Stack->GetMaximum() < h_qqMinv_min_[0]->GetMaximum()){
      qqMinv_min_Stack->SetMaximum((int)( h_qqMinv_min_[0]->GetMaximum() + 0.15 * h_qqMinv_min_[0]->GetMaximum() ));
    }
    qqMinv_min_Stack->SetMinimum(0.1);
    leg->Draw();
/*
    //cst->cd(11);
    cst->cd(3);
    gPad->SetLogy();
    ttMinv_min_Stack->Draw();
    h_ttMinv_min_[0]->Draw("same,ep");
    if(ttMinv_min_Stack->GetMaximum() < h_ttMinv_min_[0]->GetMaximum()){
      ttMinv_min_Stack->SetMaximum((int)( h_ttMinv_min_[0]->GetMaximum() + 0.15 * h_ttMinv_min_[0]->GetMaximum() ));
    }
    ttMinv_min_Stack->SetMinimum(0.1);
    leg->Draw();

    cst->cd(5);
    cst->cd(2);
    gPad->SetLogy();
    h_Phi_1stb_Stack->Draw();
    h_Phi_1stb_[0]->Draw("same,ep");
    if(h_Phi_1stb_Stack->GetMaximum() < h_Phi_1stb_[0]->GetMaximum()){
      h_Phi_1stb_Stack->SetMaximum((int)( h_Phi_1stb_[0]->GetMaximum() + 0.15 * h_Phi_1stb_[0]->GetMaximum() ));
    }
    h_Phi_1stb_Stack->SetMinimum(0.1);
    leg2->Draw();

    cst->cd(6);
    //cst->cd(5);
    //gPad->SetLogy();
    h_Phi_allb_Stack->Draw();
    h_Phi_allb_[0]->Draw("same,ep");
    if(h_Phi_allb_Stack->GetMaximum() < h_Phi_allb_[0]->GetMaximum()){
      h_Phi_allb_Stack->SetMaximum((int)( h_Phi_allb_[0]->GetMaximum() + 0.15 * h_Phi_allb_[0]->GetMaximum() ));
    }
    //h_Phi_allb_Stack->SetMinimum(0.1);
    leg2->Draw();

    cst->cd(7);
    cst->cd(3);
    gPad->SetLogy();
    h_Phi_b_minpi_Stack->Draw();
    h_Phi_b_minpi_[0]->Draw("same,ep");
    if(h_Phi_b_minpi_Stack->GetMaximum() < h_Phi_b_minpi_[0]->GetMaximum()){
      h_Phi_b_minpi_Stack->SetMaximum((int)( h_Phi_b_minpi_[0]->GetMaximum() + 0.15 * h_Phi_b_minpi_[0]->GetMaximum() ));
    }
    h_Phi_b_minpi_Stack->SetMinimum(0.1);
    leg2->Draw();
  }
*/
  cst->Update();

  return cst;
}