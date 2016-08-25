#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"

#include <string>
#include <vector>

const int nDiscrCuts = 20; int nCut = 15;

void plotROC_ak4PFJetAnalyzer(TString FileName = "")
{

  TFile* f_ = TFile::Open(FileName);
  TH1F* h_CSV_[2];
  h_CSV_[0] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv1_0"));
  h_CSV_[1] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv1_1"));

  TH1F* h_CSVv2_[2];
  h_CSVv2_[0] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv2_0"));
  h_CSVv2_[1] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv2_1"));

  TH1F* h_CSV_Cs4_[2];
  h_CSV_Cs4_[0] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv1_CS40"));
  h_CSV_Cs4_[1] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv1_CS41"));

  TH1F* h_CSVv2_Cs4_[2];
  h_CSVv2_Cs4_[0] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv2_CS40"));
  h_CSVv2_Cs4_[1] = dynamic_cast<TH1F*>(f_->Get("h_discr_CSVv2_CS41"));

  Color_t colors2[2] = {2,4};
  for(int i = 0; i < 2; i++){
    //h_CSV_b[i]->Sumw2();
    //h_CSV_b[i]->Scale(1./h_CSV_b[i]->Integral() );
    h_CSV_[i]->SetLineWidth(2);
    h_CSVv2_[i]->SetLineWidth(2);
    h_CSV_Cs4_[i]->SetLineWidth(2);
    h_CSVv2_Cs4_[i]->SetLineWidth(2);
    h_CSV_[i]->SetLineColor(colors2[i]);
    h_CSVv2_[i]->SetLineColor(colors2[i]);
    h_CSV_Cs4_[i]->SetLineColor(colors2[i]);
    h_CSVv2_Cs4_[i]->SetLineColor(colors2[i]);
  }
/*
  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(h_CSV_[0],"MC b-jet","l");
  leg->AddEntry(h_CSV_[1],"MC light-jet","l");

  TCanvas* cst = new TCanvas("cst","b vs light",10,10,1400,500);
  cst->Divide(2);
  cst->cd(1);
  gStyle->SetOptStat(0); 
  h_CSV_[1]->Draw("HIST");
  h_CSV_[0]->Draw("same, HIST");
  leg->Draw();
  cst->cd(2);
  h_CSV_Cs4_[1]->Draw("HIST");
  h_CSV_Cs4_[0]->Draw("same, HIST");
  leg->Draw();
*/

  h_CSV_[0]->SetLineColor(2);
  h_CSV_[1]->SetLineColor(2);
  h_CSVv2_[0]->SetLineColor(2);
  h_CSVv2_[1]->SetLineColor(2);
  h_CSV_Cs4_[0]->SetLineColor(4);
  h_CSV_Cs4_[1]->SetLineColor(4);
  h_CSVv2_Cs4_[0]->SetLineColor(4);
  h_CSVv2_Cs4_[1]->SetLineColor(4);
  h_CSV_[0]->SetTitle("CSVv1 for MC b-jets (PbPb 0-100);discr_CSVv1;Counts");
  h_CSV_[1]->SetTitle("CSVv1 for MC light-jets (PbPb 0-100);discr_CSVv1;Counts");
  h_CSV_Cs4_[0]->SetTitle("CSVv1 for MC b-jets (PbPb 0-100);discr_CSVv1;Counts");
  h_CSV_Cs4_[1]->SetTitle("CSVv1 for MC light-jets (PbPb 0-100);discr_CSVv1;Counts");

  TLegend* leg2 = new TLegend(0.5,0.7,0.7,0.9);
  leg2->AddEntry(h_CSV_[0],"Cs2","l");
  leg2->AddEntry(h_CSV_Cs4_[0],"Pu4","l");

  TCanvas* cst2 = new TCanvas("cst2","Cs2 vs Pu4",10,10,1400,500);
  cst2->Divide(2);
  cst2->cd(1);
  gStyle->SetOptStat(0); 
  h_CSV_Cs4_[0]->Draw("HIST");
  h_CSV_[0]->Draw("same, HIST");
  leg2->Draw();
  cst2->cd(2);
  h_CSV_Cs4_[1]->Draw("HIST");
  h_CSV_[1]->Draw("same, HIST");
  leg2->Draw();


  double btageff[nDiscrCuts], lightrejeff[nDiscrCuts], CSVcuts[nDiscrCuts];
  double btageffCs4[nDiscrCuts], lightrejeffCs4[nDiscrCuts];
  double csvCut = 0.; //(till 1)
  for(int l = 0; l < nDiscrCuts; l++){
        
    double CSVcut = csvCut + l * 0.05;        
    CSVcuts[l] = CSVcut;
	
	  int binx1 = h_CSV_[0]->GetXaxis()->FindBin( CSVcut );
    int binx2 = h_CSV_[0]->GetXaxis()->GetLast();

	  btageff[l] = h_CSV_[0]->Integral(binx1, binx2)/h_CSV_[0]->GetEntries();
	  lightrejeff[l] = 1. - h_CSV_[1]->Integral(binx1, binx2)/h_CSV_[1]->GetEntries();

	  btageffCs4[l] = h_CSV_Cs4_[0]->Integral(binx1, binx2)/h_CSV_Cs4_[0]->GetEntries();
	  lightrejeffCs4[l] = 1. - h_CSV_Cs4_[1]->Integral(binx1, binx2)/h_CSV_Cs4_[1]->GetEntries();

  }

  TGraph* g_hiBinint = new TGraph(nDiscrCuts,btageff,lightrejeff);
  TGraph* g_workpint = new TGraph(1, &btageff[nCut], &lightrejeff[nCut]);
  g_hiBinint->SetNameTitle("g_hiBinint","Performance of b-tagging: CSVv1 (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  g_hiBinint->SetMarkerStyle(20);
  g_workpint->SetMarkerStyle(29);
  g_hiBinint->SetMarkerColor(2);
  g_hiBinint->SetLineColor(2);
  g_workpint->SetMarkerColor(2);
  g_workpint->SetMarkerSize(3.5);

  TGraph* g_hiBinintCs4 = new TGraph(nDiscrCuts,btageffCs4,lightrejeffCs4);
  TGraph* g_workpintCs4 = new TGraph(1, &btageffCs4[nCut], &lightrejeffCs4[nCut]);
  g_hiBinintCs4->SetNameTitle("g_hiBinint","Performance of b-tagging: CSVv1 (PbPb 0-100%);b-jet tagging efficiency; light-jet rejection efficiency");
  g_hiBinintCs4->SetMarkerStyle(20);
  g_workpintCs4->SetMarkerStyle(29);
  g_hiBinintCs4->SetMarkerColor(4);
  g_hiBinintCs4->SetLineColor(4);
  g_workpintCs4->SetMarkerColor(4);
  g_workpintCs4->SetMarkerSize(3.5);

  TMultiGraph *mg = new TMultiGraph();
  mg->SetNameTitle("mg","Performance of b-tagging: CSVv1;b-jet tagging efficiency; light-jet rejection efficiency");
  mg->Add(g_hiBinint, "lp");
  mg->Add(g_hiBinintCs4, "lp");
  mg->Add(g_workpint,"p");
  mg->Add(g_workpintCs4,"p");
  
  TLegend* leg99 = new TLegend(0.7,0.7,0.9,0.9);
  leg99->AddEntry(g_hiBinint,"Cs2","lp");
  leg99->AddEntry(g_hiBinintCs4,"Pu4","lp");
  leg99->AddEntry(g_workpint,Form("CSV = %.2f",CSVcuts[nCut]),"p");

  TCanvas *cst99 = new TCanvas("cst99","Histograms",10,10,1800,1000);
  cst99->cd();
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(0,1);
  mg->GetYaxis()->SetRangeUser(0,1);
  leg99->Draw();

}