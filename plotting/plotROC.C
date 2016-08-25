#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "anaMuJetsSkimTree.h"

#include <string>
#include <vector>

int nJets = 4;
int nbJets = 1;

const int nCSVcuts = 21;
const int nFiles = 4;


void plotROC(std::string outFileName = "", bool CSVv2 = false, double pTLow = 0., double pTHigh = 999999.){

  anaMuJetsSkimTree ana;

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/030816/MC_tt/MCttV2.root", 0);            //NoLepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/030816/MC_W/MCWV2.root", 1);              //NoLepIso cut done
  //ana.OpenFiles("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/030816/MC_DY/MCDYV2.root", 2);            //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCttV2.root", 0);            //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCWV2.root", 1);            //NoLepIso cut done
  ana.OpenFiles("~/Documents/MCDYV2.root", 2);          //NoLepIso cut done
  ana.OpenFiles("~/Documents/DataMultijetsV2.root", 3);          //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  TH1F* h_highestCSVv1[3];
  for(int i = 0; i < 3; i++){
    if(CSVv2 == true) h_highestCSVv1[i] = new TH1F(Form("h_highestCSVv1_%d",i),"Highest CSVv2 value of jets that passed cuts;CSVv2;Counts",20,0,1);
    else h_highestCSVv1[i] = new TH1F(Form("h_highestCSVv1_%d",i),"Highest CSVv1 value of jets that passed cuts;CSVv1;Counts",20,0,1);
  }

  cout << "Making ROC plots in the following pT interval." << endl;
  cout << "     ptLow = " << pTLow << " , pTHigh = " << pTHigh << endl;
  cout << "     Using CSVv2 = " << std::boolalpha << CSVv2 << endl;
  	

  int nPassedCSVCut010[nFiles][nCSVcuts], nPassedCSVCut1030[nFiles][nCSVcuts], nPassedCSVCut30100[nFiles][nCSVcuts];
  int nDidntPassCSVCut010[nFiles][nCSVcuts], nDidntPassCSVCut1030[nFiles][nCSVcuts], nDidntPassCSVCut30100[nFiles][nCSVcuts];
  double CSVcuts[nCSVcuts];
  for(int i = 0; i < nFiles; i++){
    for(int j = 0; j < nCSVcuts; j++){
      nPassedCSVCut010[i][j] = nPassedCSVCut1030[i][j] = nPassedCSVCut30100[i][j] = 0;
      nDidntPassCSVCut010[i][j] = nDidntPassCSVCut1030[i][j] = nDidntPassCSVCut30100[i][j] = 0;
    }
  }


  for(int i = 0; i < nFiles; i++){

    ana.SetAddressBranches(i);

    for(int entry = 0; entry < (int)ana.tr_[i]->GetEntries(); entry++){
      
      ana.tr_[i]->GetEntry(entry);

      //Choose leading muon when there are 2 or more
      int indexMuon;
      if(i == 3)  indexMuon = ana.FindLeadingMuonAfterCuts(4); //4=Reversed LepIsolation option
      else        indexMuon = ana.FindLeadingMuonAfterCuts(1); //1=LepIsolation option
      if(indexMuon == -999) continue;  

      //Store the jets that pass the basic cuts
      ana.JetSmearing = false;
      std::vector<int> indexJets = ana.FindJetsAfterCuts(indexMuon); //JetSmearing turned off
      std::vector<int> srt_indexJets = ana.BuildpTVectorLeadingbJets(indexJets);
      if((int)indexJets.size() < nJets) continue;

      ana.SetCSVCut(1.05);
      for(int l = 0; l < nCSVcuts; l++){
        ana.SetCSVCut( ana.GetCSVCut() - 0.05 );
        CSVcuts[l] = ana.GetCSVCut();
  
        //std::vector<int> indexbJets = ana.FindbJets(indexJets);
        std::vector<int> indexbJets;
        if(CSVv2 == true) indexbJets = ana.FindbJetsCSVv2(indexJets);
        else              indexbJets = ana.FindbJets(indexJets);

        std::vector<int> srt_indexbJets = ana.BuildpTVectorLeadingbJets(indexbJets);
        std::vector<int> srt_CSV_indexbJets;
        if(CSVv2 == true) srt_CSV_indexbJets = ana.BuildCSVVectorLeadingbJetsCSVv2(indexbJets);
        else              srt_CSV_indexbJets = ana.BuildCSVVectorLeadingbJets(indexbJets);

        if(ana.GethiBin() < 20){

          if((int)indexbJets.size() < nbJets){
            if(ana.GetjetPt(indexJets[srt_indexJets[0]]) >= pTLow && ana.GetjetPt(indexJets[srt_indexJets[0]]) < pTHigh){
              nDidntPassCSVCut010[i][l]++;
            }
          } else {
            if(ana.GetjetPt(indexbJets[srt_indexbJets[0]]) >= pTLow && ana.GetjetPt(indexbJets[srt_indexbJets[0]]) < pTHigh){
              nPassedCSVCut010[i][l]++;
            }

            if(i==0 && l==20){ 
              if(CSVv2 == true) h_highestCSVv1[0]->Fill( ana.GetjetCSVv2( indexbJets[srt_CSV_indexbJets[0]] ) );
              else              h_highestCSVv1[0]->Fill( ana.GetjetCSVv1( indexbJets[srt_CSV_indexbJets[0]] ) );
            }
          }

        } else if(ana.GethiBin() >= 20 && ana.GethiBin() < 60){

          if((int)indexbJets.size() < nbJets){
            if(ana.GetjetPt(indexJets[srt_indexJets[0]]) >= pTLow && ana.GetjetPt(indexJets[srt_indexJets[0]]) < pTHigh){
              nDidntPassCSVCut1030[i][l]++;
            }
          } else {
            if(ana.GetjetPt(indexbJets[srt_indexbJets[0]]) >= pTLow && ana.GetjetPt(indexbJets[srt_indexbJets[0]]) < pTHigh){
              nPassedCSVCut1030[i][l]++;
            }

            if(i==0 && l==20){ 
              if(CSVv2 == true) h_highestCSVv1[1]->Fill( ana.GetjetCSVv2( indexbJets[srt_CSV_indexbJets[0]] ) );
              else              h_highestCSVv1[1]->Fill( ana.GetjetCSVv1( indexbJets[srt_CSV_indexbJets[0]] ) );
            }
          }

        } else {

          if((int)indexbJets.size() < nbJets){
            if(ana.GetjetPt(indexJets[srt_indexJets[0]]) >= pTLow && ana.GetjetPt(indexJets[srt_indexJets[0]]) < pTHigh){
              nDidntPassCSVCut30100[i][l]++;
            }
          } else {
            if(ana.GetjetPt(indexbJets[srt_indexbJets[0]]) >= pTLow && ana.GetjetPt(indexbJets[srt_indexbJets[0]]) < pTHigh){
              nPassedCSVCut30100[i][l]++;
            }

            if(i==0 && l==20){ 
              if(CSVv2 == true) h_highestCSVv1[2]->Fill( ana.GetjetCSVv2( indexbJets[srt_CSV_indexbJets[0]] ) );
              else              h_highestCSVv1[2]->Fill( ana.GetjetCSVv1( indexbJets[srt_CSV_indexbJets[0]] ) );
            }
          }

        }
      }
    }
  }

  double Normaliz[4] = {1, 0.71091, 0.084058, 2696./13332.};
  for(int i = 0; i < nFiles; i++){
    cout << "Going to normalize File " << i << " using " << Normaliz[i] << endl;
    for(int j = 0; j < nCSVcuts; j++){
      nPassedCSVCut010[i][j] = Normaliz[i] * nPassedCSVCut010[i][j];
      nDidntPassCSVCut010[i][j] = Normaliz[i] * nDidntPassCSVCut010[i][j];
      nPassedCSVCut1030[i][j] = Normaliz[i] * nPassedCSVCut1030[i][j];
      nDidntPassCSVCut1030[i][j] = Normaliz[i] * nDidntPassCSVCut1030[i][j];
      nPassedCSVCut30100[i][j] = Normaliz[i] * nPassedCSVCut30100[i][j];
      nDidntPassCSVCut30100[i][j] = Normaliz[i] * nDidntPassCSVCut30100[i][j];
    }
  }

  Color_t colors[3] = {8,4,2};
  for(int i = 0; i < 3; i++){
    h_highestCSVv1[i]->Sumw2();
    h_highestCSVv1[i]->Scale(1./h_highestCSVv1[i]->Integral() );
    //h_highestCSVv1[i]->SetMarkerStyle(20);
    //h_highestCSVv1[i]->SetMarkerColor(colors[i]);
    h_highestCSVv1[i]->SetLineColor(colors[i]);
    h_highestCSVv1[i]->SetLineWidth(2);
  }

  double x010[nCSVcuts], y010[nCSVcuts], x1030[nCSVcuts], y1030[nCSVcuts], x30100[nCSVcuts], y30100[nCSVcuts];
  double ex010[nCSVcuts], ey010[nCSVcuts], ex1030[nCSVcuts], ey1030[nCSVcuts], ex30100[nCSVcuts], ey30100[nCSVcuts];
  for(int i = 0; i < nFiles; i++){
    cout << "File: " << i << endl;
    for(int j = 0; j < nCSVcuts; j++){

      cout << "  CSV cut : " << CSVcuts[j] << endl;
      cout << "     HiBin 0-10 : " << nPassedCSVCut010[i][j] << " +- " << sqrt(nPassedCSVCut010[i][j]) << " , " << nDidntPassCSVCut010[i][j] << " +- " << sqrt(nDidntPassCSVCut010[i][j])<< endl;
      cout << "     HiBin 10-30 : " << nPassedCSVCut1030[i][j] << " +- " << sqrt(nPassedCSVCut1030[i][j]) << " , " << nDidntPassCSVCut1030[i][j]<< " +- " << sqrt(nDidntPassCSVCut1030[i][j]) << endl;
      cout << "     HiBin 30-100 : " << nPassedCSVCut30100[i][j] << " +- " << sqrt(nPassedCSVCut30100[i][j]) << " , " << nDidntPassCSVCut30100[i][j] << " +- " << sqrt(nDidntPassCSVCut30100[i][j])<< endl;

      if(i==0){
        x010[j] = 100.* nPassedCSVCut010[i][j]/(nPassedCSVCut010[i][j] + nDidntPassCSVCut010[i][j]);
        x1030[j] = 100.* nPassedCSVCut1030[i][j]/(nPassedCSVCut1030[i][j] + nDidntPassCSVCut1030[i][j]);
        x30100[j] = 100.* nPassedCSVCut30100[i][j]/(nPassedCSVCut30100[i][j] + nDidntPassCSVCut30100[i][j]);

        if(nPassedCSVCut010[i][j] > 0) ex010[j] = x010[j] * sqrt( 1./(double)nPassedCSVCut010[i][j] + 1./( (double)nPassedCSVCut010[i][j] + (double)nDidntPassCSVCut010[i][j] ) );
        else ex010[j] = 0;

        if(nPassedCSVCut1030[i][j] > 0) ex1030[j] = x1030[j] * sqrt( 1./(double)nPassedCSVCut1030[i][j] + 1./( (double)nPassedCSVCut1030[i][j] + (double)nDidntPassCSVCut30100[i][j] ) );
        else ex1030[j] = -0;

        if(nPassedCSVCut30100[i][j] > 0) ex30100[j] = x30100[j] * sqrt( 1./(double)nPassedCSVCut30100[i][j] + 1./( (double)nPassedCSVCut1030[i][j] + (double)nDidntPassCSVCut30100[i][j] ) );
        else ex30100[j] = 0;
      }
      else if(i == 1){ 
        //nPassedCSVCut010[1][j] = nPassedCSVCut1030[1][j] = nPassedCSVCut30100[1][j] = 0;

        //y010[j] = 100.* (nPassedCSVCut010[i][j] + nPassedCSVCut010[2][j])/(nPassedCSVCut010[i][j] + nDidntPassCSVCut010[i][j] + nPassedCSVCut010[2][j] + nDidntPassCSVCut010[2][j]);
        //y1030[j] = 100.* (nPassedCSVCut1030[i][j] + nPassedCSVCut1030[2][j])/(nPassedCSVCut1030[i][j] + nDidntPassCSVCut1030[i][j] + nPassedCSVCut1030[2][j] + nDidntPassCSVCut1030[2][j]);
        //y30100[j] = 100.* (nPassedCSVCut30100[i][j] + nPassedCSVCut30100[2][j])/(nPassedCSVCut30100[i][j] + nDidntPassCSVCut30100[i][j] + nPassedCSVCut30100[2][j] + nDidntPassCSVCut30100[2][j]);
        
        y010[j] = 100.* (nPassedCSVCut010[i][j] + nPassedCSVCut010[2][j] + nPassedCSVCut010[3][j])/(nPassedCSVCut010[i][j] + nDidntPassCSVCut010[i][j] + nPassedCSVCut010[2][j] + nDidntPassCSVCut010[2][j] + nPassedCSVCut010[3][j] + nDidntPassCSVCut010[3][j]);
        y1030[j] = 100.* (nPassedCSVCut1030[i][j] + nPassedCSVCut1030[2][j] + nPassedCSVCut1030[3][j])/(nPassedCSVCut1030[i][j] + nDidntPassCSVCut1030[i][j] + nPassedCSVCut1030[2][j] + nDidntPassCSVCut1030[2][j] + nPassedCSVCut1030[3][j] + nDidntPassCSVCut1030[3][j]);
        y30100[j] = 100.* (nPassedCSVCut30100[i][j] + nPassedCSVCut30100[2][j] + nPassedCSVCut30100[3][j])/(nPassedCSVCut30100[i][j] + nDidntPassCSVCut30100[i][j] + nPassedCSVCut30100[2][j] + nDidntPassCSVCut30100[2][j] + nPassedCSVCut30100[3][j] + nDidntPassCSVCut30100[3][j]);

        //if(nPassedCSVCut010[i][j] > 0 || nPassedCSVCut010[2][j] > 0) ey010[j] = y010[j] * sqrt( 1./((double)nPassedCSVCut010[i][j] + (double)nPassedCSVCut010[2][j]) + 1./( (double)nPassedCSVCut010[i][j] + (double)nDidntPassCSVCut010[i][j] + (double)nPassedCSVCut010[2][j] + (double)nDidntPassCSVCut010[2][j] ) );
        //else ey010[j] = 0;
        //y010[j] = 100. - y010[j];

        if(nPassedCSVCut010[i][j] > 0 || nPassedCSVCut010[2][j] > 0 || nPassedCSVCut010[3][j] > 0) ey010[j] = y010[j] * sqrt( 1./((double)nPassedCSVCut010[i][j] + (double)nPassedCSVCut010[2][j] + (double)nPassedCSVCut010[3][j]) + 1./( (double)nPassedCSVCut010[i][j] + (double)nDidntPassCSVCut010[i][j] + (double)nPassedCSVCut010[2][j] + (double)nDidntPassCSVCut010[2][j] + (double)nPassedCSVCut010[3][j] + (double)nDidntPassCSVCut010[3][j]) );
        else ey010[j] = 0;
        y010[j] = 100. - y010[j];
        
        //if(nPassedCSVCut1030[i][j] > 0 || nPassedCSVCut1030[2][j] > 0) ey1030[j] = y1030[j] * sqrt( 1./((double)nPassedCSVCut1030[i][j] + (double)nPassedCSVCut1030[2][j]) + 1./( (double)nPassedCSVCut1030[i][j] + (double)nDidntPassCSVCut1030[i][j] + (double)nPassedCSVCut1030[2][j] + (double)nDidntPassCSVCut1030[2][j] ) );
        //else ey1030[j] = 0;
        //y1030[j] = 100. - y1030[j];

        if(nPassedCSVCut1030[i][j] > 0 || nPassedCSVCut1030[2][j] > 0 || nPassedCSVCut1030[3][j] > 0) ey1030[j] = y1030[j] * sqrt( 1./((double)nPassedCSVCut1030[i][j] + (double)nPassedCSVCut1030[2][j] + (double)nPassedCSVCut1030[3][j]) + 1./( (double)nPassedCSVCut1030[i][j] + (double)nDidntPassCSVCut1030[i][j] + (double)nPassedCSVCut1030[2][j] + (double)nDidntPassCSVCut1030[2][j] + (double)nPassedCSVCut1030[3][j] + (double)nDidntPassCSVCut1030[3][j]) );
        else ey1030[j] = 0;
        y1030[j] = 100. - y1030[j];
        
        //if(nPassedCSVCut30100[i][j] > 0 || nPassedCSVCut30100[2][j] > 0) ey30100[j] = y30100[j] * sqrt( 1./((double)nPassedCSVCut30100[i][j] + (double)nPassedCSVCut30100[2][j]) + 1./( (double)nPassedCSVCut30100[i][j] + (double)nDidntPassCSVCut30100[i][j] + (double)nPassedCSVCut30100[2][j] + (double)nDidntPassCSVCut30100[2][j] ) );
        //else ey30100[j] = 0;
        //y30100[j] = 100. - y30100[j];

        if(nPassedCSVCut30100[i][j] > 0 || nPassedCSVCut30100[2][j] > 0 || nPassedCSVCut30100[3][j] > 0) ey30100[j] = y30100[j] * sqrt( 1./((double)nPassedCSVCut30100[i][j] + (double)nPassedCSVCut30100[2][j] + (double)nPassedCSVCut30100[3][j]) + 1./( (double)nPassedCSVCut30100[i][j] + (double)nDidntPassCSVCut30100[i][j] + (double)nPassedCSVCut30100[2][j] + (double)nDidntPassCSVCut30100[2][j] + (double)nPassedCSVCut30100[3][j] + (double)nDidntPassCSVCut30100[3][j]) );
        else ey30100[j] = 0;
        y30100[j] = 100. - y30100[j];
      }
    }
  }

  //TGraph* CSV_effcut010 = new TGraph(nCSVcuts,x010,y010);
  TGraphErrors* CSV_effcut010 = new TGraphErrors(nCSVcuts,x010,y010,ex010,ey010);
  TGraph* CSV_workp010 = new TGraph(1, &x010[5], &y010[5]);
  //TGraph* CSV_effcut1030 = new TGraph(nCSVcuts,x1030,y1030);
  TGraphErrors* CSV_effcut1030 = new TGraphErrors(nCSVcuts,x1030,y1030,ex1030,ey1030);
  TGraph* CSV_workp1030 = new TGraph(1, &x1030[5], &y1030[5]);
  //TGraph* CSV_effcut30100 = new TGraph(nCSVcuts,x30100,y30100);
  TGraphErrors* CSV_effcut30100 = new TGraphErrors(nCSVcuts,x30100,y30100,ex30100,ey30100);
  TGraph* CSV_workp30100 = new TGraph(1, &x30100[5], &y30100[5]);
  TGraph* CSV_workp30100leg = new TGraph(1, &x30100[5], &y30100[5]);

  CSV_effcut010->SetNameTitle("g_CSV_effcut010",Form("Efficiency CSV cut in pT(%d, %d), Centrality 0-10;Sig eff; Back rejection eff (W + DY + MJ)",(int)pTLow, (int)pTHigh));
  CSV_effcut1030->SetNameTitle("g_CSV_effcut1030",Form("Efficiency CSV cut in pT(%d, %d), Centrality 10-30;Sig eff; Back rejection eff (W + DY + MJ)",(int)pTLow, (int)pTHigh));
  CSV_effcut30100->SetNameTitle("g_CSV_effcut30100",Form("Efficiency CSV cut in pT(%d, %d), Centrality 30-100;Sig eff; Back rejection eff (W + DY + MJ)",(int)pTLow, (int)pTHigh));

  CSV_effcut010->SetMarkerStyle(20);
  CSV_workp010->SetMarkerStyle(29);
  CSV_effcut1030->SetMarkerStyle(20);
  CSV_workp1030->SetMarkerStyle(29);
  CSV_effcut30100->SetMarkerStyle(20);
  CSV_workp30100->SetMarkerStyle(29);
  CSV_workp30100leg->SetMarkerStyle(29);

  CSV_effcut010->SetMarkerColor(8);
  CSV_effcut010->SetLineColor(8);
  CSV_workp010->SetMarkerColor(8);
  CSV_workp010->SetMarkerSize(2.5);
  CSV_effcut1030->SetMarkerColor(kBlue);
  CSV_effcut1030->SetLineColor(kBlue);
  CSV_workp1030->SetMarkerColor(kBlue);
  CSV_workp1030->SetMarkerSize(2.5);
  CSV_effcut30100->SetMarkerColor(kRed);
  CSV_effcut30100->SetLineColor(kRed);
  CSV_workp30100->SetMarkerColor(kRed);
  CSV_workp30100->SetMarkerSize(2.5);
  CSV_workp30100leg->SetMarkerSize(2);

  TMultiGraph *mg = new TMultiGraph();
  if(CSVv2 == true) mg->SetNameTitle("mg",Form("ROC: CSVv2 cut in pT range(%d, %d);Sig eff; Back rejection eff (W + DY + MJ)",(int)pTLow, (int)pTHigh));
  else              mg->SetNameTitle("mg",Form("ROC: CSVv1 cut in pT range(%d, %d);Sig eff; Back rejection eff (W + DY + MJ)",(int)pTLow, (int)pTHigh));
  mg->Add(CSV_effcut010, "lp");
  mg->Add(CSV_effcut1030, "lp");
  mg->Add(CSV_effcut30100, "lp");
  mg->Add(CSV_workp010,"p");
  mg->Add(CSV_workp1030,"p");
  mg->Add(CSV_workp30100,"p");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
  leg->AddEntry(CSV_effcut010,"PbPb 0-10%","p");
  leg->AddEntry(CSV_effcut1030,"PbPb 10-30%","p");
  leg->AddEntry(CSV_effcut30100,"PbPb 30-100%","p");
  leg->AddEntry(CSV_workp30100leg,"Working points: CSV > 0.75","p");

  TLegend* leg2 = new TLegend(0.12,0.68,0.32,0.88);
  leg2->AddEntry(h_highestCSVv1[0],"PbPb 0-10%","l");
  leg2->AddEntry(h_highestCSVv1[1],"PbPb 10-30%","l");
  leg2->AddEntry(h_highestCSVv1[2],"PbPb 30-100%","l");

  TCanvas *cst = new TCanvas("cst","Histograms",10,10,1800,1000);
  cst->cd();
  mg->Draw("A");
  mg->GetXaxis()->SetLimits(0,100);
  mg->GetYaxis()->SetRangeUser(0,100);
  leg->Draw();

  TCanvas* cst2 = new TCanvas("cst2","CSV different centrality",10,10,1800,1000);
  cst2->cd();
  gStyle->SetOptStat(0);
  h_highestCSVv1[2]->Draw("HIST");
  h_highestCSVv1[1]->Draw("same, HIST");
  h_highestCSVv1[0]->Draw("same, HIST");
  leg2->Draw();

  cout << endl << endl << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  cst->Write();
  cst2->Write();
  outFile->Write();
  delete outFile;
}