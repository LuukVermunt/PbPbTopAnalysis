#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "anaMuJetsSkimTree.h"
#include "histControlDistributions.h"

#include <string>
#include <vector>

void OpenFiles(TString FileName);

  TTree* tr_;
  TFile* f_;

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

  //Max generated particles = 5000;
  int nGen;
  int genPdg[5000];
  float genPt[5000];
  float genPhi[5000];
  float genEta[5000];
  float genChg[5000];




void findMCbjet_test(std::string outFileName = ""){

  if(!strcmp(outFileName.c_str(), "")){
    std::cout << "No output specified. return" << std::endl;
    return;
  }
  //OpenFiles("~/Documents/Universiteit/SummerSchoolProject/topskim/test_MuJetsSkim_pdgNoJtCuts.root");          //NoLepIso cut done
  OpenFiles("~/Documents/MCtt_pdg_noJtCuts.root");          //NoLepIso cut done
  TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

  TH1F* h_pT_matchedMC[4];
  TH1F* h_pT_didntmatchedMC[4];
  for(int i = 0; i < 4; i++){
    h_pT_matchedMC[i] = new TH1F(Form("h_pT_matchedMC_%d",i),Form("pT MC %d b-quark (Matched);pT;Counts",i),150,0,300);
    h_pT_didntmatchedMC[i] = new TH1F(Form("h_pT_didntmatchedMC_%d",i),Form("pT MC %d b-quark (Not Matched);pT;Counts",i),150,0,300);
  }


  cout << endl << "First try to make own MC b-tag jet script" << endl;

  //BuildHistograms();

  tr_->SetBranchAddress("run", &run_);
  tr_->SetBranchAddress("evt", &evt_);
  tr_->SetBranchAddress("lumi", &lumi_);
  tr_->SetBranchAddress("hiBin", &hiBin);
  tr_->SetBranchAddress("vz", &vz);

  tr_->SetBranchAddress("nLep", &nLep);
  tr_->SetBranchAddress("lepID", lepID);
  tr_->SetBranchAddress("lepPt", lepPt);
  tr_->SetBranchAddress("lepPhi", lepPhi);
  tr_->SetBranchAddress("lepEta", lepEta);
  tr_->SetBranchAddress("lepChg", lepChg);
  tr_->SetBranchAddress("lepIso", lepIso);

  tr_->SetBranchAddress("nJt", &nJt);
  tr_->SetBranchAddress("jtPt", jtPt);
  tr_->SetBranchAddress("jtPhi", jtPhi);
  tr_->SetBranchAddress("jtEta", jtEta);
  tr_->SetBranchAddress("jtM", jtM);
  tr_->SetBranchAddress("discr_csvV1", discr_csvV1); 
  tr_->SetBranchAddress("discr_csvV2", discr_csvV2); 

  tr_->SetBranchAddress("nGen", &nGen);
  tr_->SetBranchAddress("genPdg", genPdg);
  tr_->SetBranchAddress("genPt", genPt);
  tr_->SetBranchAddress("genPhi", genPhi);
  tr_->SetBranchAddress("genEta", genEta);
  tr_->SetBranchAddress("genChg", genChg);

  int ibj = 0;
  int missingbj1(0), missingbj2(0), missingbbarj1(0), missingbbarj2(0);
  int doubleb(0), doublebbar(0);

  for(int entry = 0; entry < (int)tr_->GetEntries(); entry++){

    tr_->GetEntry(entry);


    //Choose leading muon
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
    
    if((int)indexJets.size() < 4) continue;





    bool bquark = false;
    bool tquark = false;
    for(int igen = 0; igen < nGen; igen++){
      if( genPdg[igen] == 5) bquark = true;
      if( genPdg[igen] == 6) tquark = true;
    }
    if(bquark == false || tquark == false){
      cout << "Event " << entry << ". No b-quark = " << std::boolalpha << bquark << " or t-quark " << tquark << endl;
      cout << "Drop event........................" << endl;
      continue;
    }

    bool bbarquark = false;
    bool tbarquark = false;
    for(int igen = 0; igen < nGen; igen++){
      if( genPdg[igen] == -5) bbarquark = true;
      if( genPdg[igen] == -6) tbarquark = true;
    }
    if(bbarquark == false || tbarquark == false){
      cout << "Event " << entry << ". No bbar-quark = " << std::boolalpha << bquark << " or tbar-quark " << tquark << endl;
      cout << "Drop event........................" << endl;
      continue;
    }






    std::vector<int> indexb, indexbbar, indexbjet, indexbbarjet;

    for(int igen = 0; igen < nGen; igen++){
      if(genPdg[igen] == 5) indexb.push_back(igen);
      if(genPdg[igen] == -5) indexbbar.push_back(igen);
    }

    //cout << "Event " << entry << endl;
    //for(int i = 0; i < (int)indexbbar.size(); i++){
    //  cout << indexbbar[i] << " (Eta,Phi,Pt) = (" << genEta[indexbbar[i]] <<", "<<genPhi[indexbbar[i]]<<", "<<genPt[indexbbar[i]]<<")" <<endl;
    //}







    bool possible2b = false;
    int i = (int)indexb.size() - 2;
    if((int)indexb.size() > 1){
      if(genPt[indexb[i]] - genPt[indexb[i+1]] < -10.) possible2b = true;
      if(TMath::Abs(genPhi[indexb[i+1]] - genPhi[indexb[i]]) > 0.5) possible2b = true;
      if(TMath::Abs(genEta[indexb[i+1]] - genEta[indexb[i]]) > 0.5) possible2b = true;
    }

    bool possible2bbar = false;
    i = (int)indexbbar.size() - 2;
    if((int)indexbbar.size() > 1){
      if(genPt[indexbbar[i]] - genPt[indexbbar[i+1]] < -10.) possible2bbar = true;
      if(TMath::Abs(genPhi[indexbbar[i+1]] - genPhi[indexbbar[i]]) > 0.5) possible2bbar = true;
      if(TMath::Abs(genEta[indexbbar[i+1]] - genEta[indexbbar[i]]) > 0.5) possible2bbar = true;
    }






    int bquark1, bquark2, bbarquark1, bbarquark2;
    if(possible2b == true){ 
      //cout << "Event " << entry << ". There are possibly multiple b-quarks..." << endl;
      doubleb++;
      bquark1 = indexb[(int)indexb.size() - 1];
      bquark2 = indexb[(int)indexb.size() - 2];
    } else {
      bquark1 = indexb[(int)indexb.size() - 1];
    }
    if(possible2bbar == true){ 
      //cout << "Event " << entry << ". There are possibly multiple bbar-quarks..." << endl;
      doublebbar++;
      bbarquark1 = indexbbar[(int)indexbbar.size() - 1];
      bbarquark2 = indexbbar[(int)indexbbar.size() - 2];
    } else {
      bbarquark1 = indexbbar[(int)indexbbar.size() - 1];
    }


    /*if(genPt[bquark1] < 30.){
      if(possible2b == true){
        bquark1 = bquark2;
        possible2b = false;
        doubleb--;
      } else {
        continue;
      }
    }
    if(genPt[bbarquark1] < 30.){
      if(possible2bbar == true){
        bbarquark1 = bbarquark2;
        possible2bbar = false;
        doublebbar--;
      } else {
        continue;
      }
    }*/




    int bjet1(-999), bjet2(-999), bbarjet1(-999), bbarjet2(-999);
    double temp_dR = 999.;
    for(int ijet = 0; ijet < nJt; ijet++){

      double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bquark1] - jtPhi[ijet]) );
      double deltaeta = TMath::Abs( genEta[bquark1] - jtEta[ijet] );
      double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bquark1]-jtEta[ijet],2) );

      if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;

      if(drJetToGen < temp_dR){
        bjet1 = ijet;
        temp_dR = drJetToGen;
      }
      //cout << "Event " << entry << " bjet1 " << ijet << endl;
      //cout << "         drJetToGen = " << drJetToGen << endl;
      //cout << "         (Eta,Phi,Pt) b, jet = (" << genEta[bquark1] <<", "<<genPhi[bquark1]<<", "<<genPt[bquark1]<<"), ("<<jtEta[ijet]<<", "<<jtPhi[ijet]<<", "<<jtPt[ijet]<<")"<<endl;
    }
    temp_dR = 999.;
    if(possible2b == true){
      for(int ijet = 0; ijet < nJt; ijet++){

        double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bquark2] - jtPhi[ijet]) );
        double deltaeta = TMath::Abs( genEta[bquark2] - jtEta[ijet] );
        double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bquark2]-jtEta[ijet],2) );

        if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;
        
        if(drJetToGen < temp_dR){
          bjet2 = ijet;
          temp_dR = drJetToGen;
        }
        //cout << "Event " << entry << " bjet2 " << ijet << endl;
        //cout << "         drJetToGen = " << drJetToGen << endl;
        //cout << "         (Eta,Phi,Pt) b, jet = (" << genEta[bquark2] <<", "<<genPhi[bquark2]<<", "<<genPt[bquark2]<<"), ("<<jtEta[ijet]<<", "<<jtPhi[ijet]<<", "<<jtPt[ijet]<<")"<<endl;
      }
    }

    temp_dR = 999.;
    for(int ijet = 0; ijet < nJt; ijet++){

      double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bbarquark1] - jtPhi[ijet]) );
      double deltaeta = TMath::Abs( genEta[bbarquark1] - jtEta[ijet] );
      double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bbarquark1]-jtEta[ijet],2) );

      if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;
      
      if(drJetToGen < temp_dR){
        bbarjet1 = ijet;
        temp_dR = drJetToGen;
      }
      //cout << "Event " << entry << " bbarjet1 " << ijet << endl;
      //cout << "         drJetToGen = " << drJetToGen << endl;
      //cout << "         (Eta,Phi,Pt) b, jet = (" << genEta[bbarquark1] <<", "<<genPhi[bbarquark1]<<", "<<genPt[bbarquark1]<<"), ("<<jtEta[ijet]<<", "<<jtPhi[ijet]<<", "<<jtPt[ijet]<<")"<<endl;
    }

    temp_dR = 999.;
    if(possible2bbar == true){
      for(int ijet = 0; ijet < nJt; ijet++){

        double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(genPhi[bbarquark2] - jtPhi[ijet]) );
        double deltaeta = TMath::Abs( genEta[bbarquark2] - jtEta[ijet] );
        double drJetToGen = sqrt( pow( deltaphi ,2) + pow(genEta[bbarquark2]-jtEta[ijet],2) );

        if(drJetToGen > 0.3 || deltaeta > 0.3 || deltaphi > 0.3) continue;
       
        if(drJetToGen < temp_dR){
          bbarjet2 = ijet;
          temp_dR = drJetToGen;
        }
        //cout << "Event " << entry << " bbarjet2 " << ijet << endl;
        //cout << "         drJetToGen = " << drJetToGen << endl;
        //cout << "         (Eta,Phi,Pt) b, jet = (" << genEta[bbarquark2] <<", "<<genPhi[bbarquark2]<<", "<<genPt[bbarquark2]<<"), ("<<jtEta[ijet]<<", "<<jtPhi[ijet]<<", "<<jtPt[ijet]<<")"<<endl;
      }
    }




    if(bjet1 == -999){
      missingbj1++;
      h_pT_didntmatchedMC[0]->Fill(genPt[bquark1]);
    /*  cout << "Event " << entry << endl;
      for(int igen = 0; igen < nGen; igen++){
        if(genPdg[igen] == 5){
          cout << "         (Eta,Phi,Pt) b = (" << genEta[igen] <<", "<<genPhi[igen]<<", "<<genPt[igen]<<")"<<endl;
        }
      }
      for(int ijet = 0; ijet < nJt; ijet++){
        cout << "         (Eta,Phi,Pt) jet = (" << jtEta[ijet]<<", "<<jtPhi[ijet]<<", "<<jtPt[ijet]<<")"<<endl;
      }
      cout << endl;*/
    } else if(bjet1 != -999) {
      h_pT_matchedMC[0]->Fill(genPt[bquark1]);
    }
    if(bbarjet1 == -999){
      missingbbarj1++;
      h_pT_didntmatchedMC[2]->Fill(genPt[bbarquark1]);
    } else if(bbarjet1 != -999) {
      h_pT_matchedMC[2]->Fill(genPt[bbarquark1]);
    }
    if(possible2b == true && bjet2 == -999){
      missingbj2++;
      h_pT_didntmatchedMC[1]->Fill(genPt[bquark2]);
    } else if(possible2b == true && bjet2 != -999){
      h_pT_matchedMC[1]->Fill(genPt[bquark2]);
    }
    if(possible2bbar == true && bbarjet2 == -999){ 
      missingbbarj2++;
      h_pT_didntmatchedMC[3]->Fill(genPt[bbarquark2]);
    } else if(possible2bbar == true && bbarjet2 != -999){
      h_pT_matchedMC[3]->Fill(genPt[bbarquark2]);
    }


    //if((int)bjet1.size() > 1 || (int)bjet2.size() > 1 || (int)bbarjet1.size() > 1 || (int)bbarjet2.size() > 1){
    //if((int)bjet1.size() != 1 || (int)bbarjet1.size() != 1){
    if(bjet1 == -999 ||bbarjet1 == -999){
      //cout << ibj << " Multiple options for jets in event " << entry << " (" << (int)bjet1.size() << " " << (int)bjet2.size() << " " << (int)bbarjet1.size() << " " << (int)bbarjet2.size() <<")"<< endl;
      //cout << ibj << " Missing option in event " << entry << " (" << bjet1 << " " << bjet2 << " " << bbarjet1 << " " << bbarjet2 <<")"<< endl;
    }
    ibj++;

  }

  cout << "Total events: " << ibj << endl;
  cout << "  Missing bjet1 comb:    " << missingbj1 << endl;
  cout << "  Missing bjet2 comb:    " << missingbj2 << " (" << doubleb << ")" << endl;
  cout << "  Missing bbarjet1 comb: " << missingbbarj1 << endl;
  cout << "  Missing bbarjet2 comb: " << missingbbarj2 << " (" << doublebbar << ")" << endl;


  cout << endl << endl << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  outFile->Write();
  delete outFile;
}




void OpenFiles(TString FileName){
  f_ = TFile::Open(FileName);
  tr_ = dynamic_cast<TTree*>(f_->Get("skimTree"));

  if(f_ == 0){
    cout << "ERROR: Failed to open file " << FileName << endl;
  } else {
    cout << "Opened file: " << FileName << endl;
  }
}