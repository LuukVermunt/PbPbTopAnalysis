#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include <string>
#include <vector>

#include "anaMuJetsSkimTree.h"
#include "histControlDistributions.h"


anaMuJetsSkimTree::anaMuJetsSkimTree(void)
{
  cout << "anaMuJetsSkimTree object is created" << endl;
}

anaMuJetsSkimTree::~anaMuJetsSkimTree(void)
{
  cout << "anaMuJetsSkimTree object is deleted" << endl;
}


void anaMuJetsSkimTree::OpenFiles(TString FileName, int i){
  f_[i] = TFile::Open(FileName);
  tr_[i] = dynamic_cast<TTree*>(f_[i]->Get("skimTree"));

  if(f_[i] == 0){
  	cout << "ERROR: Failed to open file " << FileName << endl;
  } else {
  	cout << "Opened file: " << FileName << endl;
  }
}

void anaMuJetsSkimTree::SetAddressBranches(int i){

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
  tr_[i]->SetBranchAddress("jtPt", jetPt);
  tr_[i]->SetBranchAddress("jtPhi", jtPhi);
  tr_[i]->SetBranchAddress("jtEta", jtEta);
  tr_[i]->SetBranchAddress("jtM", jtM);
  tr_[i]->SetBranchAddress("discr_csvV1", discr_csvV1); 
  tr_[i]->SetBranchAddress("discr_csvV2", discr_csvV2); 
}

void anaMuJetsSkimTree::BuildHistograms(void){

  //_0 = PbPb data CMS
  //_1 = MC tt signal
  //_2 = MC W back
  //_3 = MC DY back
  //_4 = Multijets (Data, no LepIso cut)

  

  //for(int i = 0; i < 5; i++){
  	
    for(int i = 0; i < 5; i++) h_pTselectbjet_[i] = new TH1F(Form("h_pTselectbjet_%d",i),"p_{T} of tagged b-jet",30,0,100/*300*/);
    for(int i = 0; i < 5; i++) h_ttMinv_min_[i] = new TH1F(Form("h_ttMinv_min_%d",i),"M_{ttbar} (plotted min of each event)",/*100*/70,0,800);
    for(int i = 0; i < 5; i++) h_lepPt_[i] = new TH1F(Form("h_lepPt_%d",i),"p_{T} of muon",30,0,150/*300*/);
    for(int i = 0; i < 5; i++) h_lepPseu_[i] = new TH1F(Form("h_lepPseu_%d",i),"|#eta| of muon",12,0,3);
    for(int i = 0; i < 5; i++) h_jetHt_[i] = new TH1F(Form("h_jetHt_%d",i),"H_{T} of all jets",40,0,1800);
    for(int i = 0; i < 5; i++) h_jetCSV_[i] = new TH1F(Form("h_jetCSV_%d",i),"CSV of all jets",100,0,1);
    for(int i = 0; i < 5; i++) h_lepbjetMinv_min_[i] = new TH1F(Form("h_lepbjetMinv_min_%d",i),"M_{lb} (plotted min of each event)",/*100*/70,0,300);
    for(int i = 0; i < 5; i++) h_qqMinv_min_[i] = new TH1F(Form("h_qqMinv_min_%d",i),"M_{qqbar} (plotted min of each event)",/*100*/70,0,300/*800*/);
    for(int i = 0; i < 5; i++) h_Phi_1stb_[i] = new TH1F(Form("h_Phi_1stb_%d",i),"Delta-Phi of muon - highest CSVv1 tagged b-jet in opposite hemisphere",30,0.5*TMath::Pi(),TMath::Pi());
    for(int i = 0; i < 5; i++) h_Phi_allb_[i] = new TH1F(Form("h_Phi_allb_%d",i),"Delta-Phi of muon - all tagged b-jets in opposite hemisphere",30,0.5*TMath::Pi(),TMath::Pi());
    for(int i = 0; i < 5; i++) h_Phi_b_minpi_[i] = new TH1F(Form("h_Phi_b_minpi_%d",i),"Delta-Phi most close to #pi of muon - a tagged b-jet",30,0,TMath::Pi());

    for(int i = 0; i < 5; i++) h_2dCSV_[i] = new TH2F(Form("h_2dCSV_%d",i),"CSV of 1st vs 2nd b-jet",110,-0.1,1,110,-0.1,1);
    for(int i = 0; i < 5; i++) h_2dCSV2_[i] = new TH2F(Form("h_2dCSV2_%d",i),"CSV of 1st vs 3rd b-jet",110,-0.1,1,110,-0.1,1);

    for(int i = 0; i < 5; i++) h_events_[i] = new TH1F(Form("h_events_%d",i),"Number selected events",3,0,3);
    for(int i = 0; i < 5; i++) h_totevents[i] = new TH1F(Form("h_totevents%d",i),"Total SkimTree events",3,0,3);
    for(int i = 0; i < 5; i++) h_mult_btaggedjets_[i] = new TH1F(Form("h_mult_btaggedjets_%d",i),"Multiplicity b-tagged jets",20,0,20);
    h_smearing_perc_ = new TH1F("h_smearing","p_{T} smearing factor for jets",100,0,1);
  //}
}


void anaMuJetsSkimTree::NormalizeHistograms(int option, int nBjets){


  if(option == 0){
    //histogram is data, no normalization needed

  } else {
    //histogram is MC, normalization needed
    cout << "Normalizing histograms from file " << option << endl;

    CalculateNormalizationHistograms(h_lepPt_, option, nBjets);
    CalculateNormalizationHistograms(h_lepPseu_, option, nBjets);
    CalculateNormalizationHistograms(h_jetHt_, option, nBjets);
    CalculateNormalizationHistograms(h_jetCSV_, option, nBjets);
    CalculateNormalizationHistograms(h_pTselectbjet_, option, nBjets);
    CalculateNormalizationHistograms(h_lepbjetMinv_min_, option, nBjets);
    CalculateNormalizationHistograms(h_ttMinv_min_, option, nBjets);
    CalculateNormalizationHistograms(h_qqMinv_min_, option, nBjets);
    CalculateNormalizationHistograms(h_Phi_1stb_, option, nBjets);
    CalculateNormalizationHistograms(h_Phi_allb_, option, nBjets);
    CalculateNormalizationHistograms(h_Phi_b_minpi_, option, nBjets);
  }
}

void anaMuJetsSkimTree::NormalizeHistogramsToOne(int i, int nBjets){
  double scale_lepPt[3], scale_lepPseu[3], scale_jetHt[3], scale_jetCSV[3], scale_lepbjetMinv_min[3], scale_ttMinv_min[3];

  scale_lepPt[i] = 1./h_lepPt_[i+1]->Integral();
  scale_lepPseu[i] = 1./h_lepPseu_[i+1]->Integral();
  scale_jetHt[i] = 1./h_jetHt_[i+1]->Integral();
  h_lepPt_[i+1]->Scale(scale_lepPt[i]);
  h_lepPseu_[i+1]->Scale(scale_lepPseu[i]);
  h_jetHt_[i+1]->Scale(scale_jetHt[i]);

  if(nBjets == 0){ 
    scale_jetCSV[i] = 1./h_jetCSV_[i+1]->Integral();
    h_jetCSV_[i+1]->Scale(scale_jetCSV[i]);
  }
  else{
    //ADD NEW HISTOGRAMS
    scale_lepbjetMinv_min[i] = 1./h_lepbjetMinv_min_[i+1]->Integral();
    scale_ttMinv_min[i] = 1./h_ttMinv_min_[i+1]->Integral();   
    h_lepbjetMinv_min_[i+1]->Scale(scale_lepbjetMinv_min[i]);
    h_ttMinv_min_[i+1]->Scale(scale_ttMinv_min[i]);
  }
}


void anaMuJetsSkimTree::CalculateNormalizationHistograms(TH1F* h[4], int option, int nBjets){

  //float xsections_[4] = {0.404, 0.02603, 0.64, 0.4256};			                 //from MC HiForest files
  //float xsections_[5] = {0.404, 0.45*2815.476, 875572.232, 82016.04, -999};    
  float xsections_[4] = {0.404, 0.342*2985.216, 928358.912, 86960.64};          //Scaled up pp crosssection with A^2 = 208^2. Factor 0.342 is because MC forcing
  float nGenerated[5] = {1943558, 95686, 497575, 394186, -999};                //from HiForest files

  //Number of selected events after cuts on analyse niveau
  double nSelected[5] = {h_events_[0]->GetEntries(), h_events_[1]->GetEntries(), h_events_[2]->GetEntries(), h_events_[3]->GetEntries(), h_events_[4]->GetEntries()};
  double Normfactor = xsections_[option] * xsections_[0]/(nGenerated[option]);
  double nEventsExp = Normfactor * nSelected[option];

  //Follows from plotABCD.C script
  if(option == 4 && nBjets == 1) nEventsExp = 2617;//(4j1b My smear) 2954;//(3j1b My Smear) 2557;//(My no smear, 4j1b)
  else if(option == 4 && nBjets == 0) nEventsExp = 8071;//(4j0b My smear) 10397;//(3j0b My Smear)
  else if(option == 4 && nBjets > 1) nEventsExp = 711;//(4j2b My smear)
  
  cout << "Option " << option << " nEventsExp " << nEventsExp << " Normfactor " << Normfactor << endl;

  //Normalise MC histograms
  for(int i = h[option]->GetXaxis()->GetFirst(); i < h[option]->GetXaxis()->GetLast() + 1; i++){
    h[option]->SetBinContent(i, (nEventsExp * h[option]->GetBinContent(i))/( h_events_[option]->GetEntries() ) );
  }

  /*  
  //Print number of events/normalized number of events
  cout << h[option]->GetName() << endl;
  cout << " Normfactor option " << option << " " << Normfactor << endl;
  cout << "  Expected events after basic muon-cuts: " << Normfactor * h_totevents[option]->GetEntries() << " (" << h_totevents[option]->GetEntries() << ")" << endl;
  cout << "  Expected events after basic jet-cuts: " << nEventsExp << " (" << nSelected[option] << ")" << endl;

  cout << "  DATA: Expected events after basic muon-cuts: " << h_totevents[0]->GetEntries() << endl;
  cout << "  DATA: Expected events after basic jet-cuts: " << nSelected[0] << endl;
  */
}



double anaMuJetsSkimTree::CalculateDeltaPhi(int ilep, int index_b){
/*Calculates deltaphi (between 0 and pi) between a muon and a jet*/

  double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(lepPhi[ilep] - jtPhi[index_b]) );

  return deltaphi;
}


std::vector<int> anaMuJetsSkimTree::BuildCSVVectorLeadingbJets(std::vector<int> indexJets){
/*Returns a vector with the indices of the jets from highest CSV to lowest CSV*/

  const int array_size = (int) indexJets.size();

  int sortedindex_indexJets[array_size];
  double CSV_array[array_size];

  for(int i = 0; i < (int)indexJets.size(); i++){ 
    sortedindex_indexJets[i] = i;
    CSV_array[i] = discr_csvV1[indexJets[i]];
  }
  TMath::Sort(array_size, CSV_array, sortedindex_indexJets, true);

  std::vector<int> strind_indJets(sortedindex_indexJets, sortedindex_indexJets + array_size);
  return strind_indJets;
}

std::vector<int> anaMuJetsSkimTree::BuildCSVVectorLeadingbJetsCSVv2(std::vector<int> indexJets){
/*Returns a vector with the indices of the jets from highest CSV to lowest CSV*/

  const int array_size = (int) indexJets.size();

  int sortedindex_indexJets[array_size];
  double CSV_array[array_size];

  for(int i = 0; i < (int)indexJets.size(); i++){ 
    sortedindex_indexJets[i] = i;
    CSV_array[i] = discr_csvV2[indexJets[i]];
  }
  TMath::Sort(array_size, CSV_array, sortedindex_indexJets, true);

  std::vector<int> strind_indJets(sortedindex_indexJets, sortedindex_indexJets + array_size);
  return strind_indJets;
}


std::vector<int> anaMuJetsSkimTree::BuildpTVectorLeadingbJets(std::vector<int> indexbJets){
/*Returns a vector with the indices of the jets from highest pT to lowest pT*/

  const int array_size = (int) indexbJets.size();

  int sortedindex_indexbJets[array_size];
  double pT_array[array_size];

  for(int i = 0; i < (int)indexbJets.size(); i++){ 
    sortedindex_indexbJets[i] = i;
    pT_array[i] = jetPt[indexbJets[i]];
  }
  TMath::Sort(array_size, pT_array, sortedindex_indexbJets, true);

  std::vector<int> strind_indJets(sortedindex_indexbJets, sortedindex_indexbJets + array_size);
  return strind_indJets;
}


int anaMuJetsSkimTree::FindLeadingMuonAfterCuts(int option){
  
  
  //Choose leading muon
  int indexMuon = -999;
  double Ptlepfirst = 0.;
  for(int ilep = 0; ilep<nLep; ilep++){
    if(lepPt[ilep] > Ptlepfirst){
      Ptlepfirst = lepPt[ilep];
      indexMuon = ilep;
    }
  }

  if(lepPt[indexMuon] < mupTCut_) return -999;

  if(option == 4){
    if(hiBin<20 && lepIso[indexMuon]<=0.58) return -999;
    else if(hiBin>=20 && hiBin<60 && lepIso[indexMuon]<=0.45) return -999;
    else if(hiBin>=60 && hiBin<100 && lepIso[indexMuon]<=0.3) return -999;
    else if(hiBin>=100 && hiBin<140 && lepIso[indexMuon]<=0.24) return -999;
    else if(hiBin>=140 && lepIso[indexMuon]<=0.18) return -999;
  }
  else {
    if(hiBin<20 && lepIso[indexMuon]>0.58) return -999;
    else if(hiBin>=20 && hiBin<60 && lepIso[indexMuon]>0.45) return -999;
    else if(hiBin>=60 && hiBin<100 && lepIso[indexMuon]>0.3) return -999;
    else if(hiBin>=100 && hiBin<140 && lepIso[indexMuon]>0.24) return -999;
    else if(hiBin>=140 && lepIso[indexMuon]>0.18) return -999;
  }
  
  return indexMuon;
}


std::vector<int> anaMuJetsSkimTree::FindJetsAfterCuts(int indexMuon){

  std::vector<int> indexJets;
  for(int ij = 0; ij < nJt; ij++){

    if(JetSmearing == true){ 
      double smearPercentage = rand->Gaus(0.85,0.03);
      h_smearing_perc_->Fill(smearPercentage);
      //if(jetPt[ij]>=30. && smearPercentage * jetPt[ij] < 30) cout << "Smeared jet from: " << jetPt[ij] << " to " << smearPercentage * jetPt[ij] << endl;
      jetPt[ij] = smearPercentage * jetPt[ij];
    }

    if(jetPt[ij] < jtpTCut_) continue;
    if(fabs(jtEta[ij]) > jtEtaCut_) continue;

    double drJetToMuon = sqrt( pow( CalculateDeltaPhi(indexMuon,ij) ,2) + pow(jtEta[ij]-lepEta[indexMuon],2) ); 
    if(drJetToMuon < drJetToMuonCut_)  continue;

    indexJets.push_back(ij);
  }

  return indexJets;
}


std::vector<int> anaMuJetsSkimTree::FindbJets(std::vector<int> indexJets){
  
  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(discr_csvV1[indexJets[ij]] > CSVCut_)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}

std::vector<int> anaMuJetsSkimTree::FindbJetsCSVv2(std::vector<int> indexJets){
  
  std::vector<int> indexbJets;
 
  for(int ij = 0; ij < (int)indexJets.size(); ij++){
    if(discr_csvV2[indexJets[ij]] > CSVCut_)  indexbJets.push_back(indexJets[ij]);
  }

  return indexbJets;
}


void anaMuJetsSkimTree::FillHistogramsBeforeCuts(int option){

  //_0 = PbPb data CMS
  //_1 = MC tt signal
  //_2 = MC W back
  //_3 = MC DY back
  //_4 = Multijets (Data, no LepIso cut)
  
  cout << "Filling histograms from file: " << option+1 << " of 5." << endl;
  int countEvents = 0;
  SetAddressBranches(option);

  //Sum over each entry
  for(int entry = 0; entry < (int)tr_[option]->GetEntries(); entry++){

    tr_[option]->GetEntry(entry); 
    h_totevents[option]->Fill(1.);


    //Cuts on muon
    int indexMuon = FindLeadingMuonAfterCuts(option);
    if(indexMuon == -999) continue;  
    

    //Store the jets that pass the basic cuts
    std::vector<int> indexJets; 
    //if(option != 0 && option != 4)  indexJets = FindJetsAfterCuts(indexMuon);
    //else                            indexJets = FindJetsAfterCuts(indexMuon);
    if(option != 0 && option != 4){ JetSmearing = true; indexJets = FindJetsAfterCuts(indexMuon); }
    else{                           JetSmearing = false; indexJets = FindJetsAfterCuts(indexMuon);}
    std::vector<int> srt_indexJets = BuildCSVVectorLeadingbJets(indexJets);


    //There have to be at least 4 jets that passed the cuts. True? -> Fill histograms
    if((int)indexJets.size() >= numbJets_){

      countEvents++;

      h_events_[option]->Fill( 1. );
      h_lepPt_[option]->Fill( lepPt[indexMuon] );
      h_lepPseu_[option]->Fill( TMath::Abs(lepEta[indexMuon]) );

      int sum_jtPt = 0;
      for(int l = 0; l < indexJets.size(); l++){
        h_jetCSV_[option]->Fill( discr_csvV1[ indexJets[l] ] );

        sum_jtPt += jetPt[ indexJets[l] ];
      }
      h_jetHt_[option]->Fill(sum_jtPt);


      //2d plot for CSV value for 1st and 2nd b-jet
      if(discr_csvV1[indexJets[srt_indexJets[0]]] < 0){  
        //No jets in this event with a correctly calculated discr_csvV1
        h_2dCSV_[option]->Fill( -0.049, -0.049 );
      } else if(discr_csvV1[indexJets[srt_indexJets[1]]] < 0){
        //Only one jet with a correctly calculated discr_scvV1
        h_2dCSV_[option]->Fill( discr_csvV1[indexJets[srt_indexJets[0]]], -0.049 );
      }
      else{
        h_2dCSV_[option]->Fill( discr_csvV1[indexJets[srt_indexJets[0]]], discr_csvV1[indexJets[srt_indexJets[1]]] );
      } 
     
      //2d plot for CSV value for 1st and 3rd b-jet  
      if(discr_csvV1[indexJets[srt_indexJets[0]]] < 0){  
        //No jets in this event with a correctly calculated discr_csvV1
        h_2dCSV2_[option]->Fill( -0.049, -0.049 );
      } else if(discr_csvV1[indexJets[srt_indexJets[2]]] < 0){
        //Only one (or two) jet with a correctly calculated discr_scvV1
        h_2dCSV2_[option]->Fill( discr_csvV1[indexJets[srt_indexJets[0]]], -0.049 );
      }
      else{
        h_2dCSV2_[option]->Fill( discr_csvV1[indexJets[srt_indexJets[0]]], discr_csvV1[indexJets[srt_indexJets[2]]] );
      }
    }
  }
  cout << " File " << option << " , Events = " << countEvents << endl;
}


void anaMuJetsSkimTree::FillHistogramsAfterCuts(int option){

  //_0 = PbPb data CMS
  //_1 = MC tt signal
  //_2 = MC W back
  //_3 = MC DY back
  //_4 = Multijets (Data, no LepIso cut)
  
  cout << "Filling histograms from file: " << option+1 << " of 5." << endl;
  int countEvents = 0;
  SetAddressBranches(option);

  //Sum over each entry
  for(int entry = 0; entry < (int)tr_[option]->GetEntries(); entry++){

    tr_[option]->GetEntry(entry); 
    h_totevents[option]->Fill(1.);

    //Cuts on muon
    int indexMuon = FindLeadingMuonAfterCuts(option);
    if(indexMuon == -999) continue;  

    //Store the jets that pass the basic cuts
    std::vector<int> indexJets; 
    if(option != 0 && option != 4){ JetSmearing = true; indexJets = FindJetsAfterCuts(indexMuon); }
    else{                           JetSmearing = false; indexJets = FindJetsAfterCuts(indexMuon);}
    std::vector<int> indexbJets = FindbJets(indexJets);
    std::vector<int> srt_indexJets = BuildCSVVectorLeadingbJets(indexJets);


    //There have to be at least 4 jets that passed the cuts of which at least 1 has to be b-jets
    if((int)indexbJets.size() >= numbB_ && (int)indexJets.size() >= numbJets_){

      countEvents++;

      h_events_[option]->Fill(1.);
      h_mult_btaggedjets_[option]->Fill( (int)indexbJets.size() );
      h_lepPt_[option]->Fill(lepPt[indexMuon]);
      h_lepPseu_[option]->Fill( TMath::Abs(lepEta[indexMuon]) );

      int sum_jtPt = 0;
      for(int l = 0; l < indexJets.size(); l++) sum_jtPt += jetPt[ indexJets[l] ];
      h_jetHt_[option]->Fill(sum_jtPt);


      std::vector<int> srt_indexbJets = BuildCSVVectorLeadingbJets( indexbJets );
      std::vector<int> srt_indexJets = BuildCSVVectorLeadingbJets( indexJets );

      h_pTselectbjet_[option]->Fill( jetPt[indexbJets[srt_indexbJets[0]]] );

      //Make control distributions for deltaphi between muon and b-jet
      bool accept = true;
      double dphi_min = 0.;
      for(int i = 0; i < (int) indexbJets.size(); i++ ){
        double dphi = CalculateDeltaPhi(indexMuon, indexbJets[srt_indexbJets[i]]);

        if(dphi > 0.5*TMath::Pi() ){ 
          if(accept == true){
            h_Phi_1stb_[option]->Fill(dphi);  //Leading tagged b-jet in opposite hemisphere
            accept = false;
          }
          h_Phi_allb_[option]->Fill(dphi);  //All tagged b-jets in opposite hemisphere
        }

        if( TMath::Abs( TMath::Pi() - dphi) < TMath::Abs( TMath::Pi() - dphi_min) ) dphi_min = dphi;
      }
      h_Phi_b_minpi_[option]->Fill(dphi_min);  //B-jet which is closest to deltaphi = pi


      double minM_lepbjet(999999.), minM_tt(999999.), minM_qq(999999.);
      for(int j = 0; j < (int)indexbJets.size(); j++){

        int indexleadBjet = indexbJets[srt_indexbJets[j]];

        TLorentzVector LV_lep, LV_bjet, LV_lepbjet;

        LV_lep.SetPtEtaPhiM(lepPt[indexMuon], lepEta[indexMuon], lepPhi[indexMuon], muM);
        LV_bjet.SetPtEtaPhiM(jetPt[ indexleadBjet ], jtEta[ indexleadBjet ], jtPhi[ indexleadBjet ], jtM[ indexleadBjet ]);
        LV_lepbjet = LV_lep + LV_bjet;

        //Calculate minimal .M()
        bool newMinimal_lepbjet_mass = false;
        if(minM_lepbjet > LV_lepbjet.M()){
          minM_lepbjet = LV_lepbjet.M();
          newMinimal_lepbjet_mass = true;
        }
        if(newMinimal_lepbjet_mass == false) continue;
        else minM_tt = minM_qq = 999999.;

        if(numbB_ == 1 && numbJets_ > 3){

          for(int l = 0; l < (int)indexJets.size(); l++){
            for(int ll = l+1; ll < (int)indexJets.size(); ll++){
              for(int lll = ll+1; lll < (int)indexJets.size(); lll++){

                //Check is selected bjets are not the same as selected light-
                if(indexJets[l] == indexleadBjet || indexJets[ll] == indexleadBjet || indexJets[lll] == indexleadBjet) continue;
                
                TLorentzVector LV_lightjet1, LV_lightjet2, LV_lightjet3, LV_tt, LV_qq;

                LV_lightjet1.SetPtEtaPhiM( jetPt[l], jtEta[l], jtPhi[l], jtM[l]);
                LV_lightjet2.SetPtEtaPhiM( jetPt[ll], jtEta[ll], jtPhi[ll], jtM[ll]);
                LV_lightjet3.SetPtEtaPhiM( jetPt[lll], jtEta[lll], jtPhi[lll], jtM[lll]);

                //Calculate minimal .M()
                LV_tt = LV_lep + LV_bjet + LV_lightjet1 + LV_lightjet2 + LV_lightjet3;
                LV_qq = LV_lightjet1 + LV_lightjet2;
                if(minM_tt > LV_tt.M()){ 
//                if(minM_qq > LV_qq.M()){ 
                  minM_tt = LV_tt.M();
                  minM_qq = LV_qq.M();
                }
              }
            }
          }
        } else if(numbB_ == 2 && numbJets_ > 3){
          
          for(int k = 0; k < (int)indexbJets.size(); k++){
            if(indexbJets[k] == indexleadBjet) continue;
            
            for(int l = 0; l < (int)indexJets.size(); l++){
              for(int ll = l+1; ll < (int)indexJets.size(); ll++){

                int index2ndBjet = indexbJets[srt_indexbJets[k]];

                //Check is selected bjets are not the same as selected light-jets
                if(indexJets[l] == indexleadBjet || indexJets[l] == index2ndBjet) continue;
                if(indexJets[ll] == indexleadBjet || indexJets[ll] == index2ndBjet) continue;

                TLorentzVector LV_bjet2, LV_lightjet1, LV_lightjet2, LV_tt, LV_qq;

                LV_bjet2.SetPtEtaPhiM(jetPt[ index2ndBjet ], jtEta[ index2ndBjet ], jtPhi[ index2ndBjet ], jtM[ index2ndBjet ]);
                LV_lightjet1.SetPtEtaPhiM( jetPt[l], jtEta[l], jtPhi[l], jtM[l]);
                LV_lightjet2.SetPtEtaPhiM( jetPt[ll], jtEta[ll], jtPhi[ll], jtM[ll]);
                LV_tt = LV_lepbjet + LV_bjet2 + LV_lightjet1 + LV_lightjet2;
                LV_qq = LV_lightjet1 + LV_lightjet2;

                //Calculate minimal .M()
                if(minM_tt > LV_tt.M()){ 
//                if(minM_qq > LV_qq.M()){ 
                  minM_tt = LV_tt.M();
                  minM_qq = LV_qq.M();
                }
              }
            }
          }
        } else if(numbB_ == 1 && numbJets_ == 3){

          for(int l = 0; l < (int)indexJets.size(); l++){
            for(int ll = l+1; ll < (int)indexJets.size(); ll++){

              //Check is selected bjets are not the same as selected light-
              if(indexJets[l] == indexleadBjet || indexJets[ll] == indexleadBjet) continue;
                
              TLorentzVector LV_lightjet1, LV_lightjet2, LV_qq, LV_tt;

              LV_lightjet1.SetPtEtaPhiM( jetPt[l], jtEta[l], jtPhi[l], jtM[l]);
              LV_lightjet2.SetPtEtaPhiM( jetPt[ll], jtEta[ll], jtPhi[ll], jtM[ll]);

              //Calculate minimal .M()
              LV_qq = LV_lightjet1 + LV_lightjet2;
              LV_tt = LV_lepbjet + LV_qq;
//              if(minM_qq > LV_qq.M()){ 
              if(minM_tt > LV_tt.M()){
                minM_qq = LV_qq.M();
                minM_tt = LV_tt.M();
              }
            }
          }
        }
      }

      if(minM_tt == 999999. || minM_lepbjet == 999999. || minM_qq == 999999.) cout << "Problem, stil 999999. " << (int)indexbJets.size() << endl;

      h_lepbjetMinv_min_[option]->Fill( minM_lepbjet );
      h_ttMinv_min_[option]->Fill( minM_tt );
      h_qqMinv_min_[option]->Fill( minM_qq );
    }
  }
  cout << " File " << option << " , Events = " << countEvents << endl;
}
