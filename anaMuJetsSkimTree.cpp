#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TVector2.h"

#include <string>
#include <vector>

#include "anaMuJetsSkimTree.h"


anaMuJetsSkimTree::anaMuJetsSkimTree(void):
	CSVCut1_(0.75),
	CSVCut2_(0.75)
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


void anaMuJetsSkimTree::BuildHistograms(void){

  //_0 = PbPb data CMS
  //_1 = MC tt signal
  //_2 = MC W back
  //_3 = MC DY back

  for(int i = 0; i < 4; i++){
  	h_events_[i] = new TH1F(Form("h_events_%d",i),"Number selected events",3,0,3);
  	h_totevents[i] = new TH1F(Form("h_totevents%d",i),"Total SkimTree events",3,0,3);

    h_lepPt_[i] = new TH1F(Form("h_lepPt_%d",i),"pT of muon",30,0,300);
    h_lepPseu_[i] = new TH1F(Form("h_lepPseu_%d",i),"Pseudorapidity of muon",12,0,3);
    h_jetHt_[i] = new TH1F(Form("h_jetHt_%d",i),"H_T of all jets",40,0,1800);
    h_jetCSV_[i] = new TH1F(Form("h_jetCSV_%d",i),"CSV of all jets",100,0,1);
    //h_lepbjetMinv_[i] = new TH1F(Form("h_lepbjetMinv_%d",i),"M_inv of Lepton and b-tagged jets",100,0,300);
    h_lepbjetMinv_min_[i] = new TH1F(Form("h_lepbjetMinv_min_%d",i),"M_inv of Lepton and b-tagged jet (plotted min of each event)",100,0,300);
    //h_ttMinv_[i] = new TH1F(Form("h_ttMinv_%d",i),"M_inv of ttbar",100,0,800);
    h_ttMinv_min_[i] = new TH1F(Form("h_ttMinv_min_%d",i),"M_inv of ttbar (plotted min of each event)",100,0,800);
    h_Phi_1stb_[i] = new TH1F(Form("h_Phi_1stb_%d",i),"Delta-Phi of muon-leading b-jet in opposite hemisphere",30,0,3.2);
    h_Phi_allb_[i] = new TH1F(Form("h_Phi_allb_%d",i),"Delta-Phi of muon-all b-jets in opposite hemisphere",30,0,3.2);
    h_Phi_b_minpi_[i] = new TH1F(Form("h_Phi_b_minpi_%d",i),"Delta-Phi of muon-b-jet most close to pi",30,0,3.2);

    h_2dCSV_[i] = new TH2F(Form("h_2dCSV_%d",i),"CSV of 1st vs 2nd b-jet",110,-0.1,1,110,-0.1,1);
    h_2dCSV2_[i] = new TH2F(Form("h_2dCSV2_%d",i),"CSV of 1st vs 3rd b-jet",110,-0.1,1,110,-0.1,1);
  }

}

void anaMuJetsSkimTree::CalculateNormalizationHistograms(TH1F* h[4], int option){

	//float xsections_[4] = {0.404, 0.02603, 0.64, 0.4256};			           //from MC HiForest files
	float xsections_[4] = {0.404, 0.45*2815.476, 875572.232, 82016.04};    //Pedro's mail (0.45 is BR-factor (only needed for ttbar))
	float nGenerated[4] = {1943558, 95686, 497575, 394186};                //from HiForest files

  //Number of selected events after cuts on analyse niveau
	double nSelected[4] = {h_events_[0]->GetEntries(), h_events_[1]->GetEntries(), h_events_[2]->GetEntries(), h_events_[3]->GetEntries()};
  double Normfactor = xsections_[option] * xsections_[0]/(nGenerated[option]);
	double nEventsExp = Normfactor * nSelected[option];
  
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


void anaMuJetsSkimTree::NormalizeHistograms(int option, int Drawoption){


	if(option == 0){
		//histogram is data, no normalization needed

	} else {
    //histogram is MC, normalization needed

		CalculateNormalizationHistograms(h_lepPt_, option);
		CalculateNormalizationHistograms(h_lepPseu_, option);
		CalculateNormalizationHistograms(h_jetHt_, option);
		
    if(Drawoption == 1)	CalculateNormalizationHistograms(h_jetCSV_, option);

		if(Drawoption == 11){	
      //CalculateNormalizationHistograms(h_lepbjetMinv_, option);
      //CalculateNormalizationHistograms(h_ttMinv_, option);
      CalculateNormalizationHistograms(h_lepbjetMinv_min_, option);
      CalculateNormalizationHistograms(h_ttMinv_min_, option);
    }
	}
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
  
}


void anaMuJetsSkimTree::FillHistogramsBeforeCuts(int option, int i){

  //_0 = PbPb data CMS
  //_1 = MC tt signal
  //_2 = MC W back
  //_3 = MC DY back
  
  cout << "Filling histograms from file: " << i+1 << " of 21." << endl << endl;
  SetAddressBranches(i);

  //Sum over each entry
  for(int entry = 0; entry < (int)tr_[i]->GetEntries(); entry++){
    tr_[i]->GetEntry(entry);

    h_totevents[option]->Fill(1.);


    //Choose leading muon when nLep > 1
    int indexMuon;
    if(nLep > 1){
      double Ptlepfirst = 0.;
      for(int ilep = 0; ilep<nLep; ilep++){
        if(lepPt[ilep] > Ptlepfirst){
          Ptlepfirst = lepPt[ilep];
          indexMuon = ilep;
        }
      }
    } 
    else {
      indexMuon = 0;
    }


    //Cuts on muon
    if(lepPt[indexMuon] < 18) continue;


    //Store the jets that pass the basic cuts
    std::vector<int> indexJets;
    for(int ij = 0; ij < nJt; ij++){
      if(jetPt[ij]<30.) continue;
      if(fabs(jtEta[ij])>2.) continue;

      double drJetToMuon = sqrt(pow(acos(cos(jtPhi[ij]-lepPhi[indexMuon])),2)+pow(jtEta[ij]-lepEta[indexMuon],2));
      if(drJetToMuon<0.3)  continue;

      indexJets.push_back(ij);
    }


    //There have to be at least 4 jets that passed the cuts. True? -> Fill histograms
    if((int)indexJets.size() > 3){

      h_events_[option]->Fill( 1. );
      h_lepPt_[option]->Fill( lepPt[indexMuon] );
      h_lepPseu_[option]->Fill( TMath::Abs(lepEta[indexMuon]) );

      int sum_jtPt = 0;
      for(int l = 0; l < indexJets.size(); l++){
        h_jetCSV_[option]->Fill( discr_csvV1[ indexJets[l] ] );

        sum_jtPt += jetPt[ indexJets[l] ];
      }
      h_jetHt_[option]->Fill(sum_jtPt);

      std::vector<int> srt_indexJets = BuildCSVVectorLeadingbJets(indexJets);

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
}


void anaMuJetsSkimTree::FillHistogramsAfterCuts(int option, int i){

  //_0 = PbPb data CMS
  //_1 = MC tt signal
  //_2 = MC W back
  //_3 = MC DY back
  
  cout << "Filling histograms from file: " << i+1 << " of 21." << endl << endl;
  SetAddressBranches(i);

  for(int entry = 0; entry < (int)tr_[i]->GetEntries(); entry++){
    tr_[i]->GetEntry(entry);

    h_totevents[option]->Fill(1.);


    //Choose leading muon when there are more than 2
    int indexMuon;
    if(nLep > 1){
      double Ptlepfirst = 0.;
      for(int ilep = 0; ilep<nLep; ilep++){
        if(lepPt[ilep] > Ptlepfirst){
          Ptlepfirst = lepPt[ilep];
          indexMuon = ilep;
        }
      }
    } 
    else {
      indexMuon = 0;
    }


    //Cuts on muon
    if(lepPt[indexMuon] < 18) continue;


    //Store the jets that pass the cuts
    std::vector<int> indexJets, indexbJets, indexAllJets;
    for(int ij = 0; ij < nJt; ij++){
      if(jetPt[ij]<30.) continue;
      if(fabs(jtEta[ij])>2.) continue;

      double drJetToMuon = sqrt(pow(acos(cos(jtPhi[ij]-lepPhi[indexMuon])),2)+pow(jtEta[ij]-lepEta[indexMuon],2));
      if(drJetToMuon<0.3)  continue;

      //Cuts on b-jet and save index for jets that passed the cuts
      indexAllJets.push_back(ij);
      if(discr_csvV1[ij] > CSVCut1_){ 
        indexbJets.push_back(ij);
      } else {
        indexJets.push_back(ij);
      }
    }

    //There have to be at least 4 jets that passed the cuts of which at least 2 has to be b-jets
    if((int)indexbJets.size() > 0 && (int)indexJets.size() + (int)indexbJets.size() > 3){

      h_events_[option]->Fill(1.);
      h_lepPt_[option]->Fill(lepPt[indexMuon]);
      h_lepPseu_[option]->Fill( TMath::Abs(lepEta[indexMuon]) );

      int sum_jtPt = 0;
      for(int l = 0; l < indexAllJets.size(); l++) sum_jtPt += jetPt[ indexAllJets[l] ];
      h_jetHt_[option]->Fill(sum_jtPt);


      std::vector<int> srt_indexJets = BuildCSVVectorLeadingbJets(indexbJets);


      //Make control distributions for deltaphi between muon and b-jet
      bool accept = true;
      double dphi_min = 0.;
      for(int i = 0; i < (int) indexbJets.size(); i++ ){
        double dphi = CalculateDeltaPhi(indexMuon, indexbJets[srt_indexJets[i]]);

        if(dphi > 1.57079){ 

          if(accept == true){
            //Leading b-jet in opposite hemisphere
            h_Phi_1stb_[option]->Fill(dphi);
            accept = false;
          }

          //All b-jets in opposite hemisphere
          h_Phi_allb_[option]->Fill(dphi);
        }

        //B-jet which is closest to deltaphi = pi
        if(3.14159265 - dphi < 3.14159265 - dphi_min) dphi_min = dphi;
      }
      h_Phi_b_minpi_[option]->Fill(dphi_min);


      //Make invariant mass control distributions (only plot combination of muon+2bjets+2jets with minimal .M() )
      double minM_lepbjet, minM_tt;
      for(int j = 0; j < (int)indexbJets.size(); j++){

        TLorentzVector LV_lep, LV_bjet, LV_lepbjet;

        LV_lep.SetPtEtaPhiM(lepPt[indexMuon], lepEta[indexMuon], lepPhi[indexMuon], muM);
        LV_bjet.SetPtEtaPhiM(jetPt[ indexbJets[srt_indexJets[j]] ], jtEta[ indexbJets[srt_indexJets[j]] ], jtPhi[ indexbJets[srt_indexJets[j]] ], jtM[ indexbJets[srt_indexJets[j]] ]);
        
        //Calculate minimal .M()
        LV_lepbjet = LV_lep + LV_bjet;
        //h_lepbjetMinv_[option]->Fill( LV_lepbjet.M() );
        if(j == 0)  minM_lepbjet = LV_lepbjet.M();
        else{
          if(minM_lepbjet > LV_lepbjet.M()) minM_lepbjet = LV_lepbjet.M();
        }

        for(int k = j + 1; k < (int)indexbJets.size(); k++){
          for(int l = 0; l < (int)indexAllJets.size(); l++){
            for(int ll = l+1; l < (int)indexAllJets.size(); l++){

              //Check is selected bjets are not the same as selected light-jets
              if(indexAllJets[l] == indexbJets[srt_indexJets[j]] || indexAllJets[l] == indexbJets[srt_indexJets[k]]) continue;
              if(indexAllJets[ll] == indexbJets[srt_indexJets[j]] || indexAllJets[ll] == indexbJets[srt_indexJets[k]]) continue;

              TLorentzVector LV_bjet2, LV_lightjet1, LV_lightjet2, LV_tt;

              LV_bjet2.SetPtEtaPhiM(jetPt[ indexbJets[srt_indexJets[k]] ], jtEta[ indexbJets[srt_indexJets[k]] ], jtPhi[ indexbJets[srt_indexJets[k]] ], jtM[ indexbJets[srt_indexJets[k]] ]);
              LV_lightjet1.SetPtEtaPhiM( jetPt[l], jtEta[l], jtPhi[l], jtM[l]);
              LV_lightjet2.SetPtEtaPhiM( jetPt[ll], jtEta[ll], jtPhi[ll], jtM[ll]);

              //Calculate minimal .M()
              LV_tt = LV_lepbjet + LV_bjet2 + LV_lightjet1 + LV_lightjet2;
              //h_ttMinv_[option]->Fill(LV_tt.M() );
              if(k == 1)  minM_tt = LV_tt.M();
              else{
                if(minM_tt > LV_tt.M()) minM_tt = LV_tt.M();
              }
            }
          }
        }
      }
      h_lepbjetMinv_min_[option]->Fill( minM_lepbjet );
      h_ttMinv_min_[option]->Fill( minM_tt );
    }
  }
}




//Have to clean up the two functions below...

void anaMuJetsSkimTree::LayoutHistograms(int option){

  if(option == 1 || option == 11){
    //Layout for the control distributions with normalization using cross sections

    for(int i = 0; i < 4; i++){      
      NormalizeHistograms(i, option);
    }

  } else if(option == 2){
    //Layout for 2d plots for CSV value of leading b-jets
    
    h_2dCSV_[0]->SetTitle("CSV of 1st vs 2nd b-jet (Data)");
    h_2dCSV_[1]->SetTitle("CSV of 1st vs 2nd b-jet (MC tt)");
    h_2dCSV_[2]->SetTitle("CSV of 1st vs 2nd b-jet (MC W)");
    h_2dCSV_[3]->SetTitle("CSV of 1st vs 2nd b-jet (MC DY)");

    h_2dCSV2_[0]->SetTitle("CSV of 1st vs 3rd b-jet (Data)");
    h_2dCSV2_[1]->SetTitle("CSV of 1st vs 3rd b-jet (MC tt)");
    h_2dCSV2_[2]->SetTitle("CSV of 1st vs 3rd b-jet (MC W)");
    h_2dCSV2_[3]->SetTitle("CSV of 1st vs 3rd b-jet (MC DY)");

  } else if(option == 32){
    //Layout for efficiency plot of the CSV cut

    CSV_effcut->SetName("g_CSV_effcut");
    CSV_effcut->SetTitle("Efficiency CSV cut 1st b-jet");
    CSV_effcut->GetXaxis()->SetTitle("Sig eff");
    CSV_effcut->GetYaxis()->SetTitle("Back rejection eff (W + DY)");

    CSV_effcut->SetMarkerStyle(24);
    CSV_effcut->SetMarkerColor(kRed);

    CSV_effcut2->SetName("g_CSV_effcut2");
    CSV_effcut2->SetTitle("Efficiency CSV cut 2nd b-jet (Cut CSV-b1 = 0.75)");
    CSV_effcut2->GetXaxis()->SetTitle("Sig eff");
    CSV_effcut2->GetYaxis()->SetTitle("Back rejection eff (W + DY)");

    CSV_effcut2->SetMarkerStyle(24);
    CSV_effcut2->SetMarkerColor(kRed);

  } else if(option == 13){
    //Layout for deltaphi muon b-jet control distributions

    //h_Phi_b_minpi_
    //h_Phi_1stb_
    //h_Phi_allb_

  }
  else {
    //Layout for the control distributions with normalization of MC to 1 event

    double scale_lepPt[3]; 
    double scale_lepPseu[3]; 
    double scale_jetHt[3]; 
    double scale_jetCSV[3]; 
    //double scale_lepbjetMinv[3];
    //double scale_ttMinv[3];
    double scale_lepbjetMinv_min[3];
    double scale_ttMinv_min[3];

    for(int i = 0; i < 3; i++){
      scale_lepPt[i] = 1./h_lepPt_[i+1]->Integral();
      scale_lepPseu[i] = 1./h_lepPseu_[i+1]->Integral();
      scale_jetHt[i] = 1./h_jetHt_[i+1]->Integral();
      h_lepPt_[i+1]->Scale(scale_lepPt[i]);
      h_lepPseu_[i+1]->Scale(scale_lepPseu[i]);
      h_jetHt_[i+1]->Scale(scale_jetHt[i]);

      if(option < 10 || option > 19){ 
        scale_jetCSV[i] = 1./h_jetCSV_[i+1]->Integral();
        h_jetCSV_[i+1]->Scale(scale_jetCSV[i]);
      }
      else{
        //scale_lepbjetMinv[i] = 1./h_lepbjetMinv_[i+1]->Integral();
        //scale_ttMinv[i] = 1./h_ttMinv_[i+1]->Integral(); 
        scale_lepbjetMinv_min[i] = 1./h_lepbjetMinv_min_[i+1]->Integral();
        scale_ttMinv_min[i] = 1./h_ttMinv_min_[i+1]->Integral();   
        //h_lepbjetMinv_[i+1]->Scale(scale_lepbjetMinv[i]);
        //h_ttMinv_[i+1]->Scale(scale_ttMinv[i]);
        h_lepbjetMinv_min_[i+1]->Scale(scale_lepbjetMinv_min[i]);
        h_ttMinv_min_[i+1]->Scale(scale_ttMinv_min[i]);
      }
    }
  }

  h_lepPt_[0]->SetMarkerStyle(20);
  h_lepPseu_[0]->SetMarkerStyle(20);
  h_jetHt_[0]->SetMarkerStyle(20);
  h_jetCSV_[0]->SetMarkerStyle(20);
  //h_lepbjetMinv_[0]->SetMarkerStyle(20);
  //h_ttMinv_[0]->SetMarkerStyle(20);
  h_lepbjetMinv_min_[0]->SetMarkerStyle(20);
  h_ttMinv_min_[0]->SetMarkerStyle(20);

  h_lepPt_[1]->SetFillColor(kGray);
  h_lepPseu_[1]->SetFillColor(kGray);
  h_jetHt_[1]->SetFillColor(kGray);
  h_jetCSV_[1]->SetFillColor(kGray);
  //h_lepbjetMinv_[1]->SetFillColor(kGray);
  //h_ttMinv_[1]->SetFillColor(kGray);
  h_lepbjetMinv_min_[1]->SetFillColor(kGray);
  h_ttMinv_min_[1]->SetFillColor(kGray);

  h_lepPt_[2]->SetFillColor(kOrange);
  h_lepPseu_[2]->SetFillColor(kOrange);
  h_jetHt_[2]->SetFillColor(kOrange);
  h_jetCSV_[2]->SetFillColor(kOrange);
  //h_lepbjetMinv_[2]->SetFillColor(kOrange);
  //h_ttMinv_[2]->SetFillColor(kOrange);
  h_lepbjetMinv_min_[2]->SetFillColor(kOrange);
  h_ttMinv_min_[2]->SetFillColor(kOrange);

  h_lepPt_[3]->SetFillColor(42);
  h_lepPseu_[3]->SetFillColor(42);
  h_jetHt_[3]->SetFillColor(42);
  h_jetCSV_[3]->SetFillColor(42);
  //h_lepbjetMinv_[3]->SetFillColor(42);
  //h_ttMinv_[3]->SetFillColor(42);
  h_lepbjetMinv_min_[3]->SetFillColor(42);
  h_ttMinv_min_[3]->SetFillColor(42);
}


TCanvas* anaMuJetsSkimTree::PlotHistograms(int option2){

  TCanvas *cst = new TCanvas("cst","Histograms",10,10,1800,1000);

  lepPt_Stack = new THStack("lepPt_Stack","pT distributions (muon) for 1l+X (at least 4jets)");
  lepPseu_Stack = new THStack("lepPseu_Stack","Pseudorapidity distributions (muon) for 1l+X (at least 4jets)");
  jetHt_Stack = new THStack("jetHt_Stack","H_T distributions for 1l+X (at least 4jets)");
  jetCSV_Stack = new THStack("jetCSV_Stack","CSV distributions for all jets in 1l+X (at least 4jets)");
  //lepbjetMinv_Stack = new THStack("lepbjetMinv_Stack","M_inv distributions of Muon and b-tagged jets (at least 2 b-jets)");
  //ttMinv_Stack = new THStack("ttMinv_Stack","M_inv distributions of tt (at least 4jets with at least 2 b-jets)");
  lepbjetMinv_min_Stack = new THStack("lepbjetMinv_min_Stack","M_inv distributions of Muon and b-tagged jets (Plotted min value & at least 2 b-jets)");
  ttMinv_min_Stack = new THStack("ttMinv_min_Stack","M_inv distributions of tt (Plotted min value & at least 4jets with at least 2 b-jets)");

  TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);

  if(option2 == 1 || option2 == 11){

    LayoutHistograms(option2);

    for(int j = 3; j > 0; j--){
      lepPt_Stack->Add(h_lepPt_[j]);
      lepPseu_Stack->Add(h_lepPseu_[j]);
      jetHt_Stack->Add(h_jetHt_[j]);
      if(option2 == 1) jetCSV_Stack->Add(h_jetCSV_[j]);
      if(option2 == 11){ 
        //lepbjetMinv_Stack->Add(h_lepbjetMinv_[j]);
        //ttMinv_Stack->Add(h_ttMinv_[j]);
        lepbjetMinv_min_Stack->Add(h_lepbjetMinv_min_[j]);
        ttMinv_min_Stack->Add(h_ttMinv_min_[j]);
      }
    }

    leg->AddEntry(h_lepPt_[0],"Data","p");
    leg->AddEntry(h_lepPt_[1],"MC signal (ttbar)","f");
    leg->AddEntry(h_lepPt_[2],"MC back (W)","f");
    leg->AddEntry(h_lepPt_[3],"MC back (DY)","f");

    if(option2 == 1)  cst->Divide(2,2);
    if(option2 == 11)  cst->Divide(3,2);
    cst->cd(1);
    gStyle->SetOptStat(0);
    lepPt_Stack->Draw(); 
    if(lepPt_Stack->GetMaximum() < h_lepPt_[0]->GetMaximum()){
      lepPt_Stack->SetMaximum((int)( h_lepPt_[0]->GetMaximum() + 0.15 * h_lepPt_[0]->GetMaximum() ));
    }
    h_lepPt_[0]->Draw("same,ep");
    leg->Draw();
    cst->cd(2);
    lepPseu_Stack->Draw();
    if(lepPseu_Stack->GetMaximum() < h_lepPseu_[0]->GetMaximum()){
      lepPseu_Stack->SetMaximum((int)( h_lepPseu_[0]->GetMaximum() + 0.15 * h_lepPseu_[0]->GetMaximum() ));
    }
    h_lepPseu_[0]->Draw("same,ep");
    leg->Draw();
    cst->cd(3);
    jetHt_Stack->Draw();
    if(jetHt_Stack->GetMaximum() < h_jetHt_[0]->GetMaximum()){
      jetHt_Stack->SetMaximum((int)( h_jetHt_[0]->GetMaximum() + 0.15 * h_jetHt_[0]->GetMaximum() ));
    }
    h_jetHt_[0]->Draw("same,ep");
    leg->Draw();
    cst->cd(4);
    if(option2 == 1){
      jetCSV_Stack->Draw();
      if(jetCSV_Stack->GetMaximum() < h_jetCSV_[0]->GetMaximum()){
        jetCSV_Stack->SetMaximum((int)( h_jetCSV_[0]->GetMaximum() + 0.15 * h_jetCSV_[0]->GetMaximum() ));
      }
      h_jetCSV_[0]->Draw("same,ep");
      leg->Draw();
    } 
    else {
      lepbjetMinv_min_Stack->Draw();
      if(lepbjetMinv_min_Stack->GetMaximum() < h_lepbjetMinv_min_[0]->GetMaximum()){
        lepbjetMinv_min_Stack->SetMaximum((int)( h_lepbjetMinv_min_[0]->GetMaximum() + 0.15 * h_lepbjetMinv_min_[0]->GetMaximum() ));
      }
      h_lepbjetMinv_min_[0]->Draw("same,ep");
      leg->Draw();
      cst->cd(5);
      ttMinv_min_Stack->Draw();
      if(ttMinv_min_Stack->GetMaximum() < h_ttMinv_min_[0]->GetMaximum()){
        ttMinv_min_Stack->SetMaximum((int)( h_ttMinv_min_[0]->GetMaximum() + 0.15 * h_ttMinv_min_[0]->GetMaximum() ));
      }
      h_ttMinv_min_[0]->Draw("same,ep");
      leg->Draw();
    }
    cst->Update();

  } 
  else if(option2 == 2 || option2 == 32){

    if(option2 == 32){
      const int n = 21;
      double x[n], x2[n], y[n], y2[n];

      //cout << "{";
      for(int c = 0; c < n; c++){

        //if(c > 0) cout << " , ";

        double CSV_cut_eff = 0 + c * 0.05;
        //cout << "CSV_cut_eff = " << CSV_cut_eff << endl;
        int binx1_eff = h_2dCSV_[0]->GetXaxis()->FindBin(CSV_cut_eff);
        int binx2_eff = h_2dCSV_[0]->GetXaxis()->FindBin(1);
        int biny1_eff = h_2dCSV_[0]->GetYaxis()->GetFirst();
        int biny2_eff = h_2dCSV_[0]->GetYaxis()->FindBin(1);

        double EntriesCut[3];
        double TotEntries[3];
        for(int i = 1; i < 4; i++){
          TotEntries[i-1] = h_2dCSV_[i]->GetEntries();
          EntriesCut[i-1] = h_2dCSV_[i]->Integral(binx1_eff, binx2_eff, biny1_eff, biny2_eff);
        }

        x[c] = 100.*EntriesCut[0]/TotEntries[0];
        y[c] = 100. - 100.*(EntriesCut[1] + EntriesCut[2])/(TotEntries[1] + TotEntries[2]);
        //y[c] = 100. - 100.*(EntriesCut[1])/(TotEntries[1]);
        //y[c] = 100. - 100.*(EntriesCut[2])/(TotEntries[2]);

        //cout << "  E_sig = " << 100.*EntriesCut[0]/TotEntries[0] << " E_back,rej = " << 100. - 100.*(EntriesCut[1] + EntriesCut[2])/(TotEntries[1] + TotEntries[2]) << endl;  
        //cout << "{"<<100.*EntriesCut[0]/TotEntries[0]<<", "<<100. - 100.*(EntriesCut[1] + EntriesCut[2])/(TotEntries[1] + TotEntries[2])<<"}";
      }
      //cout << "}" << endl;

      CSV_effcut = new TGraph(n,x,y);

      for(int c = 0; c < n; c++){

        double CSV_cut_eff = 0 + c * 0.05;
        int binx1_eff = h_2dCSV_[0]->GetXaxis()->FindBin(0.75);
        int binx2_eff = h_2dCSV_[0]->GetXaxis()->FindBin(1);
        int biny1_eff = h_2dCSV_[0]->GetYaxis()->FindBin(CSV_cut_eff);
        int biny2_eff = h_2dCSV_[0]->GetYaxis()->FindBin(1);

        double EntriesCut[3];
        double TotEntries[3];
        for(int i = 1; i < 4; i++){
         TotEntries[i-1] = h_2dCSV_[i]->GetEntries();
         EntriesCut[i-1] = h_2dCSV_[i]->Integral(binx1_eff, binx2_eff, biny1_eff, biny2_eff);
        }

        x2[c] = 100.*EntriesCut[0]/TotEntries[0];
        y2[c] = 100. - 100.*(EntriesCut[1] + EntriesCut[2])/(TotEntries[1] + TotEntries[2]);
      }
      CSV_effcut2 = new TGraph(n,x2,y2);
    }

    LayoutHistograms(option2);
    
    TLine *line1 = new TLine(CSVCut2_,CSVCut1_,CSVCut2_,1);
    TLine *line2 = new TLine(CSVCut2_,CSVCut1_,1,CSVCut1_);
    line1->SetLineColor(kRed);
    line2->SetLineColor(kRed);

    int binx1 = h_2dCSV_[0]->GetXaxis()->FindBin(CSVCut1_);
    int binx2 = h_2dCSV_[0]->GetXaxis()->FindBin(1);
    int biny1 = h_2dCSV_[0]->GetYaxis()->FindBin(CSVCut2_);
    int biny2 = h_2dCSV_[0]->GetYaxis()->FindBin(1);

    double perc_cutCSV[4], perc_cutCSV2[4];
    int TotEntries[4], TotEntries2[4];
    for(int i = 0; i < 4; i++){
      TotEntries[i] = (int)h_2dCSV_[i]->GetEntries();
      TotEntries2[i] = (int)h_2dCSV2_[i]->GetEntries();
      perc_cutCSV[i] = 100.*h_2dCSV_[i]->Integral(binx1, binx2, biny1, biny2)/h_2dCSV_[i]->GetEntries();
      perc_cutCSV2[i] = 100.*h_2dCSV2_[i]->Integral(binx1, binx2, biny1, biny2)/h_2dCSV2_[i]->GetEntries();
    }

    TText *t = new TText();
    t->SetTextSize(0.03);

    if(option2 == 2){
      cst->Divide(4,2);
      cst->cd(1);
      gStyle->SetOptStat(0);
      h_2dCSV_[0]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV[0]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries[0]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV_[0]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->cd(2);
      h_2dCSV_[1]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV[1]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries[1]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV_[1]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->cd(3);
      h_2dCSV_[2]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV[2]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries[2]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV_[2]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->cd(4);
      h_2dCSV_[3]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV[3]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries[3]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV_[3]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->cd(5);
      h_2dCSV2_[0]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV2[0]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries2[0]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV2_[0]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->cd(6);
      h_2dCSV2_[1]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV2[1]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries2[1]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV2_[1]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->cd(7);
      h_2dCSV2_[2]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV2[2]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries2[2]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV2_[2]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->cd(8);
      h_2dCSV2_[3]->Draw("COLZ");
      t->DrawText(.1,.75,Form("Percentage in CSV-cut = %f",perc_cutCSV2[3]));
      t->DrawText(.1,.85,Form("Total entries = %d",TotEntries2[3]));
      t->DrawText(.1,.8,Form("Entries in CSV-cut = %d",(int)h_2dCSV2_[3]->Integral(binx1, binx2, biny1, biny2) ));
      line1->Draw();
      line2->Draw();
      cst->Update();
    } else {
      cst->Divide(2);
      cst->cd(1);
      CSV_effcut->Draw("ALP");
      //t->DrawText(1,105,"CSVv1_cuts = 0.0, 0.05, 0.1, ... 1.0");
      t->DrawText(1,101,"CSVv1_cut = 1.0");
      t->DrawText(73,78,"0.75");
      t->DrawText(80,50,"0.5");
      t->DrawText(90,14,"0.25");
      cst->cd(2);
      CSV_effcut2->Draw("ALP");
      t->DrawText(1,100.2,"CSVv1_cut = 1.0");
      t->DrawText(32,96.5,"0.75");
      t->DrawText(44.5,88,"0.5");
      t->DrawText(56.5,81,"0.25");
      cst->Update();
    }
  } 
  else if(option2 == 13){

    //LayoutHistograms(option2);
    
    cst->Divide(4,3);
    cst->cd(1);
    h_Phi_1stb_[0]->Draw();
    cst->cd(2);
    h_Phi_1stb_[1]->Draw();
    cst->cd(3);
    h_Phi_1stb_[2]->Draw();
    cst->cd(4);
    h_Phi_1stb_[3]->Draw();

    cst->cd(5);
    h_Phi_allb_[0]->Draw();
    cst->cd(6);
    h_Phi_allb_[1]->Draw();
    cst->cd(7);
    h_Phi_allb_[2]->Draw();
    cst->cd(8);
    h_Phi_allb_[3]->Draw();

    cst->cd(9);
    h_Phi_b_minpi_[0]->Draw();
    cst->cd(10);
    h_Phi_b_minpi_[1]->Draw();
    cst->cd(11);
    h_Phi_b_minpi_[2]->Draw();
    cst->cd(12);
    h_Phi_b_minpi_[3]->Draw();

    cst->Update();
  }
  else {

    LayoutHistograms(option2);

    for(int j = 3; j > 0; j--){
      lepPt_Stack->Add(h_lepPt_[j]);
      lepPseu_Stack->Add(h_lepPseu_[j]);
      jetHt_Stack->Add(h_jetHt_[j]);
      if(option2 < 10 || option2 > 19){ 
        jetCSV_Stack->Add(h_jetCSV_[j]);
      }
      else{
        //lepbjetMinv_Stack->Add(h_lepbjetMinv_[j]);
        //ttMinv_Stack->Add(h_ttMinv_[j]);
        lepbjetMinv_min_Stack->Add(h_lepbjetMinv_min_[j]);
        ttMinv_min_Stack->Add(h_ttMinv_min_[j]);
      }
    }

    leg->AddEntry(h_lepPt_[1],"MC signal (ttbar)","f");
    leg->AddEntry(h_lepPt_[2],"MC back (W)","f");
    leg->AddEntry(h_lepPt_[3],"MC back (DY)","f");

    if(option2 < 10 || option2 > 19){ 
      cst->Divide(2,2);
    }
    else{
      cst->Divide(3,2);
    }
    cst->cd(1);
    gStyle->SetOptStat(0);
    lepPt_Stack->Draw(); 
    leg->Draw();
    cst->cd(2);
    lepPseu_Stack->Draw();
    leg->Draw();
    cst->cd(3);
    jetHt_Stack->Draw();
    leg->Draw();
    cst->cd(4);
    if(option2 < 10 || option2 > 19){
      jetCSV_Stack->Draw();
      leg->Draw();
    } 
    else {
      lepbjetMinv_min_Stack->Draw();
      leg->Draw();
      cst->cd(5);
      ttMinv_min_Stack->Draw();
      leg->Draw();
    }
    cst->Update();

  }

  return cst;

}