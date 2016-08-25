#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"

//#include "LepJetsSkimTree.h"
#include <string>
#include <vector>

//#include "ForestTreeHeaders/ForestMuons.h"

const bool isDebug = false;

// Jet and lepton selection
const float jetPtCut  = 30.;
const float jetEtaCut = 1.6;//2.;
const int   minNJets  = 4;

std::vector<double> CalculateDeltaEtaAndPhiCS4(int jet02, int jet04);
std::vector<double> CalculateDeltaEtaAndPhiPu4(int jet02, int jet04);


  UInt_t run_, lumi_;
	ULong64_t evt_;
	Int_t hiBin_;
	Float_t vz_;
	Float_t weight_;

  const int maxJets = 5000;
  int           nref;
  float         jtpt[maxJets];   //[nref]
  float         jteta[maxJets];   //[nref]
  float         jtphi[maxJets];   //[nref]
  float         jtm[maxJets];   //[nref]
  float         discr_csvV1[maxJets]; //[nref]
  float         discr_csvV2[maxJets];
  float         discr_tcHighEff[maxJets];
  float         discr_tcHighPur[maxJets];
  float         discr_prob[maxJets];
  float         svtxm[maxJets];
  float         svtxpt[maxJets];
  int         refparton_flavorForB[maxJets];

  int           nrefCs4;
  float         jtptCs4[maxJets];   //[nref]
  float         jtetaCs4[maxJets];   //[nref]
  float         jtphiCs4[maxJets];   //[nref]
  float         jtmCs4[maxJets];   //[nref]
  float         discr_csvV1Cs4[maxJets]; //[nref]
  float         discr_csvV2Cs4[maxJets];
  float         discr_tcHighEffCs4[maxJets];
  float         discr_tcHighPurCs4[maxJets];
  float         discr_probCs4[maxJets];
  float         svtxmCs4[maxJets];
  float         svtxptCs4[maxJets];
  int         refparton_flavorForBCs4[maxJets];

  int           nrefPu4;
  float         jtptPu4[maxJets];   //[nref]
  float         jtetaPu4[maxJets];   //[nref]
  float         jtphiPu4[maxJets];   //[nref]
  float         jtmPu4[maxJets];   //[nref]
  float         discr_csvV1Pu4[maxJets]; //[nref]
  float         discr_csvV2Pu4[maxJets];
  float         discr_tcHighEffPu4[maxJets];
  float         discr_tcHighPurPu4[maxJets];
  float         discr_probPu4[maxJets];
  float         svtxmPu4[maxJets];
  float         svtxptPu4[maxJets];
  int         refparton_flavorForBPu4[maxJets];

  //pf particles pfId, pfPt, pfEta, pfPhi
  std::vector<int>           *pfId = 0;
  std::vector<float>         *pfPt = 0;
  std::vector<float>         *pfEta = 0;
  std::vector<float>         *pfPhi = 0;
    
  //event selections
  int phfCoincFilter = 1;
  int HBHENoiseFilterResult = 1;
  int pprimaryVertexFilter = 1;
  int pcollisionEventSelection = 1;

int hiBinLow = 0;
int hiBinHigh = 200;

void checkROC_ak4PFJetAnalyzer(std::string outFileName = "delete.root", int istart = 0, int iend = 4)
{

  //if(!strcmp(inFileName.c_str(), "")){
  //  std::cout << "No inputs specified. return" << std::endl;
  //  return;
  //}

  TFile* outFile = new TFile(outFileName.c_str(), "RECREATE");

  std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
	for(int i = istart; i < iend; i++){
  	const std::string inFileName = "root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v1/merge/HiForest_"  + std::to_string(i) + ".root";
  	inFileNames_p->push_back(inFileName);
  }
  //inFileNames_p->push_back(inFileName);

//  TChain *lepTree_p = new TChain("ggHiNtuplizer/EventTree");
  TChain *jetTree_p = new TChain("akCs2PFJetAnalyzer/t");
  TChain *jetTree4_Cs_p = new TChain("akCs4PFJetAnalyzer/t");
  TChain *jetTree4_Pu_p = new TChain("akPu4PFJetAnalyzer/t");
  TChain *hiTree_p = new TChain("hiEvtAnalyzer/HiTree");
//  TChain *hltTree_p = new TChain("hltanalysis/HltTree");
  TChain *pfTree_p = new TChain("pfcandAnalyzerCS/pfTree");
  TChain *skimAnaTree_p = new TChain("skimanalysis/HltTree");
  
  TH1F* h_discr_CSVv1[2];
  TH1F* h_discr_CSVv1_CS4[2];
  TH1F* h_discr_CSVv2[2];
  TH1F* h_discr_CSVv2_CS4[2];
  for(int i = 0; i < 2; i++) h_discr_CSVv1[i] = new TH1F(Form("h_discr_CSVv1_%d",i),Form("discr_CSVv1 distribution (PbPb %d-%d);discr_CSVv1;Counts",hiBinLow/2,hiBinHigh/2),100,0,1);
 	for(int i = 0; i < 2; i++) h_discr_CSVv1_CS4[i] = new TH1F(Form("h_discr_CSVv1_CS4%d",i),Form("discr_CSVv1 (Pu4) distribution (PbPb %d-%d);discr_CSVv1;Counts",hiBinLow/2,hiBinHigh/2),100,0,1);
  for(int i = 0; i < 2; i++) h_discr_CSVv2[i] = new TH1F(Form("h_discr_CSVv2_%d",i),Form("discr_CSVv2 distribution (PbPb %d-%d);discr_CSVv2;Counts",hiBinLow/2,hiBinHigh/2),100,0,1);
  for(int i = 0; i < 2; i++) h_discr_CSVv2_CS4[i] = new TH1F(Form("h_discr_CSVv2_CS4%d",i),Form("discr_CSVv2 (Pu4) distribution (PbPb %d-%d);discr_CSVv2;Counts",hiBinLow/2,hiBinHigh/2),100,0,1);

  
  const int nFiles = (int)inFileNames_p->size();

  for(int fileIter = 0; fileIter < nFiles; fileIter++){
    std::cout << "On file: " << fileIter << "/" << nFiles << "  " << inFileNames_p->at(fileIter).c_str()  << std::endl;
//    lepTree_p->Add(inFileNames_p->at(fileIter).c_str());
    jetTree_p->Add(inFileNames_p->at(fileIter).c_str());
    jetTree4_Cs_p->Add(inFileNames_p->at(fileIter).c_str());
    jetTree4_Pu_p->Add(inFileNames_p->at(fileIter).c_str());
    hiTree_p->Add(inFileNames_p->at(fileIter).c_str());
//    hltTree_p->Add(inFileNames_p->at(fileIter).c_str());
    pfTree_p->Add(inFileNames_p->at(fileIter).c_str());
    skimAnaTree_p->Add(inFileNames_p->at(fileIter).c_str());
  }
  
//  ForestMuons fForestMu;

  //lepTree_p->SetBranchStatus("*", 0);
/*  lepTree_p->SetBranchStatus("mu*", 1);
      
  lepTree_p->SetBranchAddress("muPt", &fForestMu.muPt);
  lepTree_p->SetBranchAddress("muPhi", &fForestMu.muPhi);
  lepTree_p->SetBranchAddress("muEta", &fForestMu.muEta);
  lepTree_p->SetBranchAddress("muCharge", &fForestMu.muCharge);
  lepTree_p->SetBranchAddress("muChi2NDF", &fForestMu.muChi2NDF);
  lepTree_p->SetBranchAddress("muInnerD0", &fForestMu.muInnerD0);
  lepTree_p->SetBranchAddress("muInnerDz", &fForestMu.muInnerDz);
  lepTree_p->SetBranchAddress("muMuonHits", &fForestMu.muMuonHits);
  lepTree_p->SetBranchAddress("muStations", &fForestMu.muStations);
  lepTree_p->SetBranchAddress("muTrkLayers", &fForestMu.muTrkLayers);
  lepTree_p->SetBranchAddress("muPixelHits", &fForestMu.muPixelHits);    
*/ 
  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jtphi", 1);
  jetTree_p->SetBranchStatus("jteta", 1);
  jetTree_p->SetBranchStatus("jtm", 1);
  jetTree_p->SetBranchStatus("discr_csvV1", 1);
  jetTree_p->SetBranchStatus("discr_csvV2", 1);
  jetTree_p->SetBranchStatus("discr_tcHighEff", 1);
  jetTree_p->SetBranchStatus("discr_tcHighPur", 1);
  jetTree_p->SetBranchStatus("discr_prob", 1);
  jetTree_p->SetBranchStatus("svtxm", 1);
  jetTree_p->SetBranchStatus("svtxpt", 1);
  jetTree_p->SetBranchStatus("refparton_flavorForB", 1);
        
  jetTree_p->SetBranchAddress("nref", &nref);
  jetTree_p->SetBranchAddress("jtpt", jtpt);
  jetTree_p->SetBranchAddress("jtphi", jtphi);
  jetTree_p->SetBranchAddress("jteta", jteta);
  jetTree_p->SetBranchAddress("jtm", jtm);
  jetTree_p->SetBranchAddress("discr_csvV1", discr_csvV1);
  jetTree_p->SetBranchAddress("discr_csvV2", discr_csvV2);
  jetTree_p->SetBranchAddress("discr_tcHighEff", discr_tcHighEff);
  jetTree_p->SetBranchAddress("discr_tcHighPur", discr_tcHighPur);
  jetTree_p->SetBranchAddress("discr_prob", discr_prob);
  jetTree_p->SetBranchAddress("svtxm", svtxm);
  jetTree_p->SetBranchAddress("svtxpt", svtxpt);
  jetTree_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);


  jetTree4_Cs_p->SetBranchStatus("*", 0);
  jetTree4_Cs_p->SetBranchStatus("nref", 1);
  jetTree4_Cs_p->SetBranchStatus("jtpt", 1);
  jetTree4_Cs_p->SetBranchStatus("jtphi", 1);
  jetTree4_Cs_p->SetBranchStatus("jteta", 1);
  jetTree4_Cs_p->SetBranchStatus("jtm", 1);
  jetTree4_Cs_p->SetBranchStatus("discr_csvV1", 1);
  jetTree4_Cs_p->SetBranchStatus("discr_csvV2", 1);
  jetTree4_Cs_p->SetBranchStatus("discr_tcHighEff", 1);
  jetTree4_Cs_p->SetBranchStatus("discr_tcHighPur", 1);
  jetTree4_Cs_p->SetBranchStatus("discr_prob", 1);
  jetTree4_Cs_p->SetBranchStatus("svtxm", 1);
  jetTree4_Cs_p->SetBranchStatus("svtxpt", 1);
  jetTree4_Cs_p->SetBranchStatus("refparton_flavorForB", 1);
        
  jetTree4_Cs_p->SetBranchAddress("nref", &nrefCs4);
  jetTree4_Cs_p->SetBranchAddress("jtpt", jtptCs4);
  jetTree4_Cs_p->SetBranchAddress("jtphi", jtphiCs4);
  jetTree4_Cs_p->SetBranchAddress("jteta", jtetaCs4);
  jetTree4_Cs_p->SetBranchAddress("jtm", jtmCs4);
  jetTree4_Cs_p->SetBranchAddress("discr_csvV1", discr_csvV1Cs4);
  jetTree4_Cs_p->SetBranchAddress("discr_csvV2", discr_csvV2Cs4);
  jetTree4_Cs_p->SetBranchAddress("discr_tcHighEff", discr_tcHighEffCs4);
  jetTree4_Cs_p->SetBranchAddress("discr_tcHighPur", discr_tcHighPurCs4);
  jetTree4_Cs_p->SetBranchAddress("discr_prob", discr_probCs4);
  jetTree4_Cs_p->SetBranchAddress("svtxm", svtxmCs4);
  jetTree4_Cs_p->SetBranchAddress("svtxpt", svtxptCs4);
  jetTree4_Cs_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForBCs4);

  jetTree4_Pu_p->SetBranchStatus("*", 0);
  jetTree4_Pu_p->SetBranchStatus("nref", 1);
  jetTree4_Pu_p->SetBranchStatus("jtpt", 1);
  jetTree4_Pu_p->SetBranchStatus("jtphi", 1);
  jetTree4_Pu_p->SetBranchStatus("jteta", 1);
  jetTree4_Pu_p->SetBranchStatus("jtm", 1);
  jetTree4_Pu_p->SetBranchStatus("discr_csvV1", 1);
  jetTree4_Pu_p->SetBranchStatus("discr_csvV2", 1);
  jetTree4_Pu_p->SetBranchStatus("discr_tcHighEff", 1);
  jetTree4_Pu_p->SetBranchStatus("discr_tcHighPur", 1);
  jetTree4_Pu_p->SetBranchStatus("discr_prob", 1);
  jetTree4_Pu_p->SetBranchStatus("svtxm", 1);
  jetTree4_Pu_p->SetBranchStatus("svtxpt", 1);
  jetTree4_Pu_p->SetBranchStatus("refparton_flavorForB", 1);
        
  jetTree4_Pu_p->SetBranchAddress("nref", &nrefPu4);
  jetTree4_Pu_p->SetBranchAddress("jtpt", jtptPu4);
  jetTree4_Pu_p->SetBranchAddress("jtphi", jtphiPu4);
  jetTree4_Pu_p->SetBranchAddress("jteta", jtetaPu4);
  jetTree4_Pu_p->SetBranchAddress("jtm", jtmPu4);
  jetTree4_Pu_p->SetBranchAddress("discr_csvV1", discr_csvV1Pu4);
  jetTree4_Pu_p->SetBranchAddress("discr_csvV2", discr_csvV2Pu4);
  jetTree4_Pu_p->SetBranchAddress("discr_tcHighEff", discr_tcHighEffPu4);
  jetTree4_Pu_p->SetBranchAddress("discr_tcHighPur", discr_tcHighPurPu4);
  jetTree4_Pu_p->SetBranchAddress("discr_prob", discr_probPu4);
  jetTree4_Pu_p->SetBranchAddress("svtxm", svtxmPu4);
  jetTree4_Pu_p->SetBranchAddress("svtxpt", svtxptPu4);
  jetTree4_Pu_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForBPu4);

  hiTree_p->SetBranchStatus("*", 0);
  hiTree_p->SetBranchStatus("run", 1);
  hiTree_p->SetBranchStatus("evt", 1);
  hiTree_p->SetBranchStatus("lumi", 1);
  hiTree_p->SetBranchStatus("hiBin", 1);
  hiTree_p->SetBranchStatus("vz", 1);
  hiTree_p->SetBranchStatus("weight", 1);
    
  hiTree_p->SetBranchAddress("run", &run_);
  hiTree_p->SetBranchAddress("evt", &evt_);
  hiTree_p->SetBranchAddress("lumi", &lumi_);
  hiTree_p->SetBranchAddress("hiBin", &hiBin_);
  hiTree_p->SetBranchAddress("vz", &vz_);
  hiTree_p->SetBranchAddress("weight", &weight_);

  pfTree_p->SetBranchAddress("pfId", &pfId);
  pfTree_p->SetBranchAddress("pfPt", &pfPt);
  pfTree_p->SetBranchAddress("pfEta", &pfEta);
  pfTree_p->SetBranchAddress("pfPhi", &pfPhi);
    
//  hltTree_p->SetBranchStatus("HLT_HIL2Mu15_v2",1);
//  hltTree_p->SetBranchAddress("HLT_HIL2Mu15_v2",&trig);

  skimAnaTree_p->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
  skimAnaTree_p->SetBranchAddress("HBHENoiseFilterResult",&HBHENoiseFilterResult);
  skimAnaTree_p->SetBranchAddress("pprimaryVertexFilter",&pprimaryVertexFilter);
  skimAnaTree_p->SetBranchAddress("pcollisionEventSelection",&pcollisionEventSelection);

  int nEntries = (int)jetTree_p->GetEntries();
  //nEntries = 5000;
  int entryDiv = ((int)(nEntries/20));

  int matches(0), nomatches(0);

  for(int entry = 0; entry < nEntries; entry++){

  	if(entry%entryDiv == 0 && nEntries >= 10000) std::cout << "Entry # " << entry << "/" << nEntries << std::endl;

  	hiTree_p->GetEntry(entry);
    //lepTree_p->GetEntry(entry);
    jetTree_p->GetEntry(entry);
    jetTree4_Pu_p->GetEntry(entry);
    jetTree4_Cs_p->GetEntry(entry);
    //hltTree_p->GetEntry(entry);
    pfTree_p->GetEntry(entry);
    skimAnaTree_p->GetEntry(entry);

    if(!phfCoincFilter) continue;
    if(!HBHENoiseFilterResult) continue;
    if(!pcollisionEventSelection) continue;
    if(!pprimaryVertexFilter) continue;
    if(TMath::Abs(vz_) > 15) continue;

    int njets = 0;
    std::vector<int> indexJets;
    for(int jetIter = 0; jetIter < nref; jetIter++){
      if(jtpt[jetIter]<jetPtCut) continue;
      if(fabs(jteta[jetIter])>jetEtaCut) continue;

      indexJets.push_back(jetIter);

      ++njets;
    }
//    if(njets<minNJets) continue;

    std::vector<int> indexLinkedJets;

    for(int i02 = 0; i02 < (int)indexJets.size(); i02++){

    	bool match = false;
    	double eta_temp(0.), phi_temp(0);
    	int imatch;

    	for(int j04 = 0; j04 < nrefPu4; j04++){

    		//std::vector<double> deltaetaphi = CalculateDeltaEtaAndPhiCS4(indexJets[i02], j04);
    		std::vector<double> deltaetaphi = CalculateDeltaEtaAndPhiPu4(indexJets[i02], j04);


    		if(deltaetaphi[0] < 0.4 && deltaetaphi[1] < 0.4){
    			if(match == false){
    				match = true;
    				eta_temp = deltaetaphi[0];
    				phi_temp = deltaetaphi[1];
    				imatch = j04;
    				//cout << "Match found for jet " << indexJets[i02] << " in event " << entry << " (" << deltaetaphi[0] << ", " << deltaetaphi[1] <<")" << endl;
    			} else {
    				//cout << "Multiple matches for jet " << indexJets[i02] << " in event " << entry << endl;
    				//cout << "   ("<< deltaetaphi[0] << ", " << deltaetaphi[1] <<") vs (" << eta_temp << ", " << phi_temp <<")"<<endl;

    				if(eta_temp > deltaetaphi[0]){
    					if(phi_temp < deltaetaphi[1]){ 
    						if(deltaetaphi[1] - phi_temp < 0.5*(eta_temp - deltaetaphi[0]) ){
    							//Keep like this
    						} else {
    							eta_temp = deltaetaphi[0];
    							phi_temp = deltaetaphi[1];
    							imatch = j04;
    						}
    					} else {
    						eta_temp = deltaetaphi[0];
    						phi_temp = deltaetaphi[1];
    						imatch = j04;
    					}
    				} else {
    					//Keep like this
    				}
    			}
    		} else {
    			//cout << "No match for jet " << indexJets[i02] << " in event " << entry << " (" << deltaetaphi[0] << ", " << deltaetaphi[1] <<")" << endl;
    		}

    	}

    	if(match == false){ 
    		//cout << "No match for jet " << indexJets[i02] << " in event " << entry << endl;
    		nomatches++;
    		indexLinkedJets.push_back(-99);
    	} else {
    		matches++;
    		indexLinkedJets.push_back(imatch);
    	}


    }

    for(int i = 0; i < (int)indexJets.size(); i++){
    	if(indexLinkedJets[i] == -99){
    		if(fabs(refparton_flavorForB[indexJets[i]]) == 5){
    			h_discr_CSVv1[0]->Fill(discr_csvV1[indexJets[i]]);
    			h_discr_CSVv1_CS4[0]->Fill(-2);

    			h_discr_CSVv2[0]->Fill(discr_csvV2[indexJets[i]]);
    			h_discr_CSVv2_CS4[0]->Fill(-2);
    		} else {
    			h_discr_CSVv1[1]->Fill(discr_csvV1[indexJets[i]]);
    			h_discr_CSVv1_CS4[1]->Fill(-2);

    			h_discr_CSVv2[1]->Fill(discr_csvV2[indexJets[i]]);
    			h_discr_CSVv2_CS4[1]->Fill(-2);
    		}
    	} else {
    		if(fabs(refparton_flavorForB[indexJets[i]]) == 5){
    			h_discr_CSVv1[0]->Fill(discr_csvV1[indexJets[i]]);
    			//h_discr_CSVv1_CS4[0]->Fill(discr_csvV1Cs4[indexLinkedJets[i]]);
    			h_discr_CSVv1_CS4[0]->Fill(discr_csvV1Pu4[indexLinkedJets[i]]);

    			h_discr_CSVv2[0]->Fill(discr_csvV2[indexJets[i]]);
    			//h_discr_CSVv2_CS4[0]->Fill(discr_csvV2Cs4[indexLinkedJets[i]]);
    			h_discr_CSVv2_CS4[0]->Fill(discr_csvV2Pu4[indexLinkedJets[i]]);
    		} else {
    			h_discr_CSVv1[1]->Fill(discr_csvV1[indexJets[i]]);
    			//h_discr_CSVv1_CS4[1]->Fill(discr_csvV1Cs4[indexLinkedJets[i]]);
    			h_discr_CSVv1_CS4[1]->Fill(discr_csvV1Pu4[indexLinkedJets[i]]);

    			h_discr_CSVv2[1]->Fill(discr_csvV2[indexJets[i]]);
    			//h_discr_CSVv2_CS4[1]->Fill(discr_csvV2Cs4[indexLinkedJets[i]]);
    			h_discr_CSVv2_CS4[1]->Fill(discr_csvV2Pu4[indexLinkedJets[i]]);
    		}
    	}
    }

  }

  cout << endl << "Matches = " << matches << ", No Matches = " << nomatches << endl;


  cout << endl << endl << "Saving histogram(s) to file........" << endl;
  cout << "     Using following name: " << outFileName << endl;
  
  gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
  outFile->Write();
  delete outFile;
}


std::vector<double> CalculateDeltaEtaAndPhiCS4(int jet02, int jet04){

  //Don't take into account the events with multiple W's. Think of something

  double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(jtphi[jet02] - jtphiCs4[jet04]) );
  double deltaeta = TMath::Abs(jteta[jet02] - jtetaCs4[jet04] );

  std::vector<double> EtaAndPhi;
  EtaAndPhi.push_back(deltaeta);
  EtaAndPhi.push_back(deltaphi);

  return EtaAndPhi;

}

std::vector<double> CalculateDeltaEtaAndPhiPu4(int jet02, int jet04){

  //Don't take into account the events with multiple W's. Think of something

  double deltaphi = TMath::Abs(TVector2::Phi_mpi_pi(jtphi[jet02] - jtphiPu4[jet04]) );
  double deltaeta = TMath::Abs(jteta[jet02] - jtetaPu4[jet04] );

  std::vector<double> EtaAndPhi;
  EtaAndPhi.push_back(deltaeta);
  EtaAndPhi.push_back(deltaphi);

  return EtaAndPhi;

}
