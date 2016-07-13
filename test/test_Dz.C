#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"


#include <string>
#include <vector>

  //Global variables
  unsigned int run_, lumi_;
  ULong64_t evt_;
  int hiBin;
  float vz;
  
  int nLep;
  int lepID[4];
  float lepPt[4];
  float lepPhi[4];
  float lepEta[4];
  int lepChg[4];
  float lepIso[4];
  float lepInnerDz[4];

  int nJt;
  double sum_jtPt;
  int numb_bJets;
  float jetPt[500];
  float jtPhi[500];
  float jtEta[500];
  float jtM[500];
  float discr_csvV1[500];


const int version = 1;


std::vector<int> BuildCSVVectorLeadingbJets(std::vector<int> indexJets);


void test_Dz(const std::string outFileName = Form("test_DZ_120716_v%d.root",version) ) {

	TH1F* h_lepInnerDz_ = new TH1F("h_lepInnerDz_","InnerDz of muon",100,0,20);

	TH1F* h_dz_cent_below_ = new TH1F("h_dz_cent_below_","InnerDz of muon (CSV <= 0.75), Central",100,0,20);
	TH1F* h_dz_cent_above_ = new TH1F("h_dz_cent_above_","InnerDz of muon (CSV > 0.75), Central",100,0,20);

	TH1F* h_dz_peri_below_ = new TH1F("h_dz_peri_below_","InnerDz of muon (CSV <= 0.75), Peripheral",100,0,20);
	TH1F* h_dz_peri_above_ = new TH1F("h_dz_peri_above_","InnerDz of muon (CSV > 0.75), Peripheral",100,0,20);

	TH1F* h_dz_oth_below_ = new TH1F("h_dz_oth_below_","InnerDz of muon (CSV <= 0.75), Other",100,0,20);
	TH1F* h_dz_oth_above_ = new TH1F("h_dz_oth_above_","InnerDz of muon (CSV > 0.75), Other",100,0,20);

	TH2F* h_lepDz_CSV_cent_ = new TH2F("h_lepDz_CSV_cent_","InnerDz of muon vs highest CSV of selected jets (central);lepInnerDz;discr_scvV1",20,0,20,20,0,1);
	TH2F* h_lepDz_CSV_peri_ = new TH2F("h_lepDz_CSV_peri_","InnerDz of muon vs highest CSV of selected jets (peripheral);lepInnerDz;discr_scvV1",20,0,20,20,0,1);


	TFile* f_[4];
  	TTree* tr_[4];

  	for(int i = 0; i < 4; i++){
		TString filename;  		
		filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/MC_tt/dzproblem_v1_%d.root",i);
		f_[i] = TFile::Open(filename);
  		tr_[i] = dynamic_cast<TTree*>(f_[i]->Get("skimTree"));

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
	  	tr_[i]->SetBranchAddress("lepInnerDz", lepInnerDz);

	  	tr_[i]->SetBranchAddress("nJt", &nJt);
	  	tr_[i]->SetBranchAddress("jtPt", jetPt);
	  	tr_[i]->SetBranchAddress("jtPhi", jtPhi);
	  	tr_[i]->SetBranchAddress("jtEta", jtEta);
	  	tr_[i]->SetBranchAddress("jtM", jtM);
		tr_[i]->SetBranchAddress("discr_csvV1", discr_csvV1);

		for(int entry = 0; entry < (int)tr_[i]->GetEntries(); entry++){
  	  		
  	  		tr_[i]->GetEntry(entry);

  	  		//Choose leading muon when there are 2 or more
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
		    	h_lepInnerDz_->Fill( fabs(lepInnerDz[indexMuon]) );

	      		std::vector<int> srt_indexJets = BuildCSVVectorLeadingbJets(indexJets);

		    	if(hiBin < 40){
		    		h_lepDz_CSV_cent_->Fill(lepInnerDz[indexMuon], discr_csvV1[indexJets[srt_indexJets[0]]]);
		    		if(discr_csvV1[indexJets[srt_indexJets[0]]] > 0.75){
		    			h_dz_cent_above_->Fill(lepInnerDz[indexMuon]);
		    		}
		    		else{
		    			h_dz_cent_below_->Fill(lepInnerDz[indexMuon]);
		    		}

		    	} else if(hiBin >120){
					h_lepDz_CSV_peri_->Fill(lepInnerDz[indexMuon], discr_csvV1[indexJets[srt_indexJets[0]]]);
					if(discr_csvV1[indexJets[srt_indexJets[0]]] > 0.75){
		    			h_dz_peri_above_->Fill(lepInnerDz[indexMuon]);
		    		}
		    		else{
		    			h_dz_peri_below_->Fill(lepInnerDz[indexMuon]);
		    		}
		    	} else {
		    		if(discr_csvV1[indexJets[srt_indexJets[0]]] > 0.75){
		    			h_dz_oth_above_->Fill(lepInnerDz[indexMuon]);
		    		}
		    		else{
		    			h_dz_oth_below_->Fill(lepInnerDz[indexMuon]);
		    		}
		    	}
		    }
  	  	}
  	}


	TCanvas* cst = new TCanvas("cst","Histograms",10,10,1800,1000);
	//cst->Divide(2,3);
	cst->Divide(2);
	cst->cd(1);
	//h_dz_cent_below_->Draw();
	//h_dz_peri_below_->Draw();
	h_dz_oth_below_->Draw();
	cst->cd(2);
	//h_dz_cent_above_->Draw();
	//h_dz_peri_above_->Draw();
	h_dz_oth_above_->Draw();
	/*cst->cd(3);
	h_dz_peri_below_->Draw();
	cst->cd(4);
	h_dz_peri_above_->Draw();
	cst->cd(5);
	h_dz_oth_below_->Draw();
	cst->cd(6);
	h_dz_oth_above_->Draw();*/

	TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

	gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
	outFile->Write();
	h_lepInnerDz_->Write();
	h_lepDz_CSV_cent_->Write();
	h_lepDz_CSV_peri_->Write();
	h_dz_cent_below_->Write();
	h_dz_cent_above_->Write();
	h_dz_peri_below_->Write();
	h_dz_peri_above_->Write();
	h_dz_oth_below_->Write();
	h_dz_oth_above_->Write();
	cst->Write();
	delete outFile;
}



std::vector<int> BuildCSVVectorLeadingbJets(std::vector<int> indexJets){
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
