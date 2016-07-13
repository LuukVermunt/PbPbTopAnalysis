#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include "anaMuJetsSkimTree.h"

#include <string>
#include <vector>


const int nFiledata_ = 11;
const int nFileMCsig_ = 4;
const int nFileMCback_ = 3;
const int nFiletotal_ = 21;


/////////////////////////////
// Draw options:
// 1-9 -> Make control plots with only the basic cuts.
// 10-19 -> Make control plots after cuts
// 32 -> Make an efficiency plot of the cut
/////////////////////////////
int Drawoption = 11; 


void AnalyseMuJets(const std::string outFileName = Form("ControlPlotsmuJets_120716_%d.root",Drawoption) ) {

	//CalculateNormalizationHistograms();

	anaMuJetsSkimTree ana;

	TString filename;
	for(int i = 0; i < nFiledata_; i++){
		//filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/topskim/muJetsSkim_%d.root",i);
		filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/Data/muJetsSkim_data_v1_%d.root",i);
		ana.OpenFiles(filename,i);
	}

	for(int i = 0; i < nFileMCsig_; i++){
		//filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/topskim_MC/topskim/muJetsSkim_MCsignal_v2_%d.root",i);
		filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/MC_tt/muJetsSkim_MCtt_v1_%d.root",i);
		ana.OpenFiles(filename, i + nFiledata_);
	}

	for(int i = 0; i < nFileMCback_; i++){
		//filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/topskim_MC_backW/topskim/muJetsSkim_MCbackW_v1_%d.root",i);
		filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/MC_W/muJetsSkim_MCW_v1_%d.root",i);
		ana.OpenFiles(filename,i + nFiledata_ + nFileMCsig_);
	}

	for(int i = 0; i < nFileMCback_; i++){
		//filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/topskim_MC_backDY/topskim/muJetsSkim_MCbackDY_v1_%d.root",i);
		filename.Form("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/120716/MC_DY/muJetsSkim_MCDY_v1_%d.root",i);
		ana.OpenFiles(filename,i + nFiledata_ + nFileMCsig_ + nFileMCback_);
	}


	TFile* outFile = new TFile(outFileName.c_str(),"RECREATE");

	ana.BuildHistograms();

	if(Drawoption < 10 || Drawoption > 19){
		for(int i = 0; i < nFiletotal_; i++){
			if(i < nFiledata_) ana.FillHistogramsBeforeCuts(0,i);
			else if(i >= nFiledata_ && i < nFiledata_ + nFileMCsig_)	ana.FillHistogramsBeforeCuts(1,i);
			else if(i >= nFiledata_ + nFileMCsig_ && i < nFiledata_ + nFileMCsig_ + nFileMCback_) ana.FillHistogramsBeforeCuts(2,i);
			else ana.FillHistogramsBeforeCuts(3,i);
		}
	} else {
		for(int i = 0; i < nFiletotal_; i++){
			if(i < nFiledata_) ana.FillHistogramsAfterCuts(0,i);
			else if(i >= nFiledata_ && i < nFiledata_ + nFileMCsig_) ana.FillHistogramsAfterCuts(1,i);
			else if(i >= nFiledata_ + nFileMCsig_ && i < nFiledata_ + nFileMCsig_ + nFileMCback_) ana.FillHistogramsAfterCuts(2,i);
			else ana.FillHistogramsAfterCuts(3,i);
		}
	}

	TCanvas* cst = ana.PlotHistograms(Drawoption);

	gSystem->cd("/afs/cern.ch/user/l/lvermunt/CERNSummerProject/Analysis");
	outFile->Write();
	cst->Write();
	delete outFile;
}