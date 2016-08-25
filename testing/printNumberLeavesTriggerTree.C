#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

int printNumberLeavesTriggerTree(){

	std::ifstream input( "events_merge_test.txt" );
	std::ofstream output( "NumberLeaves.txt" );

	int i = 0;
	for( std::string line; getline( input, line ); ){

		output << line << endl;
		if(i%100 == 0) cout << i << " " << endl;

		//TFile* f = TFile::Open("root://eoscms//eos/cms/store/cmst3/group/hintt/mverweij/PbPb5TeV/data/HIEWQExo/crab_FilteredSingleMuHighPt_v3/160809_160502/0000/HiForestAOD_102.root");
		TFile* f = TFile::Open(line.c_str());
		TTree* tr = dynamic_cast<TTree*>(f->Get("hltanalysis/HltTree"));

		TObjArray* ListOfLeaves = tr->GetListOfLeaves();

		output << "   Number of Leaves = " << ListOfLeaves->GetSize() << endl;

		ListOfLeaves->Print();

		delete f;
		i++;

		//if(i > 100) break;
	}

	output.close();
	return 0;
}
