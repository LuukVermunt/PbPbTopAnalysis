#include "TFile.h"
#include "TDatime.h"
#include "TNamed.h"
#include "TMath.h"
#include "TLorentzVector.h"
#include "TChain.h"

#include <string>
#include <vector>

void checkbjetsHiForestFiles(const std::string inFileName = "")
{
  if(!strcmp(inFileName.c_str(), "")){
    std::cout << "No inputs specified. return" << std::endl;
    return;
  }

  std::vector<std::string>* inFileNames_p = new std::vector<std::string>;
  inFileNames_p->push_back(inFileName);

  TChain *jetTree_p = new TChain("akCs2PFJetAnalyzer/t");
  TChain *GenPartTree_p = new TChain("HiGenParticleAna/hi");
  
  const int nFiles = (int)inFileNames_p->size();

  for(int fileIter = 0; fileIter < nFiles; fileIter++){
    std::cout << "On file: " << fileIter << "/" << nFiles << "  " << inFileNames_p->at(fileIter).c_str()  << std::endl;
    jetTree_p->Add(inFileNames_p->at(fileIter).c_str());
    GenPartTree_p->Add(inFileNames_p->at(fileIter).c_str());
  }
        
  const int maxJets = 5000;
  int           nref;
  float         jtpt[maxJets];   //[nref]
  float         jteta[maxJets];   //[nref]
  float         jtphi[maxJets];   //[nref]
  float         jtm[maxJets];   //[nref]
  float         discr_csvV1[maxJets]; //[nref]
  int        	refparton_flavorForB[maxJets];
  
  const int maxPart = 5000;
  int 			mult;
  float			pdg[maxPart]; 	//[mult]

  jetTree_p->SetBranchStatus("*", 0);
  jetTree_p->SetBranchStatus("nref", 1);
  jetTree_p->SetBranchStatus("jtpt", 1);
  jetTree_p->SetBranchStatus("jtphi", 1);
  jetTree_p->SetBranchStatus("jteta", 1);
  jetTree_p->SetBranchStatus("jtm", 1);
  jetTree_p->SetBranchStatus("discr_csvV1", 1);
  jetTree_p->SetBranchStatus("refparton_flavorForB", 1);
        
  jetTree_p->SetBranchAddress("nref", &nref);
  jetTree_p->SetBranchAddress("jtpt", jtpt);
  jetTree_p->SetBranchAddress("jtphi", jtphi);
  jetTree_p->SetBranchAddress("jteta", jteta);
  jetTree_p->SetBranchAddress("jtm", jtm);
  jetTree_p->SetBranchAddress("discr_csvV1", discr_csvV1);
  jetTree_p->SetBranchAddress("refparton_flavorForB", refparton_flavorForB);

  GenPartTree_p->SetBranchStatus("*",0);
  GenPartTree_p->SetBranchStatus("pdg",1);
  GenPartTree_p->SetBranchStatus("mult",1);
  GenPartTree_p->SetBranchAddress("pdg",pdg);
  GenPartTree_p->SetBranchAddress("mult",&mult);

  cout << "Entries GenPartTree_p " << (int)GenPartTree_p->GetEntries() << endl;
  int count_mult = 0;
  for(int entry = 0; entry < (int)GenPartTree_p->GetEntries(); entry++){

  	GenPartTree_p->GetEntry(entry);
  	count_mult += mult;
  	
  	cout << "Test ";
  	for(int pd = 0; pd < mult; pd++){
  		cout << pdg[pd] << " ";
  	}
  	cout << endl;

  }
  cout << "Total mult: " << count_mult << endl;

  int nEntries = (int)jetTree_p->GetEntries();
  cout << "There are " << nEntries << " events in this file" << endl;

  int nobjet(0), onebjet(0), twobjet(0), morebjet(0);

  for(int entry = 0; entry < nEntries; entry++){

    jetTree_p->GetEntry(entry);
    //GenPartTree_p->GetEntry(entry);

    int bjet = 0;

    /*cout << "What is in pdg: ";
    for(int ij = 0; ij < 100; ij++){
    	if(fabs(pdg[ij]) < 2000 && fabs(pdg[ij]) != 0){
    		cout << pdg[ij] << " ";
    	}
    }
    cout << endl;*/

    bool isthereabjet = false;
    for(int ij = 0; ij < nref; ij++){
    	if(fabs(refparton_flavorForB[ij]) == 5){ 
    		isthereabjet = true;
    		bjet++;
    	}
    }

    if(isthereabjet == false){	
    	//cout << "No b-jet in this event" << endl;
    	//for(int ij = 0; ij < nref; ij++){
    	//	cout << " " << refparton_flavorForB[ij];
    	//}
    	nobjet++;
    	//	cout << endl;
    } else {
    	if(bjet == 1)		onebjet++;
    	else if(bjet == 2)	twobjet++;
    	else 				morebjet++;
    }
  }

  cout << "In " << nobjet << " events (of the " << nEntries << "), there was no b-jet" << endl;
  cout << "In " << onebjet << " events, there was one b-jet" << endl;
  cout << "In " << twobjet << " events, there were two b-jet" << endl;
  cout << "In " << morebjet << " events, there were more b-jet" << endl;

  return;
}
