#ifndef ForestMuons_h
#define ForestMuons_h

#define maxForestMuons 10

#include <iostream>
#include <vector>
#include "TBranch.h"

class ForestMuons {
public :
   ForestMuons(){};
   ~ForestMuons(){};

   // Declaration of leaf types
   // Event info
   Int_t           Run;
   Int_t           Event;
   Int_t           Lumi;
   Int_t           CentBin;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   // GenParticle info
   Int_t           Gen_nptl;
   Int_t           Gen_pid[maxForestMuons];   //[Gen_nptl]
   Int_t           Gen_mom[maxForestMuons];   //[Gen_nptl] pid mother
   Int_t           Gen_status[maxForestMuons];   //[Gen_nptl]
   Float_t         Gen_p[maxForestMuons];   //[Gen_nptl]
   Float_t         Gen_pt[maxForestMuons];   //[Gen_nptl]
   Float_t         Gen_eta[maxForestMuons];   //[Gen_nptl]
   Float_t         Gen_phi[maxForestMuons];   //[Gen_nptl]

   // RecoMuon info
   int             nMu;
   std::vector<float>   *muPt = 0;
   std::vector<float>   *muEta = 0;
   std::vector<float>   *muPhi = 0;
   std::vector<int>     *muCharge = 0;
   std::vector<int>     *muType = 0;
   std::vector<int>     *muIsGood = 0;
   std::vector<float>   *muD0 = 0;
   std::vector<float>   *muDz = 0;
   std::vector<float>   *muChi2NDF = 0;
   std::vector<float>   *muInnerD0 = 0;
   std::vector<float>   *muInnerDz = 0;
   std::vector<int>     *muTrkLayers = 0;
   std::vector<int>     *muPixelLayers = 0;
   std::vector<int>     *muPixelHits = 0;
   std::vector<int>     *muMuonHits = 0;
   std::vector<int>     *muTrkQuality = 0;
   std::vector<int>     *muStations = 0;
   std::vector<float>   *muIsoTrk = 0;
   std::vector<float>   *muPFChIso = 0;
   std::vector<float>   *muPFPhoIso = 0;
   std::vector<float>   *muPFNeuIso = 0;
   std::vector<float>   *muPFPUIso = 0;
   
   // List of branches
   // Event info
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_cbin;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   // GenParticle info
   TBranch        *b_Gen_nptl;   //!
   TBranch        *b_Gen_pid;   //!
   TBranch        *b_Gen_mom;   //!
   TBranch        *b_Gen_status;   //!
   TBranch        *b_Gen_p;   //!
   TBranch        *b_Gen_pt;   //!
   TBranch        *b_Gen_eta;   //!
   TBranch        *b_Gen_phi;   //!
   // RecoMuon info
   TBranch        *b_nMu;   //!
   TBranch        *b_muPt;   //!
   TBranch        *b_muEta;   //!
   TBranch        *b_muPhi;   //!
   TBranch        *b_muCharge;   //!
   TBranch        *b_muType;   //!
   TBranch        *b_muIsGood;   //!
   TBranch        *b_muD0;   //!
   TBranch        *b_muDz;   //!
   TBranch        *b_muChi2NDF;   //!
   TBranch        *b_muInnerD0;   //!
   TBranch        *b_muInnerDz;   //!
   TBranch        *b_muTrkLayers;   //!
   TBranch        *b_muPixelLayers;   //!
   TBranch        *b_muPixelHits;   //!
   TBranch        *b_muMuonHits;   //!
   TBranch        *b_muTrkQuality;   //!
   TBranch        *b_muStations;   //!
   TBranch        *b_muIsoTrk;   //!
   TBranch        *b_muPFChIso;   //!
   TBranch        *b_muPFPhoIso;   //!
   TBranch        *b_muPFNeuIso;   //!
   TBranch        *b_muPFPUIso;   //!
 };
#endif
