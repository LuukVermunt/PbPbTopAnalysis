#!/bin/bash
cd /afs/cern.ch/user/l/lvermunt/CMSSW_7_5_8_patch3/src
eval `scram r -sh`
cd /afs/cern.ch/user/l/lvermunt/CERNSummerProject/Files/160716/MC_tt
root -b -q "makeMuJetsSkim.C(\"MCtt_v1_0.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v1/merge/HiForest_0.root\",true)"
root -b -q "makeMuJetsSkim.C(\"MCtt_v1_1.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v1/merge/HiForest_1.root\",true)"
root -b -q "makeMuJetsSkim.C(\"MCtt_v1_2.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v1/merge/HiForest_2.root\",true)"
root -b -q "makeMuJetsSkim.C(\"MCtt_v1_3.root\",\"root://eoscms//eos/cms/store/cmst3/group/hintt/CMSSW_7_5_8_patch3/TT1l_172v5_PowhegV2_hvq/Forest/v1/merge/HiForest_3.root\",true)"
