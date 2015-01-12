#ifndef __setOutputTree__
#define __setOutputTree__

#include "TTree.h"
#include "TChain.h"

// Declaration of leaf types
extern int run;
extern int event;
extern int nJets;
extern int nVtx;
extern float met;
extern float met_px;
extern float met_py;
extern float met_pz;
extern float leptonPt;
extern float leptonEta;
extern float leptonPhi;
extern float leptonE;
extern float AK8jetPt[10];
extern float AK8jetEta[10];
extern float AK8jetPhi[10];
extern float AK8jetE[10];
extern float AK8jetPrunedMass[10];
extern float AK8jetTrimmedMass[10];
extern float AK8jetFilteredMass[10];
extern float AK8jetTau1[10];
extern float AK8jetTau2[10];
extern float AK8jetTau3[10];
extern float jetPt[10];
extern float jetEta[10];
extern float jetPhi[10];
extern float jetE[10];

// List of branches
extern TBranch *b_run;
extern TBranch *b_event;
extern TBranch *b_nJets;
extern TBranch *b_nVtx;
extern TBranch *b_met;
extern TBranch *b_met_px;
extern TBranch *b_met_py;
extern TBranch *b_met_pz;
extern TBranch *b_leptonPt;
extern TBranch *b_leptonEta;
extern TBranch *b_leptonPhi;
extern TBranch *b_leptonE;
extern TBranch *b_AK8jetPt;
extern TBranch *b_AK8jetEta;
extern TBranch *b_AK8jetPhi;
extern TBranch *b_AK8jetE;
extern TBranch *b_AK8jetPrunedMass;
extern TBranch *b_AK8jetTrimmedMass;
extern TBranch *b_AK8jetFilteredMass;
extern TBranch *b_AK8jetTau1;
extern TBranch *b_AK8jetTau2;
extern TBranch *b_AK8jetTau3;
extern TBranch *b_jetPt;
extern TBranch *b_jetEta;
extern TBranch *b_jetPhi;
extern TBranch *b_jetE;

void InitRecoTree(TTree* nt);

void init();

void SetOutTree(TTree* outTree);

#endif
