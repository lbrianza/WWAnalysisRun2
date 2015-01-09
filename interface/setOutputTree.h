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
extern float leptonPt;
extern float leptonEta;
extern float leptonPhi;
extern float leptonE;
extern float AK8jetPt[10];
extern float AK8jetEta[10];
extern float AK8jetPhi[10];
extern float AK8jetE[10];

// List of branches
extern TBranch *b_run;
extern TBranch *b_event;
extern TBranch *b_nJets;
extern TBranch *b_nVtx;
extern TBranch *b_met;
extern TBranch *b_leptonPt;
extern TBranch *b_leptonEta;
extern TBranch *b_leptonPhi;
extern TBranch *b_leptonE;
extern TBranch *b_AK8jetPt;
extern TBranch *b_AK8jetEta;
extern TBranch *b_AK8jetPhi;
extern TBranch *b_AK8jetE;

void InitRecoTree(TTree* nt);

void init();

void SetOutTree(TTree* outTree);

#endif
