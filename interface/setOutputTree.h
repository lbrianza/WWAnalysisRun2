#ifndef __setOutputTree__
#define __setOutputTree__

#include "TTree.h"
#include "TChain.h"

// Declaration of leaf types
extern int run;
extern int event;
extern float met;

// List of branches
extern TBranch *b_run;
extern TBranch *b_event;
extern TBranch *b_met;

void InitRecoTree(TTree* nt);

void init();

void SetOutTree(TTree* outTree);

#endif
