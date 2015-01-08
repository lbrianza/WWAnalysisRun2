#ifndef __initInputTree__
#define __initInputTree__

#include "TTree.h"
#include "TChain.h"

// Declaration of leaf types
extern unsigned int RunNum;
extern unsigned int EvtNum;
extern float MET;

// List of branches
extern TBranch *b_RunNum; 
extern TBranch *b_EvtNum; 
extern TBranch *b_MET; 

void InitTree(TChain* nt);

#endif
