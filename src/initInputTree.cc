#include "../interface/initInputTree.h"

// Declaration of leaf types
unsigned int RunNum;
unsigned int EvtNum;
float MET;

// List of branches
TBranch *b_RunNum; 
TBranch *b_EvtNum; 
TBranch *b_MET; 

void InitTree(TChain* nt)
{
  nt->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
  nt->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
  nt->SetBranchAddress("MET", &MET, &b_MET);
}
