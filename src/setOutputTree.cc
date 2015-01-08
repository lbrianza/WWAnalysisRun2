#include "../interface/setOutputTree.h"

int run;
int event;
float met;

// List of branches
TBranch *b_run;
TBranch *b_event;
TBranch *b_met;

void init()
{
}

void SetOutTree(TTree* outTree)
{
  outTree->Branch("run",&event,"run/I");
  outTree->Branch("event",&event,"event/I");
  outTree->Branch("met",&met,"met/F");
}

void InitRecoTree(TTree* nt)
{
  nt->SetBranchAddress("run", &run, &b_run);
  nt->SetBranchAddress("event", &event, &b_event);
  nt->SetBranchAddress("met", &met, &b_met);
}
