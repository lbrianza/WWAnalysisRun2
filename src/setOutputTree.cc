#include "../interface/setOutputTree.h"

int run;
int event;
int nJets;
int nVtx;
float met;
float leptonPt;
float leptonEta;
float leptonPhi;
float leptonE;
float AK8jetPt[10];
float AK8jetEta[10];
float AK8jetPhi[10];
float AK8jetE[10];

// List of branches
TBranch *b_run;
TBranch *b_event;
TBranch *b_nJets;
TBranch *b_nVtx;
TBranch *b_met;
TBranch *b_leptonPt;
TBranch *b_leptonEta;
TBranch *b_leptonPhi;
TBranch *b_leptonE;
TBranch *b_AK8jetPt;
TBranch *b_AK8jetEta;
TBranch *b_AK8jetPhi;
TBranch *b_AK8jetE;

void init()
{
}

void SetOutTree(TTree* outTree)
{
  outTree->Branch("run",&event,"run/I");
  outTree->Branch("event",&event,"event/I");
  outTree->Branch("nJets",&nJets,"nJets/I");
  outTree->Branch("nVtx",&nVtx,"nVtx/I");
  outTree->Branch("met",&met,"met/F");
  outTree->Branch("leptonPt",&leptonPt,"leptonPt/F");
  outTree->Branch("leptonEta",&leptonEta,"leptonEta/F");
  outTree->Branch("leptonPhi",&leptonPhi,"leptonPhi/F");
  outTree->Branch("leptonE",&leptonE,"leptonE/F");
  outTree->Branch("AK8jetPt",&AK8jetPt,"AK8jetPt[10]/F");
  outTree->Branch("AK8jetEta",&AK8jetEta,"AK8jetEta[10]/F");
  outTree->Branch("AK8jetPhi",&AK8jetPhi,"AK8jetPhi[10]/F");
  outTree->Branch("AK8jetE",&AK8jetE,"AK8jetE[10]/F");
}

void InitRecoTree(TTree* nt)
{
  nt->SetBranchAddress("run", &run, &b_run);
  nt->SetBranchAddress("event", &event, &b_event);
  nt->SetBranchAddress("met", &met, &b_met);
  nt->SetBranchAddress("leptonPt",&leptonPt,&b_leptonPt);
  nt->SetBranchAddress("leptonEta",&leptonEta,&b_leptonEta);
  nt->SetBranchAddress("leptonPhi",&leptonPhi,&b_leptonPhi);
  nt->SetBranchAddress("leptonE",&leptonE,&b_leptonE);
  nt->SetBranchAddress("AK8jetPt",&AK8jetPt,&b_AK8jetPt);
  nt->SetBranchAddress("AK8jetEta",&AK8jetEta,&b_AK8jetEta);
  nt->SetBranchAddress("AK8jetPhi",&AK8jetPhi,&b_AK8jetPhi);
  nt->SetBranchAddress("AK8jetE",&AK8jetE,&b_AK8jetE);
}
