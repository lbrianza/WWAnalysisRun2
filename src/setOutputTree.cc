#include "../interface/setOutputTree.h"

int run;
int event;
int nJets;
int nVtx;
float met;
float met_px;
float met_py;
float met_pz;
float leptonPt;
float leptonEta;
float leptonPhi;
float leptonE;
float AK8jetPt[10];
float AK8jetEta[10];
float AK8jetPhi[10];
float AK8jetE[10];
float AK8jetPrunedMass[10];
float AK8jetTrimmedMass[10];
float AK8jetFilteredMass[10];
float AK8jetTau1[10];
float AK8jetTau2[10];
float AK8jetTau3[10];
float jetPt[10];
float jetEta[10];
float jetPhi[10];
float jetE[10];
int genBosonPdgId[10];
float genBosonPt[10];
float genBosonEta[10];
float genBosonPhi[10];
float genBosonE[10];
int genLeptonPdgId[10];
float genLeptonPt[10];
float genLeptonEta[10];
float genLeptonPhi[10];
float genLeptonE[10];
float genNuPt[10];
float genNuEta[10];
float genNuPhi[10];
float genNuE[10];

// List of branches
TBranch *b_run;
TBranch *b_event;
TBranch *b_nJets;
TBranch *b_nVtx;
TBranch *b_met;
TBranch *b_met_px;
TBranch *b_met_py;
TBranch *b_met_pz;
TBranch *b_leptonPt;
TBranch *b_leptonEta;
TBranch *b_leptonPhi;
TBranch *b_leptonE;
TBranch *b_AK8jetPt;
TBranch *b_AK8jetEta;
TBranch *b_AK8jetPhi;
TBranch *b_AK8jetE;
TBranch *b_AK8jetPrunedMass;
TBranch *b_AK8jetTrimmedMass;
TBranch *b_AK8jetFilteredMass;
TBranch *b_AK8jetTau1;
TBranch *b_AK8jetTau2;
TBranch *b_AK8jetTau3;
TBranch *b_jetPt;
TBranch *b_jetEta;
TBranch *b_jetPhi;
TBranch *b_jetE;
TBranch *b_genBosonPdgId;
TBranch *b_genBosonPt;
TBranch *b_genBosonEta;
TBranch *b_genBosonPhi;
TBranch *b_genBosonE;
TBranch *b_genLeptonPdgId;
TBranch *b_genLeptonPt;
TBranch *b_genLeptonEta;
TBranch *b_genLeptonPhi;
TBranch *b_genLeptonE;
TBranch *b_genNuPt;
TBranch *b_genNuEta;
TBranch *b_genNuPhi;
TBranch *b_genNuE;


void init()
{
  run=-999;
  event=-999;
  nJets=-999;
  nVtx=-999;
  met=-999;
  met_px=-999;
  met_py=-999;
  met_pz=-999;
  leptonPt=-999;
  leptonEta=-999;
  leptonPhi=-999;
  leptonE=-999;

 for (int i=0; i<10; i++) {
   AK8jetPt[i]=-999;
   AK8jetEta[i]=-999;
   AK8jetPhi[i]=-999;
   AK8jetE[i]=-999;
   AK8jetPrunedMass[i]=-999;
   AK8jetTrimmedMass[i]=-999;
   AK8jetFilteredMass[i]=-999;
   AK8jetTau1[i]=-999;
   AK8jetTau2[i]=-999;
   AK8jetTau3[i]=-999;
   jetPt[i]=-999;
   jetEta[i]=-999;
   jetPhi[i]=-999;
   jetE[i]=-999;
   genBosonPdgId[i]=-999;
   genBosonPt[i]=-999;
   genBosonEta[i]=-999;
   genBosonPhi[i]=-999;
   genBosonE[i]=-999;
   genLeptonPdgId[i]=-999;
   genLeptonPt[i]=-999;
   genLeptonEta[i]=-999;
   genLeptonPhi[i]=-999;
   genLeptonE[i]=-999;
   genNuPt[i]=-999;
   genNuEta[i]=-999;
   genNuPhi[i]=-999;
   genNuE[i]=-999;
 }
}

void SetOutTree(TTree* outTree)
{
  outTree->Branch("run",&event,"run/I");
  outTree->Branch("event",&event,"event/I");
  outTree->Branch("nJets",&nJets,"nJets/I");
  outTree->Branch("nVtx",&nVtx,"nVtx/I");
  outTree->Branch("met",&met,"met/F");
  outTree->Branch("met_px",&met_px,"met_px/F");
  outTree->Branch("met_py",&met_py,"met_py/F");
  outTree->Branch("met_pz",&met_pz,"met_pz/F");
  outTree->Branch("leptonPt",&leptonPt,"leptonPt/F");
  outTree->Branch("leptonEta",&leptonEta,"leptonEta/F");
  outTree->Branch("leptonPhi",&leptonPhi,"leptonPhi/F");
  outTree->Branch("leptonE",&leptonE,"leptonE/F");
  outTree->Branch("AK8jetPt",&AK8jetPt,"AK8jetPt[10]/F");
  outTree->Branch("AK8jetEta",&AK8jetEta,"AK8jetEta[10]/F");
  outTree->Branch("AK8jetPhi",&AK8jetPhi,"AK8jetPhi[10]/F");
  outTree->Branch("AK8jetE",&AK8jetE,"AK8jetE[10]/F");
  outTree->Branch("AK8jetPrunedMass",&AK8jetPrunedMass,"AK8jetPrunedMass[10]/F");
  outTree->Branch("AK8jetTrimmedMass",&AK8jetTrimmedMass,"AK8jetTrimmedMass[10]/F");
  outTree->Branch("AK8jetFilteredMass",&AK8jetFilteredMass,"AK8jetFilteredMass[10]/F");
  outTree->Branch("AK8jetTau1",&AK8jetTau1,"AK8jetTau1[10]/F");
  outTree->Branch("AK8jetTau2",&AK8jetTau2,"AK8jetTau2[10]/F");
  outTree->Branch("AK8jetTau3",&AK8jetTau3,"AK8jetTau3[10]/F");
  outTree->Branch("jetPt",&jetPt,"jetPt[10]/F");
  outTree->Branch("jetEta",&jetEta,"jetEta[10]/F");
  outTree->Branch("jetPhi",&jetPhi,"jetPhi[10]/F");
  outTree->Branch("jetE",&jetE,"jetE[10]/F");
  outTree->Branch("genBosonPdgId",&genBosonPdgId,"genBosonPdgId[10]/I");
  outTree->Branch("genBosonPt",&genBosonPt,"genBosonPt[10]/F");
  outTree->Branch("genBosonEta",&genBosonEta,"genBosonEta[10]/F");
  outTree->Branch("genBosonPhi",&genBosonPhi,"genBosonPhi[10]/F");
  outTree->Branch("genBosonE",&genBosonE,"genBosonE[10]/F");
  outTree->Branch("genLeptonPdgId",&genLeptonPdgId,"genLeptonPdgId[10]/I");
  outTree->Branch("genLeptonPt",&genLeptonPt,"genLeptonPt[10]/F");
  outTree->Branch("genLeptonEta",&genLeptonEta,"genLeptonEta[10]/F");
  outTree->Branch("genLeptonPhi",&genLeptonPhi,"genLeptonPhi[10]/F");
  outTree->Branch("genLeptonE",&genLeptonE,"genLeptonE[10]/F");
  outTree->Branch("genNuPt",&genNuPt,"genNuPt[10]/F");
  outTree->Branch("genNuEta",&genNuEta,"genNuEta[10]/F");
  outTree->Branch("genNuPhi",&genNuPhi,"genNuPhi[10]/F");
  outTree->Branch("genNuE",&genNuE,"genNuE[10]/F");

}

void InitRecoTree(TTree* nt)
{
  nt->SetBranchAddress("run", &run, &b_run);
  nt->SetBranchAddress("event", &event, &b_event);
  nt->SetBranchAddress("met", &met, &b_met);
  nt->SetBranchAddress("met_px", &met_px, &b_met_px);
  nt->SetBranchAddress("met_py", &met_py, &b_met_py);
  nt->SetBranchAddress("met_pz", &met_pz, &b_met_pz);
  nt->SetBranchAddress("leptonPt",&leptonPt,&b_leptonPt);
  nt->SetBranchAddress("leptonEta",&leptonEta,&b_leptonEta);
  nt->SetBranchAddress("leptonPhi",&leptonPhi,&b_leptonPhi);
  nt->SetBranchAddress("leptonE",&leptonE,&b_leptonE);
  nt->SetBranchAddress("AK8jetPt",&AK8jetPt,&b_AK8jetPt);
  nt->SetBranchAddress("AK8jetEta",&AK8jetEta,&b_AK8jetEta);
  nt->SetBranchAddress("AK8jetPhi",&AK8jetPhi,&b_AK8jetPhi);
  nt->SetBranchAddress("AK8jetE",&AK8jetE,&b_AK8jetE);
  nt->SetBranchAddress("AK8jetPrunedMass",&AK8jetPrunedMass,&b_AK8jetPrunedMass);
  nt->SetBranchAddress("AK8jetTrimmedMass",&AK8jetTrimmedMass,&b_AK8jetTrimmedMass);
  nt->SetBranchAddress("AK8jetFilteredMass",&AK8jetFilteredMass,&b_AK8jetFilteredMass);
  nt->SetBranchAddress("AK8jetTau1",&AK8jetTau1,&b_AK8jetTau1);
  nt->SetBranchAddress("AK8jetTau2",&AK8jetTau2,&b_AK8jetTau2);
  nt->SetBranchAddress("AK8jetTau3",&AK8jetTau3,&b_AK8jetTau3);
  nt->SetBranchAddress("jetPt",&jetPt,&b_jetPt);
  nt->SetBranchAddress("jetEta",&jetEta,&b_jetEta);
  nt->SetBranchAddress("jetPhi",&jetPhi,&b_jetPhi);
  nt->SetBranchAddress("jetE",&jetE,&b_jetE);
  nt->SetBranchAddress("genBosonPdgId",&genBosonPdgId,&b_genBosonPdgId);
  nt->SetBranchAddress("genBosonPt",&genBosonPt,&b_genBosonPt);
  nt->SetBranchAddress("genBosonEta",&genBosonEta,&b_genBosonEta);
  nt->SetBranchAddress("genBosonPhi",&genBosonPhi,&b_genBosonPhi);
  nt->SetBranchAddress("genBosonE",&genBosonE,&b_genBosonE);
  nt->SetBranchAddress("genLeptonPdgId",&genLeptonPdgId,&b_genLeptonPdgId);
  nt->SetBranchAddress("genLeptonPt",&genLeptonPt,&b_genLeptonPt);
  nt->SetBranchAddress("genLeptonEta",&genLeptonEta,&b_genLeptonEta);
  nt->SetBranchAddress("genLeptonPhi",&genLeptonPhi,&b_genLeptonPhi);
  nt->SetBranchAddress("genLeptonE",&genLeptonE,&b_genLeptonE);
  nt->SetBranchAddress("genNuPt",&genNuPt,&b_genNuPt);
  nt->SetBranchAddress("genNuEta",&genNuEta,&b_genNuEta);
  nt->SetBranchAddress("genNuPhi",&genNuPhi,&b_genNuPhi);
  nt->SetBranchAddress("genNuE",&genNuE,&b_genNuE);
}
