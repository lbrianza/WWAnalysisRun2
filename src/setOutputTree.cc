#include "../interface/setOutputTree.h"

// Declaration of leaf types
int event_runNo;
int event;
int njets;
int nPV;
float pfMET;
float pfMET_Phi;
float nu_pz_type0;
float nu_pz_type2;
float l_pt;
float l_eta;
float l_phi;
float l_e;
float ungroomed_jet_pt;
float ungroomed_jet_eta;
float ungroomed_jet_phi;
float ungroomed_jet_e;
float jet_mass_pr;
float jet_mass_tr;
float jet_mass_fi;
float jet_tau2tau1;
/*int genBosonPdgId[10];
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
*/
float deltaR_lak8jet;
float deltaphi_METak8jet;
float deltaphi_Vak8jet;
float v_pt;
float v_eta;
float v_phi;
float v_mt;
float mass_lvj_type0;
float mass_lvj_type2;

// List of branches
TBranch *b_event_runNo;
TBranch *b_event;
TBranch *b_njets;
TBranch *b_nPV;
TBranch *b_pfMET;
TBranch *b_pfMET_Phi;
TBranch *b_nu_pz_type0;
TBranch *b_nu_pz_type2;
TBranch *b_l_pt;
TBranch *b_l_eta;
TBranch *b_l_phi;
TBranch *b_l_e;
TBranch *b_ungroomed_jet_pt;
TBranch *b_ungroomed_jet_eta;
TBranch *b_ungroomed_jet_phi;
TBranch *b_ungroomed_jet_e;
TBranch *b_jet_mass_pr;
TBranch *b_jet_mass_tr;
TBranch *b_jet_mass_fi;
TBranch *b_jet_tau2tau1;
/*TBranch *b_genBosonPdgId[10];
TBranch *b_genBosonPt[10];
TBranch *b_genBosonEta[10];
TBranch *b_genBosonPhi[10];
TBranch *b_genBosonE[10];
TBranch *b_genLeptonPdgId[10];
TBranch *b_genLeptonPt[10];
TBranch *b_genLeptonEta[10];
TBranch *b_genLeptonPhi[10];
TBranch *b_genLeptonE[10];
TBranch *b_genNuPt[10];
TBranch *b_genNuEta[10];
TBranch *b_genNuPhi[10];
TBranch *b_genNuE[10];
*/
TBranch *b_deltaR_lak8jet;
TBranch *b_deltaphi_METak8jet;
TBranch *b_deltaphi_Vak8jet;
TBranch *b_v_pt;
TBranch *b_v_eta;
TBranch *b_v_phi;
TBranch *b_v_mt;
TBranch *b_mass_lvj_type0;
TBranch *b_mass_lvj_type2;

// List of branches

void init()
{
  event_runNo=-999;
  event=-999;
  njets=-999;
  nPV=-999;
  pfMET=-999;
  pfMET_Phi=-999;
  nu_pz_type0=-999;
  nu_pz_type2=-999;
  l_pt=-999;
  l_eta=-999;
  l_phi=-999;
  l_e=-999;
  deltaR_lak8jet=-999;
  deltaphi_METak8jet=-999;
  deltaphi_Vak8jet=-999;
  v_pt=-999;
  v_eta=-999;
  v_phi=-999;
  v_mt=-999;
  mass_lvj_type0=-999;
  mass_lvj_type2=-999;
  ungroomed_jet_pt=-999;
  ungroomed_jet_eta=-999;
  ungroomed_jet_phi=-999;
  ungroomed_jet_e=-999;
  jet_mass_pr=-999;
  jet_mass_tr=-999;
  jet_mass_fi=-999;
  jet_tau2tau1=-999;

   /* for (int i=0; i<10; i++) {
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
   */
}

void SetOutTree(TTree* outTree)
{
  outTree->Branch("event_runNo",&event_runNo,"event_runNo/I");
  outTree->Branch("event",&event,"event/I");
  outTree->Branch("njets",&njets,"njets/I");
  outTree->Branch("nPV",&nPV,"nPV/I");
  outTree->Branch("pfMET",&pfMET,"pfMET/F");
  outTree->Branch("pfMET_Phi",&pfMET_Phi,"pfMET_Phi/F");
  outTree->Branch("nu_pz_type0",&nu_pz_type0,"nu_pz_type0/F");
  outTree->Branch("nu_pz_type2",&nu_pz_type2,"nu_pz_type2/F");
  outTree->Branch("l_pt",&l_pt,"l_pt/F");
  outTree->Branch("l_eta",&l_eta,"l_eta/F");
  outTree->Branch("l_phi",&l_phi,"l_phi/F");
  outTree->Branch("l_e",&l_e,"l_e/F");
  outTree->Branch("ungroomed_jet_pt",&ungroomed_jet_pt,"ungroomed_jet_pt/F");
  outTree->Branch("ungroomed_jet_eta",&ungroomed_jet_eta,"ungroomed_jet_eta/F");
  outTree->Branch("ungroomed_jet_phi",&ungroomed_jet_phi,"ungroomed_jet_phi/F");
  outTree->Branch("ungroomed_jet_e",&ungroomed_jet_e,"ungroomed_jet_e/F");
  outTree->Branch("jet_mass_pr",&jet_mass_pr,"jet_mass_pr");
  outTree->Branch("jet_mass_tr",&jet_mass_tr,"jet_mass_tr");
  outTree->Branch("jet_mass_fi",&jet_mass_fi,"jet_mass_fi");
  outTree->Branch("jet_tau2tau1",&jet_tau2tau1,"jet_tau2tau1");
  /*  outTree->Branch("jetPt",&jetPt,"jetPt[10]/F");
  outTree->Branch("jetEta",&jetEta,"jetEta[10]/F");
  outTree->Branch("jetPhi",&jetPhi,"jetPhi[10]/F");
  outTree->Branch("jetE",&jetE,"jetE[10]/F");
  outTree->Branch("jet_bDiscr",&jet_bDiscr,"jet_bDiscr[10]/F");
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
  */
  outTree->Branch("deltaR_lak8jet",&deltaR_lak8jet,"deltaR_lak8jet/F");
  outTree->Branch("deltaphi_METak8jet",&deltaphi_METak8jet,"deltaphi_METak8jet/F");
  outTree->Branch("deltaphi_Vak8jet",&deltaphi_Vak8jet,"deltaphi_Vak8jet/F");
  outTree->Branch("v_pt",&v_pt,"v_pt/F");
  outTree->Branch("v_eta",&v_eta,"v_eta/F");
  outTree->Branch("v_phi",&v_phi,"v_phi/F");
  outTree->Branch("v_mt",&v_mt,"v_mt/F");
  outTree->Branch("mass_lvj_type0",&mass_lvj_type0,"mass_lvj_type0/F");
  outTree->Branch("mass_lvj_type2",&mass_lvj_type2,"mass_lvj_type2/F");
}

/*
void InitRecoTree(TTree* nt)
{
  nt->SetBranchAddress("run", &run, &b_run);
  nt->SetBranchAddress("event", &event, &b_event);
  nt->SetBranchAddress("met", &met, &b_met);
  nt->SetBranchAddress("met_px", &met_px, &b_met_px);
  nt->SetBranchAddress("met_py", &met_py, &b_met_py);
  nt->SetBranchAddress("met_pz_type0", &met_pz_type0, &b_met_pz_type0);
  nt->SetBranchAddress("met_pz_type2", &met_pz_type2, &b_met_pz_type2);
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
  nt->SetBranchAddress("jet_bDiscr",&jet_bDiscr,&b_jet_bDiscr);
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
  nt->SetBranchAddress("deltaR_lak8jet",&deltaR_lak8jet,&b_deltaR_lak8jet);
  nt->SetBranchAddress("deltaphi_METak8jet",&deltaphi_METak8jet,&b_deltaphi_METak8jet);
  nt->SetBranchAddress("deltaphi_Vak8jet",&deltaphi_Vak8jet,&b_deltaphi_Vak8jet);
  nt->SetBranchAddress("W_pt",&W_pt,&b_W_pt);
  nt->SetBranchAddress("W_eta",&W_eta,&b_W_eta);
  nt->SetBranchAddress("W_phi",&W_phi,&b_W_phi);
  nt->SetBranchAddress("W_E",&W_E,&b_W_E);
  nt->SetBranchAddress("W_mt",&W_mt,&b_W_mt);
  nt->SetBranchAddress("boosted_lvj_m_type0",&boosted_lvj_m_type0,&b_boosted_lvj_m_type0);  
  nt->SetBranchAddress("boosted_lvj_m_type2",&boosted_lvj_m_type2,&b_boosted_lvj_m_type2);  
}
*/
