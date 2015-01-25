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
float vbf_maxpt_j1_pt;
float vbf_maxpt_j1_eta;
float vbf_maxpt_j1_phi;
float vbf_maxpt_j1_e;
float vbf_maxpt_j1_bDiscriminatorCSV;
float vbf_maxpt_j2_pt;
float vbf_maxpt_j2_eta;
float vbf_maxpt_j2_phi;
float vbf_maxpt_j2_e;
float vbf_maxpt_j2_bDiscriminatorCSV;
float vbf_maxpt_jj_pt;
float vbf_maxpt_jj_eta;
float vbf_maxpt_jj_phi;
float vbf_maxpt_jj_m;


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
  vbf_maxpt_j1_pt=-999;
  vbf_maxpt_j1_eta=-999;
  vbf_maxpt_j1_phi=-999;
  vbf_maxpt_j1_e=-999;
  vbf_maxpt_j1_bDiscriminatorCSV=-999;
  vbf_maxpt_j2_pt=-999;
  vbf_maxpt_j2_eta=-999;
  vbf_maxpt_j2_phi=-999;
  vbf_maxpt_j2_e=-999;
  vbf_maxpt_j2_bDiscriminatorCSV=-999;
  vbf_maxpt_jj_pt=-999;
  vbf_maxpt_jj_eta=-999;
  vbf_maxpt_jj_phi=-999;
  vbf_maxpt_jj_m=-999;
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
  /*  
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
  outTree->Branch("vbf_maxpt_j1_pt",&vbf_maxpt_j1_pt,"vbf_maxpt_j1_pt/F");
  outTree->Branch("vbf_maxpt_j1_eta",&vbf_maxpt_j1_eta,"vbf_maxpt_j1_eta/F");
  outTree->Branch("vbf_maxpt_j1_phi",&vbf_maxpt_j1_phi,"vbf_maxpt_j1_phi/F");
  outTree->Branch("vbf_maxpt_j1_e",&vbf_maxpt_j1_e,"vbf_maxpt_j1_e/F");
  outTree->Branch("vbf_maxpt_j1_bDiscriminatorCSV",&vbf_maxpt_j1_bDiscriminatorCSV,"vbf_maxpt_j1_bDiscriminatorCSV/F");
  outTree->Branch("vbf_maxpt_j2_pt",&vbf_maxpt_j2_pt,"vbf_maxpt_j2_pt/F");
  outTree->Branch("vbf_maxpt_j2_eta",&vbf_maxpt_j2_eta,"vbf_maxpt_j2_eta/F");
  outTree->Branch("vbf_maxpt_j2_phi",&vbf_maxpt_j2_phi,"vbf_maxpt_j2_phi/F");
  outTree->Branch("vbf_maxpt_j2_e",&vbf_maxpt_j2_e,"vbf_maxpt_j2_e/F");
  outTree->Branch("vbf_maxpt_j2_bDiscriminatorCSV",&vbf_maxpt_j2_bDiscriminatorCSV,"vbf_maxpt_j2_bDiscriminatorCSV/F");
  outTree->Branch("vbf_maxpt_jj_pt",&vbf_maxpt_jj_pt,"vbf_maxpt_jj_pt/F");
  outTree->Branch("vbf_maxpt_jj_eta",&vbf_maxpt_jj_eta,"vbf_maxpt_jj_eta/F");
  outTree->Branch("vbf_maxpt_jj_phi",&vbf_maxpt_jj_phi,"vbf_maxpt_jj_phi/F");
  outTree->Branch("vbf_maxpt_jj_m",&vbf_maxpt_jj_m,"vbf_maxpt_jj_m/F");
}

