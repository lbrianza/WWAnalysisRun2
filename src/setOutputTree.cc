#include "../interface/setOutputTree.h"

setOutputTree::setOutputTree(TTree* outTree){
  fTree = outTree;
  setBranches();
}

setOutputTree::~setOutputTree(){
  delete fTree;
}

void setOutputTree::initializeVariables()
{
  event_runNo=-999;
  event=-999;
  njets=-999;
  nPV=-999;
  issignal=-999;
  wSampleWeight=-999;
  totalEventWeight=-999;
  eff_and_pu_Weight=-999;
  pfMET=-999;
  pfMET_Phi=-999;
  nu_pz_type0=-999;
  nu_pz_type2=-999;
  l_pt=-999;
  l_eta=-999;
  l_phi=-999;
  l_e=-999;
  W_pt_gen=-999;
  W_pz_gen=-999;
  genGravMass=-999;
  nu_pz_gen=-999;
  deltaR_lak8jet=-999;
  deltaphi_METak8jet=-999;
  deltaphi_Vak8jet=-999;
  v_pt=-999;
  v_eta=-999;
  v_phi=-999;
  v_mt=-999;
  mass_lvj_type0=-999;
  mass_lvj_type2=-999;
  nBTagJet_loose=-999;
  nBTagJet_medium=-999;
  nBTagJet_tight=-999;
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
}

void setOutputTree::setBranches()
{
  fTree->Branch("event_runNo",&event_runNo,"event_runNo/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("njets",&njets,"njets/I");
  fTree->Branch("nPV",&nPV,"nPV/I");
  fTree->Branch("issignal",&issignal,"issignal/I");
  fTree->Branch("wSampleWeight",&wSampleWeight,"wSampleWeight/F");
  fTree->Branch("totalEventWeight",&totalEventWeight,"totalEventWeight/F");
  fTree->Branch("eff_and_pu_Weight",&eff_and_pu_Weight,"eff_and_pu_Weight/F");
  fTree->Branch("pfMET",&pfMET,"pfMET/F");
  fTree->Branch("pfMET_Phi",&pfMET_Phi,"pfMET_Phi/F");
  fTree->Branch("nu_pz_type0",&nu_pz_type0,"nu_pz_type0/F");
  fTree->Branch("nu_pz_type2",&nu_pz_type2,"nu_pz_type2/F");
  fTree->Branch("l_pt",&l_pt,"l_pt/F");
  fTree->Branch("l_eta",&l_eta,"l_eta/F");
  fTree->Branch("l_phi",&l_phi,"l_phi/F");
  fTree->Branch("l_e",&l_e,"l_e/F");
  fTree->Branch("ungroomed_jet_pt",&ungroomed_jet_pt,"ungroomed_jet_pt/F");
  fTree->Branch("ungroomed_jet_eta",&ungroomed_jet_eta,"ungroomed_jet_eta/F");
  fTree->Branch("ungroomed_jet_phi",&ungroomed_jet_phi,"ungroomed_jet_phi/F");
  fTree->Branch("ungroomed_jet_e",&ungroomed_jet_e,"ungroomed_jet_e/F");
  fTree->Branch("jet_mass_pr",&jet_mass_pr,"jet_mass_pr");
  fTree->Branch("jet_mass_tr",&jet_mass_tr,"jet_mass_tr");
  fTree->Branch("jet_mass_fi",&jet_mass_fi,"jet_mass_fi");
  fTree->Branch("jet_tau2tau1",&jet_tau2tau1,"jet_tau2tau1");
  fTree->Branch("W_pt_gen",&W_pt_gen,"W_pt_gen");
  fTree->Branch("W_pz_gen",&W_pz_gen,"W_pz_gen");
  fTree->Branch("genGravMass",&genGravMass,"genGravMass");
  fTree->Branch("nu_pz_gen",&nu_pz_gen,"nu_pz_gen");
  fTree->Branch("deltaR_lak8jet",&deltaR_lak8jet,"deltaR_lak8jet/F");
  fTree->Branch("deltaphi_METak8jet",&deltaphi_METak8jet,"deltaphi_METak8jet/F");
  fTree->Branch("deltaphi_Vak8jet",&deltaphi_Vak8jet,"deltaphi_Vak8jet/F");
  fTree->Branch("v_pt",&v_pt,"v_pt/F");
  fTree->Branch("v_eta",&v_eta,"v_eta/F");
  fTree->Branch("v_phi",&v_phi,"v_phi/F");
  fTree->Branch("v_mt",&v_mt,"v_mt/F");
  fTree->Branch("mass_lvj_type0",&mass_lvj_type0,"mass_lvj_type0/F");
  fTree->Branch("mass_lvj_type2",&mass_lvj_type2,"mass_lvj_type2/F");
  fTree->Branch("nBTagJet_loose",&nBTagJet_loose,"nBTagJet_loose/I");
  fTree->Branch("nBTagJet_medium",&nBTagJet_medium,"nBTagJet_medium/I");
  fTree->Branch("nBTagJet_tight",&nBTagJet_tight,"nBTagJet_tight/I");
  fTree->Branch("vbf_maxpt_j1_pt",&vbf_maxpt_j1_pt,"vbf_maxpt_j1_pt/F");
  fTree->Branch("vbf_maxpt_j1_eta",&vbf_maxpt_j1_eta,"vbf_maxpt_j1_eta/F");
  fTree->Branch("vbf_maxpt_j1_phi",&vbf_maxpt_j1_phi,"vbf_maxpt_j1_phi/F");
  fTree->Branch("vbf_maxpt_j1_e",&vbf_maxpt_j1_e,"vbf_maxpt_j1_e/F");
  fTree->Branch("vbf_maxpt_j1_bDiscriminatorCSV",&vbf_maxpt_j1_bDiscriminatorCSV,"vbf_maxpt_j1_bDiscriminatorCSV/F");
  fTree->Branch("vbf_maxpt_j2_pt",&vbf_maxpt_j2_pt,"vbf_maxpt_j2_pt/F");
  fTree->Branch("vbf_maxpt_j2_eta",&vbf_maxpt_j2_eta,"vbf_maxpt_j2_eta/F");
  fTree->Branch("vbf_maxpt_j2_phi",&vbf_maxpt_j2_phi,"vbf_maxpt_j2_phi/F");
  fTree->Branch("vbf_maxpt_j2_e",&vbf_maxpt_j2_e,"vbf_maxpt_j2_e/F");
  fTree->Branch("vbf_maxpt_j2_bDiscriminatorCSV",&vbf_maxpt_j2_bDiscriminatorCSV,"vbf_maxpt_j2_bDiscriminatorCSV/F");
  fTree->Branch("vbf_maxpt_jj_pt",&vbf_maxpt_jj_pt,"vbf_maxpt_jj_pt/F");
  fTree->Branch("vbf_maxpt_jj_eta",&vbf_maxpt_jj_eta,"vbf_maxpt_jj_eta/F");
  fTree->Branch("vbf_maxpt_jj_phi",&vbf_maxpt_jj_phi,"vbf_maxpt_jj_phi/F");
  fTree->Branch("vbf_maxpt_jj_m",&vbf_maxpt_jj_m,"vbf_maxpt_jj_m/F");
}

