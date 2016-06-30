#include "../interface/setOutputTree.h"

setOutputTree::setOutputTree(TTree* outTree){
  fTree = outTree;
  fTree->SetAutoSave(500000000);
  setBranches();
}

setOutputTree::~setOutputTree(){
  delete fTree;
}

void setOutputTree::initializeVariables()
{
  run=-999;
  event=-999;
  nEvents=0;
  nNegEvents=0;
  lumi=-999;
  nPV=-999;
  issignal=0;
  issignal_PuppiAK8=0;
  issignal_AK4jetjet=0;
  issignal_PuppiAK4jetjet=0;
  wSampleWeight=-999;
  genWeight=1;
  top1_NNLO_Weight=1.;
  top2_NNLO_Weight=1.;
  trig_eff_Weight=1.;
  id_eff_Weight=1.;
  gen_top1_pt=-999;
  gen_top2_pt=-999;
  totalEventWeight=-999;
  totalEventWeight_2=-999;
  totalEventWeight_3=-999;
  eff_and_pu_Weight=-999;
  eff_and_pu_Weight_2=-999;
  eff_and_pu_Weight_3=-999;
  pfMET=-999;
  pfMET_jes_up=-999;
  pfMET_jes_dn=-999;
  pfMET_jer=-999;
  pfMET_jer_up=-999;
  pfMET_jer_dn=-999;
  pfMET_Phi=-999;
  pfMETpuppi=-999;
  pfMETpuppi_jes_up=-999;
  pfMETpuppi_jes_dn=-999;
  pfMETpuppi_jer=-999;
  pfMETpuppi_jer_up=-999;
  pfMETpuppi_jer_dn=-999;
  pfMETpuppi_Phi=-999;
  nu_pz_type0=-999;
  nu_pz_type2=-999;
  nu_pz_run2=-999;
  nu_pz_run2_oth=-999;
  nu_pz_run2_type=-999;
  nu_pz_isre=1;
  l_pt=-999;
  l_eta=-999;
  l_phi=-999;
  l_e=-999;
  W_pt_gen=-999;
  W_pz_gen=-999;
  W_rap_gen=-999;
  genGravMass=-999;
  nu_pz_gen=-999;
  nu_pt_gen=-999;
  nu_phi_gen=-999;
  nu_eta_gen=-999;
  hadW_pt_gen=-999;
  hadW_eta_gen=-999;
  hadW_phi_gen=-999;
  hadW_e_gen=-999;
  hadW_m_gen=-999;
  lepW_pt_gen=-999;
  lepW_eta_gen=-999;
  lepW_phi_gen=-999;
  lepW_e_gen=-999;
  lepW_m_gen=-999;
  AK8_pt_gen=-999;
  AK8_eta_gen=-999;
  AK8_phi_gen=-999;
  AK8_e_gen=-999;
  AK8_mass_gen=-999;
  AK8_pruned_mass_gen=-999;
  AK8_softdrop_mass_gen=-999;
  AK8_softdrop_pt_gen=-999;
  deltaR_lak8jet=-999;
  deltaphi_METak8jet=-999;
  deltaphi_Vak8jet=-999;
  deltaR_lPuppiak8jet=-999;
  deltaphi_METPuppiak8jet=-999;
  deltaphi_VPuppiak8jet=-999;
  deltaR_lak4jetjet=-999;
  deltaphi_METak4jetjet=-999;
  deltaphi_Vak4jetjet=-999;
  deltaR_lPuppiak4jetjet=-999;
  deltaphi_METPuppiak4jetjet=-999;
  deltaphi_VPuppiak4jetjet=-999;
  v_pt=-999;
  v_eta=-999;
  v_phi=-999;
  v_mt=-999;
  v_puppi_pt=-999;
  v_puppi_eta=-999;
  v_puppi_phi=-999;
  v_puppi_mt=-999;
  mass_lvj_type0=-999;
  mass_lvj_type0_met_jes_up=-999;
  mass_lvj_type0_met_jes_dn=-999;
  mass_lvj_type0_met_jer=-999;
  mass_lvj_type0_met_jer_up=-999;
  mass_lvj_type0_met_jer_dn=-999;
  mass_lvj_type0_met_PuppiAK8_jes_up=-999;
  mass_lvj_type0_met_PuppiAK8_jes_dn=-999;
  mass_lvj_type2=-999;
  mass_lvj_run2=-999;
  mass_lvj_type0_PuppiAK8=-999;
  mass_lvj_type2_PuppiAK8=-999;
  mass_lvj_run2_PuppiAK8=-999;
  mass_lvjj_type0_AK4=-999;
  mass_lvjj_type0_met_jes_up_AK4=-999;
  mass_lvjj_type0_met_jes_dn_AK4=-999;
  mass_lvjj_type2_AK4=-999;
  mass_lvjj_run2_AK4=-999;
  mass_lvjj_type0_PuppiAK4=-999;
  mass_lvjj_type0_met_jes_up_PuppiAK4=-999;
  mass_lvjj_type0_met_jes_dn_PuppiAK4=-999;
  mass_lvjj_type2_PuppiAK4=-999;
  mass_lvjj_run2_PuppiAK4=-999;
  njets=-999;
  njetsPuppi=-999;
  njets_unmerged=-999;
  njetsPuppi_unmerged=-999;
  nBTagJet_loose=-999;
  nBTagJet_medium=-999;
  nBTagJet_tight=-999;
  nBTagJetPuppi_loose=-999;
  nBTagJetPuppi_medium=-999;
  nBTagJetPuppi_tight=-999;
  nBTagJet_loose_unmerged=-999;
  nBTagJet_medium_unmerged=-999;
  nBTagJet_tight_unmerged=-999;
  nBTagJetPuppi_loose_unmerged=-999;
  nBTagJetPuppi_medium_unmerged=-999;
  nBTagJetPuppi_tight_unmerged=-999;
  ungroomed_jet_pt=-999;
  ungroomed_jet_pt_jes_up=-999;
  ungroomed_jet_pt_jes_dn=-999;
  ungroomed_jet_pt_jer=-999;
  ungroomed_jet_pt_jer_up=-999;
  ungroomed_jet_pt_jer_dn=-999;
  ungroomed_jet_eta=-999;
  ungroomed_jet_phi=-999;
  ungroomed_jet_e=-999;
  jet_mass_pr=-999;
  jet_mass_pr_jes_up=-999;
  jet_mass_pr_jes_dn=-999;
  jet_mass_pr_jer=-999;
  jet_mass_pr_jer_up=-999;
  jet_mass_pr_jer_dn=-999;
  jet_mass_so=-999;
  jet_pt_so=-999;
  jet_mass_tr=-999;
  jet_mass_fi=-999;
  jet_mass=-999;
  jet_tau2tau1=-999;
  AK4_jetjet_pt=-999;
  AK4_jetjet_mass=-999;
  AK4_jetjet_deltaeta=-999;
  AK4_jetjet_deltaphi=-999;
  AK4_jetjet_deltar=-999;
  PuppiAK4_jetjet_pt=-999;
  PuppiAK4_jetjet_mass=-999;
  PuppiAK4_jetjet_deltaeta=-999;
  PuppiAK4_jetjet_deltaphi=-999;
  PuppiAK4_jetjet_deltar=-999;
  ungroomed_PuppiAK8_jet_pt=-999;
  ungroomed_PuppiAK8_jet_pt_jes_up=-999;
  ungroomed_PuppiAK8_jet_pt_jes_dn=-999;
  ungroomed_PuppiAK8_jet_pt_jer=-999;
  ungroomed_PuppiAK8_jet_pt_jer_up=-999;
  ungroomed_PuppiAK8_jet_pt_jer_dn=-999;
  ungroomed_PuppiAK8_jet_eta=-999;
  ungroomed_PuppiAK8_jet_phi=-999;
  ungroomed_PuppiAK8_jet_e=-999;
  PuppiAK8_jet_mass_pr=-999;
  PuppiAK8_jet_mass_pr_jes_up=-999;
  PuppiAK8_jet_mass_pr_jes_dn=-999;
  PuppiAK8_jet_mass_pr_jer=-999;
  PuppiAK8_jet_mass_pr_jer_up=-999;
  PuppiAK8_jet_mass_pr_jer_dn=-999;
  PuppiAK8_jet_mass_so=-999;
  PuppiAK8_jet_pt_so=-999;
  PuppiAK8_jet_mass_tr=-999;
  PuppiAK8_jet_mass_fi=-999;
  PuppiAK8_jet_mass=-999;
  PuppiAK8_jet_tau2tau1=-999;
  AK4_jet1_pt=-999;
  AK4_jet1_pt_jes_up=-999;
  AK4_jet1_pt_jes_dn=-999;
  AK4_jet1_pt_jer=-999;
  AK4_jet1_pt_jer_up=-999;
  AK4_jet1_pt_jer_dn=-999;
  AK4_jet1_eta=-999;
  AK4_jet1_phi=-999;
  AK4_jet1_e=-999;
  AK4_jet2_pt=-999;
  AK4_jet2_pt_jes_up=-999;
  AK4_jet2_pt_jes_dn=-999;
  AK4_jet2_pt_jer=-999;
  AK4_jet2_pt_jer_up=-999;
  AK4_jet2_pt_jer_dn=-999;
  AK4_jet2_eta=-999;
  AK4_jet2_phi=-999;
  AK4_jet2_e=-999;
  PuppiAK4_jet1_pt=-999;
  PuppiAK4_jet1_pt_jes_up=-999;
  PuppiAK4_jet1_pt_jes_dn=-999;
  PuppiAK4_jet1_pt_jer=-999;
  PuppiAK4_jet1_pt_jer_up=-999;
  PuppiAK4_jet1_pt_jer_dn=-999;
  PuppiAK4_jet1_eta=-999;
  PuppiAK4_jet1_phi=-999;
  PuppiAK4_jet1_e=-999;
  PuppiAK4_jet2_pt=-999;
  PuppiAK4_jet2_pt_jes_up=-999;
  PuppiAK4_jet2_pt_jes_dn=-999;
  PuppiAK4_jet2_pt_jer=-999;
  PuppiAK4_jet2_pt_jer_up=-999;
  PuppiAK4_jet2_pt_jer_dn=-999;
  PuppiAK4_jet2_eta=-999;
  PuppiAK4_jet2_phi=-999;
  PuppiAK4_jet2_e=-999;
  ttb_ungroomed_jet_pt=-999;
  ttb_ungroomed_jet_eta=-999;
  ttb_ungroomed_jet_phi=-999;
  ttb_ungroomed_jet_e=-999;
  ttb_jet_mass_pr=-999;
  ttb_jet_mass_so=-999;
  ttb_jet_pt_so=-999;
  ttb_jet_mass_tr=-999;
  ttb_jet_mass_fi=-999;
  ttb_jet_tau2tau1=-999;
  ttb_deltaeta_lak8jet=-999;
  mass_leptonic_closerjet=-999;
  mass_ungroomedjet_closerjet=-999;
  AK8_closerjet_pt=-999;
  AK8_closerjet_eta=-999;
  AK8_closerjet_phi=-999;
  AK8_closerjet_e=-999;
  vbf_maxpt_j1_pt=-999; 
  vbf_maxpt_j1_pt_jes_up=-999;
  vbf_maxpt_j1_pt_jes_dn=-999;
  vbf_maxpt_j1_pt_jer=-999;
  vbf_maxpt_j1_pt_jer_up=-999;
  vbf_maxpt_j1_pt_jer_dn=-999;
  vbf_maxpt_j1_eta=-999;
  vbf_maxpt_j1_eta_jes_up=-999;
  vbf_maxpt_j1_eta_jes_dn=-999;
  vbf_maxpt_j1_eta_jer=-999;
  vbf_maxpt_j1_eta_jer_up=-999;
  vbf_maxpt_j1_eta_jer_dn=-999;
  vbf_maxpt_j1_phi=-999;
  vbf_maxpt_j1_e=-999;
  vbf_maxpt_j1_bDiscriminatorCSV=-999;
  vbf_maxpt_j2_pt=-999;
  vbf_maxpt_j2_pt_jes_up=-999;
  vbf_maxpt_j2_pt_jes_dn=-999;
  vbf_maxpt_j2_pt_jer=-999;
  vbf_maxpt_j2_pt_jer_up=-999;
  vbf_maxpt_j2_pt_jer_dn=-999;
  vbf_maxpt_j2_eta=-999;
  vbf_maxpt_j2_eta_jes_up=-999;
  vbf_maxpt_j2_eta_jes_dn=-999;
  vbf_maxpt_j2_eta_jer=-999;
  vbf_maxpt_j2_eta_jer_up=-999;
  vbf_maxpt_j2_eta_jer_dn=-999;
  vbf_maxpt_j2_phi=-999;
  vbf_maxpt_j2_e=-999;
  vbf_maxpt_j2_bDiscriminatorCSV=-999;
  vbf_maxpt_jj_pt=-999;
  vbf_maxpt_jj_eta=-999;
  vbf_maxpt_jj_phi=-999;
  vbf_maxpt_jj_m=-999;
  jet2_pt=0;
  jet2_eta=0;
  jet2_phi=0;
  jet2_e=0;
  jet2_btag=0;
  jet3_pt=0;
  jet3_eta=0;
  jet3_phi=0;
  jet3_e=0;
  jet3_btag=0;
  deltaR_AK8_closestBtagJet=0;
  deltaR_AK8_closestBtagJet_loose=0;
  vbf_maxpt_deltaR =-999;
  AK4_1_pt_gen=-999;
  AK4_1_eta_gen=-999;
  AK4_1_phi_gen=-999;
  AK4_1_e_gen=-999;
  AK4_1_mass_gen=-999;
  AK4_2_pt_gen=-999;
  AK4_2_eta_gen=-999;
  AK4_2_phi_gen=-999;
  AK4_2_e_gen=-999;
  AK4_2_mass_gen=-999;
  AK4_BIG_gen_mass=-999;
  deltaR_AK4=-999;

}

void setOutputTree::setBranches()
{
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("nEvents",&nEvents,"nEvents/I");
  fTree->Branch("nNegEvents",&nNegEvents,"nNegEvents/I");
  fTree->Branch("lumi",&lumi,"lumi/I");
  fTree->Branch("nPV",&nPV,"nPV/I");
  fTree->Branch("issignal",&issignal,"issignal/I");
  fTree->Branch("issignal_PuppiAK8",&issignal_PuppiAK8,"issignal_PuppiAK8/I");
  fTree->Branch("issignal_AK4jetjet",&issignal_AK4jetjet,"issignal_AK4jetjet/I");
  fTree->Branch("issignal_PuppiAK4jetjet",&issignal_PuppiAK4jetjet,"issignal_PuppiAK4jetjet/I");
  fTree->Branch("wSampleWeight",&wSampleWeight,"wSampleWeight/F");
  fTree->Branch("genWeight",&genWeight,"genWeight/F");
  fTree->Branch("top1_NNLO_Weight",&top1_NNLO_Weight,"top1_NNLO_Weight/F");
  fTree->Branch("top2_NNLO_Weight",&top2_NNLO_Weight,"top2_NNLO_Weight/F");
  fTree->Branch("gen_top1_pt",&gen_top1_pt,"gen_top1_pt/F");
  fTree->Branch("gen_top2_pt",&gen_top2_pt,"gen_top2_pt/F");
  fTree->Branch("trig_eff_Weight",&trig_eff_Weight,"trig_eff_Weight/F");
  fTree->Branch("id_eff_Weight",&id_eff_Weight,"id_eff_Weight/F");
  fTree->Branch("totalEventWeight",&totalEventWeight,"totalEventWeight/F");
  fTree->Branch("eff_and_pu_Weight",&eff_and_pu_Weight,"eff_and_pu_Weight/F");
  fTree->Branch("totalEventWeight_2",&totalEventWeight_2,"totalEventWeight_2/F");
  fTree->Branch("eff_and_pu_Weight_2",&eff_and_pu_Weight_2,"eff_and_pu_Weight_2/F");
  fTree->Branch("totalEventWeight_3",&totalEventWeight_3,"totalEventWeight_3/F");
  fTree->Branch("eff_and_pu_Weight_3",&eff_and_pu_Weight_3,"eff_and_pu_Weight_3/F");
  fTree->Branch("pfMET",&pfMET,"pfMET/F");
  fTree->Branch("pfMET_jes_up",&pfMET_jes_up,"pfMET_jes_up/F");
  fTree->Branch("pfMET_jes_dn",&pfMET_jes_dn,"pfMET_jes_dn/F");
  fTree->Branch("pfMET_jer",&pfMET_jer,"pfMET_jer/F");
  fTree->Branch("pfMET_jer_up",&pfMET_jer_up,"pfMET_jer_up/F");
  fTree->Branch("pfMET_jer_dn",&pfMET_jer_dn,"pfMET_jer_dn/F");
  fTree->Branch("pfMET_Phi",&pfMET_Phi,"pfMET_Phi/F");
  fTree->Branch("pfMETpuppi",&pfMETpuppi,"pfMETpuppi/F");
  fTree->Branch("pfMETpuppi_jes_up",&pfMETpuppi_jes_up,"pfMETpuppi_jes_up/F");
  fTree->Branch("pfMETpuppi_jes_dn",&pfMETpuppi_jes_dn,"pfMETpuppi_jes_dn/F");
  fTree->Branch("pfMETpuppi_jer",&pfMETpuppi_jer,"pfMETpuppi_jer/F");
  fTree->Branch("pfMETpuppi_jer_up",&pfMETpuppi_jer_up,"pfMETpuppi_jer_up/F");
  fTree->Branch("pfMETpuppi_jer_dn",&pfMETpuppi_jer_dn,"pfMETpuppi_jer_dn/F");
  fTree->Branch("pfMETpuppi_Phi",&pfMETpuppi_Phi,"pfMETpuppi_Phi/F");
  fTree->Branch("nu_pz_type0",&nu_pz_type0,"nu_pz_type0/F");
  fTree->Branch("nu_pz_type2",&nu_pz_type2,"nu_pz_type2/F");
  fTree->Branch("nu_pz_run2",&nu_pz_run2,"nu_pz_run2/F");
  fTree->Branch("nu_pz_run2_oth",&nu_pz_run2_oth,"nu_pz_run2_oth/F");
  fTree->Branch("nu_pz_run2_type",&nu_pz_run2_type,"nu_pz_run2_type/I");
  fTree->Branch("nu_pz_isre",&nu_pz_isre,"nu_pz_isre/I");
  fTree->Branch("l_pt",&l_pt,"l_pt/F");
  fTree->Branch("l_eta",&l_eta,"l_eta/F");
  fTree->Branch("l_phi",&l_phi,"l_phi/F");
  fTree->Branch("l_e",&l_e,"l_e/F");
  fTree->Branch("ungroomed_jet_pt",&ungroomed_jet_pt,"ungroomed_jet_pt/F");
  fTree->Branch("ungroomed_jet_pt_jes_up",&ungroomed_jet_pt_jes_up,"ungroomed_jet_pt_jes_up/F");
  fTree->Branch("ungroomed_jet_pt_jes_dn",&ungroomed_jet_pt_jes_dn,"ungroomed_jet_pt_jes_dn/F");
  fTree->Branch("ungroomed_jet_pt_jer",&ungroomed_jet_pt_jer,"ungroomed_jet_pt_jer/F");
  fTree->Branch("ungroomed_jet_pt_jer_up",&ungroomed_jet_pt_jer_up,"ungroomed_jet_pt_jer_up/F");
  fTree->Branch("ungroomed_jet_pt_jer_dn",&ungroomed_jet_pt_jer_dn,"ungroomed_jet_pt_jer_dn/F");
  fTree->Branch("ungroomed_jet_eta",&ungroomed_jet_eta,"ungroomed_jet_eta/F");
  fTree->Branch("ungroomed_jet_phi",&ungroomed_jet_phi,"ungroomed_jet_phi/F");
  fTree->Branch("ungroomed_jet_e",&ungroomed_jet_e,"ungroomed_jet_e/F");
  fTree->Branch("jet_mass_pr",&jet_mass_pr,"jet_mass_pr");
  fTree->Branch("jet_mass_pr_jes_up",&jet_mass_pr_jes_up,"jet_mass_pr_jes_up/F");
  fTree->Branch("jet_mass_pr_jes_dn",&jet_mass_pr_jes_dn,"jet_mass_pr_jes_dn/F");
  fTree->Branch("jet_mass_pr_jer",&jet_mass_pr_jer,"jet_mass_pr_jer/F");
  fTree->Branch("jet_mass_pr_jer_up",&jet_mass_pr_jer_up,"jet_mass_pr_jer_up/F");
  fTree->Branch("jet_mass_pr_jer_dn",&jet_mass_pr_jer_dn,"jet_mass_pr_jer_dn/F");
  fTree->Branch("jet_mass_so",&jet_mass_so,"jet_mass_so");
  fTree->Branch("jet_pt_so",&jet_pt_so,"jet_pt_so");
  fTree->Branch("jet_mass_tr",&jet_mass_tr,"jet_mass_tr");
  fTree->Branch("jet_mass_fi",&jet_mass_fi,"jet_mass_fi");
  fTree->Branch("jet_mass",&jet_mass,"jet_mass");
  fTree->Branch("jet_tau2tau1",&jet_tau2tau1,"jet_tau2tau1");
  fTree->Branch("AK4_jetjet_pt",&AK4_jetjet_pt,"AK4_jetjet_pt");
  fTree->Branch("AK4_jetjet_mass",&AK4_jetjet_mass,"AK4_jetjet_mass");
  fTree->Branch("AK4_jetjet_deltaeta",&AK4_jetjet_deltaeta,"AK4_jetjet_deltaeta");
  fTree->Branch("AK4_jetjet_deltaphi",&AK4_jetjet_deltaphi,"AK4_jetjet_deltaphi");
  fTree->Branch("AK4_jetjet_deltar",&AK4_jetjet_deltar,"AK4_jetjet_deltar");
  fTree->Branch("PuppiAK4_jetjet_pt",&PuppiAK4_jetjet_pt,"PuppiAK4_jetjet_pt");
  fTree->Branch("PuppiAK4_jetjet_mass",&PuppiAK4_jetjet_mass,"PuppiAK4_jetjet_mass");
  fTree->Branch("PuppiAK4_jetjet_deltaeta",&PuppiAK4_jetjet_deltaeta,"PuppiAK4_jetjet_deltaeta");
  fTree->Branch("PuppiAK4_jetjet_deltaphi",&PuppiAK4_jetjet_deltaphi,"PuppiAK4_jetjet_deltaphi");
  fTree->Branch("PuppiAK4_jetjet_deltar",&PuppiAK4_jetjet_deltar,"PuppiAK4_jetjet_deltar");
  fTree->Branch("ungroomed_PuppiAK8_jet_pt",&ungroomed_PuppiAK8_jet_pt,"ungroomed_PuppiAK8_jet_pt/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_pt_jes_up",&ungroomed_PuppiAK8_jet_pt_jes_up,"ungroomed_PuppiAK8_jet_pt_jes_up/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_pt_jes_dn",&ungroomed_PuppiAK8_jet_pt_jes_dn,"ungroomed_PuppiAK8_jet_pt_jes_dn/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_pt_jer",&ungroomed_PuppiAK8_jet_pt_jer,"ungroomed_PuppiAK8_jet_pt_jer/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_pt_jer_up",&ungroomed_PuppiAK8_jet_pt_jer_up,"ungroomed_PuppiAK8_jet_pt_jer_up/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_pt_jer_dn",&ungroomed_PuppiAK8_jet_pt_jer_dn,"ungroomed_PuppiAK8_jet_pt_jer_dn/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_eta",&ungroomed_PuppiAK8_jet_eta,"ungroomed_PuppiAK8_jet_eta/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_phi",&ungroomed_PuppiAK8_jet_phi,"ungroomed_PuppiAK8_jet_phi/F");
  fTree->Branch("ungroomed_PuppiAK8_jet_e",&ungroomed_PuppiAK8_jet_e,"ungroomed_PuppiAK8_jet_e/F");
  fTree->Branch("PuppiAK8_jet_mass_pr",&PuppiAK8_jet_mass_pr,"PuppiAK8_jet_mass_pr");
  fTree->Branch("PuppiAK8_jet_mass_pr_jes_up",&PuppiAK8_jet_mass_pr_jes_up,"PuppiAK8_jet_mass_pr_jes_up/F");
  fTree->Branch("PuppiAK8_jet_mass_pr_jes_dn",&PuppiAK8_jet_mass_pr_jes_dn,"PuppiAK8_jet_mass_pr_jes_dn/F");
  fTree->Branch("PuppiAK8_jet_mass_pr_jer",&PuppiAK8_jet_mass_pr_jer,"PuppiAK8_jet_mass_pr_jer/F");
  fTree->Branch("PuppiAK8_jet_mass_pr_jer_up",&PuppiAK8_jet_mass_pr_jer_up,"PuppiAK8_jet_mass_pr_jer_up/F");
  fTree->Branch("PuppiAK8_jet_mass_pr_jer_dn",&PuppiAK8_jet_mass_pr_jer_dn,"PuppiAK8_jet_mass_pr_jer_dn/F");
  fTree->Branch("PuppiAK8_jet_mass_so",&PuppiAK8_jet_mass_so,"PuppiAK8_jet_mass_so");
  fTree->Branch("PuppiAK8_jet_pt_so",&PuppiAK8_jet_pt_so,"PuppiAK8_jet_pt_so");
  fTree->Branch("PuppiAK8_jet_mass_tr",&PuppiAK8_jet_mass_tr,"PuppiAK8_jet_mass_tr");
  fTree->Branch("PuppiAK8_jet_mass_fi",&PuppiAK8_jet_mass_fi,"PuppiAK8_jet_mass_fi");
  fTree->Branch("PuppiAK8_jet_mass",&PuppiAK8_jet_mass,"PuppiAK8_jet_mass");
  fTree->Branch("PuppiAK8_jet_tau2tau1",&PuppiAK8_jet_tau2tau1,"PuppiAK8_jet_tau2tau1");
  fTree->Branch("AK4_jet1_pt",&AK4_jet1_pt,"AK4_jet1_pt/F");
  fTree->Branch("AK4_jet1_pt_jes_up",&AK4_jet1_pt_jes_up,"AK4_jet1_pt_jes_up/F");
  fTree->Branch("AK4_jet1_pt_jes_dn",&AK4_jet1_pt_jes_dn,"AK4_jet1_pt_jes_dn/F");
  fTree->Branch("AK4_jet1_pt_jer",&AK4_jet1_pt_jer,"AK4_jet1_pt_jer/F");
  fTree->Branch("AK4_jet1_pt_jer_up",&AK4_jet1_pt_jer_up,"AK4_jet1_pt_jer_up/F");
  fTree->Branch("AK4_jet1_pt_jer_dn",&AK4_jet1_pt_jer_dn,"AK4_jet1_pt_jer_dn/F");
  fTree->Branch("AK4_jet1_eta",&AK4_jet1_eta,"AK4_jet1_eta/F");
  fTree->Branch("AK4_jet1_phi",&AK4_jet1_phi,"AK4_jet1_phi/F");
  fTree->Branch("AK4_jet1_e",&AK4_jet1_e,"AK4_jet1_e/F");
  fTree->Branch("AK4_jet2_pt",&AK4_jet2_pt,"AK4_jet2_pt/F");
  fTree->Branch("AK4_jet2_pt_jes_up",&AK4_jet2_pt_jes_up,"AK4_jet2_pt_jes_up/F");
  fTree->Branch("AK4_jet2_pt_jes_dn",&AK4_jet2_pt_jes_dn,"AK4_jet2_pt_jes_dn/F");
  fTree->Branch("AK4_jet2_pt_jer",&AK4_jet2_pt_jer,"AK4_jet2_pt_jer/F");
  fTree->Branch("AK4_jet2_pt_jer_up",&AK4_jet2_pt_jer_up,"AK4_jet2_pt_jer_up/F");
  fTree->Branch("AK4_jet2_pt_jer_dn",&AK4_jet2_pt_jer_dn,"AK4_jet2_pt_jer_dn/F");
  fTree->Branch("AK4_jet2_eta",&AK4_jet2_eta,"AK4_jet2_eta/F");
  fTree->Branch("AK4_jet2_phi",&AK4_jet2_phi,"AK4_jet2_phi/F");
  fTree->Branch("AK4_jet2_e",&AK4_jet2_e,"AK4_jet2_e/F");
  fTree->Branch("PuppiAK4_jet1_pt",&PuppiAK4_jet1_pt,"PuppiAK4_jet1_pt/F");
  fTree->Branch("PuppiAK4_jet1_pt_jes_up",&PuppiAK4_jet1_pt_jes_up,"PuppiAK4_jet1_pt_jes_up/F");
  fTree->Branch("PuppiAK4_jet1_pt_jes_dn",&PuppiAK4_jet1_pt_jes_dn,"PuppiAK4_jet1_pt_jes_dn/F");
  fTree->Branch("PuppiAK4_jet1_pt_jer",&PuppiAK4_jet1_pt_jer,"PuppiAK4_jet1_pt_jer/F");
  fTree->Branch("PuppiAK4_jet1_pt_jer_up",&PuppiAK4_jet1_pt_jer_up,"PuppiAK4_jet1_pt_jer_up/F");
  fTree->Branch("PuppiAK4_jet1_pt_jer_dn",&PuppiAK4_jet1_pt_jer_dn,"PuppiAK4_jet1_pt_jer_dn/F");
  fTree->Branch("PuppiAK4_jet1_eta",&PuppiAK4_jet1_eta,"PuppiAK4_jet1_eta/F");
  fTree->Branch("PuppiAK4_jet1_phi",&PuppiAK4_jet1_phi,"PuppiAK4_jet1_phi/F");
  fTree->Branch("PuppiAK4_jet1_e",&PuppiAK4_jet1_e,"PuppiAK4_jet1_e/F");
  fTree->Branch("PuppiAK4_jet2_pt",&PuppiAK4_jet2_pt,"PuppiAK4_jet2_pt/F");
  fTree->Branch("PuppiAK4_jet2_pt_jes_up",&PuppiAK4_jet2_pt_jes_up,"PuppiAK4_jet2_pt_jes_up/F");
  fTree->Branch("PuppiAK4_jet2_pt_jes_dn",&PuppiAK4_jet2_pt_jes_dn,"PuppiAK4_jet2_pt_jes_dn/F");
  fTree->Branch("PuppiAK4_jet2_pt_jer",&PuppiAK4_jet2_pt_jer,"PuppiAK4_jet2_pt_jer/F");
  fTree->Branch("PuppiAK4_jet2_pt_jer_up",&PuppiAK4_jet2_pt_jer_up,"PuppiAK4_jet2_pt_jer_up/F");
  fTree->Branch("PuppiAK4_jet2_pt_jer_dn",&PuppiAK4_jet2_pt_jer_dn,"PuppiAK4_jet2_pt_jer_dn/F");
  fTree->Branch("PuppiAK4_jet2_eta",&PuppiAK4_jet2_eta,"PuppiAK4_jet2_eta/F");
  fTree->Branch("PuppiAK4_jet2_phi",&PuppiAK4_jet2_phi,"PuppiAK4_jet2_phi/F");
  fTree->Branch("PuppiAK4_jet2_e",&PuppiAK4_jet2_e,"PuppiAK4_jet2_e/F");
  fTree->Branch("ttb_ungroomed_jet_pt",&ttb_ungroomed_jet_pt,"ttb_ungroomed_jet_pt/F");
  fTree->Branch("ttb_ungroomed_jet_eta",&ttb_ungroomed_jet_eta,"ttb_ungroomed_jet_eta/F");
  fTree->Branch("ttb_ungroomed_jet_phi",&ttb_ungroomed_jet_phi,"ttb_ungroomed_jet_phi/F");
  fTree->Branch("ttb_ungroomed_jet_e",&ttb_ungroomed_jet_e,"ttb_ungroomed_jet_e/F");
  fTree->Branch("ttb_jet_mass_pr",&ttb_jet_mass_pr,"ttb_jet_mass_pr");
  fTree->Branch("ttb_jet_mass_so",&ttb_jet_mass_so,"ttb_jet_mass_so");
  fTree->Branch("ttb_jet_pt_so",&ttb_jet_pt_so,"ttb_jet_pt_so");
  fTree->Branch("ttb_jet_mass_tr",&ttb_jet_mass_tr,"ttb_jet_mass_tr");
  fTree->Branch("ttb_jet_mass_fi",&ttb_jet_mass_fi,"ttb_jet_mass_fi");
  fTree->Branch("ttb_jet_tau2tau1",&ttb_jet_tau2tau1,"ttb_jet_tau2tau1");
  fTree->Branch("ttb_deltaeta_lak8jet",&ttb_deltaeta_lak8jet,"ttb_deltaeta_lak8jet/F");
  fTree->Branch("W_pt_gen",&W_pt_gen,"W_pt_gen");
  fTree->Branch("W_pz_gen",&W_pz_gen,"W_pz_gen");
  fTree->Branch("W_rap_gen",&W_rap_gen,"W_rap_gen");
  fTree->Branch("genGravMass",&genGravMass,"genGravMass");
  fTree->Branch("nu_pz_gen",&nu_pz_gen,"nu_pz_gen");
  fTree->Branch("nu_pt_gen",&nu_pt_gen,"nu_pt_gen");
  fTree->Branch("nu_phi_gen",&nu_phi_gen,"nu_phi_gen");
  fTree->Branch("nu_eta_gen",&nu_eta_gen,"nu_eta_gen");
  fTree->Branch("hadW_pt_gen",&hadW_pt_gen,"hadW_pt_gen");
  fTree->Branch("hadW_eta_gen",&hadW_eta_gen,"hadW_eta_gen");
  fTree->Branch("hadW_phi_gen",&hadW_phi_gen,"hadW_phi_gen");
  fTree->Branch("hadW_e_gen",&hadW_e_gen,"hadW_e_gen");
  fTree->Branch("hadW_m_gen",&hadW_m_gen,"hadW_m_gen");
  fTree->Branch("lepW_pt_gen",&lepW_pt_gen,"lepW_pt_gen");
  fTree->Branch("lepW_eta_gen",&lepW_eta_gen,"lepW_eta_gen");
  fTree->Branch("lepW_phi_gen",&lepW_phi_gen,"lepW_phi_gen");
  fTree->Branch("lepW_e_gen",&lepW_e_gen,"lepW_e_gen");
  fTree->Branch("lepW_m_gen",&lepW_m_gen,"lepW_m_gen");
  fTree->Branch("AK8_pt_gen",&AK8_pt_gen,"AK8_pt_gen");
  fTree->Branch("AK8_eta_gen",&AK8_eta_gen,"AK8_eta_gen");
  fTree->Branch("AK8_phi_gen",&AK8_phi_gen,"AK8_phi_gen");
  fTree->Branch("AK8_e_gen",&AK8_e_gen,"AK8_e_gen");
  fTree->Branch("AK8_mass_gen",&AK8_mass_gen,"AK8_mass_gen");
  fTree->Branch("AK8_pruned_mass_gen",&AK8_pruned_mass_gen,"AK8_pruned_mass_gen");
  fTree->Branch("AK8_softdrop_mass_gen",&AK8_softdrop_mass_gen,"AK8_softdrop_mass_gen");
  fTree->Branch("AK8_softdrop_pt_gen",&AK8_softdrop_pt_gen,"AK8_softdrop_pt_gen");
  fTree->Branch("deltaR_lak8jet",&deltaR_lak8jet,"deltaR_lak8jet/F");
  fTree->Branch("deltaphi_METak8jet",&deltaphi_METak8jet,"deltaphi_METak8jet/F");
  fTree->Branch("deltaphi_Vak8jet",&deltaphi_Vak8jet,"deltaphi_Vak8jet/F");
  fTree->Branch("deltaR_lPuppiak8jet",&deltaR_lPuppiak8jet,"deltaR_lPuppiak8jet/F");
  fTree->Branch("deltaphi_METPuppiak8jet",&deltaphi_METPuppiak8jet,"deltaphi_METPuppiak8jet/F");
  fTree->Branch("deltaphi_VPuppiak8jet",&deltaphi_VPuppiak8jet,"deltaphi_VPuppiak8jet/F");
  fTree->Branch("deltaR_lak4jetjet",&deltaR_lak4jetjet,"deltaR_lak4jetjet/F");
  fTree->Branch("deltaphi_METak4jetjet",&deltaphi_METak4jetjet,"deltaphi_METak4jetjet/F");
  fTree->Branch("deltaphi_Vak4jetjet",&deltaphi_Vak4jetjet,"deltaphi_Vak4jetjet/F");
  fTree->Branch("deltaR_lPuppiak4jetjet",&deltaR_lPuppiak4jetjet,"deltaR_lPuppiak4jetjet/F");
  fTree->Branch("deltaphi_METPuppiak4jetjet",&deltaphi_METPuppiak4jetjet,"deltaphi_METPuppiak4jetjet/F");
  fTree->Branch("deltaphi_VPuppiak4jetjet",&deltaphi_VPuppiak4jetjet,"deltaphi_VPuppiak4jetjet/F");
  fTree->Branch("v_pt",&v_pt,"v_pt/F");
  fTree->Branch("v_eta",&v_eta,"v_eta/F");
  fTree->Branch("v_phi",&v_phi,"v_phi/F");
  fTree->Branch("v_mt",&v_mt,"v_mt/F");
  fTree->Branch("v_puppi_pt",&v_puppi_pt,"v_puppi_pt/F");
  fTree->Branch("v_puppi_eta",&v_puppi_eta,"v_puppi_eta/F");
  fTree->Branch("v_puppi_phi",&v_puppi_phi,"v_puppi_phi/F");
  fTree->Branch("v_puppi_mt",&v_puppi_mt,"v_puppi_mt/F");
  fTree->Branch("mass_lvj_type0",&mass_lvj_type0,"mass_lvj_type0/F");
  fTree->Branch("mass_lvj_type0_met_jes_up",&mass_lvj_type0_met_jes_up,"mass_lvj_type0_met_jes_up/F");
  fTree->Branch("mass_lvj_type0_met_jes_dn",&mass_lvj_type0_met_jes_dn,"mass_lvj_type0_met_jes_dn/F");
  fTree->Branch("mass_lvj_type0_met_jer",&mass_lvj_type0_met_jer,"mass_lvj_type0_met_jer/F");
  fTree->Branch("mass_lvj_type0_met_jer_up",&mass_lvj_type0_met_jer_up,"mass_lvj_type0_met_jer_up/F");
  fTree->Branch("mass_lvj_type0_met_jer_dn",&mass_lvj_type0_met_jer_dn,"mass_lvj_type0_met_jer_dn/F");
  fTree->Branch("mass_lvj_type0_met_PuppiAK8_jes_up",&mass_lvj_type0_met_PuppiAK8_jes_up,"mass_lvj_type0_met_PuppiAK8_jes_up/F");
  fTree->Branch("mass_lvj_type0_met_PuppiAK8_jes_dn",&mass_lvj_type0_met_PuppiAK8_jes_dn,"mass_lvj_type0_met_PuppiAK8_jes_dn/F");
  fTree->Branch("mass_lvj_type2",&mass_lvj_type2,"mass_lvj_type2/F");
  fTree->Branch("mass_lvj_run2",&mass_lvj_run2,"mass_lvj_run2/F");
  fTree->Branch("mass_lvj_type0_PuppiAK8",&mass_lvj_type0_PuppiAK8,"mass_lvj_type0_PuppiAK8/F");
  fTree->Branch("mass_lvj_type2_PuppiAK8",&mass_lvj_type2_PuppiAK8,"mass_lvj_type2_PuppiAK8/F");
  fTree->Branch("mass_lvj_run2_PuppiAK8",&mass_lvj_run2_PuppiAK8,"mass_lvj_run2_PuppiAK8/F");
  fTree->Branch("mass_lvjj_type0_AK4",&mass_lvjj_type0_AK4,"mass_lvjj_type0_AK4/F");
  fTree->Branch("mass_lvjj_type0_met_jes_up_AK4",&mass_lvjj_type0_met_jes_up_AK4,"mass_lvjj_type0_met_jes_up_AK4/F");
  fTree->Branch("mass_lvjj_type0_met_jes_dn_AK4",&mass_lvjj_type0_met_jes_dn_AK4,"mass_lvjj_type0_met_jes_dn_AK4/F");
  fTree->Branch("mass_lvjj_type2_AK4",&mass_lvjj_type2_AK4,"mass_lvjj_type2_AK4/F");
  fTree->Branch("mass_lvjj_run2_AK4",&mass_lvjj_run2_AK4,"mass_lvjj_run2_AK4/F");
  fTree->Branch("mass_lvjj_type0_PuppiAK4",&mass_lvjj_type0_PuppiAK4,"mass_lvjj_type0_PuppiAK4/F");
  fTree->Branch("mass_lvjj_type0_met_jes_up_PuppiAK4",&mass_lvjj_type0_met_jes_up_PuppiAK4,"mass_lvjj_type0_met_jes_up_PuppiAK4/F");
  fTree->Branch("mass_lvjj_type0_met_jes_dn_PuppiAK4",&mass_lvjj_type0_met_jes_dn_PuppiAK4,"mass_lvjj_type0_met_jes_dn_PuppiAK4/F");
  fTree->Branch("mass_lvjj_type2_PuppiAK4",&mass_lvjj_type2_PuppiAK4,"mass_lvjj_type2_PuppiAK4/F");
  fTree->Branch("mass_lvjj_run2_PuppiAK4",&mass_lvjj_run2_PuppiAK4,"mass_lvjj_run2_PuppiAK4/F");
  fTree->Branch("njets",&njets,"njets/I");
  fTree->Branch("njetsPuppi",&njetsPuppi,"njetsPuppi/I");
  fTree->Branch("njets_unmerged",&njets_unmerged,"njets_unmerged/I");
  fTree->Branch("njetsPuppi_unmerged",&njetsPuppi_unmerged,"njetsPuppi_unmerged/I");
  fTree->Branch("nBTagJet_loose",&nBTagJet_loose,"nBTagJet_loose/I");
  fTree->Branch("nBTagJet_medium",&nBTagJet_medium,"nBTagJet_medium/I");
  fTree->Branch("nBTagJet_tight",&nBTagJet_tight,"nBTagJet_tight/I");
  fTree->Branch("nBTagJetPuppi_loose",&nBTagJetPuppi_loose,"nBTagJetPuppi_loose/I");
  fTree->Branch("nBTagJetPuppi_medium",&nBTagJetPuppi_medium,"nBTagJetPuppi_medium/I");
  fTree->Branch("nBTagJetPuppi_tight",&nBTagJetPuppi_tight,"nBTagJetPuppi_tight/I");
  fTree->Branch("nBTagJet_loose_unmerged",&nBTagJet_loose_unmerged,"nBTagJet_loose_unmerged/I");
  fTree->Branch("nBTagJet_medium_unmerged",&nBTagJet_medium_unmerged,"nBTagJet_medium_unmerged/I");
  fTree->Branch("nBTagJet_tight_unmerged",&nBTagJet_tight_unmerged,"nBTagJet_tight_unmerged/I");
  fTree->Branch("nBTagJetPuppi_loose_unmerged",&nBTagJetPuppi_loose_unmerged,"nBTagJetPuppi_loose_unmerged/I");
  fTree->Branch("nBTagJetPuppi_medium_unmerged",&nBTagJetPuppi_medium_unmerged,"nBTagJetPuppi_medium_unmerged/I");
  fTree->Branch("nBTagJetPuppi_tight_unmerged",&nBTagJetPuppi_tight_unmerged,"nBTagJetPuppi_tight_unmerged/I");
  fTree->Branch("mass_leptonic_closerjet",&mass_leptonic_closerjet,"mass_leptonic_closerjet/F");
  fTree->Branch("mass_ungroomedjet_closerjet",&mass_ungroomedjet_closerjet,"mass_ungroomedjet_closerjet/F");
  fTree->Branch("AK8_closerjet_pt",&AK8_closerjet_pt,"AK8_closerjet_pt/F");
  fTree->Branch("AK8_closerjet_eta",&AK8_closerjet_eta,"AK8_closerjet_eta/F");
  fTree->Branch("AK8_closerjet_phi",&AK8_closerjet_phi,"AK8_closerjet_phi/F");
  fTree->Branch("AK8_closerjet_e",&AK8_closerjet_e,"AK8_closerjet_e/F");
  fTree->Branch("vbf_maxpt_j1_pt",&vbf_maxpt_j1_pt,"vbf_maxpt_j1_pt/F");
  fTree->Branch("vbf_maxpt_j1_pt_jes_up",&vbf_maxpt_j1_pt_jes_up,"vbf_maxpt_j1_pt_jes_up/F");
  fTree->Branch("vbf_maxpt_j1_pt_jes_dn",&vbf_maxpt_j1_pt_jes_dn,"vbf_maxpt_j1_pt_jes_dn/F");
  fTree->Branch("vbf_maxpt_j1_pt_jer",&vbf_maxpt_j1_pt_jer,"vbf_maxpt_j1_pt_jer/F");
  fTree->Branch("vbf_maxpt_j1_pt_jer_up",&vbf_maxpt_j1_pt_jer_up,"vbf_maxpt_j1_pt_jer_up/F");
  fTree->Branch("vbf_maxpt_j1_pt_jer_dn",&vbf_maxpt_j1_pt_jer_dn,"vbf_maxpt_j1_pt_jer_dn/F");
  fTree->Branch("vbf_maxpt_j1_eta",&vbf_maxpt_j1_eta,"vbf_maxpt_j1_eta/F");
  fTree->Branch("vbf_maxpt_j1_eta_jes_up",&vbf_maxpt_j1_eta_jes_up,"vbf_maxpt_j1_eta_jes_up/F");
  fTree->Branch("vbf_maxpt_j1_eta_jes_dn",&vbf_maxpt_j1_eta_jes_dn,"vbf_maxpt_j1_eta_jes_dn/F");
  fTree->Branch("vbf_maxpt_j1_eta_jer",&vbf_maxpt_j1_eta_jer,"vbf_maxpt_j1_eta_jer/F");
  fTree->Branch("vbf_maxpt_j1_eta_jer_up",&vbf_maxpt_j1_eta_jer_up,"vbf_maxpt_j1_eta_jer_up/F");
  fTree->Branch("vbf_maxpt_j1_eta_jer_dn",&vbf_maxpt_j1_eta_jer_dn,"vbf_maxpt_j1_eta_jer_dn/F");
  fTree->Branch("vbf_maxpt_j1_phi",&vbf_maxpt_j1_phi,"vbf_maxpt_j1_phi/F");
  fTree->Branch("vbf_maxpt_j1_e",&vbf_maxpt_j1_e,"vbf_maxpt_j1_e/F");
  fTree->Branch("vbf_maxpt_j1_bDiscriminatorCSV",&vbf_maxpt_j1_bDiscriminatorCSV,"vbf_maxpt_j1_bDiscriminatorCSV/F");
  fTree->Branch("vbf_maxpt_j2_pt",&vbf_maxpt_j2_pt,"vbf_maxpt_j2_pt/F");
  fTree->Branch("vbf_maxpt_j2_pt_jes_up",&vbf_maxpt_j2_pt_jes_up,"vbf_maxpt_j2_pt_jes_up/F");
  fTree->Branch("vbf_maxpt_j2_pt_jes_dn",&vbf_maxpt_j2_pt_jes_dn,"vbf_maxpt_j2_pt_jes_dn/F");
  fTree->Branch("vbf_maxpt_j2_pt_jer",&vbf_maxpt_j2_pt_jer,"vbf_maxpt_j2_pt_jer/F");
  fTree->Branch("vbf_maxpt_j2_pt_jer_up",&vbf_maxpt_j2_pt_jer_up,"vbf_maxpt_j2_pt_jer_up/F");
  fTree->Branch("vbf_maxpt_j2_pt_jer_dn",&vbf_maxpt_j2_pt_jer_dn,"vbf_maxpt_j2_pt_jer_dn/F");
  fTree->Branch("vbf_maxpt_j2_eta",&vbf_maxpt_j2_eta,"vbf_maxpt_j2_eta/F");
  fTree->Branch("vbf_maxpt_j2_eta_jes_up",&vbf_maxpt_j2_eta_jes_up,"vbf_maxpt_j2_eta_jes_up/F");
  fTree->Branch("vbf_maxpt_j2_eta_jes_dn",&vbf_maxpt_j2_eta_jes_dn,"vbf_maxpt_j2_eta_jes_dn/F");
  fTree->Branch("vbf_maxpt_j2_eta_jer",&vbf_maxpt_j2_eta_jer,"vbf_maxpt_j2_eta_jer/F");
  fTree->Branch("vbf_maxpt_j2_eta_jer_up",&vbf_maxpt_j2_eta_jer_up,"vbf_maxpt_j2_eta_jer_up/F");
  fTree->Branch("vbf_maxpt_j2_eta_jer_dn",&vbf_maxpt_j2_eta_jer_dn,"vbf_maxpt_j2_eta_jer_dn/F");
  fTree->Branch("vbf_maxpt_j2_phi",&vbf_maxpt_j2_phi,"vbf_maxpt_j2_phi/F");
  fTree->Branch("vbf_maxpt_j2_e",&vbf_maxpt_j2_e,"vbf_maxpt_j2_e/F");
  fTree->Branch("vbf_maxpt_j2_bDiscriminatorCSV",&vbf_maxpt_j2_bDiscriminatorCSV,"vbf_maxpt_j2_bDiscriminatorCSV/F");
  fTree->Branch("vbf_maxpt_jj_pt",&vbf_maxpt_jj_pt,"vbf_maxpt_jj_pt/F");
  fTree->Branch("vbf_maxpt_jj_eta",&vbf_maxpt_jj_eta,"vbf_maxpt_jj_eta/F");
  fTree->Branch("vbf_maxpt_jj_phi",&vbf_maxpt_jj_phi,"vbf_maxpt_jj_phi/F");
  fTree->Branch("vbf_maxpt_jj_m",&vbf_maxpt_jj_m,"vbf_maxpt_jj_m/F");
  fTree->Branch("jet2_pt",&jet2_pt,"jet2_pt/F");
  fTree->Branch("jet2_eta",&jet2_eta,"jet2_eta/F");
  fTree->Branch("jet2_phi",&jet2_phi,"jet2_phi/F");
  fTree->Branch("jet2_e",&jet2_e,"jet2_e/F");
  fTree->Branch("jet2_btag",&jet2_btag,"jet2_btag/F");
  fTree->Branch("jet3_pt",&jet3_pt,"jet3_pt/F");
  fTree->Branch("jet3_eta",&jet3_eta,"jet3_eta/F");
  fTree->Branch("jet3_phi",&jet3_phi,"jet3_phi/F");
  fTree->Branch("jet3_e",&jet3_e,"jet3_e/F");
  fTree->Branch("jet3_btag",&jet3_btag,"jet3_btag/F");
  fTree->Branch("deltaR_AK8_closestBtagJet",&deltaR_AK8_closestBtagJet,"deltaR_AK8_closestBtagJet/F");
  fTree->Branch("deltaR_AK8_closestBtagJet_loose",&deltaR_AK8_closestBtagJet_loose,"deltaR_AK8_closestBtagJet_loose/F");
  fTree->Branch("vbf_maxpt_deltaR",&vbf_maxpt_deltaR,"vbf_maxpt_deltaR/F");

  fTree->Branch("AK4_1_pt_gen",&AK4_1_pt_gen,"AK4_1_pt_gen/F");
  fTree->Branch("AK4_1_eta_gen",&AK4_1_eta_gen,"AK4_1_eta_gen/F");
  fTree->Branch("AK4_1_phi_gen",&AK4_1_phi_gen,"AK4_1_phi_gen/F");
  fTree->Branch("AK4_1_e_gen",&AK4_1_e_gen,"AK4_1_e_gen/F");
  fTree->Branch("AK4_1_mass_gen",&AK4_1_mass_gen,"AK4_1_mass_gen/F");
  fTree->Branch("AK4_2_pt_gen",&AK4_2_pt_gen,"AK4_2_pt_gen/F");
  fTree->Branch("AK4_2_eta_gen",&AK4_2_eta_gen,"AK4_2_eta_gen/F");
  fTree->Branch("AK4_2_phi_gen",&AK4_2_phi_gen,"AK4_2_phi_gen/F");
  fTree->Branch("AK4_2_e_gen",&AK4_2_e_gen,"AK4_2_e_gen/F");
  fTree->Branch("AK4_2_mass_gen",&AK4_2_mass_gen,"AK4_2_mass_gen/F");
  fTree->Branch("AK4_BIG_gen_mass",&AK4_BIG_gen_mass,"AK4_BIG_gen_mass/F");
  fTree->Branch("deltaR_AK4",&deltaR_AK4,"deltaR_AK4/F");
}
