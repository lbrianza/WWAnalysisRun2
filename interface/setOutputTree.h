#ifndef __setOutputTree__
#define __setOutputTree__

#include "TTree.h"
#include "TChain.h"

class setOutputTree {

 public:

  TTree* fTree;
  
  int run;
  int event;
  int lumi;
  int njets;
  int nPV;
  int issignal;
  float wSampleWeight;
  float genWeight;
  float totalEventWeight;
  float eff_and_pu_Weight;
  float pfMET;
  float pfMET_Phi;
  float nu_pz_type0;
  float nu_pz_type2;
  float nu_pz_run2;
  float nu_pz_run2_oth;
  int nu_pz_run2_type;
  int nu_pz_isre;
  float l_pt;
  float l_eta;
  float l_phi;
  float l_e;
  float ungroomed_jet_pt;
  float ungroomed_jet_eta;
  float ungroomed_jet_phi;
  float ungroomed_jet_e;
  float jet_mass_pr;
  float jet_mass_so;
  float jet_mass_tr;
  float jet_mass_fi;
  float jet_tau2tau1;
  float ttb_ungroomed_jet_pt;
  float ttb_ungroomed_jet_eta;
  float ttb_ungroomed_jet_phi;
  float ttb_ungroomed_jet_e;
  float ttb_jet_mass_pr;
  float ttb_jet_mass_so;
  float ttb_jet_mass_tr;
  float ttb_jet_mass_fi;
  float ttb_jet_tau2tau1;
  float ttb_deltaeta_lak8jet;
  float W_pt_gen;
  float W_pz_gen;
  float W_rap_gen;
  float genGravMass;
  float nu_pz_gen;
  float nu_pt_gen;
  float nu_phi_gen;
  float nu_eta_gen;
  float hadW_pt_gen;
  float hadW_eta_gen;
  float hadW_phi_gen;
  float hadW_e_gen;
  float hadW_m_gen;
  float lepW_pt_gen;
  float lepW_eta_gen;
  float lepW_phi_gen;
  float lepW_e_gen;
  float lepW_m_gen;
  float AK8_pt_gen;
  float AK8_eta_gen;
  float AK8_phi_gen;
  float AK8_e_gen;
  float AK8_pruned_mass_gen;
  float AK8_softdrop_mass_gen;
  float deltaR_lak8jet;
  float deltaphi_METak8jet;
  float deltaphi_Vak8jet;
  float v_pt;
  float v_eta;
  float v_phi;
  float v_mt;
  float mass_lvj_type0;
  float mass_lvj_type2;
  float mass_lvj_run2;
  float mass_leptonic_closerjet;
  float mass_ungroomedjet_closerjet;
  float AK8_closerjet_pt;
  float AK8_closerjet_eta;
  float AK8_closerjet_phi;
  float AK8_closerjet_e;
  int nBTagJet_loose;
  int nBTagJet_medium;
  int nBTagJet_tight;
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
  float jet2_pt;
  float jet2_btag;
  float jet3_pt;
  float jet3_btag;

  setOutputTree(TTree* outputTree);
  //  setOutputTree(TTree *outputTree=0);
  //  setOutputTree(TFile *outputFile=0, std::string outputTreeName="WWTree");
  ~setOutputTree();

  void initializeVariables();
  
  void setBranches();

};

#endif
