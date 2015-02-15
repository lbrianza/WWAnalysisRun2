#ifndef __setOutputTree__
#define __setOutputTree__

#include "TTree.h"
#include "TChain.h"

class setOutputTree {

 public:

  TTree* fTree;
  
  int event_runNo;
  int event;
  int njets;
  int nPV;
  int issignal;
  float wSampleWeight;
  float totalEventWeight;
  float eff_and_pu_Weight;
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
  float W_pt_gen;
  float W_pz_gen;
  float genGravMass;
  float nu_pz_gen;
  float deltaR_lak8jet;
  float deltaphi_METak8jet;
  float deltaphi_Vak8jet;
  float v_pt;
  float v_eta;
  float v_phi;
  float v_mt;
  float mass_lvj_type0;
  float mass_lvj_type2;
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

  setOutputTree(TTree* outputTree);
  //  setOutputTree(TTree *outputTree=0);
  //  setOutputTree(TFile *outputFile=0, std::string outputTreeName="WWTree");
  ~setOutputTree();

  void initializeVariables();
  
  void setBranches();

};

#endif
