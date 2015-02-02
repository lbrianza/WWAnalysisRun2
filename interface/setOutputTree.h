#ifndef __setOutputTree__
#define __setOutputTree__

#include "TTree.h"
#include "TChain.h"

// Declaration of leaf types
extern int event_runNo;
extern int event;
extern int njets;
extern int nPV;
extern float pfMET;
extern float pfMET_Phi;
extern float nu_pz_type0;
extern float nu_pz_type2;
extern float l_pt;
extern float l_eta;
extern float l_phi;
extern float l_e;
extern float ungroomed_jet_pt;
extern float ungroomed_jet_eta;
extern float ungroomed_jet_phi;
extern float ungroomed_jet_e;
extern float jet_mass_pr;
extern float jet_mass_tr;
extern float jet_mass_fi;
extern float jet_tau2tau1;
extern float W_pt_gen;
extern float W_pz_gen;
extern float gen_GravMass;
extern float nu_pz_gen;
/*extern int genBosonPdgId[10];
extern float genBosonPt[10];
extern float genBosonEta[10];
extern float genBosonPhi[10];
extern float genBosonE[10];
extern int genLeptonPdgId[10];
extern float genLeptonPt[10];
extern float genLeptonEta[10];
extern float genLeptonPhi[10];
extern float genLeptonE[10];
extern float genNuPt[10];
extern float genNuEta[10];
extern float genNuPhi[10];
extern float genNuE[10];
*/
extern float deltaR_lak8jet;
extern float deltaphi_METak8jet;
extern float deltaphi_Vak8jet;
extern float v_pt;
extern float v_eta;
extern float v_phi;
extern float v_mt;
extern float mass_lvj_type0;
extern float mass_lvj_type2;
extern float vbf_maxpt_j1_pt;
extern float vbf_maxpt_j1_eta;
extern float vbf_maxpt_j1_phi;
extern float vbf_maxpt_j1_e;
extern float vbf_maxpt_j1_bDiscriminatorCSV;
extern float vbf_maxpt_j2_pt;
extern float vbf_maxpt_j2_eta;
extern float vbf_maxpt_j2_phi;
extern float vbf_maxpt_j2_e;
extern float vbf_maxpt_j2_bDiscriminatorCSV;
extern float vbf_maxpt_jj_pt;
extern float vbf_maxpt_jj_eta;
extern float vbf_maxpt_jj_phi;
extern float vbf_maxpt_jj_m;

void init();

void SetOutTree(TTree* outTree);

#endif
