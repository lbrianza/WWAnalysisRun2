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

// List of branches
extern TBranch *b_event_runNo;
extern TBranch *b_event;
extern TBranch *b_njets;
extern TBranch *b_nPV;
extern TBranch *b_pfMET;
extern TBranch *b_pfMET_Phi;
extern TBranch *b_nu_pz_type0;
extern TBranch *b_nu_pz_type2;
extern TBranch *b_l_pt;
extern TBranch *b_l_eta;
extern TBranch *b_l_phi;
extern TBranch *b_l_e;
extern TBranch *b_ungroomed_jet_pt;
extern TBranch *b_ungroomed_jet_eta;
extern TBranch *b_ungroomed_jet_phi;
extern TBranch *b_ungroomed_jet_e;
extern TBranch *b_jet_mass_pr;
extern TBranch *b_jet_mass_tr;
extern TBranch *b_jet_mass_fi;
extern TBranch *b_jet_tau2tau1;
/*extern TBranch *b_genBosonPdgId[10];
extern TBranch *b_genBosonPt[10];
extern TBranch *b_genBosonEta[10];
extern TBranch *b_genBosonPhi[10];
extern TBranch *b_genBosonE[10];
extern TBranch *b_genLeptonPdgId[10];
extern TBranch *b_genLeptonPt[10];
extern TBranch *b_genLeptonEta[10];
extern TBranch *b_genLeptonPhi[10];
extern TBranch *b_genLeptonE[10];
extern TBranch *b_genNuPt[10];
extern TBranch *b_genNuEta[10];
extern TBranch *b_genNuPhi[10];
extern TBranch *b_genNuE[10];
*/
extern TBranch *b_deltaR_lak8jet;
extern TBranch *b_deltaphi_METak8jet;
extern TBranch *b_deltaphi_Vak8jet;
extern TBranch *b_v_pt;
extern TBranch *b_v_eta;
extern TBranch *b_v_phi;
extern TBranch *b_v_mt;
extern TBranch *b_mass_lvj_type0;
extern TBranch *b_mass_lvj_type2;

void InitRecoTree(TTree* nt);

void init();

void SetOutTree(TTree* outTree);

#endif
