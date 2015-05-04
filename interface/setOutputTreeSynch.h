#ifndef __setOutputTreeSynch__
#define __setOutputTreeSynch__

#include "TTree.h"
#include "TChain.h"

class setOutputTreeSynch {

 public:

  TTree* fTree;

  int run;
  int event;
  int lumi;
  int nPV;

  float pfMET;
  float pfMETPhi;

  int nLooseEle; //number of electrons with looseID
  int nLooseMu; //number of electrons with looseID

  //SELECTED LEPTON - the most energetic one satisfying HEEP_ID/HighPtMuon_ID :
  float l_pt;
  float l_eta;
  float l_phi;

  //FAT JET: the most energetic AK8 jet satisfying loosejetID && cleaned from the all HEEP/highPtMuon leptons:
  float jet_pt;
  float jet_eta;
  float jet_phi;
  float jet_mass_pruned;
  float jet_mass_softdrop;
  float jet_tau2tau1;

  //W boson:
  float W_pt;
  float W_eta;
  float W_phi;

  //lvj MASS:
  float m_lvj;

  //AK4 JETS collection: - cleaned from the all HEEP/highPtMuon leptons && dR>=1.0 from the fat jet
  int njets; //AK4 jets
  int nbtag; //number of AK4 jets b-tagged with iCSVM
  float jet2_pt; //1st most energetic AK4
  float jet2_eta; //1st most energetic AK4
  float jet2_phi; //1st most energetic AK4
  float jet2_btag; //1st most energetic AK4
  float jet3_pt; //2nd most energetic AK4
  float jet3_eta; //2nd most energetic AK4
  float jet3_phi; //2nd most energetic AK4
  float jet3_btag; //2nd most energetic AK4 
  
  int issignal;
  float wSampleWeight;
  float totalEventWeight;
  float eff_and_pu_Weight;
  float nu_pz_type0;
  float nu_pz_type2;
  float nu_pz_run2;
  float nu_pz_run2_oth;
  int nu_pz_run2_type;
  int nu_pz_isre;
  float l_e;
  float jet_e;
  float jet_mass_tr;
  float jet_mass_fi;
  float W_pt_gen;
  float W_pz_gen;
  float W_rap_gen;
  float genGravMass;
  float nu_pz_gen;
  float nu_pt_gen;
  float nu_phi_gen;
  float nu_eta_gen;
  float deltaR_lak8jet;
  float deltaphi_METak8jet;
  float deltaphi_Vak8jet;
  float W_mt;
  float mass_lvj_type0;
  float mass_lvj_run2;
  float mass_leptonic_closerjet;
  float mass_ungroomedjet_closerjet;
  int nBTagJet_loose;
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



  setOutputTreeSynch(TTree* outputTree);
  //  setOutputTreeSynch(TTree *outputTree=0);
  //  setOutputTreeSynch(TFile *outputFile=0, std::string outputTreeName="WWTree");
  ~setOutputTreeSynch();

  void initializeVariables();
  
  void setBranches();

};

#endif
