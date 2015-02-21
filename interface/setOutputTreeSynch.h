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
  

  setOutputTreeSynch(TTree* outputTree);
  //  setOutputTreeSynch(TTree *outputTree=0);
  //  setOutputTreeSynch(TFile *outputFile=0, std::string outputTreeName="WWTree");
  ~setOutputTreeSynch();

  void initializeVariables();
  
  void setBranches();

};

#endif
