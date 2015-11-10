//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Jun  4 12:56:41 2015 by ROOT version 5.34/03
// from TTree otree/otree
// found on file: /gwteray/users/brianza/WWNtupleRun2/WWTree_mu/WWTree_pseudodata.root
//////////////////////////////////////////////////////////

#ifndef setInputPseudodata_h
#define setInputPseudodata_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class setInputPseudodata {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Int_t           njets;
   Int_t           nPV;
   Int_t           issignal;
   Float_t         wSampleWeight;
   Float_t         totalEventWeight;
   Float_t         eff_and_pu_Weight;
   Float_t         pfMET;
   Float_t         pfMET_Phi;
   Float_t         nu_pz_type0;
   Float_t         nu_pz_type2;
   Float_t         nu_pz_run2;
   Float_t         nu_pz_run2_oth;
   Int_t           nu_pz_run2_type;
   Int_t           nu_pz_isre;
   Float_t         l_pt;
   Float_t         l_eta;
   Float_t         l_phi;
   Float_t         l_e;
   Float_t         ungroomed_jet_pt;
   Float_t         ungroomed_jet_eta;
   Float_t         ungroomed_jet_phi;
   Float_t         ungroomed_jet_e;
   Float_t         jet_mass_pr;
   Float_t         jet_mass_so;
   Float_t         jet_pt_so;
   Float_t         jet_mass_tr;
   Float_t         jet_mass_fi;
   Float_t         jet_tau2tau1;
   Float_t         W_pt_gen;
   Float_t         W_pz_gen;
   Float_t         W_rap_gen;
   Float_t         genGravMass;
   Float_t         nu_pz_gen;
   Float_t         nu_pt_gen;
   Float_t         nu_phi_gen;
   Float_t         nu_eta_gen;
   Float_t         deltaR_lak8jet;
   Float_t         deltaphi_METak8jet;
   Float_t         deltaphi_Vak8jet;
   Float_t         v_pt;
   Float_t         v_eta;
   Float_t         v_phi;
   Float_t         v_mt;
   Float_t         mass_lvj_type0;
   Float_t         mass_lvj_type2;
   Float_t         mass_lvj_run2;
   Int_t           nBTagJet_loose;
   Int_t           nBTagJet_medium;
   Int_t           nBTagJet_tight;
   Float_t         mass_leptonic_closerjet;
   Float_t         mass_ungroomedjet_closerjet;
   Float_t         vbf_maxpt_j1_pt;
   Float_t         vbf_maxpt_j1_eta;
   Float_t         vbf_maxpt_j1_phi;
   Float_t         vbf_maxpt_j1_e;
   Float_t         vbf_maxpt_j1_bDiscriminatorCSV;
   Float_t         vbf_maxpt_j2_pt;
   Float_t         vbf_maxpt_j2_eta;
   Float_t         vbf_maxpt_j2_phi;
   Float_t         vbf_maxpt_j2_e;
   Float_t         vbf_maxpt_j2_bDiscriminatorCSV;
   Float_t         vbf_maxpt_jj_pt;
   Float_t         vbf_maxpt_jj_eta;
   Float_t         vbf_maxpt_jj_phi;
   Float_t         vbf_maxpt_jj_m;

   // List of branches
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_issignal;   //!
   TBranch        *b_wSampleWeight;   //!
   TBranch        *b_totalEventWeight;   //!
   TBranch        *b_eff_and_pu_Weight;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMET_Phi;   //!
   TBranch        *b_nu_pz_type0;   //!
   TBranch        *b_nu_pz_type2;   //!
   TBranch        *b_nu_pz_run2;   //!
   TBranch        *b_nu_pz_run2_oth;   //!
   TBranch        *b_nu_pz_run2_type;   //!
   TBranch        *b_nu_pz_isre;   //!
   TBranch        *b_l_pt;   //!
   TBranch        *b_l_eta;   //!
   TBranch        *b_l_phi;   //!
   TBranch        *b_l_e;   //!
   TBranch        *b_ungroomed_jet_pt;   //!
   TBranch        *b_ungroomed_jet_eta;   //!
   TBranch        *b_ungroomed_jet_phi;   //!
   TBranch        *b_ungroomed_jet_e;   //!
   TBranch        *b_jet_mass_pr;   //!
   TBranch        *b_jet_mass_so;   //!
   TBranch        *b_jet_pt_so;   //!
   TBranch        *b_jet_mass_tr;   //!
   TBranch        *b_jet_mass_fi;   //!
   TBranch        *b_jet_tau2tau1;   //!
   TBranch        *b_W_pt_gen;   //!
   TBranch        *b_W_pz_gen;   //!
   TBranch        *b_W_rap_gen;   //!
   TBranch        *b_genGravMass;   //!
   TBranch        *b_nu_pz_gen;   //!
   TBranch        *b_nu_pt_gen;   //!
   TBranch        *b_nu_phi_gen;   //!
   TBranch        *b_nu_eta_gen;   //!
   TBranch        *b_deltaR_lak8jet;   //!
   TBranch        *b_deltaphi_METak8jet;   //!
   TBranch        *b_deltaphi_Vak8jet;   //!
   TBranch        *b_v_pt;   //!
   TBranch        *b_v_eta;   //!
   TBranch        *b_v_phi;   //!
   TBranch        *b_v_mt;   //!
   TBranch        *b_mass_lvj_type0;   //!
   TBranch        *b_mass_lvj_type2;   //!
   TBranch        *b_mass_lvj_run2;   //!
   TBranch        *b_nBTagJet_loose;   //!
   TBranch        *b_nBTagJet_medium;   //!
   TBranch        *b_nBTagJet_tight;   //!
   TBranch        *b_mass_leptonic_closerjet;   //!
   TBranch        *b_mass_ungroomedjet_closerjet;   //!
   TBranch        *b_vbf_maxpt_j1_pt;   //!
   TBranch        *b_vbf_maxpt_j1_eta;   //!
   TBranch        *b_vbf_maxpt_j1_phi;   //!
   TBranch        *b_vbf_maxpt_j1_e;   //!
   TBranch        *b_vbf_maxpt_j1_bDiscriminatorCSV;   //!
   TBranch        *b_vbf_maxpt_j2_pt;   //!
   TBranch        *b_vbf_maxpt_j2_eta;   //!
   TBranch        *b_vbf_maxpt_j2_phi;   //!
   TBranch        *b_vbf_maxpt_j2_e;   //!
   TBranch        *b_vbf_maxpt_j2_bDiscriminatorCSV;   //!
   TBranch        *b_vbf_maxpt_jj_pt;   //!
   TBranch        *b_vbf_maxpt_jj_eta;   //!
   TBranch        *b_vbf_maxpt_jj_phi;   //!
   TBranch        *b_vbf_maxpt_jj_m;   //!

   setInputPseudodata(TFile* inputFile, std::string inputTree);
   virtual ~setInputPseudodata();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   void     Init();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

