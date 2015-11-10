//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri May 29 09:20:13 2015 by ROOT version 5.34/03
// from TTree PKUTree/PKUTree
// found on file: pseudoData/mu_PKUTree_pdata.root
//////////////////////////////////////////////////////////

#ifndef pseudodataNtuple_h
#define pseudodataNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class pseudodataNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           CategoryID;
   Int_t           vTagID;
   Double_t        massVhad;
   Double_t        tau21;
   Double_t        triggerWeight;
   Double_t        lumiWeight;
   Double_t        pileupWeight;
   Int_t           run;
   Int_t           event;
   Int_t           lumi;
   Int_t           nPV;
   Float_t         pfMET;
   Float_t         pfMETPhi;
   Double_t        weight;
   Int_t           nLooseEle;
   Int_t           nLooseMu;
   Float_t         l_pt;
   Float_t         l_eta;
   Float_t         l_phi;
   Float_t         jet_pt;
   Float_t         jet_eta;
   Float_t         jet_phi;
   Float_t         jet_mass_pruned;
   Float_t         jet_mass_softdrop;
   Float_t         jet_pt_softdrop;
   Float_t         jet_tau2tau1;
   Float_t         W_pt;
   Float_t         W_eta;
   Float_t         W_phi;
   Float_t         m_lvj;
   Int_t           njets;
   Int_t           nbtag;
   Float_t         jet2_pt;
   Float_t         jet2_btag;
   Float_t         jet3_pt;
   Float_t         jet3_btag;
   Double_t        jetAK8_pt;
   Double_t        jetAK8_mass;
   Double_t        jetAK8_jec;
   Double_t        METraw_et;
   Double_t        METraw_phi;
   Double_t        METraw_sumEt;
   Double_t        MET_et;
   Double_t        MET_phi;
   Double_t        MET_sumEt;
   Double_t        MET_corrPx;
   Double_t        MET_corrPy;

   // List of branches
   TBranch        *b_CategoryID;   //!
   TBranch        *b_vTagID;   //!
   TBranch        *b_massVhad;   //!
   TBranch        *b_tau21;   //!
   TBranch        *b_triggerWeight;   //!
   TBranch        *b_lumiWeight;   //!
   TBranch        *b_pileupWeight;   //!
   TBranch        *b_run;   //!
   TBranch        *b_event;   //!
   TBranch        *b_lumi;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_pfMET;   //!
   TBranch        *b_pfMETPhi;   //!
   TBranch        *b_weight;   //!
   TBranch        *b_nLooseEle;   //!
   TBranch        *b_nLooseMu;   //!
   TBranch        *b_l_pt;   //!
   TBranch        *b_l_eta;   //!
   TBranch        *b_l_phi;   //!
   TBranch        *b_jet_pt;   //!
   TBranch        *b_jet_eta;   //!
   TBranch        *b_jet_phi;   //!
   TBranch        *b_jet_mass_pruned;   //!
   TBranch        *b_jet_mass_softdrop;   //!
   TBranch        *b_jet_pt_softdrop;   //!
   TBranch        *b_jet_tau2tau1;   //!
   TBranch        *b_W_pt;   //!
   TBranch        *b_W_eta;   //!
   TBranch        *b_W_phi;   //!
   TBranch        *b_m_lvj;   //!
   TBranch        *b_njets;   //!
   TBranch        *b_nbtag;   //!
   TBranch        *b_jet2_pt;   //!
   TBranch        *b_jet2_btag;   //!
   TBranch        *b_jet3_pt;   //!
   TBranch        *b_jet3_btag;   //!
   TBranch        *b_jetAK8_pt;   //!
   TBranch        *b_jetAK8_mass;   //!
   TBranch        *b_jetAK8_jec;   //!
   TBranch        *b_METraw_et;   //!
   TBranch        *b_METraw_phi;   //!
   TBranch        *b_METraw_sumEt;   //!
   TBranch        *b_MET_et;   //!
   TBranch        *b_MET_phi;   //!
   TBranch        *b_MET_sumEt;   //!
   TBranch        *b_MET_corrPx;   //!
   TBranch        *b_MET_corrPy;   //!

   //   pseudodataNtuple(TTree *tree=0);
   pseudodataNtuple(TFile* inputFile, std::string inputTree);
   virtual ~pseudodataNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init();
   //   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif
