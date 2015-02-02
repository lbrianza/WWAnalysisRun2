//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  2 14:07:00 2015 by ROOT version 5.34/14
// from TTree PreSelection/PreSelection
// found on file: /gwteray/users/brianza/WWNtupleRun2/ReducedTree/ReducedSelection_TTbar.root
//////////////////////////////////////////////////////////

#ifndef setInputTree_h
#define setInputTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include <TTree.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class setInputTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          RunNum;
   UInt_t          LumiBlockNum;
   UInt_t          EvtNum;
   Int_t           NJets;
   Int_t           BTags;
   Int_t           Leptons;
   Int_t           NVtx;
   Float_t         Weight;
   Float_t         MHT;
   Float_t         METPt;
   Float_t         METPhi;
   Float_t         HT;
   Float_t         DeltaPhi1;
   Float_t         DeltaPhi2;
   Float_t         DeltaPhi3;
   UShort_t        selectedIDIsoMuonsNum;
   Float_t         selectedIDIsoMuonsPt[4];   //[selectedIDIsoMuonsNum]
   Float_t         selectedIDIsoMuonsEta[4];   //[selectedIDIsoMuonsNum]
   Float_t         selectedIDIsoMuonsPhi[4];   //[selectedIDIsoMuonsNum]
   Float_t         selectedIDIsoMuonsE[4];   //[selectedIDIsoMuonsNum]
   Float_t         selectedIDIsoMuonsTLorentzVector[4];   //[selectedIDIsoMuonsNum]
   UShort_t        selectedIDIsoElectronsNum;
   Float_t         selectedIDIsoElectronsPt[5];   //[selectedIDIsoElectronsNum]
   Float_t         selectedIDIsoElectronsEta[5];   //[selectedIDIsoElectronsNum]
   Float_t         selectedIDIsoElectronsPhi[5];   //[selectedIDIsoElectronsNum]
   Float_t         selectedIDIsoElectronsE[5];   //[selectedIDIsoElectronsNum]
   Float_t         selectedIDIsoElectronsTLorentzVector[5];   //[selectedIDIsoElectronsNum]
   UShort_t        IsolatedTracksNum;
   Float_t         IsolatedTracksPt[5];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksEta[5];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksPhi[5];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksE[5];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksTLorentzVector[5];   //[IsolatedTracksNum]
   UShort_t        selectedIDMuonsNum;
   Float_t         selectedIDMuonsPt[6];   //[selectedIDMuonsNum]
   Float_t         selectedIDMuonsEta[6];   //[selectedIDMuonsNum]
   Float_t         selectedIDMuonsPhi[6];   //[selectedIDMuonsNum]
   Float_t         selectedIDMuonsE[6];   //[selectedIDMuonsNum]
   Float_t         selectedIDMuonsTLorentzVector[6];   //[selectedIDMuonsNum]
   UShort_t        selectedIDElectronsNum;
   Float_t         selectedIDElectronsPt[8];   //[selectedIDElectronsNum]
   Float_t         selectedIDElectronsEta[8];   //[selectedIDElectronsNum]
   Float_t         selectedIDElectronsPhi[8];   //[selectedIDElectronsNum]
   Float_t         selectedIDElectronsE[8];   //[selectedIDElectronsNum]
   Float_t         selectedIDElectronsTLorentzVector[8];   //[selectedIDElectronsNum]
   UShort_t        GenBosonNum;
   Float_t         GenBosonPt[3];   //[GenBosonNum]
   Float_t         GenBosonEta[3];   //[GenBosonNum]
   Float_t         GenBosonPhi[3];   //[GenBosonNum]
   Float_t         GenBosonE[3];   //[GenBosonNum]
   Float_t         GenBosonTLorentzVector[3];   //[GenBosonNum]
   Int_t           GenBoson_GenBosonPDGId[3];   //[GenBosonNum]
   UShort_t        GenMuNum;
   Float_t         GenMuPt[3];   //[GenMuNum]
   Float_t         GenMuEta[3];   //[GenMuNum]
   Float_t         GenMuPhi[3];   //[GenMuNum]
   Float_t         GenMuE[3];   //[GenMuNum]
   Float_t         GenMuTLorentzVector[3];   //[GenMuNum]
   Int_t           GenMu_GenMuFromTau[3];   //[GenMuNum]
   UShort_t        GenElecNum;
   Float_t         GenElecPt[3];   //[GenElecNum]
   Float_t         GenElecEta[3];   //[GenElecNum]
   Float_t         GenElecPhi[3];   //[GenElecNum]
   Float_t         GenElecE[3];   //[GenElecNum]
   Float_t         GenElecTLorentzVector[3];   //[GenElecNum]
   Int_t           GenElec_GenElecFromTau[3];   //[GenElecNum]
   UShort_t        GenTauNum;
   Float_t         GenTauPt[3];   //[GenTauNum]
   Float_t         GenTauEta[3];   //[GenTauNum]
   Float_t         GenTauPhi[3];   //[GenTauNum]
   Float_t         GenTauE[3];   //[GenTauNum]
   Float_t         GenTauTLorentzVector[3];   //[GenTauNum]
   Int_t           GenTau_GenTauHad[3];   //[GenTauNum]
   UShort_t        GenNuNum;
   Float_t         GenNuPt[3];   //[GenNuNum]
   Float_t         GenNuEta[3];   //[GenNuNum]
   Float_t         GenNuPhi[3];   //[GenNuNum]
   Float_t         GenNuE[3];   //[GenNuNum]
   Float_t         GenNuTLorentzVector[3];   //[GenNuNum]
   UShort_t        JetsNum;
   Float_t         JetsPt[36];   //[JetsNum]
   Float_t         JetsEta[36];   //[JetsNum]
   Float_t         JetsPhi[36];   //[JetsNum]
   Float_t         JetsE[36];   //[JetsNum]
   Float_t         JetsTLorentzVector[36];   //[JetsNum]
   Float_t         Jets_bDiscriminator[36];   //[JetsNum]
   Float_t         Jets_chargedEmEnergyFraction[36];   //[JetsNum]
   Float_t         Jets_chargedHadronEnergyFraction[36];   //[JetsNum]
   Int_t           Jets_chargedHadronMultiplicity[36];   //[JetsNum]
   Int_t           Jets_electronMultiplicity[36];   //[JetsNum]
   Float_t         Jets_jetArea[36];   //[JetsNum]
   Float_t         Jets_muonEnergyFraction[36];   //[JetsNum]
   Int_t           Jets_muonMultiplicity[36];   //[JetsNum]
   Float_t         Jets_neutralEmEnergyFraction[36];   //[JetsNum]
   Int_t           Jets_neutralHadronMultiplicity[36];   //[JetsNum]
   Float_t         Jets_photonEnergyFraction[36];   //[JetsNum]
   Int_t           Jets_photonMultiplicity[36];   //[JetsNum]
   UShort_t        AK8JetsNum;
   Float_t         AK8JetsPt[9];   //[AK8JetsNum]
   Float_t         AK8JetsEta[9];   //[AK8JetsNum]
   Float_t         AK8JetsPhi[9];   //[AK8JetsNum]
   Float_t         AK8JetsE[9];   //[AK8JetsNum]
   Float_t         AK8JetsTLorentzVector[9];   //[AK8JetsNum]
   Float_t         AK8Jets_chargedEmEnergyFraction[9];   //[AK8JetsNum]
   Float_t         AK8Jets_chargedHadronEnergyFraction[9];   //[AK8JetsNum]
   Int_t           AK8Jets_chargedHadronMultiplicity[9];   //[AK8JetsNum]
   Int_t           AK8Jets_electronMultiplicity[9];   //[AK8JetsNum]
   Float_t         AK8Jets_jetArea[9];   //[AK8JetsNum]
   Float_t         AK8Jets_muonEnergyFraction[9];   //[AK8JetsNum]
   Int_t           AK8Jets_muonMultiplicity[9];   //[AK8JetsNum]
   Float_t         AK8Jets_neutralEmEnergyFraction[9];   //[AK8JetsNum]
   Int_t           AK8Jets_neutralHadronMultiplicity[9];   //[AK8JetsNum]
   Float_t         AK8Jets_photonEnergyFraction[9];   //[AK8JetsNum]
   Int_t           AK8Jets_photonMultiplicity[9];   //[AK8JetsNum]
   Float_t         AK8Jets_prunedMass[9];   //[AK8JetsNum]
   Float_t         AK8Jets_trimmedMass[9];   //[AK8JetsNum]
   Float_t         AK8Jets_filteredMass[9];   //[AK8JetsNum]
   Float_t         AK8Jets_tau1[9];   //[AK8JetsNum]
   Float_t         AK8Jets_tau2[9];   //[AK8JetsNum]
   Float_t         AK8Jets_tau3[9];   //[AK8JetsNum]

   // List of branches
   TBranch        *b_RunNum;   //!
   TBranch        *b_LumiBlockNum;   //!
   TBranch        *b_EvtNum;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_BTags;   //!
   TBranch        *b_Leptons;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_METPt;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_DeltaPhi1;   //!
   TBranch        *b_DeltaPhi2;   //!
   TBranch        *b_DeltaPhi3;   //!
   TBranch        *b_selectedIDIsoMuonsNum;   //!
   TBranch        *b_selectedIDIsoMuonsPt;   //!
   TBranch        *b_selectedIDIsoMuonsEta;   //!
   TBranch        *b_selectedIDIsoMuonsPhi;   //!
   TBranch        *b_selectedIDIsoMuonsE;   //!
   TBranch        *b_selectedIDIsoMuonsTLorentzVector;   //!
   TBranch        *b_selectedIDIsoElectronsNum;   //!
   TBranch        *b_selectedIDIsoElectronsPt;   //!
   TBranch        *b_selectedIDIsoElectronsEta;   //!
   TBranch        *b_selectedIDIsoElectronsPhi;   //!
   TBranch        *b_selectedIDIsoElectronsE;   //!
   TBranch        *b_selectedIDIsoElectronsTLorentzVector;   //!
   TBranch        *b_IsolatedTracksNum;   //!
   TBranch        *b_IsolatedTracksPt;   //!
   TBranch        *b_IsolatedTracksEta;   //!
   TBranch        *b_IsolatedTracksPhi;   //!
   TBranch        *b_IsolatedTracksE;   //!
   TBranch        *b_IsolatedTracksTLorentzVector;   //!
   TBranch        *b_selectedIDMuonsNum;   //!
   TBranch        *b_selectedIDMuonsPt;   //!
   TBranch        *b_selectedIDMuonsEta;   //!
   TBranch        *b_selectedIDMuonsPhi;   //!
   TBranch        *b_selectedIDMuonsE;   //!
   TBranch        *b_selectedIDMuonsTLorentzVector;   //!
   TBranch        *b_selectedIDElectronsNum;   //!
   TBranch        *b_selectedIDElectronsPt;   //!
   TBranch        *b_selectedIDElectronsEta;   //!
   TBranch        *b_selectedIDElectronsPhi;   //!
   TBranch        *b_selectedIDElectronsE;   //!
   TBranch        *b_selectedIDElectronsTLorentzVector;   //!
   TBranch        *b_GenBosonNum;   //!
   TBranch        *b_GenBosonPt;   //!
   TBranch        *b_GenBosonEta;   //!
   TBranch        *b_GenBosonPhi;   //!
   TBranch        *b_GenBosonE;   //!
   TBranch        *b_GenBosonTLorentzVector;   //!
   TBranch        *b_GenBoson_GenBosonPDGId;   //!
   TBranch        *b_GenMuNum;   //!
   TBranch        *b_GenMuPt;   //!
   TBranch        *b_GenMuEta;   //!
   TBranch        *b_GenMuPhi;   //!
   TBranch        *b_GenMuE;   //!
   TBranch        *b_GenMuTLorentzVector;   //!
   TBranch        *b_GenMu_GenMuFromTau;   //!
   TBranch        *b_GenElecNum;   //!
   TBranch        *b_GenElecPt;   //!
   TBranch        *b_GenElecEta;   //!
   TBranch        *b_GenElecPhi;   //!
   TBranch        *b_GenElecE;   //!
   TBranch        *b_GenElecTLorentzVector;   //!
   TBranch        *b_GenElec_GenElecFromTau;   //!
   TBranch        *b_GenTauNum;   //!
   TBranch        *b_GenTauPt;   //!
   TBranch        *b_GenTauEta;   //!
   TBranch        *b_GenTauPhi;   //!
   TBranch        *b_GenTauE;   //!
   TBranch        *b_GenTauTLorentzVector;   //!
   TBranch        *b_GenTau_GenTauHad;   //!
   TBranch        *b_GenNuNum;   //!
   TBranch        *b_GenNuPt;   //!
   TBranch        *b_GenNuEta;   //!
   TBranch        *b_GenNuPhi;   //!
   TBranch        *b_GenNuE;   //!
   TBranch        *b_GenNuTLorentzVector;   //!
   TBranch        *b_JetsNum;   //!
   TBranch        *b_JetsPt;   //!
   TBranch        *b_JetsEta;   //!
   TBranch        *b_JetsPhi;   //!
   TBranch        *b_JetsE;   //!
   TBranch        *b_JetsTLorentzVector;   //!
   TBranch        *b_Jets_bDiscriminator;   //!
   TBranch        *b_Jets_chargedEmEnergyFraction;   //!
   TBranch        *b_Jets_chargedHadronEnergyFraction;   //!
   TBranch        *b_Jets_chargedHadronMultiplicity;   //!
   TBranch        *b_Jets_electronMultiplicity;   //!
   TBranch        *b_Jets_jetArea;   //!
   TBranch        *b_Jets_muonEnergyFraction;   //!
   TBranch        *b_Jets_muonMultiplicity;   //!
   TBranch        *b_Jets_neutralEmEnergyFraction;   //!
   TBranch        *b_Jets_neutralHadronMultiplicity;   //!
   TBranch        *b_Jets_photonEnergyFraction;   //!
   TBranch        *b_Jets_photonMultiplicity;   //!
   TBranch        *b_AK8JetsNum;   //!
   TBranch        *b_AK8JetsPt;   //!
   TBranch        *b_AK8JetsEta;   //!
   TBranch        *b_AK8JetsPhi;   //!
   TBranch        *b_AK8JetsE;   //!
   TBranch        *b_AK8JetsTLorentzVector;   //!
   TBranch        *b_AK8Jets_chargedEmEnergyFraction;   //!
   TBranch        *b_AK8Jets_chargedHadronEnergyFraction;   //!
   TBranch        *b_AK8Jets_chargedHadronMultiplicity;   //!
   TBranch        *b_AK8Jets_electronMultiplicity;   //!
   TBranch        *b_AK8Jets_jetArea;   //!
   TBranch        *b_AK8Jets_muonEnergyFraction;   //!
   TBranch        *b_AK8Jets_muonMultiplicity;   //!
   TBranch        *b_AK8Jets_neutralEmEnergyFraction;   //!
   TBranch        *b_AK8Jets_neutralHadronMultiplicity;   //!
   TBranch        *b_AK8Jets_photonEnergyFraction;   //!
   TBranch        *b_AK8Jets_photonMultiplicity;   //!
   TBranch        *b_AK8Jets_prunedMass;   //!
   TBranch        *b_AK8Jets_trimmedMass;   //!
   TBranch        *b_AK8Jets_filteredMass;   //!
   TBranch        *b_AK8Jets_tau1;   //!
   TBranch        *b_AK8Jets_tau2;   //!
   TBranch        *b_AK8Jets_tau3;   //!

   setInputTree(TFile* inputFile, std::string inputTree);
   virtual ~setInputTree();
   Int_t    Cut(Long64_t entry);
   Int_t    GetEntry(Long64_t entry);
   Long64_t LoadTree(Long64_t entry);
   void     Init();
   //   virtual void     Loop();
   Bool_t   Notify();
   void     Show(Long64_t entry = -1);
};

#endif 
