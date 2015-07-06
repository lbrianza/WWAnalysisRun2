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
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   UInt_t          RunNum;
   UInt_t          LumiBlockNum;
   UInt_t          EvtNum;
   Int_t           NJets;
   Int_t           BTags;
   Int_t           NVtx;
   Float_t         Weight;
   Float_t         MHT;
   Float_t         METPt;
   Float_t         METPhi;
   Float_t         METPtRaw;
   Float_t         METPhiRaw;
   Float_t         HT;
   Float_t         DeltaPhi1;
   Float_t         DeltaPhi2;
   Float_t         DeltaPhi3;
   UShort_t        IsolatedTracksNum;
   Float_t         IsolatedTracksPt[3];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksEta[3];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksPhi[3];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksE[3];   //[IsolatedTracksNum]
   Float_t         IsolatedTracksTLorentzVector[3];   //[IsolatedTracksNum]
   UShort_t        GenBosonNum;
   Float_t         GenBosonPt[2];   //[GenBosonNum]
   Float_t         GenBosonEta[2];   //[GenBosonNum]
   Float_t         GenBosonPhi[2];   //[GenBosonNum]
   Float_t         GenBosonE[2];   //[GenBosonNum]
   Float_t         GenBosonTLorentzVector[2];   //[GenBosonNum]
   Int_t           GenBoson_GenBosonPDGId[2];   //[GenBosonNum]
   UShort_t        GenJetsNum;
   Float_t         GenJetsPt[27];   //[GenJetsNum]
   Float_t         GenJetsEta[27];   //[GenJetsNum]
   Float_t         GenJetsPhi[27];   //[GenJetsNum]
   Float_t         GenJetsE[27];   //[GenJetsNum]
   Float_t         GenJetsTLorentzVector[27];   //[GenJetsNum]
   UShort_t        GenJetsAK8Num;
   Float_t         GenJetsAK8Pt[27];   //[GenJetsNum]
   Float_t         GenJetsAK8Eta[27];   //[GenJetsNum]
   Float_t         GenJetsAK8Phi[27];   //[GenJetsNum]
   Float_t         GenJetsAK8E[27];   //[GenJetsNum]
   Float_t         GenJetsAK8TLorentzVector[27];   //[GenJetsNum]
   Float_t         GenJetsAK8_prunedMass[27];   //[GenJetsNum]
   Float_t         GenJetsAK8_softdropMass[27];   //[GenJetsNum]
   Float_t         GenJetsAK8_tau1[27];   //[GenJetsNum]
   Float_t         GenJetsAK8_tau2[27];   //[GenJetsNum]
   Float_t         GenJetsAK8_tau3[27];   //[GenJetsNum]
   UShort_t        GenMuNum;
   Float_t         GenMuPt[2];   //[GenMuNum]
   Float_t         GenMuEta[2];   //[GenMuNum]
   Float_t         GenMuPhi[2];   //[GenMuNum]
   Float_t         GenMuE[2];   //[GenMuNum]
   Float_t         GenMuTLorentzVector[2];   //[GenMuNum]
   Int_t           GenMu_GenMuFromTau[2];   //[GenMuNum]
   UShort_t        GenElecNum;
   Float_t         GenElecPt[2];   //[GenElecNum]
   Float_t         GenElecEta[2];   //[GenElecNum]
   Float_t         GenElecPhi[2];   //[GenElecNum]
   Float_t         GenElecE[2];   //[GenElecNum]
   Float_t         GenElecTLorentzVector[2];   //[GenElecNum]
   Int_t           GenElec_GenElecFromTau[2];   //[GenElecNum]
   UShort_t        GenTauNum;
   Float_t         GenTauPt[2];   //[GenTauNum]
   Float_t         GenTauEta[2];   //[GenTauNum]
   Float_t         GenTauPhi[2];   //[GenTauNum]
   Float_t         GenTauE[2];   //[GenTauNum]
   Float_t         GenTauTLorentzVector[2];   //[GenTauNum]
   Int_t           GenTau_GenTauHad[2];   //[GenTauNum]
   UShort_t        GenNuNum;
   Float_t         GenNuPt[2];   //[GenNuNum]
   Float_t         GenNuEta[2];   //[GenNuNum]
   Float_t         GenNuPhi[2];   //[GenNuNum]
   Float_t         GenNuE[2];   //[GenNuNum]
   Float_t         GenNuTLorentzVector[2];   //[GenNuNum]
   UShort_t        JetsNum;
   Float_t         JetsPt[26];   //[JetsNum]
   Float_t         JetsEta[26];   //[JetsNum]
   Float_t         JetsPhi[26];   //[JetsNum]
   Float_t         JetsE[26];   //[JetsNum]
   Float_t         JetsTLorentzVector[26];   //[JetsNum]
   Float_t         Jets_bDiscriminatorICSV[26];   //[JetsNum]
   Float_t         Jets_bDiscriminatorCSV[26];   //[JetsNum]
   Float_t         Jets_chargedEmEnergyFraction[26];   //[JetsNum]
   Float_t         Jets_chargedHadronEnergyFraction[26];   //[JetsNum]
   Int_t           Jets_chargedHadronMultiplicity[26];   //[JetsNum]
   Int_t           Jets_electronMultiplicity[26];   //[JetsNum]
   Float_t         Jets_jetArea[26];   //[JetsNum]
   Float_t         Jets_muonEnergyFraction[26];   //[JetsNum]
   Int_t           Jets_muonMultiplicity[26];   //[JetsNum]
   Float_t         Jets_neutralEmEnergyFraction[26];   //[JetsNum]
   Int_t           Jets_neutralHadronMultiplicity[26];   //[JetsNum]
   Float_t         Jets_photonEnergyFraction[26];   //[JetsNum]
   Int_t           Jets_photonMultiplicity[26];   //[JetsNum]
   UChar_t         Jets_isLooseJetId[26];   //[JetsNum]
   Float_t         Jets_PtCorr[26];   //[JetsNum]
   Float_t         Jets_EtaCorr[26];   //[JetsNum]
   Float_t         Jets_PhiCorr[26];   //[JetsNum]
   Float_t         Jets_ECorr[26];   //[JetsNum]
   UShort_t        AK8JetsNum;
   Float_t         AK8JetsPt[7];   //[AK8JetsNum]
   Float_t         AK8JetsEta[7];   //[AK8JetsNum]
   Float_t         AK8JetsPhi[7];   //[AK8JetsNum]
   Float_t         AK8JetsE[7];   //[AK8JetsNum]
   Float_t         AK8JetsTLorentzVector[7];   //[AK8JetsNum]
   Float_t         AK8Jets_bDiscriminatorICSV[26];   //[JetsNum]
   Float_t         AK8Jets_bDiscriminatorCSV[26];   //[JetsNum]
   Float_t         AK8Jets_chargedEmEnergyFraction[7];   //[AK8JetsNum]
   Float_t         AK8Jets_chargedHadronEnergyFraction[7];   //[AK8JetsNum]
   Int_t           AK8Jets_chargedHadronMultiplicity[7];   //[AK8JetsNum]
   Int_t           AK8Jets_electronMultiplicity[7];   //[AK8JetsNum]
   Float_t         AK8Jets_jetArea[7];   //[AK8JetsNum]
   Float_t         AK8Jets_muonEnergyFraction[7];   //[AK8JetsNum]
   Int_t           AK8Jets_muonMultiplicity[7];   //[AK8JetsNum]
   Float_t         AK8Jets_neutralEmEnergyFraction[7];   //[AK8JetsNum]
   Int_t           AK8Jets_neutralHadronMultiplicity[7];   //[AK8JetsNum]
   Float_t         AK8Jets_photonEnergyFraction[7];   //[AK8JetsNum]
   Int_t           AK8Jets_photonMultiplicity[7];   //[AK8JetsNum]
   Float_t         AK8Jets_prunedMass[7];   //[AK8JetsNum]
   Float_t         AK8Jets_softDropMass[7];   //[AK8JetsNum]
   Float_t         AK8Jets_trimmedMass[7];   //[AK8JetsNum]
   Float_t         AK8Jets_filteredMass[7];   //[AK8JetsNum]
   Float_t         AK8Jets_tau1[7];   //[AK8JetsNum]
   Float_t         AK8Jets_tau2[7];   //[AK8JetsNum]
   Float_t         AK8Jets_tau3[7];   //[AK8JetsNum]
   UChar_t         AK8Jets_AK8isLooseJetId[7];   //[AK8JetsNum]
   Float_t         AK8Jets_PtCorr[7];   //[AK8JetsNum]
   Float_t         AK8Jets_EtaCorr[7];   //[AK8JetsNum]
   Float_t         AK8Jets_PhiCorr[7];   //[AK8JetsNum]
   Float_t         AK8Jets_ECorr[7];   //[AK8JetsNum]
   UShort_t        ElectronsNum;
   Float_t         ElectronsPt[7];   //[ElectronsNum]
   Float_t         ElectronsEta[7];   //[ElectronsNum]
   Float_t         ElectronsPhi[7];   //[ElectronsNum]
   Float_t         ElectronsE[7];   //[ElectronsNum]
   Float_t         ElectronsTLorentzVector[7];   //[ElectronsNum]
   Int_t           Electrons_charge[7];   //[ElectronsNum]
   UChar_t         Electrons_isHEEP[7];   //[ElectronsNum]
   UChar_t         Electrons_isHEEPv50[7];   //[ElectronsNum]
   Int_t           Electrons_type[7];   //[ElectronsNum]
   Float_t         Electrons_mass[7];   //[ElectronsNum]
   Float_t         Electrons_pfDeltaCorrRelIso[7];   //[ElectronsNum]
   Float_t         Electrons_pfRhoCorrRelIso04[7];   //[ElectronsNum]
   Float_t         Electrons_pfRhoCorrRelIso03[7];   //[ElectronsNum]
   Float_t         Electrons_pfRelIso[7];   //[ElectronsNum]
   Float_t         Electrons_photonIso[7];   //[ElectronsNum]
   Float_t         Electrons_neutralHadIso[7];   //[ElectronsNum]
   Float_t         Electrons_chargedHadIso[7];   //[ElectronsNum]
   Float_t         Electrons_trackIso[7];   //[ElectronsNum]
   UChar_t         Electrons_isLoose[7];   //[ElectronsNum]
   UShort_t        MuonsNum;
   Float_t         MuonsPt[33];   //[MuonsNum]
   Float_t         MuonsEta[33];   //[MuonsNum]
   Float_t         MuonsPhi[33];   //[MuonsNum]
   Float_t         MuonsE[33];   //[MuonsNum]
   Float_t         MuonsTLorentzVector[33];   //[MuonsNum]
   Int_t           Muons_charge[33];   //[MuonsNum]
   UChar_t         Muons_isHighPt[33];   //[MuonsNum]
   Int_t           Muons_type[33];   //[MuonsNum]
   Float_t         Muons_mass[33];   //[MuonsNum]
   Float_t         Muons_pfDeltaCorrRelIso[33];   //[MuonsNum]
   Float_t         Muons_pfRelIso[33];   //[MuonsNum]
   Float_t         Muons_photonIso[33];   //[MuonsNum]
   Float_t         Muons_neutralHadIso[33];   //[MuonsNum]
   Float_t         Muons_chargedHadIso[33];   //[MuonsNum]
   Float_t         Muons_trackIso[33];   //[MuonsNum]
   UChar_t         Muons_isLoose[33];   //[MuonsNum]
   UChar_t         Muons_isPFMuon[33];   //[MuonsNum]

   // List of branches
   TBranch        *b_RunNum;   //!
   TBranch        *b_LumiBlockNum;   //!
   TBranch        *b_EvtNum;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_BTags;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_MHT;   //!
   TBranch        *b_METPt;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METPtRaw;   //!
   TBranch        *b_METPhiRaw;   //!
   TBranch        *b_HT;   //!
   TBranch        *b_DeltaPhi1;   //!
   TBranch        *b_DeltaPhi2;   //!
   TBranch        *b_DeltaPhi3;   //!
   TBranch        *b_IsolatedTracksNum;   //!
   TBranch        *b_IsolatedTracksPt;   //!
   TBranch        *b_IsolatedTracksEta;   //!
   TBranch        *b_IsolatedTracksPhi;   //!
   TBranch        *b_IsolatedTracksE;   //!
   TBranch        *b_IsolatedTracksTLorentzVector;   //!
   TBranch        *b_GenBosonNum;   //!
   TBranch        *b_GenBosonPt;   //!
   TBranch        *b_GenBosonEta;   //!
   TBranch        *b_GenBosonPhi;   //!
   TBranch        *b_GenBosonE;   //!
   TBranch        *b_GenBosonTLorentzVector;   //!
   TBranch        *b_GenBoson_GenBosonPDGId;   //!
   TBranch        *b_GenJetsNum;   //!
   TBranch        *b_GenJetsPt;   //!
   TBranch        *b_GenJetsEta;   //!
   TBranch        *b_GenJetsPhi;   //!
   TBranch        *b_GenJetsE;   //!
   TBranch        *b_GenJetsTLorentzVector;   //!
   TBranch        *b_GenJetsAK8Num;   //!
   TBranch        *b_GenJetsAK8Pt;   //!
   TBranch        *b_GenJetsAK8Eta;   //!
   TBranch        *b_GenJetsAK8Phi;   //!
   TBranch        *b_GenJetsAK8E;   //!
   TBranch        *b_GenJetsAK8TLorentzVector;   //!
   TBranch        *b_GenJetsAK8_prunedMass;   //!
   TBranch        *b_GenJetsAK8_softdropMass;   //!
   TBranch        *b_GenJetsAK8_tau1;   //!
   TBranch        *b_GenJetsAK8_tau2;   //!
   TBranch        *b_GenJetsAK8_tau3;   //!
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
   TBranch        *b_Jets_bDiscriminatorICSV;   //!
   TBranch        *b_Jets_bDiscriminatorCSV;   //!
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
   TBranch        *b_Jets_isLooseJetId;   //!
   TBranch        *b_Jets_PtCorr;   //!
   TBranch        *b_Jets_EtaCorr;   //!
   TBranch        *b_Jets_PhiCorr;   //!
   TBranch        *b_Jets_ECorr;   //!
   TBranch        *b_AK8JetsNum;   //!
   TBranch        *b_AK8JetsPt;   //!
   TBranch        *b_AK8JetsEta;   //!
   TBranch        *b_AK8JetsPhi;   //!
   TBranch        *b_AK8JetsE;   //!
   TBranch        *b_AK8JetsTLorentzVector;   //!
   TBranch        *b_AK8Jets_bDiscriminatorICSV;   //!
   TBranch        *b_AK8Jets_bDiscriminatorCSV;   //!
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
   TBranch        *b_AK8Jets_softDropMass;   //!
   TBranch        *b_AK8Jets_trimmedMass;   //!
   TBranch        *b_AK8Jets_filteredMass;   //!
   TBranch        *b_AK8Jets_tau1;   //!
   TBranch        *b_AK8Jets_tau2;   //!
   TBranch        *b_AK8Jets_tau3;   //!
   TBranch        *b_AK8Jets_AK8isLooseJetId;   //!
   TBranch        *b_AK8Jets_PtCorr;   //!
   TBranch        *b_AK8Jets_EtaCorr;   //!
   TBranch        *b_AK8Jets_PhiCorr;   //!
   TBranch        *b_AK8Jets_ECorr;   //!
   TBranch        *b_ElectronsNum;   //!
   TBranch        *b_ElectronsPt;   //!
   TBranch        *b_ElectronsEta;   //!
   TBranch        *b_ElectronsPhi;   //!
   TBranch        *b_ElectronsE;   //!
   TBranch        *b_ElectronsTLorentzVector;   //!
   TBranch        *b_Electrons_charge;   //!
   TBranch        *b_Electrons_isHEEP;   //!
   TBranch        *b_Electrons_isHEEPv50;   //!
   TBranch        *b_Electrons_type;   //!
   TBranch        *b_Electrons_mass;   //!
   TBranch        *b_Electrons_pfDeltaCorrRelIso;   //!
   TBranch        *b_Electrons_pfRhoCorrRelIso04;   //!
   TBranch        *b_Electrons_pfRhoCorrRelIso03;   //!
   TBranch        *b_Electrons_pfRelIso;   //!
   TBranch        *b_Electrons_photonIso;   //!
   TBranch        *b_Electrons_neutralHadIso;   //!
   TBranch        *b_Electrons_chargedHadIso;   //!
   TBranch        *b_Electrons_trackIso;   //!
   TBranch        *b_Electrons_isLoose;   //!
   TBranch        *b_MuonsNum;   //!
   TBranch        *b_MuonsPt;   //!
   TBranch        *b_MuonsEta;   //!
   TBranch        *b_MuonsPhi;   //!
   TBranch        *b_MuonsE;   //!
   TBranch        *b_MuonsTLorentzVector;   //!
   TBranch        *b_Muons_charge;   //!
   TBranch        *b_Muons_isHighPt;   //!
   TBranch        *b_Muons_type;   //!
   TBranch        *b_Muons_mass;   //!
   TBranch        *b_Muons_pfDeltaCorrRelIso;   //!
   TBranch        *b_Muons_pfRelIso;   //!
   TBranch        *b_Muons_photonIso;   //!
   TBranch        *b_Muons_neutralHadIso;   //!
   TBranch        *b_Muons_chargedHadIso;   //!
   TBranch        *b_Muons_trackIso;   //!
   TBranch        *b_Muons_isLoose;   //!
   TBranch        *b_Muons_isPFMuon;   //!

   setInputTree(std::string inputTree);
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
