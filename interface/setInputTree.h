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
   Int_t           NVtx;
   Int_t           npT;
   Int_t           passFilterHBHE;
   Int_t           passFilterHBHEIso;
   Int_t           passFilterCSCHalo;
   Int_t           passFilterGoodVtx;
   Int_t           passFilterEEBadSC;
   Int_t           passFilterEcalDeadCellTriggerPrimitive;
   Int_t           passFilterGlobalTightHalo2016;
   Int_t           passFilterBadChCand;
   Int_t           passFilterBadPFMuon;
   Int_t           passFilterHBHELooseRerun;
   Int_t           passFilterHBHETightRerun;
   Int_t           passFilterHBHEIsoRerun;
   Float_t         Weight;
   Double_t         METPt;
   Double_t         METPhi;
   Double_t         METPtUp;
   Double_t         METPhiUp;
   Double_t         METPtDown;
   Double_t         METPhiDown;
   Double_t         METPtRaw;
   Double_t         METPhiRaw;
   Float_t         CaloMetPt;
   Float_t         CaloMetPhi;
   Double_t         genEventWeight;
   Float_t         PUWeight;
   Double_t         METpuppiPt;
   Double_t         METpuppiPhi;
   Double_t         METpuppiPtUp;
   Double_t         METpuppiPhiUp;
   Double_t         METpuppiPtDown;
   Double_t         METpuppiPhiDown;
   Double_t         METpuppiPtRaw;
   Double_t         METpuppiPhiRaw;
   Double_t         METpuppiCaloMetPt;
   Double_t         METpuppiCaloMetPhi;
   UShort_t        GenBosonNum;
   Float_t         GenBosonPt[30];   //[GenBosonNum]
   Float_t         GenBosonEta[30];   //[GenBosonNum]
   Float_t         GenBosonPhi[30];   //[GenBosonNum]
   Float_t         GenBosonE[30];   //[GenBosonNum]
   Float_t         GenBosonTLorentzVector[30];   //[GenBosonNum]
   Int_t           GenBoson_GenBosonPDGId[30];   //[GenBosonNum]
   Int_t           GenBoson_isBosonLeptonic[30];   //[GenBosonNum]
   UShort_t        GenMuNum;
   Float_t         GenMuPt[30];   //[GenMuNum]
   Float_t         GenMuEta[30];   //[GenMuNum]
   Float_t         GenMuPhi[30];   //[GenMuNum]
   Float_t         GenMuE[30];   //[GenMuNum]
   Float_t         GenMuTLorentzVector[30];   //[GenMuNum]
   Int_t           GenMu_GenMuFromTau[30];   //[GenMuNum]
   UShort_t        GenElecNum;
   Float_t         GenElecPt[30];   //[GenElecNum]
   Float_t         GenElecEta[30];   //[GenElecNum]
   Float_t         GenElecPhi[30];   //[GenElecNum]
   Float_t         GenElecE[30];   //[GenElecNum]
   Float_t         GenElecTLorentzVector[30];   //[GenElecNum]
   Int_t           GenElec_GenElecFromTau[30];   //[GenElecNum]
   UShort_t        GenTauNum;
   Float_t         GenTauPt[30];   //[GenTauNum]
   Float_t         GenTauEta[30];   //[GenTauNum]
   Float_t         GenTauPhi[30];   //[GenTauNum]
   Float_t         GenTauE[30];   //[GenTauNum]
   Float_t         GenTauTLorentzVector[30];   //[GenTauNum]
   Int_t           GenTau_GenTauHad[30];   //[GenTauNum]
   UShort_t        GenNuNum;
   Float_t         GenNuPt[30];   //[GenNuNum]
   Float_t         GenNuEta[30];   //[GenNuNum]
   Float_t         GenNuPhi[30];   //[GenNuNum]
   Float_t         GenNuE[30];   //[GenNuNum]
   Float_t         GenNuTLorentzVector[30];   //[GenNuNum]
   UShort_t        GenTopNum;
   Float_t         GenTopPt[30];   //[GenTopNum]
   Float_t         GenTopEta[30];   //[GenTopNum]
   Float_t         GenTopPhi[30];   //[GenTopNum]
   Float_t         GenTopE[30];   //[GenTopNum]
   Float_t         GenTopTLorentzVector[30];   //[GenTopNum]
   Int_t           GenTop_GenTopPDGId[30];   //[GenTopNum]
   UShort_t        GenJetsNum;
   Float_t         GenJetsPt[60];   //[GenJetsNum]
   Float_t         GenJetsEta[60];   //[GenJetsNum]
   Float_t         GenJetsPhi[60];   //[GenJetsNum]
   Float_t         GenJetsE[60];   //[GenJetsNum]
   Float_t         GenJetsTLorentzVector[60];   //[GenJetsNum]
   UShort_t        GenJetsAK8Num;
   Float_t         GenJetsAK8Pt[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8Eta[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8Phi[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8E[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8TLorentzVector[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8_prunedMass[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8_softdropMass[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8_softdropPt[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8_tau1[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8_tau2[30];   //[GenJetsAK8Num]
   Float_t         GenJetsAK8_tau3[30];   //[GenJetsAK8Num]
   UShort_t        JetsNum;
   Float_t         JetsPt[25];   //[JetsNum]
   Float_t         JetsEta[25];   //[JetsNum]
   Float_t         JetsPhi[25];   //[JetsNum]
   Float_t         JetsE[25];   //[JetsNum]
   Float_t         JetsTLorentzVector[25];   //[JetsNum]
   Float_t         Jets_bDiscriminatorCSV[25];   //[JetsNum]
   Float_t         Jets_bDiscriminatorICSV[25];   //[JetsNum]
   Float_t         Jets_chargedEmEnergyFraction[25];   //[JetsNum]
   Float_t         Jets_chargedHadronEnergyFraction[25];   //[JetsNum]
   Int_t           Jets_chargedHadronMultiplicity[25];   //[JetsNum]
   Int_t           Jets_electronMultiplicity[25];   //[JetsNum]
   Float_t         Jets_jetArea[25];   //[JetsNum]
   Float_t         Jets_muonEnergyFraction[25];   //[JetsNum]
   Int_t           Jets_muonMultiplicity[25];   //[JetsNum]
   Float_t         Jets_neutralEmEnergyFraction[25];   //[JetsNum]
   Int_t           Jets_neutralHadronMultiplicity[25];   //[JetsNum]
   Float_t         Jets_photonEnergyFraction[25];   //[JetsNum]
   Int_t           Jets_photonMultiplicity[25];   //[JetsNum]
   UChar_t         Jets_isLooseJetId[25];   //[JetsNum]
   UChar_t         Jets_isTightJetId[25];   //[JetsNum]
   Float_t         Jets_PtCorr[25];   //[JetsNum]
   Float_t         Jets_EtaCorr[25];   //[JetsNum]
   Float_t         Jets_PhiCorr[25];   //[JetsNum]
   Float_t         Jets_ECorr[25];   //[JetsNum]
   Float_t         Jets_AK4correction[25];   //[JetsNum]
   Float_t         Jets_AK4correctionUp[25];   //[JetsNum]
   Float_t         Jets_AK4correctionDown[25];   //[JetsNum]
   UShort_t        JetsPuppiNum;
   Float_t         JetsPuppiPt[25];   //[JetsPuppiNum]
   Float_t         JetsPuppiEta[25];   //[JetsPuppiNum]
   Float_t         JetsPuppiPhi[25];   //[JetsPuppiNum]
   Float_t         JetsPuppiE[25];   //[JetsPuppiNum]
   Float_t         JetsPuppiTLorentzVector[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_bDiscriminatorCSV[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_bDiscriminatorICSV[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_chargedEmEnergyFraction[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_chargedHadronEnergyFraction[25];   //[JetsPuppiNum]
   Int_t           JetsPuppi_chargedHadronMultiplicity[25];   //[JetsPuppiNum]
   Int_t           JetsPuppi_electronMultiplicity[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_jetArea[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_muonEnergyFraction[25];   //[JetsPuppiNum]
   Int_t           JetsPuppi_muonMultiplicity[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_neutralEmEnergyFraction[25];   //[JetsPuppiNum]
   Int_t           JetsPuppi_neutralHadronMultiplicity[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_photonEnergyFraction[25];   //[JetsPuppiNum]
   Int_t           JetsPuppi_photonMultiplicity[25];   //[JetsPuppiNum]
   UChar_t         JetsPuppi_isLooseJetId[25];   //[JetsPuppiNum]
   UChar_t         JetsPuppi_isTightJetId[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_PtCorr[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_EtaCorr[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_PhiCorr[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_ECorr[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_AK4correction[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_AK4correctionUp[25];   //[JetsPuppiNum]
   Float_t         JetsPuppi_AK4correctionDown[25];   //[JetsPuppiNum]
   UShort_t        AK8JetsNum;
   Float_t         AK8JetsPt[30];   //[AK8JetsNum]
   Float_t         AK8JetsEta[30];   //[AK8JetsNum]
   Float_t         AK8JetsPhi[30];   //[AK8JetsNum]
   Float_t         AK8JetsE[30];   //[AK8JetsNum]
   Float_t         AK8JetsTLorentzVector[30];   //[AK8JetsNum]
   Float_t         AK8Jets_bDiscriminatorCSV[30];   //[AK8JetsNum]
   Float_t         AK8Jets_bDiscriminatorICSV[30];   //[AK8JetsNum]
   Float_t         AK8Jets_chargedEmEnergyFraction[30];   //[AK8JetsNum]
   Float_t         AK8Jets_chargedHadronEnergyFraction[30];   //[AK8JetsNum]
   Int_t           AK8Jets_chargedHadronMultiplicity[30];   //[AK8JetsNum]
   Int_t           AK8Jets_electronMultiplicity[30];   //[AK8JetsNum]
   Float_t         AK8Jets_jetArea[30];   //[AK8JetsNum]
   Float_t         AK8Jets_muonEnergyFraction[30];   //[AK8JetsNum]
   Int_t           AK8Jets_muonMultiplicity[30];   //[AK8JetsNum]
   Float_t         AK8Jets_neutralEmEnergyFraction[30];   //[AK8JetsNum]
   Int_t           AK8Jets_neutralHadronMultiplicity[30];   //[AK8JetsNum]
   Float_t         AK8Jets_photonEnergyFraction[30];   //[AK8JetsNum]
   Int_t           AK8Jets_photonMultiplicity[30];   //[AK8JetsNum]
   Float_t         AK8Jets_prunedMass[30];   //[AK8JetsNum]
   Float_t         AK8Jets_softDropMass[30];   //[AK8JetsNum]
   Float_t         AK8Jets_softDropPt[30];   //[AK8JetsNum]
   Float_t         AK8Jets_trimmedMass[30];   //[AK8JetsNum]
   Float_t         AK8Jets_filteredMass[30];   //[AK8JetsNum]
   Float_t         AK8Jets_tau1[30];   //[AK8JetsNum]
   Float_t         AK8Jets_tau2[30];   //[AK8JetsNum]
   Float_t         AK8Jets_tau3[30];   //[AK8JetsNum]
   UChar_t         AK8Jets_AK8isLooseJetId[30];   //[AK8JetsNum]
   UChar_t         AK8Jets_AK8isTightJetId[30];   //[AK8JetsNum]
   Float_t         AK8Jets_PtCorr[30];   //[AK8JetsNum]
   Float_t         AK8Jets_EtaCorr[30];   //[AK8JetsNum]
   Float_t         AK8Jets_PhiCorr[30];   //[AK8JetsNum]
   Float_t         AK8Jets_ECorr[30];   //[AK8JetsNum]
   Float_t         AK8Jets_mass[30];   //[AK8JetsNum]
   Float_t         AK8Jets_AK8correction[30];   //[AK8JetsNum]
   Float_t         AK8Jets_AK8correctionUp[30];   //[AK8JetsNum]
   Float_t         AK8Jets_AK8correctionDown[30];   //[AK8JetsNum]
   Float_t         AK8Jets_AK8massCorrection[30];   //[AK8JetsNum]
   Float_t         AK8Jets_AK8massCorrectionUp[30];   //[AK8JetsNum]
   Float_t         AK8Jets_AK8massCorrectionDown[30];   //[AK8JetsNum]
   UShort_t        PuppiAK8JetsNum;
   Float_t         PuppiAK8JetsPt[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8JetsEta[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8JetsPhi[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8JetsE[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8JetsTLorentzVector[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_bDiscriminatorCSV[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_bDiscriminatorICSV[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_chargedEmEnergyFraction[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_chargedHadronEnergyFraction[30];   //[PuppiAK8JetsNum]
   Int_t           PuppiAK8Jets_chargedHadronMultiplicity[30];   //[PuppiAK8JetsNum]
   Int_t           PuppiAK8Jets_electronMultiplicity[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_jetArea[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_muonEnergyFraction[30];   //[PuppiAK8JetsNum]
   Int_t           PuppiAK8Jets_muonMultiplicity[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_neutralEmEnergyFraction[30];   //[PuppiAK8JetsNum]
   Int_t           PuppiAK8Jets_neutralHadronMultiplicity[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_photonEnergyFraction[30];   //[PuppiAK8JetsNum]
   Int_t           PuppiAK8Jets_photonMultiplicity[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_prunedMass[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_softDropMass[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_softDropPt[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_trimmedMass[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_filteredMass[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_tau1[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_tau2[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_tau3[30];   //[PuppiAK8JetsNum]
   UChar_t         PuppiAK8Jets_PuppiAK8isLooseJetId[30];   //[PuppiAK8JetsNum]
   UChar_t         PuppiAK8Jets_PuppiAK8isTightJetId[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PtCorr[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_EtaCorr[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PhiCorr[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_ECorr[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_mass[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PuppiAK8correction[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PuppiAK8correctionUp[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PuppiAK8correctionDown[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PuppiAK8massCorrection[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PuppiAK8massCorrectionUp[30];   //[PuppiAK8JetsNum]
   Float_t         PuppiAK8Jets_PuppiAK8massCorrectionDown[30];   //[PuppiAK8JetsNum]
   UShort_t        ElectronsNum;
   Float_t         ElectronsPt[30];   //[ElectronsNum]
   Float_t         ElectronsEta[30];   //[ElectronsNum]
   Float_t         ElectronsPhi[30];   //[ElectronsNum]
   Float_t         ElectronsE[30];   //[ElectronsNum]
   Float_t         ElectronsTLorentzVector[30];   //[ElectronsNum]
   Int_t           Electrons_charge[30];   //[ElectronsNum]
   UChar_t         Electrons_isHEEP[30];   //[ElectronsNum]
   Int_t           Electrons_type[30];   //[ElectronsNum]
   Float_t         Electrons_mass[30];   //[ElectronsNum]
   Float_t         Electrons_pfDeltaCorrRelIso[30];   //[ElectronsNum]
   Float_t         Electrons_pfRhoCorrRelIso04[30];   //[ElectronsNum]
   Float_t         Electrons_pfRhoCorrRelIso03[30];   //[ElectronsNum]
   Float_t         Electrons_pfRelIso[30];   //[ElectronsNum]
   Float_t         Electrons_photonIso[30];   //[ElectronsNum]
   Float_t         Electrons_neutralHadIso[30];   //[ElectronsNum]
   Float_t         Electrons_chargedHadIso[30];   //[ElectronsNum]
   Float_t         Electrons_trackIso[30];   //[ElectronsNum]
   UChar_t         Electrons_isLoose[30];   //[ElectronsNum]
   UChar_t         Electrons_isMedium[30];   //[ElectronsNum]
   UChar_t         Electrons_isTight[30];   //[ElectronsNum]
   Float_t         Electrons_SCEnergy[30];   //[ElectronsNum]
   Float_t         Electrons_deltaEtaSCTracker[30];   //[ElectronsNum]
   Float_t         Electrons_deltaPhiSCTracker[30];   //[ElectronsNum]
   Float_t         Electrons_sigmaIetaIeta[30];   //[ElectronsNum]
   Float_t         Electrons_sigmaIphiIphi[30];   //[ElectronsNum]
   UShort_t        MuonsNum;
   Float_t         MuonsPt[70];   //[MuonsNum]
   Float_t         MuonsEta[70];   //[MuonsNum]
   Float_t         MuonsPhi[70];   //[MuonsNum]
   Float_t         MuonsE[70];   //[MuonsNum]
   Float_t         MuonsTLorentzVector[70];   //[MuonsNum]
   Int_t           Muons_charge[70];   //[MuonsNum]
   UChar_t         Muons_isHighPt[70];   //[MuonsNum]
   Int_t           Muons_type[70];   //[MuonsNum]
   Float_t         Muons_mass[70];   //[MuonsNum]
   Float_t         Muons_pfDeltaCorrRelIso[70];   //[MuonsNum]
   Float_t         Muons_pfRelIso[70];   //[MuonsNum]
   Float_t         Muons_photonIso[70];   //[MuonsNum]
   Float_t         Muons_neutralHadIso[70];   //[MuonsNum]
   Float_t         Muons_chargedHadIso[70];   //[MuonsNum]
   Float_t         Muons_trackIso[70];   //[MuonsNum]
   UChar_t         Muons_isLoose[70];   //[MuonsNum]
   UChar_t         Muons_isMedium[70];   //[MuonsNum]
   UChar_t         Muons_isTight[70];   //[MuonsNum]
   UChar_t         Muons_isPFMuon[70];   //[MuonsNum]
   std::vector<int>     *TriggerProducerTriggerPrescales;
   std::vector<unsigned int> *TriggerProducerTriggerPass;
   std::vector<std::string>  *TriggerProducerTriggerNames;


   // List of branches
   TBranch        *b_RunNum;   //!
   TBranch        *b_LumiBlockNum;   //!
   TBranch        *b_EvtNum;   //!
   TBranch        *b_NJets;   //!
   TBranch        *b_NVtx;   //!
   TBranch        *b_npT;   //!
   TBranch        *b_passFilterHBHE;   //!
   TBranch        *b_passFilterHBHEIso;   //!
   TBranch        *b_passFilterCSCHalo;   //!
   TBranch        *b_passFilterGoodVtx;   //!
   TBranch        *b_passFilterEEBadSC;   //!
   TBranch        *b_passFilterEcalDeadCellTriggerPrimitive;
   TBranch        *b_passFilterGlobalTightHalo2016;
   TBranch        *b_passFilterBadChCand;
   TBranch        *b_passFilterBadPFMuon;
   TBranch        *b_passFilterHBHELooseRerun;   //!
   TBranch        *b_passFilterHBHETightRerun;   //!
   TBranch        *b_passFilterHBHEIsoRerun;   //!
   TBranch        *b_Weight;   //!
   TBranch        *b_METPt;   //!
   TBranch        *b_METPhi;   //!
   TBranch        *b_METPtUp;   //!
   TBranch        *b_METPhiUp;   //!
   TBranch        *b_METPtDown;   //!
   TBranch        *b_METPhiDown;   //!
   TBranch        *b_METPtRaw;   //!
   TBranch        *b_METPhiRaw;   //!
   TBranch        *b_CaloMetPt;   //!
   TBranch        *b_CaloMetPhi;   //!
   TBranch        *b_genEventWeight;   //!
   TBranch        *b_PUWeight;   //!
   TBranch        *b_METpuppiPt;   //!
   TBranch        *b_METpuppiPhi;   //!
   TBranch        *b_METpuppiPtUp;   //!
   TBranch        *b_METpuppiPhiUp;   //!
   TBranch        *b_METpuppiPtDown;   //!
   TBranch        *b_METpuppiPhiDown;   //!
   TBranch        *b_METpuppiPtRaw;   //!
   TBranch        *b_METpuppiPhiRaw;   //!
   TBranch        *b_METpuppiCaloMetPt;   //!
   TBranch        *b_METpuppiCaloMetPhi;   //!
   TBranch        *b_GenBosonNum;   //!
   TBranch        *b_GenBosonPt;   //!
   TBranch        *b_GenBosonEta;   //!
   TBranch        *b_GenBosonPhi;   //!
   TBranch        *b_GenBosonE;   //!
   TBranch        *b_GenBosonTLorentzVector;   //!
   TBranch        *b_GenBoson_GenBosonPDGId;   //!
   TBranch        *b_GenBoson_isBosonLeptonic;   //!
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
   TBranch        *b_GenTopNum;   //!
   TBranch        *b_GenTopPt;   //!
   TBranch        *b_GenTopEta;   //!
   TBranch        *b_GenTopPhi;   //!
   TBranch        *b_GenTopE;   //!
   TBranch        *b_GenTopTLorentzVector;   //!
   TBranch        *b_GenTop_GenTopPDGId;   //!
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
   TBranch        *b_GenJetsAK8_softdropPt;   //!
   TBranch        *b_GenJetsAK8_tau1;   //!
   TBranch        *b_GenJetsAK8_tau2;   //!
   TBranch        *b_GenJetsAK8_tau3;   //!
   TBranch        *b_JetsNum;   //!
   TBranch        *b_JetsPt;   //!
   TBranch        *b_JetsEta;   //!
   TBranch        *b_JetsPhi;   //!
   TBranch        *b_JetsE;   //!
   TBranch        *b_JetsTLorentzVector;   //!
   TBranch        *b_Jets_bDiscriminatorCSV;   //!
   TBranch        *b_Jets_bDiscriminatorICSV;   //!
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
   TBranch        *b_Jets_isTightJetId;   //!
   TBranch        *b_Jets_PtCorr;   //!
   TBranch        *b_Jets_EtaCorr;   //!
   TBranch        *b_Jets_PhiCorr;   //!
   TBranch        *b_Jets_ECorr;   //!
   TBranch        *b_Jets_AK4correction;   //!
   TBranch        *b_Jets_AK4correctionUp;   //!
   TBranch        *b_Jets_AK4correctionDown;   //!
   TBranch        *b_AK8JetsNum;   //!
   TBranch        *b_AK8JetsPt;   //!
   TBranch        *b_AK8JetsEta;   //!
   TBranch        *b_AK8JetsPhi;   //!
   TBranch        *b_AK8JetsE;   //!
   TBranch        *b_AK8JetsTLorentzVector;   //!
   TBranch        *b_AK8Jets_bDiscriminatorCSV;   //!
   TBranch        *b_AK8Jets_bDiscriminatorICSV;   //!
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
   TBranch        *b_AK8Jets_softDropPt;   //!
   TBranch        *b_AK8Jets_trimmedMass;   //!
   TBranch        *b_AK8Jets_filteredMass;   //!
   TBranch        *b_AK8Jets_tau1;   //!
   TBranch        *b_AK8Jets_tau2;   //!
   TBranch        *b_AK8Jets_tau3;   //!
   TBranch        *b_AK8Jets_AK8isLooseJetId;   //!
   TBranch        *b_AK8Jets_AK8isTightJetId;   //!
   TBranch        *b_AK8Jets_PtCorr;   //!
   TBranch        *b_AK8Jets_EtaCorr;   //!
   TBranch        *b_AK8Jets_PhiCorr;   //!
   TBranch        *b_AK8Jets_ECorr;   //!
   TBranch        *b_AK8Jets_mass;   //!
   TBranch        *b_AK8Jets_AK8correction;   //!
   TBranch        *b_AK8Jets_AK8correctionUp;   //!
   TBranch        *b_AK8Jets_AK8correctionDown;   //!
   TBranch        *b_AK8Jets_AK8massCorrection;   //!
   TBranch        *b_AK8Jets_AK8massCorrectionUp;   //!
   TBranch        *b_AK8Jets_AK8massCorrectionDown;   //!
   TBranch        *b_PuppiAK8JetsNum;   //!
   TBranch        *b_PuppiAK8JetsPt;   //!
   TBranch        *b_PuppiAK8JetsEta;   //!
   TBranch        *b_PuppiAK8JetsPhi;   //!
   TBranch        *b_PuppiAK8JetsE;   //!
   TBranch        *b_PuppiAK8JetsTLorentzVector;   //!
   TBranch        *b_PuppiAK8Jets_bDiscriminatorCSV;   //!
   TBranch        *b_PuppiAK8Jets_bDiscriminatorICSV;   //!
   TBranch        *b_PuppiAK8Jets_chargedEmEnergyFraction;   //!
   TBranch        *b_PuppiAK8Jets_chargedHadronEnergyFraction;   //!
   TBranch        *b_PuppiAK8Jets_chargedHadronMultiplicity;   //!
   TBranch        *b_PuppiAK8Jets_electronMultiplicity;   //!
   TBranch        *b_PuppiAK8Jets_jetArea;   //!
   TBranch        *b_PuppiAK8Jets_muonEnergyFraction;   //!
   TBranch        *b_PuppiAK8Jets_muonMultiplicity;   //!
   TBranch        *b_PuppiAK8Jets_neutralEmEnergyFraction;   //!
   TBranch        *b_PuppiAK8Jets_neutralHadronMultiplicity;   //!
   TBranch        *b_PuppiAK8Jets_photonEnergyFraction;   //!
   TBranch        *b_PuppiAK8Jets_photonMultiplicity;   //!
   TBranch        *b_PuppiAK8Jets_prunedMass;   //!
   TBranch        *b_PuppiAK8Jets_softDropMass;   //!
   TBranch        *b_PuppiAK8Jets_softDropPt;   //!
   TBranch        *b_PuppiAK8Jets_trimmedMass;   //!
   TBranch        *b_PuppiAK8Jets_filteredMass;   //!
   TBranch        *b_PuppiAK8Jets_tau1;   //!
   TBranch        *b_PuppiAK8Jets_tau2;   //!
   TBranch        *b_PuppiAK8Jets_tau3;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8isLooseJetId;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8isTightJetId;   //!
   TBranch        *b_PuppiAK8Jets_PtCorr;   //!
   TBranch        *b_PuppiAK8Jets_EtaCorr;   //!
   TBranch        *b_PuppiAK8Jets_PhiCorr;   //!
   TBranch        *b_PuppiAK8Jets_ECorr;   //!
   TBranch        *b_PuppiAK8Jets_mass;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8correction;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8correctionUp;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8correctionDown;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8massCorrection;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8massCorrectionUp;   //!
   TBranch        *b_PuppiAK8Jets_PuppiAK8massCorrectionDown;   //!
   TBranch        *b_JetsPuppiNum;   //!
   TBranch        *b_JetsPuppiPt;   //!
   TBranch        *b_JetsPuppiEta;   //!
   TBranch        *b_JetsPuppiPhi;   //!
   TBranch        *b_JetsPuppiE;   //!
   TBranch        *b_JetsPuppiTLorentzVector;   //!
   TBranch        *b_JetsPuppi_bDiscriminatorCSV;   //!
   TBranch        *b_JetsPuppi_bDiscriminatorICSV;   //!
   TBranch        *b_JetsPuppi_chargedEmEnergyFraction;   //!
   TBranch        *b_JetsPuppi_chargedHadronEnergyFraction;   //!
   TBranch        *b_JetsPuppi_chargedHadronMultiplicity;   //!
   TBranch        *b_JetsPuppi_electronMultiplicity;   //!
   TBranch        *b_JetsPuppi_jetArea;   //!
   TBranch        *b_JetsPuppi_muonEnergyFraction;   //!
   TBranch        *b_JetsPuppi_muonMultiplicity;   //!
   TBranch        *b_JetsPuppi_neutralEmEnergyFraction;   //!
   TBranch        *b_JetsPuppi_neutralHadronMultiplicity;   //!
   TBranch        *b_JetsPuppi_photonEnergyFraction;   //!
   TBranch        *b_JetsPuppi_photonMultiplicity;   //!
   TBranch        *b_JetsPuppi_isLooseJetId;   //!
   TBranch        *b_JetsPuppi_isTightJetId;   //!
   TBranch        *b_JetsPuppi_PtCorr;   //!
   TBranch        *b_JetsPuppi_EtaCorr;   //!
   TBranch        *b_JetsPuppi_PhiCorr;   //!
   TBranch        *b_JetsPuppi_ECorr;   //!
   TBranch        *b_JetsPuppi_AK4correction;   //!
   TBranch        *b_JetsPuppi_AK4correctionUp;   //!
   TBranch        *b_JetsPuppi_AK4correctionDown;   //!
   TBranch        *b_ElectronsNum;   //!
   TBranch        *b_ElectronsPt;   //!
   TBranch        *b_ElectronsEta;   //!
   TBranch        *b_ElectronsPhi;   //!
   TBranch        *b_ElectronsE;   //!
   TBranch        *b_ElectronsTLorentzVector;   //!
   TBranch        *b_Electrons_charge;   //!
   TBranch        *b_Electrons_isHEEP;   //!
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
   TBranch        *b_Electrons_isMedium;   //!
   TBranch        *b_Electrons_isTight;   //!
   TBranch        *b_Electrons_SCEnergy;   //!
   TBranch        *b_Electrons_deltaEtaSCTracker;   //!
   TBranch        *b_Electrons_deltaPhiSCTracker;   //!
   TBranch        *b_Electrons_sigmaIetaIeta;   //!
   TBranch        *b_Electrons_sigmaIphiIphi;   //!
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
   TBranch        *b_Muons_isMedium;   //!
   TBranch        *b_Muons_isTight;   //!
   TBranch        *b_Muons_isPFMuon;   //!
   TBranch        *b_TriggerProducerTriggerPrescales;   //!
   TBranch        *b_TriggerProducerTriggerPass;   //!
   TBranch        *b_TriggerProducerTriggerNames;   //!


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
