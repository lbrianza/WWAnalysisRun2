#include "../interface/setInputTree.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

setInputTree::setInputTree(std::string inputTreeName)
{
  /*  if(inputFile == 0) {
    TFile* f = new TFile("/gwteray/users/brianza/WWNtupleRun2/ReducedTree/ReducedSelection_TTbar.root");
    fChain = (TTree*) f -> Get("WJet");
  }
  else fChain = (TTree*) inputFile -> Get(inputTreeName.c_str());
  */
  fChain = new TChain(inputTreeName.c_str());
  Init();
}

setInputTree::~setInputTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t setInputTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
  if (!fChain)  return 0; 
   return fChain->GetEntry(entry);
}
Long64_t setInputTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void setInputTree::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
  //   if (!tree) return;
  //   fChain = tree;
  TriggerProducerTriggerPrescales=0;
  TriggerProducerTriggerPass=0;
  TriggerProducerTriggerNames=0;

   fCurrent = -1;
   fChain->SetMakeClass(1);

  fChain->SetBranchAddress("RunNum", &RunNum, &b_RunNum);
   fChain->SetBranchAddress("LumiBlockNum", &LumiBlockNum, &b_LumiBlockNum);
   fChain->SetBranchAddress("EvtNum", &EvtNum, &b_EvtNum);
   fChain->SetBranchAddress("NJets", &NJets, &b_NJets);
   fChain->SetBranchAddress("NVtx", &NVtx, &b_NVtx);
   fChain->SetBranchAddress("npT", &npT, &b_npT);
   fChain->SetBranchAddress("passFilterHBHE", &passFilterHBHE, &b_passFilterHBHE);
   fChain->SetBranchAddress("passFilterHBHEIso", &passFilterHBHEIso, &b_passFilterHBHEIso);
   fChain->SetBranchAddress("passFilterCSCHalo", &passFilterCSCHalo, &b_passFilterCSCHalo);
   fChain->SetBranchAddress("passFilterGoodVtx", &passFilterGoodVtx, &b_passFilterGoodVtx);
   fChain->SetBranchAddress("passFilterEEBadSC", &passFilterEEBadSC, &b_passFilterEEBadSC);
   fChain->SetBranchAddress("passFilterEcalDeadCellTriggerPrimitive", &passFilterEcalDeadCellTriggerPrimitive, &b_passFilterEcalDeadCellTriggerPrimitive);
   fChain->SetBranchAddress("passFilterGlobalTightHalo2016", &passFilterGlobalTightHalo2016, &b_passFilterGlobalTightHalo2016);
   fChain->SetBranchAddress("passFilterBadChCand", &passFilterBadChCand, &b_passFilterBadChCand);
   fChain->SetBranchAddress("passFilterBadPFMuon", &passFilterBadPFMuon, &b_passFilterBadPFMuon);
   fChain->SetBranchAddress("passFilterHBHELooseRerun", &passFilterHBHELooseRerun, &b_passFilterHBHELooseRerun);
   fChain->SetBranchAddress("passFilterHBHETightRerun", &passFilterHBHETightRerun, &b_passFilterHBHETightRerun);
   fChain->SetBranchAddress("passFilterHBHEIsoRerun", &passFilterHBHEIsoRerun, &b_passFilterHBHEIsoRerun);
   fChain->SetBranchAddress("Weight", &Weight, &b_Weight);
   fChain->SetBranchAddress("METPt", &METPt, &b_METPt);
   fChain->SetBranchAddress("METPhi", &METPhi, &b_METPhi);
   fChain->SetBranchAddress("METPtUp", &METPtUp, &b_METPtUp);
   fChain->SetBranchAddress("METPhiUp", &METPhiUp, &b_METPhiUp);
   fChain->SetBranchAddress("METPtDown", &METPtDown, &b_METPtDown);
   fChain->SetBranchAddress("METPhiDown", &METPhiDown, &b_METPhiDown);
   fChain->SetBranchAddress("METPtRaw", &METPtRaw, &b_METPtRaw);
   fChain->SetBranchAddress("METPhiRaw", &METPhiRaw, &b_METPhiRaw);
   fChain->SetBranchAddress("CaloMetPt", &CaloMetPt, &b_CaloMetPt);
   fChain->SetBranchAddress("CaloMetPhi", &CaloMetPhi, &b_CaloMetPhi);
   fChain->SetBranchAddress("genEventWeight", &genEventWeight, &b_genEventWeight);
   fChain->SetBranchAddress("PUWeight", &PUWeight, &b_PUWeight);
   fChain->SetBranchAddress("METpuppiPt", &METpuppiPt, &b_METpuppiPt);
   fChain->SetBranchAddress("METpuppiPhi", &METpuppiPhi, &b_METpuppiPhi);
   fChain->SetBranchAddress("METpuppiPtUp", &METpuppiPtUp, &b_METpuppiPtUp);
   fChain->SetBranchAddress("METpuppiPhiUp", &METpuppiPhiUp, &b_METpuppiPhiUp);
   fChain->SetBranchAddress("METpuppiPtDown", &METpuppiPtDown, &b_METpuppiPtDown);
   fChain->SetBranchAddress("METpuppiPhiDown", &METpuppiPhiDown, &b_METpuppiPhiDown);
   fChain->SetBranchAddress("METpuppiPtRaw", &METpuppiPtRaw, &b_METpuppiPtRaw);
   fChain->SetBranchAddress("METpuppiPhiRaw", &METpuppiPhiRaw, &b_METpuppiPhiRaw);
   fChain->SetBranchAddress("METpuppiCaloMetPt", &METpuppiCaloMetPt, &b_METpuppiCaloMetPt);
   fChain->SetBranchAddress("METpuppiCaloMetPhi", &METpuppiCaloMetPhi, &b_METpuppiCaloMetPhi);
   fChain->SetBranchAddress("GenBosonNum", &GenBosonNum, &b_GenBosonNum);
   fChain->SetBranchAddress("GenBosonPt", GenBosonPt, &b_GenBosonPt);
   fChain->SetBranchAddress("GenBosonEta", GenBosonEta, &b_GenBosonEta);
   fChain->SetBranchAddress("GenBosonPhi", GenBosonPhi, &b_GenBosonPhi);
   fChain->SetBranchAddress("GenBosonE", GenBosonE, &b_GenBosonE);
   fChain->SetBranchAddress("GenBosonTLorentzVector", GenBosonTLorentzVector, &b_GenBosonTLorentzVector);
   fChain->SetBranchAddress("GenBoson_GenBosonPDGId", GenBoson_GenBosonPDGId, &b_GenBoson_GenBosonPDGId);
   fChain->SetBranchAddress("GenBoson_isBosonLeptonic", GenBoson_isBosonLeptonic, &b_GenBoson_isBosonLeptonic);
   fChain->SetBranchAddress("GenMuNum", &GenMuNum, &b_GenMuNum);
   fChain->SetBranchAddress("GenMuPt", GenMuPt, &b_GenMuPt);
   fChain->SetBranchAddress("GenMuEta", GenMuEta, &b_GenMuEta);
   fChain->SetBranchAddress("GenMuPhi", GenMuPhi, &b_GenMuPhi);
   fChain->SetBranchAddress("GenMuE", GenMuE, &b_GenMuE);
   fChain->SetBranchAddress("GenMuTLorentzVector", GenMuTLorentzVector, &b_GenMuTLorentzVector);
   fChain->SetBranchAddress("GenMu_GenMuFromTau", GenMu_GenMuFromTau, &b_GenMu_GenMuFromTau);
   fChain->SetBranchAddress("GenElecNum", &GenElecNum, &b_GenElecNum);
   fChain->SetBranchAddress("GenElecPt", GenElecPt, &b_GenElecPt);
   fChain->SetBranchAddress("GenElecEta", GenElecEta, &b_GenElecEta);
   fChain->SetBranchAddress("GenElecPhi", GenElecPhi, &b_GenElecPhi);
   fChain->SetBranchAddress("GenElecE", GenElecE, &b_GenElecE);
   fChain->SetBranchAddress("GenElecTLorentzVector", GenElecTLorentzVector, &b_GenElecTLorentzVector);
   fChain->SetBranchAddress("GenElec_GenElecFromTau", GenElec_GenElecFromTau, &b_GenElec_GenElecFromTau);
   fChain->SetBranchAddress("GenTauNum", &GenTauNum, &b_GenTauNum);
   fChain->SetBranchAddress("GenTauPt", &GenTauPt, &b_GenTauPt);
   fChain->SetBranchAddress("GenTauEta", &GenTauEta, &b_GenTauEta);
   fChain->SetBranchAddress("GenTauPhi", &GenTauPhi, &b_GenTauPhi);
   fChain->SetBranchAddress("GenTauE", &GenTauE, &b_GenTauE);
   fChain->SetBranchAddress("GenTauTLorentzVector", &GenTauTLorentzVector, &b_GenTauTLorentzVector);
   fChain->SetBranchAddress("GenTau_GenTauHad", &GenTau_GenTauHad, &b_GenTau_GenTauHad);
   fChain->SetBranchAddress("GenNuNum", &GenNuNum, &b_GenNuNum);
   fChain->SetBranchAddress("GenNuPt", GenNuPt, &b_GenNuPt);
   fChain->SetBranchAddress("GenNuEta", GenNuEta, &b_GenNuEta);
   fChain->SetBranchAddress("GenNuPhi", GenNuPhi, &b_GenNuPhi);
   fChain->SetBranchAddress("GenNuE", GenNuE, &b_GenNuE);
   fChain->SetBranchAddress("GenNuTLorentzVector", GenNuTLorentzVector, &b_GenNuTLorentzVector);
   fChain->SetBranchAddress("GenTopNum", &GenTopNum, &b_GenTopNum);
   fChain->SetBranchAddress("GenTopPt", &GenTopPt, &b_GenTopPt);
   fChain->SetBranchAddress("GenTopEta", &GenTopEta, &b_GenTopEta);
   fChain->SetBranchAddress("GenTopPhi", &GenTopPhi, &b_GenTopPhi);
   fChain->SetBranchAddress("GenTopE", &GenTopE, &b_GenTopE);
   fChain->SetBranchAddress("GenTopTLorentzVector", &GenTopTLorentzVector, &b_GenTopTLorentzVector);
   fChain->SetBranchAddress("GenTop_GenTopPDGId", &GenTop_GenTopPDGId, &b_GenTop_GenTopPDGId);
   fChain->SetBranchAddress("GenJetsNum", &GenJetsNum, &b_GenJetsNum);
   fChain->SetBranchAddress("GenJetsPt", GenJetsPt, &b_GenJetsPt);
   fChain->SetBranchAddress("GenJetsEta", GenJetsEta, &b_GenJetsEta);
   fChain->SetBranchAddress("GenJetsPhi", GenJetsPhi, &b_GenJetsPhi);
   fChain->SetBranchAddress("GenJetsE", GenJetsE, &b_GenJetsE);
   fChain->SetBranchAddress("GenJetsTLorentzVector", GenJetsTLorentzVector, &b_GenJetsTLorentzVector);
   fChain->SetBranchAddress("GenJetsAK8Num", &GenJetsAK8Num, &b_GenJetsAK8Num);
   fChain->SetBranchAddress("GenJetsAK8Pt", GenJetsAK8Pt, &b_GenJetsAK8Pt);
   fChain->SetBranchAddress("GenJetsAK8Eta", GenJetsAK8Eta, &b_GenJetsAK8Eta);
   fChain->SetBranchAddress("GenJetsAK8Phi", GenJetsAK8Phi, &b_GenJetsAK8Phi);
   fChain->SetBranchAddress("GenJetsAK8E", GenJetsAK8E, &b_GenJetsAK8E);
   fChain->SetBranchAddress("GenJetsAK8TLorentzVector", GenJetsAK8TLorentzVector, &b_GenJetsAK8TLorentzVector);
   fChain->SetBranchAddress("GenJetsAK8_prunedMass", GenJetsAK8_prunedMass, &b_GenJetsAK8_prunedMass);
   fChain->SetBranchAddress("GenJetsAK8_softdropMass", GenJetsAK8_softdropMass, &b_GenJetsAK8_softdropMass);
   fChain->SetBranchAddress("GenJetsAK8_softdropPt", GenJetsAK8_softdropPt, &b_GenJetsAK8_softdropPt);
   fChain->SetBranchAddress("GenJetsAK8_tau1", GenJetsAK8_tau1, &b_GenJetsAK8_tau1);
   fChain->SetBranchAddress("GenJetsAK8_tau2", GenJetsAK8_tau2, &b_GenJetsAK8_tau2);
   fChain->SetBranchAddress("GenJetsAK8_tau3", GenJetsAK8_tau3, &b_GenJetsAK8_tau3);
   fChain->SetBranchAddress("JetsNum", &JetsNum, &b_JetsNum);
   fChain->SetBranchAddress("JetsPt", JetsPt, &b_JetsPt);
   fChain->SetBranchAddress("JetsEta", JetsEta, &b_JetsEta);
   fChain->SetBranchAddress("JetsPhi", JetsPhi, &b_JetsPhi);
   fChain->SetBranchAddress("JetsE", JetsE, &b_JetsE);
   fChain->SetBranchAddress("JetsTLorentzVector", JetsTLorentzVector, &b_JetsTLorentzVector);
   fChain->SetBranchAddress("Jets_bDiscriminatorCSV", Jets_bDiscriminatorCSV, &b_Jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("Jets_bDiscriminatorICSV", Jets_bDiscriminatorICSV, &b_Jets_bDiscriminatorICSV);
   fChain->SetBranchAddress("Jets_chargedEmEnergyFraction", Jets_chargedEmEnergyFraction, &b_Jets_chargedEmEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronEnergyFraction", Jets_chargedHadronEnergyFraction, &b_Jets_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("Jets_chargedHadronMultiplicity", Jets_chargedHadronMultiplicity, &b_Jets_chargedHadronMultiplicity);
   fChain->SetBranchAddress("Jets_electronMultiplicity", Jets_electronMultiplicity, &b_Jets_electronMultiplicity);
   fChain->SetBranchAddress("Jets_jetArea", Jets_jetArea, &b_Jets_jetArea);
   fChain->SetBranchAddress("Jets_muonEnergyFraction", Jets_muonEnergyFraction, &b_Jets_muonEnergyFraction);
   fChain->SetBranchAddress("Jets_muonMultiplicity", Jets_muonMultiplicity, &b_Jets_muonMultiplicity);
   fChain->SetBranchAddress("Jets_neutralEmEnergyFraction", Jets_neutralEmEnergyFraction, &b_Jets_neutralEmEnergyFraction);
   fChain->SetBranchAddress("Jets_neutralHadronMultiplicity", Jets_neutralHadronMultiplicity, &b_Jets_neutralHadronMultiplicity);
   fChain->SetBranchAddress("Jets_photonEnergyFraction", Jets_photonEnergyFraction, &b_Jets_photonEnergyFraction);
   fChain->SetBranchAddress("Jets_photonMultiplicity", Jets_photonMultiplicity, &b_Jets_photonMultiplicity);
   fChain->SetBranchAddress("Jets_isLooseJetId", Jets_isLooseJetId, &b_Jets_isLooseJetId);
   fChain->SetBranchAddress("Jets_isTightJetId", Jets_isTightJetId, &b_Jets_isTightJetId);
   fChain->SetBranchAddress("Jets_PtCorr", Jets_PtCorr, &b_Jets_PtCorr);
   fChain->SetBranchAddress("Jets_EtaCorr", Jets_EtaCorr, &b_Jets_EtaCorr);
   fChain->SetBranchAddress("Jets_PhiCorr", Jets_PhiCorr, &b_Jets_PhiCorr);
   fChain->SetBranchAddress("Jets_ECorr", Jets_ECorr, &b_Jets_ECorr);
   fChain->SetBranchAddress("Jets_AK4correction", Jets_AK4correction, &b_Jets_AK4correction);
   fChain->SetBranchAddress("Jets_AK4correctionUp", Jets_AK4correctionUp, &b_Jets_AK4correctionUp);
   fChain->SetBranchAddress("Jets_AK4correctionDown", Jets_AK4correctionDown, &b_Jets_AK4correctionDown);
   fChain->SetBranchAddress("AK8JetsNum", &AK8JetsNum, &b_AK8JetsNum);
   fChain->SetBranchAddress("AK8JetsPt", AK8JetsPt, &b_AK8JetsPt);
   fChain->SetBranchAddress("AK8JetsEta", AK8JetsEta, &b_AK8JetsEta);
   fChain->SetBranchAddress("AK8JetsPhi", AK8JetsPhi, &b_AK8JetsPhi);
   fChain->SetBranchAddress("AK8JetsE", AK8JetsE, &b_AK8JetsE);
   fChain->SetBranchAddress("AK8JetsTLorentzVector", AK8JetsTLorentzVector, &b_AK8JetsTLorentzVector);
   fChain->SetBranchAddress("AK8Jets_bDiscriminatorCSV", AK8Jets_bDiscriminatorCSV, &b_AK8Jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("AK8Jets_bDiscriminatorICSV", AK8Jets_bDiscriminatorICSV, &b_AK8Jets_bDiscriminatorICSV);
   fChain->SetBranchAddress("AK8Jets_chargedEmEnergyFraction", AK8Jets_chargedEmEnergyFraction, &b_AK8Jets_chargedEmEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_chargedHadronEnergyFraction", AK8Jets_chargedHadronEnergyFraction, &b_AK8Jets_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_chargedHadronMultiplicity", AK8Jets_chargedHadronMultiplicity, &b_AK8Jets_chargedHadronMultiplicity);
   fChain->SetBranchAddress("AK8Jets_electronMultiplicity", AK8Jets_electronMultiplicity, &b_AK8Jets_electronMultiplicity);
   fChain->SetBranchAddress("AK8Jets_jetArea", AK8Jets_jetArea, &b_AK8Jets_jetArea);
   fChain->SetBranchAddress("AK8Jets_muonEnergyFraction", AK8Jets_muonEnergyFraction, &b_AK8Jets_muonEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_muonMultiplicity", AK8Jets_muonMultiplicity, &b_AK8Jets_muonMultiplicity);
   fChain->SetBranchAddress("AK8Jets_neutralEmEnergyFraction", AK8Jets_neutralEmEnergyFraction, &b_AK8Jets_neutralEmEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_neutralHadronMultiplicity", AK8Jets_neutralHadronMultiplicity, &b_AK8Jets_neutralHadronMultiplicity);
   fChain->SetBranchAddress("AK8Jets_photonEnergyFraction", AK8Jets_photonEnergyFraction, &b_AK8Jets_photonEnergyFraction);
   fChain->SetBranchAddress("AK8Jets_photonMultiplicity", AK8Jets_photonMultiplicity, &b_AK8Jets_photonMultiplicity);
   fChain->SetBranchAddress("AK8Jets_prunedMass", AK8Jets_prunedMass, &b_AK8Jets_prunedMass);
   fChain->SetBranchAddress("AK8Jets_softDropMass", AK8Jets_softDropMass, &b_AK8Jets_softDropMass);
   fChain->SetBranchAddress("AK8Jets_softDropPt", AK8Jets_softDropPt, &b_AK8Jets_softDropPt);
   fChain->SetBranchAddress("AK8Jets_trimmedMass", AK8Jets_trimmedMass, &b_AK8Jets_trimmedMass);
   fChain->SetBranchAddress("AK8Jets_filteredMass", AK8Jets_filteredMass, &b_AK8Jets_filteredMass);
   fChain->SetBranchAddress("AK8Jets_tau1", AK8Jets_tau1, &b_AK8Jets_tau1);
   fChain->SetBranchAddress("AK8Jets_tau2", AK8Jets_tau2, &b_AK8Jets_tau2);
   fChain->SetBranchAddress("AK8Jets_tau3", AK8Jets_tau3, &b_AK8Jets_tau3);
   fChain->SetBranchAddress("AK8Jets_AK8isLooseJetId", AK8Jets_AK8isLooseJetId, &b_AK8Jets_AK8isLooseJetId);
   fChain->SetBranchAddress("AK8Jets_AK8isTightJetId", AK8Jets_AK8isTightJetId, &b_AK8Jets_AK8isTightJetId);
   fChain->SetBranchAddress("AK8Jets_PtCorr", AK8Jets_PtCorr, &b_AK8Jets_PtCorr);
   fChain->SetBranchAddress("AK8Jets_EtaCorr", AK8Jets_EtaCorr, &b_AK8Jets_EtaCorr);
   fChain->SetBranchAddress("AK8Jets_PhiCorr", AK8Jets_PhiCorr, &b_AK8Jets_PhiCorr);
   fChain->SetBranchAddress("AK8Jets_ECorr", AK8Jets_ECorr, &b_AK8Jets_ECorr);
   fChain->SetBranchAddress("AK8Jets_mass", AK8Jets_mass, &b_AK8Jets_mass);
   fChain->SetBranchAddress("AK8Jets_AK8correction", AK8Jets_AK8correction, &b_AK8Jets_AK8correction);
   fChain->SetBranchAddress("AK8Jets_AK8correctionUp", AK8Jets_AK8correctionUp, &b_AK8Jets_AK8correctionUp);
   fChain->SetBranchAddress("AK8Jets_AK8correctionDown", AK8Jets_AK8correctionDown, &b_AK8Jets_AK8correctionDown);
   fChain->SetBranchAddress("AK8Jets_AK8massCorrection", AK8Jets_AK8massCorrection, &b_AK8Jets_AK8massCorrection);
   fChain->SetBranchAddress("AK8Jets_AK8massCorrectionUp", AK8Jets_AK8massCorrectionUp, &b_AK8Jets_AK8massCorrectionUp);
   fChain->SetBranchAddress("AK8Jets_AK8massCorrectionDown", AK8Jets_AK8massCorrectionDown, &b_AK8Jets_AK8massCorrectionDown);
   fChain->SetBranchAddress("PuppiAK8JetsNum", &PuppiAK8JetsNum, &b_PuppiAK8JetsNum);
   fChain->SetBranchAddress("PuppiAK8JetsPt", PuppiAK8JetsPt, &b_PuppiAK8JetsPt);
   fChain->SetBranchAddress("PuppiAK8JetsEta", PuppiAK8JetsEta, &b_PuppiAK8JetsEta);
   fChain->SetBranchAddress("PuppiAK8JetsPhi", PuppiAK8JetsPhi, &b_PuppiAK8JetsPhi);
   fChain->SetBranchAddress("PuppiAK8JetsE", PuppiAK8JetsE, &b_PuppiAK8JetsE);
   fChain->SetBranchAddress("PuppiAK8JetsTLorentzVector", PuppiAK8JetsTLorentzVector, &b_PuppiAK8JetsTLorentzVector);
   fChain->SetBranchAddress("PuppiAK8Jets_bDiscriminatorCSV", PuppiAK8Jets_bDiscriminatorCSV, &b_PuppiAK8Jets_bDiscriminatorCSV);
   fChain->SetBranchAddress("PuppiAK8Jets_bDiscriminatorICSV", PuppiAK8Jets_bDiscriminatorICSV, &b_PuppiAK8Jets_bDiscriminatorICSV);
   fChain->SetBranchAddress("PuppiAK8Jets_chargedEmEnergyFraction", PuppiAK8Jets_chargedEmEnergyFraction, &b_PuppiAK8Jets_chargedEmEnergyFraction);
   fChain->SetBranchAddress("PuppiAK8Jets_chargedHadronEnergyFraction", PuppiAK8Jets_chargedHadronEnergyFraction, &b_PuppiAK8Jets_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("PuppiAK8Jets_chargedHadronMultiplicity", PuppiAK8Jets_chargedHadronMultiplicity, &b_PuppiAK8Jets_chargedHadronMultiplicity);
   fChain->SetBranchAddress("PuppiAK8Jets_electronMultiplicity", PuppiAK8Jets_electronMultiplicity, &b_PuppiAK8Jets_electronMultiplicity);
   fChain->SetBranchAddress("PuppiAK8Jets_jetArea", PuppiAK8Jets_jetArea, &b_PuppiAK8Jets_jetArea);
   fChain->SetBranchAddress("PuppiAK8Jets_muonEnergyFraction", PuppiAK8Jets_muonEnergyFraction, &b_PuppiAK8Jets_muonEnergyFraction);
   fChain->SetBranchAddress("PuppiAK8Jets_muonMultiplicity", PuppiAK8Jets_muonMultiplicity, &b_PuppiAK8Jets_muonMultiplicity);
   fChain->SetBranchAddress("PuppiAK8Jets_neutralEmEnergyFraction", PuppiAK8Jets_neutralEmEnergyFraction, &b_PuppiAK8Jets_neutralEmEnergyFraction);
   fChain->SetBranchAddress("PuppiAK8Jets_neutralHadronMultiplicity", PuppiAK8Jets_neutralHadronMultiplicity, &b_PuppiAK8Jets_neutralHadronMultiplicity);
   fChain->SetBranchAddress("PuppiAK8Jets_photonEnergyFraction", PuppiAK8Jets_photonEnergyFraction, &b_PuppiAK8Jets_photonEnergyFraction);
   fChain->SetBranchAddress("PuppiAK8Jets_photonMultiplicity", PuppiAK8Jets_photonMultiplicity, &b_PuppiAK8Jets_photonMultiplicity);
   fChain->SetBranchAddress("PuppiAK8Jets_prunedMass", PuppiAK8Jets_prunedMass, &b_PuppiAK8Jets_prunedMass);
   fChain->SetBranchAddress("PuppiAK8Jets_softDropMass", PuppiAK8Jets_softDropMass, &b_PuppiAK8Jets_softDropMass);
   fChain->SetBranchAddress("PuppiAK8Jets_softDropPt", PuppiAK8Jets_softDropPt, &b_PuppiAK8Jets_softDropPt);
   fChain->SetBranchAddress("PuppiAK8Jets_trimmedMass", PuppiAK8Jets_trimmedMass, &b_PuppiAK8Jets_trimmedMass);
   fChain->SetBranchAddress("PuppiAK8Jets_filteredMass", PuppiAK8Jets_filteredMass, &b_PuppiAK8Jets_filteredMass);
   fChain->SetBranchAddress("PuppiAK8Jets_tau1", PuppiAK8Jets_tau1, &b_PuppiAK8Jets_tau1);
   fChain->SetBranchAddress("PuppiAK8Jets_tau2", PuppiAK8Jets_tau2, &b_PuppiAK8Jets_tau2);
   fChain->SetBranchAddress("PuppiAK8Jets_tau3", PuppiAK8Jets_tau3, &b_PuppiAK8Jets_tau3);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8isLooseJetId", PuppiAK8Jets_PuppiAK8isLooseJetId, &b_PuppiAK8Jets_PuppiAK8isLooseJetId);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8isTightJetId", PuppiAK8Jets_PuppiAK8isTightJetId, &b_PuppiAK8Jets_PuppiAK8isTightJetId);
   fChain->SetBranchAddress("PuppiAK8Jets_PtCorr", PuppiAK8Jets_PtCorr, &b_PuppiAK8Jets_PtCorr);
   fChain->SetBranchAddress("PuppiAK8Jets_EtaCorr", PuppiAK8Jets_EtaCorr, &b_PuppiAK8Jets_EtaCorr);
   fChain->SetBranchAddress("PuppiAK8Jets_PhiCorr", PuppiAK8Jets_PhiCorr, &b_PuppiAK8Jets_PhiCorr);
   fChain->SetBranchAddress("PuppiAK8Jets_ECorr", PuppiAK8Jets_ECorr, &b_PuppiAK8Jets_ECorr);
   fChain->SetBranchAddress("PuppiAK8Jets_mass", PuppiAK8Jets_mass, &b_PuppiAK8Jets_mass);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8correction", PuppiAK8Jets_PuppiAK8correction, &b_PuppiAK8Jets_PuppiAK8correction);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8correctionUp", PuppiAK8Jets_PuppiAK8correctionUp, &b_PuppiAK8Jets_PuppiAK8correctionUp);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8correctionDown", PuppiAK8Jets_PuppiAK8correctionDown, &b_PuppiAK8Jets_PuppiAK8correctionDown);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8massCorrection", PuppiAK8Jets_PuppiAK8massCorrection, &b_PuppiAK8Jets_PuppiAK8massCorrection);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8massCorrectionUp", PuppiAK8Jets_PuppiAK8massCorrectionUp, &b_PuppiAK8Jets_PuppiAK8massCorrectionUp);
   fChain->SetBranchAddress("PuppiAK8Jets_PuppiAK8massCorrectionDown", PuppiAK8Jets_PuppiAK8massCorrectionDown, &b_PuppiAK8Jets_PuppiAK8massCorrectionDown);
   fChain->SetBranchAddress("JetsPuppiNum", &JetsPuppiNum, &b_JetsPuppiNum);
   fChain->SetBranchAddress("JetsPuppiPt", JetsPuppiPt, &b_JetsPuppiPt);
   fChain->SetBranchAddress("JetsPuppiEta", JetsPuppiEta, &b_JetsPuppiEta);
   fChain->SetBranchAddress("JetsPuppiPhi", JetsPuppiPhi, &b_JetsPuppiPhi);
   fChain->SetBranchAddress("JetsPuppiE", JetsPuppiE, &b_JetsPuppiE);
   fChain->SetBranchAddress("JetsPuppiTLorentzVector", JetsPuppiTLorentzVector, &b_JetsPuppiTLorentzVector);
   fChain->SetBranchAddress("JetsPuppi_bDiscriminatorCSV", JetsPuppi_bDiscriminatorCSV, &b_JetsPuppi_bDiscriminatorCSV);
   fChain->SetBranchAddress("JetsPuppi_bDiscriminatorICSV", JetsPuppi_bDiscriminatorICSV, &b_JetsPuppi_bDiscriminatorICSV);
   fChain->SetBranchAddress("JetsPuppi_chargedEmEnergyFraction", JetsPuppi_chargedEmEnergyFraction, &b_JetsPuppi_chargedEmEnergyFraction);
   fChain->SetBranchAddress("JetsPuppi_chargedHadronEnergyFraction", JetsPuppi_chargedHadronEnergyFraction, &b_JetsPuppi_chargedHadronEnergyFraction);
   fChain->SetBranchAddress("JetsPuppi_chargedHadronMultiplicity", JetsPuppi_chargedHadronMultiplicity, &b_JetsPuppi_chargedHadronMultiplicity);
   fChain->SetBranchAddress("JetsPuppi_electronMultiplicity", JetsPuppi_electronMultiplicity, &b_JetsPuppi_electronMultiplicity);
   fChain->SetBranchAddress("JetsPuppi_jetArea", JetsPuppi_jetArea, &b_JetsPuppi_jetArea);
   fChain->SetBranchAddress("JetsPuppi_muonEnergyFraction", JetsPuppi_muonEnergyFraction, &b_JetsPuppi_muonEnergyFraction);
   fChain->SetBranchAddress("JetsPuppi_muonMultiplicity", JetsPuppi_muonMultiplicity, &b_JetsPuppi_muonMultiplicity);
   fChain->SetBranchAddress("JetsPuppi_neutralEmEnergyFraction", JetsPuppi_neutralEmEnergyFraction, &b_JetsPuppi_neutralEmEnergyFraction);
   fChain->SetBranchAddress("JetsPuppi_neutralHadronMultiplicity", JetsPuppi_neutralHadronMultiplicity, &b_JetsPuppi_neutralHadronMultiplicity);
   fChain->SetBranchAddress("JetsPuppi_photonEnergyFraction", JetsPuppi_photonEnergyFraction, &b_JetsPuppi_photonEnergyFraction);
   fChain->SetBranchAddress("JetsPuppi_photonMultiplicity", JetsPuppi_photonMultiplicity, &b_JetsPuppi_photonMultiplicity);
   fChain->SetBranchAddress("JetsPuppi_isLooseJetId", JetsPuppi_isLooseJetId, &b_JetsPuppi_isLooseJetId);
   fChain->SetBranchAddress("JetsPuppi_isTightJetId", JetsPuppi_isTightJetId, &b_JetsPuppi_isTightJetId);
   fChain->SetBranchAddress("JetsPuppi_PtCorr", JetsPuppi_PtCorr, &b_JetsPuppi_PtCorr);
   fChain->SetBranchAddress("JetsPuppi_EtaCorr", JetsPuppi_EtaCorr, &b_JetsPuppi_EtaCorr);
   fChain->SetBranchAddress("JetsPuppi_PhiCorr", JetsPuppi_PhiCorr, &b_JetsPuppi_PhiCorr);
   fChain->SetBranchAddress("JetsPuppi_ECorr", JetsPuppi_ECorr, &b_JetsPuppi_ECorr);
   fChain->SetBranchAddress("JetsPuppi_AK4correction", JetsPuppi_AK4correction, &b_JetsPuppi_AK4correction);
   fChain->SetBranchAddress("JetsPuppi_AK4correctionUp", JetsPuppi_AK4correctionUp, &b_JetsPuppi_AK4correctionUp);
   fChain->SetBranchAddress("JetsPuppi_AK4correctionDown", JetsPuppi_AK4correctionDown, &b_JetsPuppi_AK4correctionDown);
   fChain->SetBranchAddress("ElectronsNum", &ElectronsNum, &b_ElectronsNum);
   fChain->SetBranchAddress("ElectronsPt", ElectronsPt, &b_ElectronsPt);
   fChain->SetBranchAddress("ElectronsEta", ElectronsEta, &b_ElectronsEta);
   fChain->SetBranchAddress("ElectronsPhi", ElectronsPhi, &b_ElectronsPhi);
   fChain->SetBranchAddress("ElectronsE", ElectronsE, &b_ElectronsE);
   fChain->SetBranchAddress("ElectronsTLorentzVector", ElectronsTLorentzVector, &b_ElectronsTLorentzVector);
   fChain->SetBranchAddress("Electrons_charge", Electrons_charge, &b_Electrons_charge);
   fChain->SetBranchAddress("Electrons_isHEEP", Electrons_isHEEP, &b_Electrons_isHEEP);
   fChain->SetBranchAddress("Electrons_type", Electrons_type, &b_Electrons_type);
   fChain->SetBranchAddress("Electrons_mass", Electrons_mass, &b_Electrons_mass);
   fChain->SetBranchAddress("Electrons_pfDeltaCorrRelIso", Electrons_pfDeltaCorrRelIso, &b_Electrons_pfDeltaCorrRelIso);
   fChain->SetBranchAddress("Electrons_pfRhoCorrRelIso04", Electrons_pfRhoCorrRelIso04, &b_Electrons_pfRhoCorrRelIso04);
   fChain->SetBranchAddress("Electrons_pfRhoCorrRelIso03", Electrons_pfRhoCorrRelIso03, &b_Electrons_pfRhoCorrRelIso03);
   fChain->SetBranchAddress("Electrons_pfRelIso", Electrons_pfRelIso, &b_Electrons_pfRelIso);
   fChain->SetBranchAddress("Electrons_photonIso", Electrons_photonIso, &b_Electrons_photonIso);
   fChain->SetBranchAddress("Electrons_neutralHadIso", Electrons_neutralHadIso, &b_Electrons_neutralHadIso);
   fChain->SetBranchAddress("Electrons_chargedHadIso", Electrons_chargedHadIso, &b_Electrons_chargedHadIso);
   fChain->SetBranchAddress("Electrons_trackIso", Electrons_trackIso, &b_Electrons_trackIso);
   fChain->SetBranchAddress("Electrons_isLoose", Electrons_isLoose, &b_Electrons_isLoose);
   fChain->SetBranchAddress("Electrons_isMedium", Electrons_isMedium, &b_Electrons_isMedium);
   fChain->SetBranchAddress("Electrons_isTight", Electrons_isTight, &b_Electrons_isTight);
   fChain->SetBranchAddress("Electrons_SCEnergy", Electrons_SCEnergy, &b_Electrons_SCEnergy);
   fChain->SetBranchAddress("Electrons_deltaEtaSCTracker", Electrons_deltaEtaSCTracker, &b_Electrons_deltaEtaSCTracker);
   fChain->SetBranchAddress("Electrons_deltaPhiSCTracker", Electrons_deltaPhiSCTracker, &b_Electrons_deltaPhiSCTracker);
   fChain->SetBranchAddress("Electrons_sigmaIetaIeta", Electrons_sigmaIetaIeta, &b_Electrons_sigmaIetaIeta);
   fChain->SetBranchAddress("Electrons_sigmaIphiIphi", Electrons_sigmaIphiIphi, &b_Electrons_sigmaIphiIphi);
   fChain->SetBranchAddress("MuonsNum", &MuonsNum, &b_MuonsNum);
   fChain->SetBranchAddress("MuonsPt", MuonsPt, &b_MuonsPt);
   fChain->SetBranchAddress("MuonsEta", MuonsEta, &b_MuonsEta);
   fChain->SetBranchAddress("MuonsPhi", MuonsPhi, &b_MuonsPhi);
   fChain->SetBranchAddress("MuonsE", MuonsE, &b_MuonsE);
   fChain->SetBranchAddress("MuonsTLorentzVector", MuonsTLorentzVector, &b_MuonsTLorentzVector);
   fChain->SetBranchAddress("Muons_charge", Muons_charge, &b_Muons_charge);
   fChain->SetBranchAddress("Muons_isHighPt", Muons_isHighPt, &b_Muons_isHighPt);
   fChain->SetBranchAddress("Muons_type", Muons_type, &b_Muons_type);
   fChain->SetBranchAddress("Muons_mass", Muons_mass, &b_Muons_mass);
   fChain->SetBranchAddress("Muons_pfDeltaCorrRelIso", Muons_pfDeltaCorrRelIso, &b_Muons_pfDeltaCorrRelIso);
   fChain->SetBranchAddress("Muons_pfRelIso", Muons_pfRelIso, &b_Muons_pfRelIso);
   fChain->SetBranchAddress("Muons_photonIso", Muons_photonIso, &b_Muons_photonIso);
   fChain->SetBranchAddress("Muons_neutralHadIso", Muons_neutralHadIso, &b_Muons_neutralHadIso);
   fChain->SetBranchAddress("Muons_chargedHadIso", Muons_chargedHadIso, &b_Muons_chargedHadIso);
   fChain->SetBranchAddress("Muons_trackIso", Muons_trackIso, &b_Muons_trackIso);
   fChain->SetBranchAddress("Muons_isLoose", Muons_isLoose, &b_Muons_isLoose);
   fChain->SetBranchAddress("Muons_isMedium", Muons_isMedium, &b_Muons_isMedium);
   fChain->SetBranchAddress("Muons_isTight", Muons_isTight, &b_Muons_isTight);
   fChain->SetBranchAddress("Muons_isPFMuon", Muons_isPFMuon, &b_Muons_isPFMuon);
   fChain->SetBranchAddress("TriggerProducerTriggerPrescales", &TriggerProducerTriggerPrescales, &b_TriggerProducerTriggerPrescales);
   fChain->SetBranchAddress("TriggerProducerTriggerPass", &TriggerProducerTriggerPass, &b_TriggerProducerTriggerPass);
   fChain->SetBranchAddress("TriggerProducerTriggerNames", &TriggerProducerTriggerNames, &b_TriggerProducerTriggerNames);
   Notify();
  //   Notify();
}

Bool_t setInputTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void setInputTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t setInputTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
