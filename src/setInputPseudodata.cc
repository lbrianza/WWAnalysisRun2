#include "../interface/setInputPseudodata.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>


setInputPseudodata::setInputPseudodata(TFile* inputFile, std::string inputTreeName)
{
  if(inputFile == 0) {
    TFile* f = new TFile("/gwteray/users/brianza/WWNtupleRun2/ReducedTree/ReducedSelection_TTbar.root");
    fChain = (TTree*) f -> Get("WJet");
  }
  else fChain = (TTree*) inputFile -> Get(inputTreeName.c_str());
  Init();
}

setInputPseudodata::~setInputPseudodata()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t setInputPseudodata::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t setInputPseudodata::LoadTree(Long64_t entry)
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

void setInputPseudodata::Init()
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers

   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("issignal", &issignal, &b_issignal);
   fChain->SetBranchAddress("wSampleWeight", &wSampleWeight, &b_wSampleWeight);
   fChain->SetBranchAddress("totalEventWeight", &totalEventWeight, &b_totalEventWeight);
   fChain->SetBranchAddress("eff_and_pu_Weight", &eff_and_pu_Weight, &b_eff_and_pu_Weight);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMET_Phi", &pfMET_Phi, &b_pfMET_Phi);
   fChain->SetBranchAddress("nu_pz_type0", &nu_pz_type0, &b_nu_pz_type0);
   fChain->SetBranchAddress("nu_pz_type2", &nu_pz_type2, &b_nu_pz_type2);
   fChain->SetBranchAddress("nu_pz_run2", &nu_pz_run2, &b_nu_pz_run2);
   fChain->SetBranchAddress("nu_pz_run2_oth", &nu_pz_run2_oth, &b_nu_pz_run2_oth);
   fChain->SetBranchAddress("nu_pz_run2_type", &nu_pz_run2_type, &b_nu_pz_run2_type);
   fChain->SetBranchAddress("nu_pz_isre", &nu_pz_isre, &b_nu_pz_isre);
   fChain->SetBranchAddress("l_pt", &l_pt, &b_l_pt);
   fChain->SetBranchAddress("l_eta", &l_eta, &b_l_eta);
   fChain->SetBranchAddress("l_phi", &l_phi, &b_l_phi);
   fChain->SetBranchAddress("l_e", &l_e, &b_l_e);
   fChain->SetBranchAddress("ungroomed_jet_pt", &ungroomed_jet_pt, &b_ungroomed_jet_pt);
   fChain->SetBranchAddress("ungroomed_jet_eta", &ungroomed_jet_eta, &b_ungroomed_jet_eta);
   fChain->SetBranchAddress("ungroomed_jet_phi", &ungroomed_jet_phi, &b_ungroomed_jet_phi);
   fChain->SetBranchAddress("ungroomed_jet_e", &ungroomed_jet_e, &b_ungroomed_jet_e);
   fChain->SetBranchAddress("jet_mass_pr", &jet_mass_pr, &b_jet_mass_pr);
   fChain->SetBranchAddress("jet_mass_so", &jet_mass_so, &b_jet_mass_so);
   fChain->SetBranchAddress("jet_pt_so", &jet_pt_so, &b_jet_pt_so);
   fChain->SetBranchAddress("jet_mass_tr", &jet_mass_tr, &b_jet_mass_tr);
   fChain->SetBranchAddress("jet_mass_fi", &jet_mass_fi, &b_jet_mass_fi);
   fChain->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1, &b_jet_tau2tau1);
   fChain->SetBranchAddress("W_pt_gen", &W_pt_gen, &b_W_pt_gen);
   fChain->SetBranchAddress("W_pz_gen", &W_pz_gen, &b_W_pz_gen);
   fChain->SetBranchAddress("W_rap_gen", &W_rap_gen, &b_W_rap_gen);
   fChain->SetBranchAddress("genGravMass", &genGravMass, &b_genGravMass);
   fChain->SetBranchAddress("nu_pz_gen", &nu_pz_gen, &b_nu_pz_gen);
   fChain->SetBranchAddress("nu_pt_gen", &nu_pt_gen, &b_nu_pt_gen);
   fChain->SetBranchAddress("nu_phi_gen", &nu_phi_gen, &b_nu_phi_gen);
   fChain->SetBranchAddress("nu_eta_gen", &nu_eta_gen, &b_nu_eta_gen);
   fChain->SetBranchAddress("deltaR_lak8jet", &deltaR_lak8jet, &b_deltaR_lak8jet);
   fChain->SetBranchAddress("deltaphi_METak8jet", &deltaphi_METak8jet, &b_deltaphi_METak8jet);
   fChain->SetBranchAddress("deltaphi_Vak8jet", &deltaphi_Vak8jet, &b_deltaphi_Vak8jet);
   fChain->SetBranchAddress("v_pt", &v_pt, &b_v_pt);
   fChain->SetBranchAddress("v_eta", &v_eta, &b_v_eta);
   fChain->SetBranchAddress("v_phi", &v_phi, &b_v_phi);
   fChain->SetBranchAddress("v_mt", &v_mt, &b_v_mt);
   fChain->SetBranchAddress("mass_lvj_type0", &mass_lvj_type0, &b_mass_lvj_type0);
   fChain->SetBranchAddress("mass_lvj_type2", &mass_lvj_type2, &b_mass_lvj_type2);
   fChain->SetBranchAddress("mass_lvj_run2", &mass_lvj_run2, &b_mass_lvj_run2);
   fChain->SetBranchAddress("nBTagJet_loose", &nBTagJet_loose, &b_nBTagJet_loose);
   fChain->SetBranchAddress("nBTagJet_medium", &nBTagJet_medium, &b_nBTagJet_medium);
   fChain->SetBranchAddress("nBTagJet_tight", &nBTagJet_tight, &b_nBTagJet_tight);
   fChain->SetBranchAddress("mass_leptonic_closerjet", &mass_leptonic_closerjet, &b_mass_leptonic_closerjet);
   fChain->SetBranchAddress("mass_ungroomedjet_closerjet", &mass_ungroomedjet_closerjet, &b_mass_ungroomedjet_closerjet);
   fChain->SetBranchAddress("vbf_maxpt_j1_pt", &vbf_maxpt_j1_pt, &b_vbf_maxpt_j1_pt);
   fChain->SetBranchAddress("vbf_maxpt_j1_eta", &vbf_maxpt_j1_eta, &b_vbf_maxpt_j1_eta);
   fChain->SetBranchAddress("vbf_maxpt_j1_phi", &vbf_maxpt_j1_phi, &b_vbf_maxpt_j1_phi);
   fChain->SetBranchAddress("vbf_maxpt_j1_e", &vbf_maxpt_j1_e, &b_vbf_maxpt_j1_e);
   fChain->SetBranchAddress("vbf_maxpt_j1_bDiscriminatorCSV", &vbf_maxpt_j1_bDiscriminatorCSV, &b_vbf_maxpt_j1_bDiscriminatorCSV);
   fChain->SetBranchAddress("vbf_maxpt_j2_pt", &vbf_maxpt_j2_pt, &b_vbf_maxpt_j2_pt);
   fChain->SetBranchAddress("vbf_maxpt_j2_eta", &vbf_maxpt_j2_eta, &b_vbf_maxpt_j2_eta);
   fChain->SetBranchAddress("vbf_maxpt_j2_phi", &vbf_maxpt_j2_phi, &b_vbf_maxpt_j2_phi);
   fChain->SetBranchAddress("vbf_maxpt_j2_e", &vbf_maxpt_j2_e, &b_vbf_maxpt_j2_e);
   fChain->SetBranchAddress("vbf_maxpt_j2_bDiscriminatorCSV", &vbf_maxpt_j2_bDiscriminatorCSV, &b_vbf_maxpt_j2_bDiscriminatorCSV);
   fChain->SetBranchAddress("vbf_maxpt_jj_pt", &vbf_maxpt_jj_pt, &b_vbf_maxpt_jj_pt);
   fChain->SetBranchAddress("vbf_maxpt_jj_eta", &vbf_maxpt_jj_eta, &b_vbf_maxpt_jj_eta);
   fChain->SetBranchAddress("vbf_maxpt_jj_phi", &vbf_maxpt_jj_phi, &b_vbf_maxpt_jj_phi);
   fChain->SetBranchAddress("vbf_maxpt_jj_m", &vbf_maxpt_jj_m, &b_vbf_maxpt_jj_m);
   Notify();
}

Bool_t setInputPseudodata::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void setInputPseudodata::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t setInputPseudodata::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}


