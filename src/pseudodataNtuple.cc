#include "../interface/pseudodataNtuple.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

pseudodataNtuple::pseudodataNtuple(TFile* inputFile, std::string inputTreeName)
{
  if(inputFile == 0) {
    TFile* f = new TFile("/gwteray/users/brianza/WWNtupleRun2/ReducedTree/ReducedSelection_TTbar.root");
    fChain = (TTree*) f -> Get("WJet");
  }
  else fChain = (TTree*) inputFile -> Get(inputTreeName.c_str());
  Init();
}

/*pseudodataNtuple::pseudodataNtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("pseudoData/mu_PKUTree_pdata.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("pseudoData/mu_PKUTree_pdata.root");
      }
      f->GetObject("PKUTree",tree);

   }
   Init(tree);
   }*/

pseudodataNtuple::~pseudodataNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pseudodataNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pseudodataNtuple::LoadTree(Long64_t entry)
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

void pseudodataNtuple::Init()
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

   fChain->SetBranchAddress("CategoryID", &CategoryID, &b_CategoryID);
   fChain->SetBranchAddress("vTagID", &vTagID, &b_vTagID);
   fChain->SetBranchAddress("massVhad", &massVhad, &b_massVhad);
   fChain->SetBranchAddress("tau21", &tau21, &b_tau21);
   fChain->SetBranchAddress("triggerWeight", &triggerWeight, &b_triggerWeight);
   fChain->SetBranchAddress("lumiWeight", &lumiWeight, &b_lumiWeight);
   fChain->SetBranchAddress("pileupWeight", &pileupWeight, &b_pileupWeight);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("lumi", &lumi, &b_lumi);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("pfMET", &pfMET, &b_pfMET);
   fChain->SetBranchAddress("pfMETPhi", &pfMETPhi, &b_pfMETPhi);
   fChain->SetBranchAddress("weight", &weight, &b_weight);
   fChain->SetBranchAddress("nLooseEle", &nLooseEle, &b_nLooseEle);
   fChain->SetBranchAddress("nLooseMu", &nLooseMu, &b_nLooseMu);
   fChain->SetBranchAddress("l_pt", &l_pt, &b_l_pt);
   fChain->SetBranchAddress("l_eta", &l_eta, &b_l_eta);
   fChain->SetBranchAddress("l_phi", &l_phi, &b_l_phi);
   fChain->SetBranchAddress("jet_pt", &jet_pt, &b_jet_pt);
   fChain->SetBranchAddress("jet_eta", &jet_eta, &b_jet_eta);
   fChain->SetBranchAddress("jet_phi", &jet_phi, &b_jet_phi);
   fChain->SetBranchAddress("jet_mass_pruned", &jet_mass_pruned, &b_jet_mass_pruned);
   fChain->SetBranchAddress("jet_mass_softdrop", &jet_mass_softdrop, &b_jet_mass_softdrop);
   fChain->SetBranchAddress("jet_pt_softdrop", &jet_pt_softdrop, &b_jet_pt_softdrop);
   fChain->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1, &b_jet_tau2tau1);
   fChain->SetBranchAddress("W_pt", &W_pt, &b_W_pt);
   fChain->SetBranchAddress("W_eta", &W_eta, &b_W_eta);
   fChain->SetBranchAddress("W_phi", &W_phi, &b_W_phi);
   fChain->SetBranchAddress("m_lvj", &m_lvj, &b_m_lvj);
   fChain->SetBranchAddress("njets", &njets, &b_njets);
   fChain->SetBranchAddress("nbtag", &nbtag, &b_nbtag);
   fChain->SetBranchAddress("jet2_pt", &jet2_pt, &b_jet2_pt);
   fChain->SetBranchAddress("jet2_btag", &jet2_btag, &b_jet2_btag);
   fChain->SetBranchAddress("jet3_pt", &jet3_pt, &b_jet3_pt);
   fChain->SetBranchAddress("jet3_btag", &jet3_btag, &b_jet3_btag);
   fChain->SetBranchAddress("jetAK8_pt", &jetAK8_pt, &b_jetAK8_pt);
   fChain->SetBranchAddress("jetAK8_mass", &jetAK8_mass, &b_jetAK8_mass);
   fChain->SetBranchAddress("jetAK8_jec", &jetAK8_jec, &b_jetAK8_jec);
   fChain->SetBranchAddress("METraw_et", &METraw_et, &b_METraw_et);
   fChain->SetBranchAddress("METraw_phi", &METraw_phi, &b_METraw_phi);
   fChain->SetBranchAddress("METraw_sumEt", &METraw_sumEt, &b_METraw_sumEt);
   fChain->SetBranchAddress("MET_et", &MET_et, &b_MET_et);
   fChain->SetBranchAddress("MET_phi", &MET_phi, &b_MET_phi);
   fChain->SetBranchAddress("MET_sumEt", &MET_sumEt, &b_MET_sumEt);
   fChain->SetBranchAddress("MET_corrPx", &MET_corrPx, &b_MET_corrPx);
   fChain->SetBranchAddress("MET_corrPy", &MET_corrPy, &b_MET_corrPy);
   Notify();
}

Bool_t pseudodataNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pseudodataNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pseudodataNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
