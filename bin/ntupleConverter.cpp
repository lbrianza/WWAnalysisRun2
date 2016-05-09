#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <stdint.h>
#include <iomanip>
#include <ctime>
#include <map>
#include <math.h>

#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TChain.h"
#include "TProfile.h"
#include "TMath.h"
#include "TString.h"
#include "TClass.h"
#include "TApplication.h"
#include "TLorentzVector.h"

#include "../interface/setOutputTree.h"
#include "../interface/analysisUtils.h"
#include "../interface/pseudodataNtuple.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  std::string leptonName = argv[3];
  std::string inputTreeName = argv[4];
  std::string inputFile = argv[5];


  std::cout<<"file: "<<(inputFolder+inputFile).c_str()<<std::endl;
  //  TFile *MyFile = new TFile((inputFolder+inputFile).c_str(),"READ");
  TFile *MyFile = TFile::Open((inputFolder+inputFile).c_str());
  pseudodataNtuple *ReducedTree = new pseudodataNtuple (MyFile, inputTreeName.c_str());
  if (ReducedTree->fChain == 0) return (-1);
  ReducedTree->Init();
  //  TTree *chain = (TTree *) MyFile->Get(inputTreeName.c_str());
  //  setInputTree(chain);

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_")+leptonName+std::string("/")+outputFile).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTree *WWTree = new setOutputTree(outTree);

  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    ReducedTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                                                                                                                           

    if(iEntry % 1000 == 0)    
      cout << "read entry: " << iEntry << endl;

    WWTree->initializeVariables(); //initialize all variables
    
    WWTree->issignal = 0;

    WWTree->wSampleWeight = 1.; //xsec/numberOfEntries
    WWTree->totalEventWeight = 1.; //temporary value
    WWTree->eff_and_pu_Weight = 1.; //temporary value
    
    //save event variables
    WWTree->run   = ReducedTree->run;
    WWTree->event = ReducedTree->event;
    WWTree->lumi = ReducedTree->lumi;
   // WWTree->njets = ReducedTree->NJets;
    WWTree->nPV  = ReducedTree->nPV;
    
    WWTree->l_pt  = ReducedTree->l_pt;
    WWTree->l_eta = ReducedTree->l_eta;
    WWTree->l_phi = ReducedTree->l_phi;	

    WWTree->pfMET   = ReducedTree->pfMET;
    WWTree->pfMET_Phi = ReducedTree->pfMETPhi;
    
    WWTree->v_pt = ReducedTree->W_pt;
    WWTree->v_eta = ReducedTree->W_eta;
    WWTree->v_phi = ReducedTree->W_phi;

    WWTree->ungroomed_jet_pt  = ReducedTree->jet_pt;
    WWTree->ungroomed_jet_eta = ReducedTree->jet_eta;
    WWTree->ungroomed_jet_phi = ReducedTree->jet_phi;
    //    WWTree->jet_mass_pr   = ReducedTree->jet_mass_pruned;
    WWTree->jet_mass_pr   = ReducedTree->massVhad;
    WWTree->jet_mass_so   = ReducedTree->jet_mass_softdrop;
    WWTree->jet_pt_so   = ReducedTree->jet_pt_softdrop;
    WWTree->jet_tau2tau1   = ReducedTree->jet_tau2tau1;

    WWTree->deltaR_lak8jet = deltaR(WWTree->l_eta, WWTree->l_phi, WWTree->ungroomed_jet_eta, WWTree->ungroomed_jet_phi);
    WWTree->deltaphi_METak8jet = getDeltaPhi(WWTree->pfMET_Phi,WWTree->ungroomed_jet_phi);
    WWTree->deltaphi_Vak8jet = getDeltaPhi(WWTree->v_phi,WWTree->ungroomed_jet_phi);
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak8jet)>2.0 && fabs(WWTree->deltaphi_Vak8jet)>2.0)
      WWTree->issignal=1;

    WWTree->njets=ReducedTree->njets;
    WWTree->nBTagJet_medium=ReducedTree->nbtag;
    WWTree->mass_lvj_type0 = ReducedTree->m_lvj;
    WWTree->mass_lvj_type2 = ReducedTree->m_lvj;

    //fill the tree
    outTree->Fill();
  }


  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
