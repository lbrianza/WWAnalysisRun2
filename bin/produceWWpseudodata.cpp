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
#include "TRandom3.h"

#include "../interface/setInputPseudodata.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/METzCalculator_Run2.h"
#include "../interface/analysisUtils.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFile = argv[1];

  std::cout<<"file: "<<(inputFile).c_str()<<std::endl;
  //  TFile *MyFile = new TFile((inputFolder+inputFile).c_str(),"READ");
  TFile *MyFile = TFile::Open((inputFile).c_str());
  setInputPseudodata *ReducedTree = new setInputPseudodata (MyFile, "otree");
  if (ReducedTree->fChain == 0) return (-1);
  ReducedTree->Init();

  //---------output tree----------------
  TFile* outROOT = TFile::Open("pseudodata.root","recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTree *WWTree = new setOutputTree(outTree);

  int j=0;
    TRandom3* random1 = new TRandom3;
    TRandom3* random2 = new TRandom3;
      random1->SetSeed(3000);
      random2->SetSeed(3001);

  //---------start loop on events------------
    while (j<1703) {
    //  while (j<100) {
    double rand_no1 = random1->Rndm();
    double rand_no2 = random2->Rndm();

    int iEntry = rand_no1*ReducedTree->fChain->GetEntries();
    iEntry = (int)iEntry;
    ReducedTree->fChain->GetEntry(iEntry);   

    double iW = rand_no2*0.00036;

    if (iW > ReducedTree->wSampleWeight) continue;
    if (ReducedTree->issignal==1 && ReducedTree->v_pt>200 && ReducedTree->pfMET>40 && ReducedTree->l_pt>50 && ReducedTree->ungroomed_jet_pt>200 && ReducedTree->nBTagJet_medium <1 && ReducedTree->jet_mass_pr > 40 && ReducedTree->jet_mass_pr < 130 && ReducedTree->jet_tau2tau1<0.5) {

    std::cout<<"random n: "<<rand_no1<<" "<<rand_no2<<" "<<iEntry<<" "<<iW<<std::endl;

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
    WWTree->pfMET_Phi = ReducedTree->pfMET_Phi;
    
    WWTree->v_pt = ReducedTree->v_pt;
    WWTree->v_eta = ReducedTree->v_eta;
    WWTree->v_phi = ReducedTree->v_phi;

    WWTree->ungroomed_jet_pt  = ReducedTree->ungroomed_jet_pt;
    WWTree->ungroomed_jet_eta = ReducedTree->ungroomed_jet_eta;
    WWTree->ungroomed_jet_phi = ReducedTree->ungroomed_jet_phi;
    WWTree->jet_mass_pr   = ReducedTree->jet_mass_pr;
    WWTree->jet_mass_so   = ReducedTree->jet_mass_so;
    WWTree->jet_pt_so   = ReducedTree->jet_pt_so;
    WWTree->jet_tau2tau1   = ReducedTree->jet_tau2tau1;

    WWTree->deltaR_lak8jet = deltaR(WWTree->l_eta, WWTree->l_phi, WWTree->ungroomed_jet_eta, WWTree->ungroomed_jet_phi);
    WWTree->deltaphi_METak8jet = getDeltaPhi(WWTree->pfMET_Phi,WWTree->ungroomed_jet_phi);
    WWTree->deltaphi_Vak8jet = getDeltaPhi(WWTree->v_phi,WWTree->ungroomed_jet_phi);
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && fabs(WWTree->deltaphi_METak8jet)>2.0 && fabs(WWTree->deltaphi_Vak8jet)>2.0)
      WWTree->issignal=1;

    WWTree->njets=ReducedTree->njets;
    WWTree->nBTagJet_medium=ReducedTree->nBTagJet_medium;
    WWTree->mass_lvj_type0 = ReducedTree->mass_lvj_type0;
    WWTree->mass_lvj_type2 = ReducedTree->mass_lvj_type2;

    if(j % 10 == 0)    
      cout << "fill entry: " << j << endl;
    
    j++;

    outTree->Fill();
    }
  }

  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
