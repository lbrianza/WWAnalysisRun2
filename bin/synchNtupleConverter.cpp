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

#include "../interface/analysisUtils.h"
#include "../interface/setOutputTreeSynch.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  std::string leptonName = argv[3];
  std::string inputFile = argv[4];

  std::cout<<"file: "<<(inputFile).c_str()<<std::endl;
  TFile * fS = new TFile((inputFolder+leptonName+"/WWTree_"+inputFile+".root").c_str());
  TTree * inputTree = (TTree *)fS->Get("otree");

  int run;
  int event;
  int lumi;
  int njets;
  int nPV;
  int issignal;
  float pfMET;
  float pfMET_Phi;
  float l_pt;
  float l_eta;
  float l_phi;
  float l_e;
  float ungroomed_jet_pt;
  float ungroomed_jet_eta;
  float ungroomed_jet_phi;
  float ungroomed_jet_e;
  float jet_mass_pr;
  float jet_mass_so;
  float jet_tau2tau1;
  float v_pt;
  float v_eta;
  float v_phi;
  float v_mt;
  float mass_lvj_type0;
  int nBTagJet_medium;
  float jet2_pt;
  float jet2_btag;
  float jet3_pt;
  float jet3_btag;

  inputTree->SetBranchAddress("run", &run);
  inputTree->SetBranchAddress("event", &event);
  inputTree->SetBranchAddress("lumi", &lumi);
  inputTree->SetBranchAddress("njets", &njets);
  inputTree->SetBranchAddress("nPV", &nPV);
  inputTree->SetBranchAddress("issignal", &issignal);
  inputTree->SetBranchAddress("pfMET", &pfMET);
  inputTree->SetBranchAddress("pfMET_Phi", &pfMET_Phi);
  inputTree->SetBranchAddress("l_pt", &l_pt);
  inputTree->SetBranchAddress("l_eta", &l_eta);
  inputTree->SetBranchAddress("l_phi", &l_phi);
  inputTree->SetBranchAddress("l_e", &l_e);
  inputTree->SetBranchAddress("ungroomed_jet_pt", &ungroomed_jet_pt);
  inputTree->SetBranchAddress("ungroomed_jet_eta", &ungroomed_jet_eta);
  inputTree->SetBranchAddress("ungroomed_jet_phi", &ungroomed_jet_phi);
  inputTree->SetBranchAddress("ungroomed_jet_e", &ungroomed_jet_e);
  inputTree->SetBranchAddress("jet_mass_pr", &jet_mass_pr);
  inputTree->SetBranchAddress("jet_mass_so", &jet_mass_so);
  inputTree->SetBranchAddress("jet_tau2tau1", &jet_tau2tau1);
  inputTree->SetBranchAddress("v_pt", &v_pt);
  inputTree->SetBranchAddress("v_eta", &v_eta);
  inputTree->SetBranchAddress("v_phi", &v_phi);
  inputTree->SetBranchAddress("v_mt", &v_mt);
  inputTree->SetBranchAddress("mass_lvj_type2", &mass_lvj_type0);
  inputTree->SetBranchAddress("nBTagJet_medium", &nBTagJet_medium);
  inputTree->SetBranchAddress("jet2_pt", &jet2_pt);
  inputTree->SetBranchAddress("jet2_btag", &jet2_btag);
  inputTree->SetBranchAddress("jet3_pt", &jet3_pt);
  inputTree->SetBranchAddress("jet3_btag", &jet3_btag);

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_synch_")+leptonName+std::string("/")+std::string("WWTree_")+outputFile+std::string(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);
  setOutputTreeSynch *WWTree = new setOutputTreeSynch(outTree);
  
  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<inputTree->GetEntries();jentry++) {

    inputTree->GetEntry(jentry);

    WWTree->initializeVariables(); //initialize all variables

    if(jentry % 1000 == 0)    
      cout << "read entry: " << jentry << endl;

    WWTree->issignal = issignal;

    //save event variables
    WWTree->run   = run;
    WWTree->event = event;
    WWTree->lumi = lumi;
    WWTree->njets = njets;
    WWTree->nPV  = nPV;
    
    WWTree->l_pt  = l_pt;
    WWTree->l_eta = l_eta;
    WWTree->l_phi = l_phi;	

    WWTree->pfMET   = pfMET;
    WWTree->pfMETPhi = pfMET_Phi;
    
    WWTree->W_pt = v_pt;
    WWTree->W_eta = v_eta;
    WWTree->W_phi = v_phi;

    WWTree->jet_pt  = ungroomed_jet_pt;
    WWTree->jet_eta = ungroomed_jet_eta;
    WWTree->jet_phi = ungroomed_jet_phi;
    WWTree->jet_mass_pruned   = jet_mass_pr;
    WWTree->jet_mass_softdrop   = jet_mass_so;
    WWTree->jet_tau2tau1   = jet_tau2tau1;

    WWTree->nbtag=nBTagJet_medium;
    WWTree->m_lvj = mass_lvj_type0;

    WWTree->jet2_pt = jet2_pt;
    WWTree->jet2_btag = jet2_btag;
    WWTree->jet3_pt = jet3_pt;
    WWTree->jet3_btag = jet3_btag;

    //fill the tree
    if(strcmp(leptonName.c_str(),"mu")==0 && WWTree->issignal==1 && WWTree->W_pt>200 && WWTree->pfMET>40 && WWTree->l_pt>53 && WWTree->jet_pt>200 && WWTree->nbtag <1 && ((WWTree->jet_mass_pruned > 40 && WWTree->jet_mass_pruned<65) || (WWTree->jet_mass_pruned > 135 && WWTree->jet_mass_pruned<150))) {
      //WWTree->jet_mass_pruned > 40 && WWTree->jet_mass_pruned < 130) {//&& WWTree->jet_tau2tau1<0.5) {
      outTree->Fill();
    }
    if(strcmp(leptonName.c_str(),"el")==0 && WWTree->issignal==1 && WWTree->W_pt>200 && WWTree->pfMET>80 && WWTree->l_pt>120 && WWTree->jet_pt>200 && WWTree->nbtag <1 && ((WWTree->jet_mass_pruned > 40 && WWTree->jet_mass_pruned<65) || (WWTree->jet_mass_pruned > 135 && WWTree->jet_mass_pruned<150))) {
      //WWTree->jet_mass_pruned > 40 && WWTree->jet_mass_pruned < 130) {//&& WWTree->jet_tau2tau1<0.5) {
      outTree->Fill();
    }
  }

  //--------close everything-------------
  outTree->Write();
  outROOT->Close();

  return(0);
}
