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
  float AK10_jet_mass_pr;
  float AK10_jet_mass_so;
  float AK12_jet_mass_pr;
  float AK12_jet_mass_so;
  float PuppiAK8_jet_mass_pr;
  float PuppiAK8_jet_mass_so;
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
  float AK8_pruned_mass_gen;
  float AK8_softdrop_mass_gen;
  float AK10_pruned_mass_gen;
  float AK10_softdrop_mass_gen;
  float AK12_pruned_mass_gen;
  float AK12_softdrop_mass_gen;

  inputTree->SetBranchAddress("AK8_pruned_mass_gen", &AK8_pruned_mass_gen);
  inputTree->SetBranchAddress("AK8_softdrop_mass_gen", &AK8_softdrop_mass_gen);
  inputTree->SetBranchAddress("AK10_pruned_mass_gen", &AK10_pruned_mass_gen);
  inputTree->SetBranchAddress("AK10_softdrop_mass_gen", &AK10_softdrop_mass_gen);
  inputTree->SetBranchAddress("AK12_pruned_mass_gen", &AK12_pruned_mass_gen);
  inputTree->SetBranchAddress("AK12_softdrop_mass_gen", &AK12_softdrop_mass_gen);
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
  inputTree->SetBranchAddress("AK10_jet_mass_pr", &AK10_jet_mass_pr);
  inputTree->SetBranchAddress("AK10_jet_mass_so", &AK10_jet_mass_so);
  inputTree->SetBranchAddress("AK12_jet_mass_pr", &AK12_jet_mass_pr);
  inputTree->SetBranchAddress("AK12_jet_mass_so", &AK12_jet_mass_so);
  inputTree->SetBranchAddress("PuppiAK8_jet_mass_pr", &PuppiAK8_jet_mass_pr);
  inputTree->SetBranchAddress("PuppiAK8_jet_mass_so", &PuppiAK8_jet_mass_so);
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

  TH1F *h_AK8_pruned_mass = new TH1F("h_AK8_pruned_mass","h_AK8_pruned_mass",200,-30,30);
  TH1F *h_AK10_pruned_mass = new TH1F("h_AK10_pruned_mass","h_AK10_pruned_mass",200,-30,30);
  TH1F *h_AK12_pruned_mass = new TH1F("h_AK12_pruned_mass","h_AK12_pruned_mass",200,-30,30);
  TH1F *h_PuppiAK8_pruned_mass = new TH1F("h_PuppiAK8_pruned_mass","h_PuppiAK8_pruned_mass",200,-30,30);
  TH1F *h_AK8_softdrop_mass = new TH1F("h_AK8_softdrop_mass","h_AK8_softdrop_mass",200,-30,30);
  TH1F *h_AK10_softdrop_mass = new TH1F("h_AK10_softdrop_mass","h_AK10_softdrop_mass",200,-30,30);
  TH1F *h_AK12_softdrop_mass = new TH1F("h_AK12_softdrop_mass","h_AK12_softdrop_mass",200,-30,30);
  TH1F *h_PuppiAK8_softdrop_mass = new TH1F("h_PuppiAK8_softdrop_mass","h_PuppiAK8_softdrop_mass",200,-30,30);

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("hadronicWstudies_")+leptonName+outputFile+std::string(".root")).c_str(),"recreate");
  outROOT->cd();
  //  outTree->SetDirectory(0);
  
  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<inputTree->GetEntries();jentry++) {

    inputTree->GetEntry(jentry);

    if(jentry % 1000 == 0)    
      cout << "read entry: " << jentry << endl;

    h_AK8_pruned_mass->Fill(jet_mass_pr-AK8_pruned_mass_gen);
    h_AK10_pruned_mass->Fill(AK10_jet_mass_pr-AK10_pruned_mass_gen);
    h_AK12_pruned_mass->Fill(AK12_jet_mass_pr-AK12_pruned_mass_gen);
    h_PuppiAK8_pruned_mass->Fill(PuppiAK8_jet_mass_pr-AK8_pruned_mass_gen);
    h_AK8_softdrop_mass->Fill(jet_mass_pr-AK8_softdrop_mass_gen);
    h_AK10_softdrop_mass->Fill(AK10_jet_mass_pr-AK10_softdrop_mass_gen);
    h_AK12_softdrop_mass->Fill(AK12_jet_mass_pr-AK12_softdrop_mass_gen);
    h_PuppiAK8_softdrop_mass->Fill(PuppiAK8_jet_mass_pr-AK8_softdrop_mass_gen);

  }

  h_AK8_pruned_mass->Write();
  h_AK10_pruned_mass->Write();
  h_AK12_pruned_mass->Write();
  h_PuppiAK8_pruned_mass->Write();
  h_AK8_softdrop_mass->Write();
  h_AK10_softdrop_mass->Write();
  h_AK12_softdrop_mass->Write();
  h_PuppiAK8_softdrop_mass->Write();

  //--------close everything-------------
  //  outTree->Write();
  outROOT->Close();

  return(0);
}
