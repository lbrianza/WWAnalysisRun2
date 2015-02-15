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

#include "../interface/setInputTree.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"
#include "../interface/analysisUtils.h"

using namespace std;

//*******MAIN*******************************************************************

int main (int argc, char** argv)
{ 
  std::string inputFolder = argv[1];
  std::string outputFile = argv[2];
  bool isMC = argv[3];
  std::string leptonName = argv[4];
  std::string inputTreeName = argv[5];
  std::string inputFile = argv[6];
  std::string xSecWeight = argv[7];
  std::string numberOfEntries = argv[8];
  float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }
    
  //--------input tree-----------
  //  TChain* chain = new TChain("TreeMaker2/PreSelection");
  //  TChain* chain = new TChain(inputTreeName.c_str());
  //  InitTree(chain);
  //  chain->Add((inputFolder+inputFile).c_str());

  TFile *MyFile = new TFile((inputFolder+inputFile).c_str(),"READ");
  setInputTree *ReducedTree = new setInputTree (MyFile, inputTreeName.c_str());
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
    int nb = ReducedTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                                                                                                                           

    if(iEntry % 1000 == 0)    
      cout << "read entry: " << iEntry << endl;

    WWTree->initializeVariables(); //initialize all variables

    WWTree->wSampleWeight = weight;
    WWTree->totalEventWeight = 1.; //temporary value
    WWTree->eff_and_pu_Weight = 1.; //temporary value

    //require exactly one lepton
    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->selectedIDIsoElectronsNum!=1) continue; 
    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->selectedIDIsoMuonsNum!=1) continue;      
    
    if (ReducedTree->AK8JetsNum < 1) continue; //at least one jet

    if (ReducedTree->AK8JetsPt[0] < 150) continue; 
    if (ReducedTree->METPt < 50) continue;

    //lepton Pt selection
    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->selectedIDIsoElectronsPt[0]<30) continue; 
    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->selectedIDIsoMuonsPt[0]<30) continue; 
    
    //save event variables
    WWTree->event_runNo   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->njets = ReducedTree->NJets;
    WWTree->nPV  = ReducedTree->NVtx;
    
    /////////////////LEPTON
    if (strcmp(leptonName.c_str(),"el")==0) {
	WWTree->l_pt  = ReducedTree->selectedIDIsoElectronsPt[0];
	WWTree->l_eta = ReducedTree->selectedIDIsoElectronsEta[0];
	WWTree->l_phi = ReducedTree->selectedIDIsoElectronsPhi[0];	
	WWTree->l_e = ReducedTree->selectedIDIsoElectronsE[0];	
      }
    else if (strcmp(leptonName.c_str(),"mu")==0) {
	WWTree->l_pt  = ReducedTree->selectedIDIsoMuonsPt[0];
	WWTree->l_eta = ReducedTree->selectedIDIsoMuonsEta[0];
	WWTree->l_phi = ReducedTree->selectedIDIsoMuonsPhi[0];
	WWTree->l_e = ReducedTree->selectedIDIsoMuonsE[0];
      }

    //////////////MET

    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_mu, W_Met;

    W_mu.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    W_Met.SetPxPyPzE(ReducedTree->METPt * TMath::Cos(ReducedTree->METPhi), ReducedTree->METPt * TMath::Sin(ReducedTree->METPhi), 0., sqrt(ReducedTree->METPt*ReducedTree->METPt));

    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }


    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    double pz2_type0 = NeutrinoPz_type0.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type0_met; 
    W_neutrino_type0_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    //    W_mass_type0_met = (W_neutrino_type0_met+W_mu).M();
    //    W_pz_type0_met = (W_neutrino_type0_met+W_mu).Pz();
    //    W_nu1_pz_type0_met = pz1_type0;
    //    W_nu2_pz_type0_met = pz2_type0;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type0; 
    W_neutrino_type0.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type0,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type0*pz1_type0));

    if (NeutrinoPz_type0.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type0.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type0.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type0 = W_neutrino_1;
      else W_neutrino_type0 = W_neutrino_2;
    }

    //    W_mass_type0 = (W_mu+W_neutrino_type0).M();
    //    W_pz_type0 = (W_mu+W_neutrino_type0).Pz();
    //    W_nu1_pz_type0 = pz1_type0;
    //    W_nu2_pz_type0 = pz2_type0;

    // type2 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type2;
    NeutrinoPz_type2.SetMET(W_Met);
    NeutrinoPz_type2.SetLepton(W_mu);
    NeutrinoPz_type2.SetLeptonType(leptonName.c_str());

    double pz1_type2 = NeutrinoPz_type2.Calculate(2); // Default one -> according to type2
    double pz2_type2 = NeutrinoPz_type2.getOther(); // Default one

    // don't touch the neutrino pT
    TLorentzVector W_neutrino_type2_met; 
    W_neutrino_type2_met.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));
    //    W_mass_type2_met = (W_neutrino_type2_met+W_mu).M();
    //    W_pz_type2_met = (W_neutrino_type2_met+W_mu).Pz();
    //    W_nu1_pz_type2_met = pz1_type2;
    //    W_nu2_pz_type2_met = pz2_type2;

    // chenge the neutrino pT in case of complex solution in order to make it real
    TLorentzVector W_neutrino_type2; 
    W_neutrino_type2.SetPxPyPzE(W_Met.Px(),W_Met.Py(),pz1_type2,sqrt(W_Met.Pt()*W_Met.Pt()+pz1_type2*pz1_type2));

    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = NeutrinoPz_type2.getPtneutrino(1);
      double nu_pt2 = NeutrinoPz_type2.getPtneutrino(2);
      TLorentzVector W_neutrino_1;
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt1 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(ReducedTree->METPhi),
			      nu_pt2 * TMath::Sin(ReducedTree->METPhi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }

    //    W_mass_type2 = (W_mu+W_neutrino_type2).M();
    //    W_pz_type2 = (W_mu+W_neutrino_type2).Pz();
    //    W_nu1_pz_type2 = pz1_type2;
    //    W_nu2_pz_type2 = pz2_type2;

    WWTree->pfMET   = sqrt(ReducedTree->METPt*ReducedTree->METPt);
    WWTree->pfMET_Phi = ReducedTree->METPhi;
    WWTree->nu_pz_type0 = pz1_type0;
    WWTree->nu_pz_type2 = pz1_type2;


    /////////////////LEPTONIC W

    TLorentzVector *W = new TLorentzVector();
    TLorentzVector *LEP = new TLorentzVector();
    TLorentzVector *NU0  = new TLorentzVector();
    TLorentzVector *NU2  = new TLorentzVector();
    
    LEP->SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,WWTree->l_e);
    NU0->SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type0,WWTree->pfMET);
    NU2->SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),WWTree->nu_pz_type2,WWTree->pfMET);
    *W = *LEP + *NU0;
    
    WWTree->v_pt = W->Pt();
    WWTree->v_eta = W->Eta();
    WWTree->v_phi = W->Phi();
    WWTree->v_mt = TMath::Sqrt(2*LEP->Et()*NU0->Et()*(1-TMath::Cos(LEP->DeltaPhi(*NU0))));
    //    W_mt = W->Mt();

    //////////////////ANGULAR VARIABLES

    TLorentzVector *JET = new TLorentzVector();
    JET->SetPtEtaPhiE(ReducedTree->AK8JetsPt[0],ReducedTree->AK8JetsEta[0],ReducedTree->AK8JetsPhi[0],ReducedTree->AK8JetsE[0]);
    WWTree->deltaR_lak8jet = JET->DeltaR(*LEP);
    WWTree->deltaphi_METak8jet = JET->DeltaPhi(*NU0);
    WWTree->deltaphi_Vak8jet = JET->DeltaPhi(*W);
    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && WWTree->deltaphi_METak8jet>2.0 && WWTree->deltaphi_Vak8jet>2.0)
      WWTree->issignal=1;

    //FOUR-BODY INVARIANT MASS
    WWTree->mass_lvj_type0 = (*LEP + *NU0 + *JET).M();
    WWTree->mass_lvj_type2 = (*LEP + *NU2 + *JET).M();

    //delete all the TLorentzVector before a new selection
    delete W;
    delete LEP;
    delete NU0;
    delete NU2;
    delete JET;

    if (WWTree->v_pt < 150) continue;
    if (WWTree->deltaR_lak8jet < (TMath::Pi()/2.0))   continue;

    ///////////JETS
    float tempPt=0.;
    for (unsigned int i=0; i<ReducedTree->AK8JetsNum; i++)
      {
	if (ReducedTree->AK8JetsPt[i]<30 || ReducedTree->AK8JetsEta[i]>4.7)  continue;
	if (ReducedTree->AK8JetsPt[i]<=tempPt) continue; //to save the jet with largest pt
	WWTree->ungroomed_jet_pt  = ReducedTree->AK8JetsPt[i];
	WWTree->ungroomed_jet_eta = ReducedTree->AK8JetsEta[i];
	WWTree->ungroomed_jet_phi = ReducedTree->AK8JetsPhi[i];
	WWTree->ungroomed_jet_e   = ReducedTree->AK8JetsE[i];
	WWTree->jet_mass_pr   = ReducedTree->AK8Jets_prunedMass[i];
	WWTree->jet_mass_tr   = ReducedTree->AK8Jets_trimmedMass[i];
	WWTree->jet_mass_fi   = ReducedTree->AK8Jets_filteredMass[i];
	WWTree->jet_tau2tau1   = ReducedTree->AK8Jets_tau2[i]/ReducedTree->AK8Jets_tau1[i];
	tempPt = WWTree->ungroomed_jet_pt;
      }

    /////////VBF and b-tag section
    bool fillVBF = true;

    TLorentzVector *HADW = new TLorentzVector();
    HADW->SetPtEtaPhiE(WWTree->ungroomed_jet_pt,WWTree->ungroomed_jet_eta,WWTree->ungroomed_jet_phi,WWTree->ungroomed_jet_e); //AK8 fat jet (hadronic W)
    TLorentzVector *AK4 = new TLorentzVector();
    std::vector<int> indexGoodJets;
    indexGoodJets.clear();
    if (indexGoodJets.size()!=0)  fillVBF=false;

    WWTree->nBTagJet_loose=0;
    WWTree->nBTagJet_medium=0;
    WWTree->nBTagJet_tight=0;
   
    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {
	//fill B-Tag info
	if (ReducedTree->Jets_bDiscriminator[i]>0.244)   WWTree->nBTagJet_loose++;
	if (ReducedTree->Jets_bDiscriminator[i]>0.679)   WWTree->nBTagJet_medium++;
	if (ReducedTree->Jets_bDiscriminator[i]>0.898)   WWTree->nBTagJet_tight++;

	if (ReducedTree->JetsPt[i]<30 || ReducedTree->JetsEta[i]>4.7)  continue;
	AK4->SetPtEtaPhiE(ReducedTree->JetsPt[i],ReducedTree->JetsEta[i],ReducedTree->JetsPhi[i],ReducedTree->JetsE[i]);
	float deltaR = HADW->DeltaR(*AK4);
	if (deltaR<0.8) continue; //the vbf jets must be outside the had W cone
	indexGoodJets.push_back(i); //save index of the "good" vbf jets candidate
      }
    if (indexGoodJets.size()<2)  fillVBF=false; //check if at least 2 jets are inside the collection
    delete HADW;
    delete AK4;

    if (fillVBF) 
      {
	TLorentzVector *VBF1 = new TLorentzVector();
	TLorentzVector *VBF2 = new TLorentzVector();
	TLorentzVector *TOT = new TLorentzVector();
	float tempPtMax=0.;
	int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
    
	for (unsigned int i=0; i<indexGoodJets.size()-1; i++) {
	  for (unsigned int ii=i+1; ii<indexGoodJets.size(); ii++) {
	    VBF1->SetPtEtaPhiE(ReducedTree->JetsPt[indexGoodJets.at(i)],ReducedTree->JetsEta[indexGoodJets.at(i)],ReducedTree->JetsPhi[indexGoodJets.at(i)],ReducedTree->JetsE[indexGoodJets.at(i)]);
	    VBF2->SetPtEtaPhiE(ReducedTree->JetsPt[indexGoodJets.at(ii)],ReducedTree->JetsEta[indexGoodJets.at(ii)],ReducedTree->JetsPhi[indexGoodJets.at(ii)],ReducedTree->JetsE[indexGoodJets.at(ii)]);
	    *TOT = *VBF1 + *VBF2;
	    if (TOT->Pt() < tempPtMax) continue;
	    tempPtMax = TOT->Pt(); //take the jet pair with largest Pt
	    nVBF1 = indexGoodJets.at(i); //save position of the 1st vbf jet
	    nVBF2 = indexGoodJets.at(ii); //save position of the 2nd vbf jet
	  }
	}
	
	if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
	  {
	    VBF1->SetPtEtaPhiE(ReducedTree->JetsPt[nVBF1],ReducedTree->JetsEta[nVBF1],ReducedTree->JetsPhi[nVBF1],ReducedTree->JetsE[nVBF1]);
	    VBF2->SetPtEtaPhiE(ReducedTree->JetsPt[nVBF2],ReducedTree->JetsEta[nVBF2],ReducedTree->JetsPhi[nVBF2],ReducedTree->JetsE[nVBF2]);
	    *TOT = *VBF1 + *VBF2;
	    
	    WWTree->vbf_maxpt_j1_pt = ReducedTree->JetsPt[nVBF1];
	    WWTree->vbf_maxpt_j1_eta = ReducedTree->JetsEta[nVBF1];
	    WWTree->vbf_maxpt_j1_phi = ReducedTree->JetsPhi[nVBF1];
	    WWTree->vbf_maxpt_j1_e = ReducedTree->JetsE[nVBF1];
	    WWTree->vbf_maxpt_j1_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminator[nVBF1];
	    WWTree->vbf_maxpt_j2_pt = ReducedTree->JetsPt[nVBF2];
	    WWTree->vbf_maxpt_j2_eta = ReducedTree->JetsEta[nVBF2];
	    WWTree->vbf_maxpt_j2_phi = ReducedTree->JetsPhi[nVBF2];
	    WWTree->vbf_maxpt_j2_e = ReducedTree->JetsE[nVBF2];
	    WWTree->vbf_maxpt_j2_bDiscriminatorCSV = ReducedTree->Jets_bDiscriminator[nVBF2];
	    WWTree->vbf_maxpt_jj_pt = TOT->Pt();
	    WWTree->vbf_maxpt_jj_eta = TOT->Eta();
	    WWTree->vbf_maxpt_jj_phi = TOT->Phi();
	    WWTree->vbf_maxpt_jj_m = TOT->M();	
	  }
	
	delete VBF1;
	delete VBF2;
	delete TOT;
      }

    /////////////////MC Infos
    if (isMC)
      {
	TLorentzVector temp, temp2;
	//	std::cout<<"entry: "<<iEntry<<" "<<GenNuNum<<std::endl;
	double deltaPhiOld=100.;
	for (int i=0; i<ReducedTree->GenBosonNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenBosonPhi[i],WWTree->v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;
	  //	  std::cout<<"bosone: "<<i<<" "<<ReducedTree->GenBosonPhi[i]<<" "<<v_phi<<std::endl;
	  temp.SetPtEtaPhiE(ReducedTree->GenBosonPt[i],ReducedTree->GenBosonEta[i],ReducedTree->GenBosonPhi[i],ReducedTree->GenBosonE[i]);
	  WWTree->W_pt_gen = ReducedTree->GenBosonPt[i];
	  WWTree->W_pz_gen = temp.Pz();
	  deltaPhiOld = deltaPhi;
	}	
	if (ReducedTree->GenBosonNum==2) {
	  temp.SetPtEtaPhiE(ReducedTree->GenBosonPt[0],ReducedTree->GenBosonEta[0],ReducedTree->GenBosonPhi[0],ReducedTree->GenBosonE[0]);
	  temp2.SetPtEtaPhiE(ReducedTree->GenBosonPt[1],ReducedTree->GenBosonEta[1],ReducedTree->GenBosonPhi[1],ReducedTree->GenBosonE[1]);
	  WWTree->genGravMass=(temp+temp2).M();	
	}

	deltaPhiOld=100.;
       	for (int i=0; i<ReducedTree->GenNuNum; i++) {
	  double deltaPhi = getDeltaPhi(ReducedTree->GenNuPhi[i],WWTree->v_phi);
	  if (abs(deltaPhi)>abs(deltaPhiOld))   continue;	  
	  temp.SetPtEtaPhiE(ReducedTree->GenNuPt[i],ReducedTree->GenNuEta[i],ReducedTree->GenNuPhi[i],ReducedTree->GenNuE[i]);
	  WWTree->nu_pz_gen=temp.Pz();	  
	  deltaPhiOld = deltaPhi;
	}		
      }
    
    //fill the tree
    outTree->Fill();
  }

  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
