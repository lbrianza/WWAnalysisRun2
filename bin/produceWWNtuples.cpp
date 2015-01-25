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

#include "../interface/initInputTree.h"
#include "../interface/setOutputTree.h"
#include "../interface/METzCalculator.h"

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
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }
    
  //--------input tree-----------
  //  TChain* chain = new TChain("TreeMaker2/PreSelection");
  TChain* chain = new TChain(inputTreeName.c_str());
  InitTree(chain);
  chain->Add((inputFolder+inputFile).c_str());

  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);
  SetOutTree(outTree);

  //---------start loop on events------------
  for(int iEntry=0; iEntry<chain->GetEntries(); iEntry++){
    if(iEntry % 1000 == 0)      cout << "read entry: " << iEntry << endl;

    //get entry
    chain->GetEntry(iEntry);
    init(); //initialize all variables

    //require exactly one lepton
    if ( strcmp(leptonName.c_str(),"el")==0 && selectedIDIsoElectronsNum!=1) continue; 
    if ( strcmp(leptonName.c_str(),"mu")==0 && selectedIDIsoMuonsNum!=1) continue;      
    
    if (AK8JetsNum < 1) continue; //at least one jet

    if (AK8JetsPt[0] < 150) continue; 
    if (METPt < 50) continue;

    //lepton Pt selection
    if ( strcmp(leptonName.c_str(),"el")==0 && selectedIDIsoElectronsPt[0]<35) continue; 
    if ( strcmp(leptonName.c_str(),"mu")==0 && selectedIDIsoMuonsPt[0]<30) continue; 
    
    //save event variables
    event_runNo   = RunNum;
    event = EvtNum;
    njets = NJets;
    nPV  = NVtx;
    
    /////////////////LEPTON
    if (strcmp(leptonName.c_str(),"el")==0) {
	l_pt  = selectedIDIsoElectronsPt[0];
	l_eta = selectedIDIsoElectronsEta[0];
	l_phi = selectedIDIsoElectronsPhi[0];	
	l_e = selectedIDIsoElectronsE[0];	
      }
    else if (strcmp(leptonName.c_str(),"mu")==0) {
	l_pt  = selectedIDIsoMuonsPt[0];
	l_eta = selectedIDIsoMuonsEta[0];
	l_phi = selectedIDIsoMuonsPhi[0];
	l_e = selectedIDIsoMuonsE[0];
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

    W_mu.SetPtEtaPhiE(l_pt,l_eta,l_phi,l_e);
    W_Met.SetPxPyPzE(METPt * TMath::Cos(METPhi), METPt * TMath::Sin(METPhi), 0., sqrt(METPt*METPt));

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
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(METPhi),
			      nu_pt1 * TMath::Sin(METPhi), pz1_type0, sqrt(nu_pt1*nu_pt1 + pz1_type0*pz1_type0) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(METPhi),
			      nu_pt2 * TMath::Sin(METPhi), pz1_type0, sqrt(nu_pt2*nu_pt2 + pz1_type0*pz1_type0) );

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
      W_neutrino_1.SetPxPyPzE(nu_pt1 * TMath::Cos(METPhi),
			      nu_pt1 * TMath::Sin(METPhi), pz1_type2, sqrt(nu_pt1*nu_pt1 + pz1_type2*pz1_type2) );
      TLorentzVector W_neutrino_2;
      W_neutrino_2.SetPxPyPzE(nu_pt2 * TMath::Cos(METPhi),
			      nu_pt2 * TMath::Sin(METPhi), pz1_type2, sqrt(nu_pt2*nu_pt2 + pz1_type2*pz1_type2) );

      if ( fabs((W_mu+W_neutrino_1).M()-Wmass) < fabs((W_mu+W_neutrino_2).M()-Wmass) ) W_neutrino_type2 = W_neutrino_1;
      else W_neutrino_type2 = W_neutrino_2;
    }

    //    W_mass_type2 = (W_mu+W_neutrino_type2).M();
    //    W_pz_type2 = (W_mu+W_neutrino_type2).Pz();
    //    W_nu1_pz_type2 = pz1_type2;
    //    W_nu2_pz_type2 = pz2_type2;

    pfMET   = sqrt(METPt*METPt);
    pfMET_Phi = METPhi;
    nu_pz_type0 = pz1_type0;
    nu_pz_type2 = pz1_type2;


    /////////////////LEPTONIC W

    TLorentzVector *W = new TLorentzVector();
    TLorentzVector *LEP = new TLorentzVector();
    TLorentzVector *NU0  = new TLorentzVector();
    TLorentzVector *NU2  = new TLorentzVector();
    
    LEP->SetPtEtaPhiE(l_pt,l_eta,l_phi,l_e);
    NU0->SetPxPyPzE(METPt*TMath::Cos(METPhi),METPt*TMath::Sin(METPhi),nu_pz_type0,pfMET);
    NU2->SetPxPyPzE(METPt*TMath::Cos(METPhi),METPt*TMath::Sin(METPhi),nu_pz_type2,pfMET);
    *W = *LEP + *NU0;
    
    v_pt = W->Pt();
    v_eta = W->Eta();
    v_phi = W->Phi();
    v_mt = TMath::Sqrt(2*LEP->Et()*NU0->Et()*(1-TMath::Cos(LEP->DeltaPhi(*NU0))));
    //    W_mt = W->Mt();

    //////////////////ANGULAR VARIABLES

    TLorentzVector *JET = new TLorentzVector();
    JET->SetPtEtaPhiE(AK8JetsPt[0],AK8JetsEta[0],AK8JetsPhi[0],AK8JetsE[0]);
    deltaR_lak8jet = JET->DeltaR(*LEP);
    deltaphi_METak8jet = JET->DeltaPhi(*NU0);
    deltaphi_Vak8jet = JET->DeltaPhi(*W);

    //FOUR-BODY INVARIANT MASS
    mass_lvj_type0 = (*LEP + *NU0 + *JET).M();
    mass_lvj_type2 = (*LEP + *NU2 + *JET).M();

    //delete all the TLorentzVector before a new selection
    delete W;
    delete LEP;
    delete NU0;
    delete NU2;
    delete JET;

    if (v_pt < 150) continue;
    if (deltaR_lak8jet < (TMath::Pi()/2.0))   continue;

    ///////////JETS
    float tempPt=0.;
    for (unsigned int i=0; i<AK8JetsNum; i++)
      {
	if (AK8JetsPt[i]<30 || AK8JetsEta[i]>4.7)  continue;
	if (AK8JetsPt[i]<=tempPt) continue; //to save the jet with largest pt
	ungroomed_jet_pt  = AK8JetsPt[i];
	ungroomed_jet_eta = AK8JetsEta[i];
	ungroomed_jet_phi = AK8JetsPhi[i];
	ungroomed_jet_e   = AK8JetsE[i];
	jet_mass_pr   = AK8Jets_prunedMass[i];
	jet_mass_tr   = AK8Jets_trimmedMass[i];
	jet_mass_fi   = AK8Jets_filteredMass[i];
	jet_tau2tau1   = AK8Jets_tau2[i]/AK8Jets_tau1[i];
	tempPt = ungroomed_jet_pt;
      }

    /////////VBF PART
    TLorentzVector *HADW = new TLorentzVector();
    HADW->SetPtEtaPhiE(ungroomed_jet_pt,ungroomed_jet_eta,ungroomed_jet_phi,ungroomed_jet_e); //AK8 fat jet (hadronic W)
    TLorentzVector *AK4 = new TLorentzVector();
    std::vector<int> indexGoodJets;

    for (unsigned int i=0; i<JetsNum; i++) //loop on AK4 jet
      {
	if (JetsPt[i]<30 || JetsEta[i]>4.7)  continue;
	AK4->SetPtEtaPhiE(JetsPt[i],JetsEta[i],JetsPhi[i],JetsE[i]);
	float deltaR = HADW->DeltaR(*AK4);
	if (deltaR<0.8) continue; //the vbf jets must be outside the had W cone
	indexGoodJets.push_back(i); //save index of the "good" vbf jets candidate
      }

    delete HADW;
    delete AK4;

    TLorentzVector *VBF1 = new TLorentzVector();
    TLorentzVector *VBF2 = new TLorentzVector();
    TLorentzVector *TOT = new TLorentzVector();
    float tempPtMax=0.;
    int nVBF1=-1, nVBF2=-1; //position of the two vbf jets
    
    for (unsigned int i=0; i<indexGoodJets.size()-1; i++) {
      for (unsigned int ii=i+1; ii<indexGoodJets.size(); ii++) {
	VBF1->SetPtEtaPhiE(JetsPt[indexGoodJets.at(i)],JetsEta[indexGoodJets.at(i)],JetsPhi[indexGoodJets.at(i)],JetsE[indexGoodJets.at(i)]);
	VBF2->SetPtEtaPhiE(JetsPt[indexGoodJets.at(ii)],JetsEta[indexGoodJets.at(ii)],JetsPhi[indexGoodJets.at(ii)],JetsE[indexGoodJets.at(ii)]);
	*TOT = *VBF1 + *VBF2;
	if (TOT->Pt() < tempPtMax) continue;
	tempPtMax = TOT->Pt(); //take the jet pair with largest Pt
	nVBF1 = indexGoodJets.at(i); //save position of the 1st vbf jet
	nVBF2 = indexGoodJets.at(ii); //save position of the 2nd vbf jet
      }
    }

    if (nVBF1!=-1 && nVBF2!=-1) //save infos for vbf jet pair
      {
	VBF1->SetPtEtaPhiE(JetsPt[nVBF1],JetsEta[nVBF1],JetsPhi[nVBF1],JetsE[nVBF1]);
	VBF2->SetPtEtaPhiE(JetsPt[nVBF2],JetsEta[nVBF2],JetsPhi[nVBF2],JetsE[nVBF2]);
	*TOT = *VBF1 + *VBF2;

	vbf_maxpt_j1_pt = JetsPt[nVBF1];
	vbf_maxpt_j1_eta = JetsEta[nVBF1];
	vbf_maxpt_j1_phi = JetsPhi[nVBF1];
	vbf_maxpt_j1_e = JetsE[nVBF1];
	vbf_maxpt_j1_bDiscriminatorCSV = Jets_bDiscriminator[nVBF1];
	vbf_maxpt_j2_pt = JetsPt[nVBF2];
	vbf_maxpt_j2_eta = JetsEta[nVBF2];
	vbf_maxpt_j2_phi = JetsPhi[nVBF2];
	vbf_maxpt_j2_e = JetsE[nVBF2];
	vbf_maxpt_j2_bDiscriminatorCSV = Jets_bDiscriminator[nVBF2];
	vbf_maxpt_jj_pt = TOT->Pt();
	vbf_maxpt_jj_eta = TOT->Eta();
	vbf_maxpt_jj_phi = TOT->Phi();
	vbf_maxpt_jj_m = TOT->M();	
      }

    delete VBF1;
    delete VBF2;
    delete TOT;

    /////////////////MC Infos
    /*    if (isMC)
      {
	for (int i=0; i<GenBosonNum; i++) {
	  genBosonPdgId[i] = GenBoson_GenBosonPDGId[i];
	  genBosonPt[i]    = GenBosonPt[i];
	  genBosonEta[i]   = GenBosonEta[i];
	  genBosonPhi[i]   = GenBosonPhi[i];
	  genBosonE[i]     = GenBosonE[i];
	}	

	int start=0;
	for (int i=0; i<GenElecNum; i++) {
	  genLeptonPdgId[i] = 11;
	  genLeptonPt[i]    = GenElecPt[i];
	  genLeptonEta[i]   = GenElecEta[i];
	  genLeptonPhi[i]   = GenElecPhi[i];
	  genLeptonE[i]     = GenElecE[i];
	  start++;
	}	

	for (int i=0; i<GenMuNum; i++) {
	  genLeptonPdgId[start+i] = 13;
	  genLeptonPt[start+i]    = GenMuPt[i];
	  genLeptonEta[start+i]   = GenMuEta[i];
	  genLeptonPhi[start+i]   = GenMuPhi[i];
	  genLeptonE[start+i]     = GenMuE[i];
	}	

	for (int i=0; i<GenNuNum; i++) {
	  genLeptonPt[i]    = GenNuPt[i];
	  genLeptonEta[i]   = GenNuEta[i];
	  genLeptonPhi[i]   = GenNuPhi[i];
	  genLeptonE[i]     = GenNuE[i];
	}	

      }
    */
    //fill the tree
    outTree->Fill();
  }

  //--------close everything-------------
  chain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
