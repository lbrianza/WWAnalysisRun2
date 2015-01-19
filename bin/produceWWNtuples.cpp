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
  std::string inputFile = argv[1];
  std::string outputFile = argv[2];
  bool isMC = argv[3];

  //--------input tree-----------
  TChain* chain = new TChain("TreeMaker2/PreSelection");
  InitTree(chain);
  chain->Add(inputFile.c_str());

  //---------output tree----------------
  TFile* outROOT = TFile::Open((outputFile+".root").c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("WWTree", "WWTree");
  outTree->SetDirectory(0);
  SetOutTree(outTree);

  //---------start loop on events------------
  for(int iEntry=0; iEntry<chain->GetEntries(); iEntry++){
    if(iEntry % 1000 == 0)      cout << "read entry: " << iEntry << endl;

    //get entry
    chain->GetEntry(iEntry);
    init(); //initialize all variables

    if ( (selectedIDIsoElectronsNum+selectedIDIsoMuonsNum)!=1)  continue;      //require exactly one lepton
    if (JetsNum < 1 || AK8JetsNum < 1) continue; //at least one jet

    if (AK8JetsPt[0] < 150) continue; 
    if (METPt < 50) continue;
    if (selectedIDIsoElectronsPt[0]<35 || selectedIDIsoMuonsPt[0]<30) continue; //lepton pt selection

    //save variables
    run   = RunNum;
    event = EvtNum;
    nJets = NJets;
    nVtx  = NVtx;
    
    /////////////////LEPTON
    std::string leptonName="";
    if (selectedIDIsoElectronsNum==1)   
      {
	leptonName="electron";
	leptonPt  = selectedIDIsoElectronsPt[0];
	leptonEta = selectedIDIsoElectronsEta[0];
	leptonPhi = selectedIDIsoElectronsPhi[0];
	leptonE   = selectedIDIsoElectronsE[0];
      }
    else if (selectedIDIsoMuonsNum==1)      
      {
	leptonName="muon";
	leptonPt  = selectedIDIsoMuonsPt[0];
	leptonEta = selectedIDIsoMuonsEta[0];
	leptonPhi = selectedIDIsoMuonsPhi[0];
	leptonE   = selectedIDIsoMuonsE[0];
      }
    else  {
      cout<<"Error!! No leptons. "<<endl;
      continue;
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

    W_mu.SetPtEtaPhiE(leptonPt,leptonEta,leptonPhi,leptonE);
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

    met   = sqrt(METPt*METPt);
    met_px = METPt*TMath::Cos(METPhi);
    met_py = METPt*TMath::Sin(METPhi);
    met_pz_type0 = pz1_type0;
    met_pz_type2 = pz1_type2;


    /////////////////LEPTONIC W

    TLorentzVector *W = new TLorentzVector();
    TLorentzVector *LEP = new TLorentzVector();
    TLorentzVector *NU  = new TLorentzVector();
    
    LEP->SetPtEtaPhiE(leptonPt,leptonEta,leptonPhi,leptonE);
    NU->SetPxPyPzE(met_px,met_py,met_pz_type0,met);
    *W = *LEP + *NU;
    
    W_pt = W->Pt();
    W_eta = W->Eta();
    W_phi = W->Phi();
    W_E = W->E();
    W_mt = W->Mt();

    //////////////////ANGULAR VARIABLES

    TLorentzVector *JET = new TLorentzVector();
    JET->SetPtEtaPhiE(AK8JetsPt[0],AK8JetsEta[0],AK8JetsPhi[0],AK8JetsE[0]);
    deltaR_lak8jet = JET->DeltaR(*LEP);
    deltaphi_METak8jet = JET->DeltaPhi(*NU);
    deltaphi_Vak8jet = JET->DeltaPhi(*W);

    //delete all the TLorentzVector before a new selection
    delete W;
    delete LEP;
    delete NU;
    delete JET;

    if (W_pt < 150) continue;
    if (deltaR_lak8jet < (TMath::Pi()/2.0))   continue;


    ///////////JETS
    for (unsigned int i=0; i<AK8JetsNum; i++)
      {
	if (AK8JetsPt[i]<30 || AK8JetsEta[i]>4.7)  continue;
	AK8jetPt[i]  = AK8JetsPt[i];
	AK8jetEta[i] = AK8JetsEta[i];
	AK8jetPhi[i] = AK8JetsPhi[i];
	AK8jetE[i]   = AK8JetsE[i];
	AK8jetPrunedMass[i]   = AK8Jets_prunedMass[i];
	AK8jetTrimmedMass[i]   = AK8Jets_trimmedMass[i];
	AK8jetFilteredMass[i]   = AK8Jets_filteredMass[i];
	AK8jetTau1[i]   = AK8Jets_tau1[i];
	AK8jetTau2[i]   = AK8Jets_tau2[i];
	AK8jetTau3[i]   = AK8Jets_tau3[i];
      }

    for (unsigned int i=0; i<JetsNum; i++)
      {
	if (JetsPt[i]<30 || JetsEta[i]>4.7)  continue;
	jetPt[i]  = JetsPt[i];
	jetEta[i] = JetsEta[i];
	jetPhi[i] = JetsPhi[i];
	jetE[i]   = JetsE[i];
	jet_bDiscr[i] = Jets_bDiscriminator[i];
      }

    /////////////////MC Infos
    if (isMC)
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

    //fill the tree
    outTree->Fill();
  }
  
  //--------close everything-------------
  chain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
