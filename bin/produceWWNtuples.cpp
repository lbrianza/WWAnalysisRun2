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
    if(iEntry % 100 == 0)      cout << "read entry: " << iEntry << endl;

    //get entry
    chain->GetEntry(iEntry);

    if (selectedIDIsoElectronsNum!=1 && selectedIDIsoMuonsNum!=1)  continue;      

    //    cout<<"qui"<<JetsPt->size()<<endl;
    if (JetsNum < 1 || AK8JetsNum < 1) continue;
    if (AK8JetsPt[0] < 200) continue;
    if (METPt < 50) continue;

    //save variables
    run   = RunNum;
    event = EvtNum;
    nJets = NJets;
    nVtx  = NVtx;
    
    //Lepton info
    std::string leptonName="";
    if      (selectedIDIsoElectronsNum==1)   
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

    //MET

    // Fill lepton information
    int isReal_type0;
    TLorentzVector mup;
    mup.SetPtEtaPhiE(leptonPt, leptonEta, leptonPhi, leptonE );

    TLorentzVector b_metpt; 
    b_metpt.SetPxPyPzE(METPt * cos(METPhi), METPt * sin(METPhi), 0, sqrt(METPt*METPt) );

    METzCalculator b_metpz_type0;
    b_metpz_type0.SetMET(b_metpt);
    b_metpz_type0.SetLepton(mup);
    b_metpz_type0.SetLeptonType(leptonName);

    double b_nvpz1_type0 = b_metpz_type0.Calculate(0); // Default one
    double b_nvpz2_type0 = b_metpz_type0.getOther() ;

    if(!b_metpz_type0.IsComplex()) isReal_type0=1;

    TLorentzVector b_nvp_type0_met;
    b_nvp_type0_met.SetPxPyPzE(b_metpt.Px(), b_metpt.Py(), b_nvpz1_type0, sqrt(b_metpt.Px()*b_metpt.Px() + b_metpt.Py()*b_metpt.Py() + b_nvpz1_type0*b_nvpz1_type0) );
    TLorentzVector b_nvp_type0;
    b_nvp_type0.SetPxPyPzE(b_metpt.Px(), b_metpt.Py(), b_nvpz1_type0, sqrt(b_metpt.Px()*b_metpt.Px() + b_metpt.Py()*b_metpt.Py() + b_nvpz1_type0*b_nvpz1_type0) );
    double W_mass_type0_met = (mup+b_nvp_type0_met).M(); 
    double W_pz_type0_met = (mup+b_nvp_type0_met).Pz(); 
    double W_nu1_pz_type0_met = b_nvpz1_type0; 
    double W_nu2_pz_type0_met = b_nvpz2_type0;
    //std::cout<<" type0 : pz1 "<<W_nu1_pz_type0<<" pz2 : "<<W_nu2_pz_type0<<" W_mass "<<W_mass<<" W_mass new "<<W_mass_type0_met<<std::endl;

    if (b_metpz_type0.IsComplex()) {// if this is a complix, change MET
      double nu_pt1 = b_metpz_type0.getPtneutrino(1);
      double nu_pt2 = b_metpz_type0.getPtneutrino(2);
      TLorentzVector tmpp1_type0;
      tmpp1_type0.SetPxPyPzE(nu_pt1 * cos(METPhi), nu_pt1 * sin(METPhi), b_nvpz1_type0, sqrt(nu_pt1*nu_pt1 + b_nvpz1_type0*b_nvpz1_type0) );
      TLorentzVector tmpp2_type0;
      tmpp2_type0.SetPxPyPzE(nu_pt2 * cos(METPhi), nu_pt2 * sin(METPhi), b_nvpz1_type0, sqrt(nu_pt2*nu_pt2 + b_nvpz1_type0*b_nvpz1_type0) );
      b_nvp_type0 = tmpp1_type0; if ( fabs((mup+tmpp1_type0).M()-80.4) > fabs((mup+tmpp2_type0).M()-80.4) ) b_nvp_type0 = tmpp2_type0;
    }
    double W_mass_type0 = (mup+b_nvp_type0).M(); 
    double W_pz_type0 = (mup+b_nvp_type0).Pz(); 
    double W_nu1_pz_type0 = b_nvpz1_type0; 
    double W_nu2_pz_type0 = b_nvpz2_type0;
    //std::cout<<" type0 : pz1 "<<W_nu1_pz_type0<<" pz2 : "<<W_nu2_pz_type0<<" W_mass "<<W_mass<<" W_mass new "<<W_mass_type0<<std::endl;

    met   = METPt;
    met_px = METPt*TMath::Cos(METPhi);
    met_py = METPt*TMath::Sin(METPhi);
    met_pz = b_nvpz1_type0;

    //jet infos
    for (unsigned int i=0; i<AK8JetsNum; i++)
      {
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
	jetPt[i]  = JetsPt[i];
	jetEta[i] = JetsEta[i];
	jetPhi[i] = JetsPhi[i];
	jetE[i]   = JetsE[i];
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
