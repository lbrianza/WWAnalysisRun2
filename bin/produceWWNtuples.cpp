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

#include "../interface/initInputTree.h"
#include "../interface/setOutputTree.h"

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
    met   = METPt;
    met_px = METPt*TMath::Cos(METPhi);
    met_py = METPt*TMath::Sin(METPhi);
    met_pz = -999;  //TO BE FIXED
    nJets = NJets;
    nVtx  = NVtx;

    if      (selectedIDIsoElectronsNum==1)   
      {
	leptonPt  = selectedIDIsoElectronsPt[0];
	leptonEta = selectedIDIsoElectronsEta[0];
	leptonPhi = selectedIDIsoElectronsPhi[0];
	leptonE   = selectedIDIsoElectronsE[0];
      }
    else if (selectedIDIsoMuonsNum==1)      
      {
	leptonPt  = selectedIDIsoMuonsPt[0];
	leptonEta = selectedIDIsoMuonsEta[0];
	leptonPhi = selectedIDIsoMuonsPhi[0];
	leptonE   = selectedIDIsoMuonsE[0];
      }
    else  cout<<"Error!! No leptons. "<<endl;

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
