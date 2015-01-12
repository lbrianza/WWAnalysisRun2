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

    //    cout<<"qui"<<MHTJetsPt->size()<<endl;
    if (MHTJetsNum < 1) continue;
    if (MHTJetsPt[0] < 200) continue;
    if (MET < 50) continue;

    //save variables
    run   = RunNum;
    event = EvtNum;
    met   = MET;
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

    for (int i=0; i<MHTJetsNum; i++)
      {
	AK8jetPt[i]  = MHTJetsPt[i];
	AK8jetEta[i] = MHTJetsEta[i];
	AK8jetPhi[i] = MHTJetsPhi[i];
	AK8jetE[i]   = MHTJetsE[i];
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
