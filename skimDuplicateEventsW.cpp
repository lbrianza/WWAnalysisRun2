//c++ -o skimDuplicateEventsW skimDuplicateEventsW.cpp `root-config --cflags --glibs`

#include <iostream>
#include <fstream>
#include <cstdio>
#include <cstdlib>
#include "TROOT.h"
#include "TSystem.h"
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TString.h"
#include "TCut.h"
#include "TMath.h"
#include "TApplication.h"
#include "TError.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TPad.h"
#include "TStyle.h"
#include "TLorentzVector.h"
#include "TChain.h"

#include "TMath.h"
#include <cmath>

using namespace std;
int main(int argc, char** argv)
{
  TFile  myFile("WTOT.root"); // Open the file
  TTree* myTree = (TTree*) myFile.Get("Wtree"); // Get the Muon tree from the file

  int run,event;
  float invMass, ele1_eta, ele1_SCE, ele1_ET, caloMet;

  myTree->SetBranchAddress("run", &run);
  myTree->SetBranchAddress("event", &event);

  myTree->SetBranchAddress("caloMet",&caloMet);
  myTree->SetBranchAddress("ele1_eta",&ele1_eta);
  myTree->SetBranchAddress("ele1_SCE",&ele1_SCE);
  myTree->SetBranchAddress("ele1_ET",&ele1_ET);

  TFile *fileOut = new TFile("Wtree_final.root","RECREATE");
  TTree *treeOut  = new TTree("Wtree","Wtree");
  int Event, Run;
  float caloMET, ele1eta, ele1SCE, ele1ET;

  treeOut->Branch("Run",&Run,"Run/I");
  treeOut->Branch("Event",&Event,"Event/I");
  treeOut->Branch("caloMET",&caloMET,"caloMET/F");
  treeOut->Branch("ele1eta",&ele1eta,"ele1eta/F");
  treeOut->Branch("ele1SCE",&ele1SCE,"ele1SCE/F");
  treeOut->Branch("ele1ET",&ele1ET,"ele1ET/F");

  std::vector<int> runList;
  std::vector<int> eventList;

  std::cout<<"Events: "<<myTree->GetEntries()<<std::endl<<"Press any key to continue"<<std::endl;
  //  getchar();

  bool skipEvent=false;

  for (long int i=0; i<myTree->GetEntries(); i++)
    {
      myTree->GetEntry(i);

      if(i % 1000 == 0)
	cout << "read entry: " << i << "/" << myTree->GetEntries()<<endl;

      for (int j=0; j<runList.size(); j++) {
	if (run==runList.at(j) && event==eventList.at(j))
	  skipEvent=true;
      }

      if (skipEvent)  continue;

      runList.push_back(run);
      eventList.push_back(event);

      Event = event;
      Run   = run;

      ele1SCE = ele1_SCE;
      ele1eta = ele1_eta;
      ele1ET = ele1_ET;
      caloMET = caloMet;

      treeOut->Fill();

    }    

      fileOut->Write();
      fileOut->Print();

  return 0;
}
