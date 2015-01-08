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

    chain->GetEntry(iEntry);
    run   = RunNum;
    event = EvtNum;
    met   = MET;

    outTree->Fill();
  }
  
  //--------close everything-------------
  chain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
