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

#include "../interface/setOutputTree.h"
#include "../interface/analysisUtils.h"
#include "../interface/pseudodataNtuple.h"

using namespace std;

//*******MAIN*******************************************************************                                                            
int main (int argc, char** argv)
{

TFile *temp = new TFile("PU.root", "RECREATE");

 // float inix= 0;
 // float finx= 60;
 float nbinx= 60.0; 
 /* 
TFile * fD = new TFile("EOS/cms/store/caf/user/lbrianza/WWReducedTree_run2/data_forPU.root");
TTree * tree = (TTree *)fD->Get("TreeMaker2/PreSelection");
int nVtx;
 float METPt;
 float ElectronsPt[10];
 float MuonsPt[10];
tree->SetBranchAddress("NVtx", &nVtx);
tree->SetBranchAddress("METPt", &METPt);
tree->SetBranchAddress("ElectronsPt", &ElectronsPt);
tree->SetBranchAddress("MuonsPt", &MuonsPt);

Long64_t nentries = tree->GetEntriesFast();

TH1F *hD= new TH1F("hD","hD",nbinx,inix,finx);
TH1F *hRatio= new TH1F("hRatio","hRatio",nbinx,inix,finx);
 hD->Sumw2();
 hRatio->Sumw2();
 //for(int i2=0; i2<= nentries; i2++)
  //{
  //    tree->GetEntry(i2);
    //    if ((ElectronsPt[0]>90 || MuonsPt[0]>40) && METPt>20) {
      tree->Draw("NVtx>>hD","(ElectronsPt[0]>90 || MuonsPt[0]>40) && METPt>20 && AK8JetsPt[0]>100");
      tree->Draw("NVtx>>hRatio","(ElectronsPt[0]>90 || MuonsPt[0]>40) && METPt>20 && AK8JetsPt[0]>100");

      //       hD->Fill(nVtx); 
      //       hRatio->Fill(nVtx);
      //    }
      //}
      */
 TFile * fD = new TFile("/afs/cern.ch/user/l/lbrianza/work/public/WWAnalysis/pileupDATA.root");
 TH1F* hD = (TH1F*)fD->Get("pileup");
 TH1F* hRatio = (TH1F*)fD->Get("pileup");

 hD->Sumw2();
 hRatio->Sumw2();

 double nDError;
 double nDCount = hD->IntegralAndError(1,nbinx-1,nDError);
 cout<< "nDCount: (data)" << nDCount<<"nDError: "<<nDError<<endl;

 nDCount = hD->IntegralAndError(1,nbinx,nDError);
 cout<< "nDCount: " << nDCount<<"nDError: "<<nDError<<endl;

 TFile * fMC = new TFile("/afs/cern.ch/user/l/lbrianza/work/public/WWAnalysis/pileupMC.root");
 TH1F* h1 = (TH1F*)fMC->Get("pileup");

 h1->Sumw2();

double n1Error;
double n1Count = h1->IntegralAndError(1,nbinx-1,n1Error);
cout<<"n1Count: (MC)"<<n1Count<<"n1Error: "<<n1Error<<endl;

hRatio->Scale(n1Count/nDCount);


 for (int iBin=0; iBin<60; iBin++) {
   //   std::cout<<iBin<<" "<<hRatio->GetBinContent(hRatio->FindBin(iBin))/(float)h1->GetBinContent(h1->FindBin(iBin))<<" "<<hD->GetBinContent(iBin)<<" "<<h1->GetBinContent(iBin)<<" "<<std::endl;

   if (h1->GetBinContent(iBin)==0 || hD->GetBinContent(iBin)==0)
     hRatio->SetBinContent(iBin,0.);
   else {
     hRatio->SetBinContent(iBin,(float)hRatio->GetBinContent(hRatio->FindBin(iBin))/(float)h1->GetBinContent(h1->FindBin(iBin)));
     // std::cout<<iBin<<" "<<(float)hRatio->GetBinContent(hRatio->FindBin(iBin))/(float)h1->GetBinContent(h1->FindBin(iBin))<<std::endl;
   }
   //   hRatio->Divide(h1);  
 } 

 for (int iBin=0; iBin<60; iBin++) {
   std::cout<<iBin<<" "<<hRatio->GetBinContent(iBin)<<std::endl;
 }
temp->cd();
hRatio->Write();
temp->Close();
delete temp;

 return 0;
}


