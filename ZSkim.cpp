//c++ -o ZSkim ZSkim.cpp `root-config --cflags --glibs`

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
  TH1F h1 ("h1","h1",100,60,120);
  //  TH1F h1 ("h1","h1",1000,0,1000);

  h1.SetDirectory(0);

  TChain* myTree = new TChain("TreeMaker2/PreSelection");

  TLorentzVector ele1 = TLorentzVector();
  TLorentzVector ele2 = TLorentzVector();

  float SCEle[2], etaEle[2], PhiEle[2], PtEle[2], deltaEtaSCTracker[2], deltaPhiSCTracker[2], sigmaIetaIeta[2];
  UShort_t ElectronsNum;
  UInt_t EvtNum, RunNum;

  float  ElectronsPt[200];        
  float  ElectronsEta[200];       
  float  ElectronsPhi[200];       
  float  Electrons_SCEnergy[200];       
  float Electrons_deltaEtaSCTracker[200];
  float Electrons_deltaPhiSCTracker[200];
  float Electrons_sigmaIetaIeta[200];
  float Electrons_neutralHadIso[200];
  float Electrons_chargedHadIso[200];
  float Electrons_photonIso[200];

  myTree->SetBranchAddress("EvtNum", &EvtNum);
  myTree->SetBranchAddress("RunNum", &RunNum);

  myTree->SetBranchAddress("ElectronsNum", &ElectronsNum);
  myTree->SetBranchAddress("ElectronsPt",        &ElectronsPt[0]);
  myTree->SetBranchAddress("ElectronsEta",        &ElectronsEta[0]);
  myTree->SetBranchAddress("ElectronsPhi",       &ElectronsPhi[0]);
  myTree->SetBranchAddress("Electrons_SCEnergy",       &Electrons_SCEnergy[0]);
  myTree->SetBranchAddress("Electrons_deltaEtaSCTracker", &Electrons_deltaEtaSCTracker[0]);
  myTree->SetBranchAddress("Electrons_deltaPhiSCTracker", &Electrons_deltaPhiSCTracker[0]);
  myTree->SetBranchAddress("Electrons_sigmaIetaIeta", &Electrons_sigmaIetaIeta[0]);

  myTree->SetBranchAddress("Electrons_neutralHadIso", &Electrons_neutralHadIso[0]);
  myTree->SetBranchAddress("Electrons_chargedHadIso", &Electrons_chargedHadIso[0]);
  myTree->SetBranchAddress("Electrons_photonIso", &Electrons_photonIso[0]);

  ifstream rootList ("ntupleList6.txt");

  std::string string1 = "ciao";
  std::string string2 = "";

  while (!rootList.eof())
    {
      char iRun_tW[700];
    
      rootList >> iRun_tW;
      string1=string(iRun_tW);

      if (string1 == string2) continue;
      string2=string(iRun_tW);
      
      std::cout<<iRun_tW<<std::endl;
      TChain* tTemp = new TChain("TreeMaker2/PreSelection");
      tTemp->Add(iRun_tW);

      myTree->Add(iRun_tW);

      tTemp->Delete();
    }

  //  TFile  myFile("EOS/cms/store/group/dpg_ecal/alca_ecalcalib/ecalMIBI/lbrianza/ntuple_doubleEG.root"); // Open the file
  //  TTree* myTree = (TTree*) myFile.Get("TreeMaker2/PreSelection"); // Get the Muon tree from the file

  /*
  t->SetBranchAddress("Electrons_SCEnergy", SCEle);
  t->SetBranchAddress("ElectronsEta", etaEle);
  t->SetBranchAddress("ElectronsPhi", PhiEle);
  t->SetBranchAddress("ElectronsPt", PtEle);


  t->SetBranchAddress("EvtNum", &EvtNum);
  t->SetBranchAddress("ElectronsNum", &ElectronsNum);
  */
  TFile *fileOut = new TFile("Ztree6.root","RECREATE");
  TTree *treeOut  = new TTree("Ztree","Ztree");
  int event, run;
  float invMass, ele1_eta, ele2_eta, ele1_SCE, ele2_SCE;

  treeOut->Branch("run",&run,"run/I");
  treeOut->Branch("event",&event,"event/I");
  treeOut->Branch("invMass",&invMass,"invMass/F");
  treeOut->Branch("ele1_eta",&ele1_eta,"ele1_eta/F");
  treeOut->Branch("ele2_eta",&ele2_eta,"ele2_eta/F");
  treeOut->Branch("ele1_SCE",&ele1_SCE,"ele1_SCE/F");
  treeOut->Branch("ele2_SCE",&ele2_SCE,"ele2_SCE/F");

  std::cout<<"Events: "<<myTree->GetEntries()<<std::endl<<"Press any key to continue"<<std::endl;
  //  getchar();

  for (long int i=0; i<myTree->GetEntries(); i++)
  //  for (long int i=0; i<300000; i++)
    {
      myTree->GetEntry(i);
      if(i % 1000 == 0)
	cout << "read entry: " << i << "/" << myTree->GetEntries()<<endl;

      if (ElectronsNum<2) continue;

      if (fabs(Electrons_deltaEtaSCTracker[0])>0.04 || fabs(Electrons_deltaEtaSCTracker[1])>0.04 || fabs(Electrons_deltaPhiSCTracker[0])>0.04 || fabs(Electrons_deltaPhiSCTracker[1])>0.04 || fabs(Electrons_sigmaIetaIeta[0])>0.015 || fabs(Electrons_sigmaIetaIeta[1])>0.05) continue;

	float theta1 = TMath::Sin(2*atan(exp(-1*ElectronsEta[0])));
	float theta2 = TMath::Sin(2*atan(exp(-1*ElectronsEta[1])));

	//	std::cout<<ElectronsPt[0]<<" "<<theta1<<std::endl;

	float iso1 = (Electrons_chargedHadIso[0]+Electrons_neutralHadIso[0]+Electrons_photonIso[0])/(Electrons_SCEnergy[0]*theta1);
	float iso2 = (Electrons_chargedHadIso[1]+Electrons_neutralHadIso[1]+Electrons_photonIso[1])/(Electrons_SCEnergy[1]*theta2);
	//		std::cout<<iso1<<" "<<iso2<<std::endl;

		if (Electrons_SCEnergy[0]*theta1<20 || Electrons_SCEnergy[1]*theta2<20) continue;
		if (iso1>0.3 || iso2>0.3) continue;

      //	if (SCEle[0]<15 || SCEle[1]<15) continue;
      //      std::cout<<"entry "<<EvtNum<<std::endl;
      //      if (i==17877) continue;

      // std::cout<<etaEle[0]<<" "<<etaEle[1]<<std::endl;


	ele1.SetPtEtaPhiE(Electrons_SCEnergy[0]*theta1,ElectronsEta[0],ElectronsPhi[0],Electrons_SCEnergy[0]);
	ele2.SetPtEtaPhiE(Electrons_SCEnergy[1]*theta2,ElectronsEta[1],ElectronsPhi[1],Electrons_SCEnergy[1]);

	//	ele1.SetPtEtaPhiE(ElectronsPt[0]*theta1,ElectronsEta[0],ElectronsPhi[0],Electrons_SCEnergy[0]);
	//      ele2.SetPtEtaPhiE(ElectronsPt[1]*theta2,ElectronsEta[1],ElectronsPhi[1],Electrons_SCEnergy[1]);

      //      ele2.SetPtEtaPhiE(PtEle[1],etaEle[1],PhiEle[1],SCEle[1]);

      //      std::cerr<<"ciao2"<<std::endl;

      invMass = (ele1+ele2).M();
      if (invMass<60 || invMass>120) continue;
      //      std::cout<<(ele1+ele2).M()<<std::endl;
      //      std::cout<<TMath::Cos(1)<<std::endl;

      //      std::cout<<Electrons_deltaPhiSCTracker[1]<<std::endl;
      //      std::cerr<<invMass<<std::endl;
      //float invMass = sqrt(2*Electrons_SCEnergy1*Electrons_SCEnergy2);//*(1-Math::Cos((ElectronsEta[0]-ElectronsEta[1])))); 
      //      std::cout<<iso1<<std::endl;
      h1.Fill(invMass);

      event = EvtNum;
      run = RunNum;

      ele1_SCE = Electrons_SCEnergy[0];
      ele2_SCE = Electrons_SCEnergy[1];
      ele1_eta = ElectronsEta[0];
      ele2_eta = ElectronsEta[1];

      treeOut->Fill();
      //        std::cout<<h1.GetEntries()<<std::endl;
      //      std::cerr<<"ciao"<<std::endl;
      //      std::cout<<i<<" "<<myTree->GetEntries()-1<<std::endl;

    }    
      std::cout<<h1.GetEntries()<<std::endl;
      h1.SaveAs("Zpeak6.root","root"); 	

      fileOut->Write();
      fileOut->Print();

  return 0;
}
