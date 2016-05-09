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
#include "../interface/setOutputTreeSynch.h"
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
  // float weight = std::atof(xSecWeight.c_str())/std::atof(numberOfEntries.c_str());
  if (strcmp(leptonName.c_str(),"el")!=0 && strcmp(leptonName.c_str(),"mu")!=0) {
    std::cout<<"Error: wrong lepton category"<<std::endl;
    return(-1);
  }

  setInputTree *ReducedTree = new setInputTree (inputTreeName.c_str());
  ReducedTree->Init();


  char command1[3000];
  sprintf(command1, "xrd eoscms dirlist %s/%s/  | awk '{print \"root://eoscms.cern.ch/\"$5}' > listTemp_%s.txt", (inputFolder).c_str(), (inputFile).c_str(), outputFile.c_str());
  std::cout<<command1<<std::endl;
  system(command1);
  char list1[2000];
  sprintf (list1, "listTemp_%s.txt", outputFile.c_str());
  ifstream rootList (list1);

  while (!rootList.eof())
    {
      char iRun_tW[700];
      rootList >> iRun_tW;
      ReducedTree->fChain->Add(iRun_tW);
    }

  char command3[300];
  sprintf(command3, "rm listTemp_%s.txt", outputFile.c_str());
  system(command3);

  // int cutEff[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  //---------output tree----------------
  TFile* outROOT = TFile::Open((std::string("output/output_")+leptonName+std::string("/")+outputFile+(".root")).c_str(),"recreate");
  outROOT->cd();
  TTree* outTree = new TTree("otree", "otree");
  outTree->SetDirectory(0);

  setOutputTreeSynch *WWTree = new setOutputTreeSynch(outTree);

  //---------start loop on events------------
  for (Long64_t jentry=0; jentry<ReducedTree->fChain->GetEntries();jentry++) {

    Long64_t iEntry = ReducedTree->LoadTree(jentry);
    if (iEntry < 0) break;
    ReducedTree->fChain->GetEntry(jentry);   
    // if (Cut(ientry) < 0) continue;                                                                                                                           

    if(iEntry % 1000 == 0)    
      cout << "read entry: " << iEntry << endl;

    WWTree->initializeVariables(); //initialize all variables
    
    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->ElectronsNum==0) continue; 
    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->MuonsNum==0) continue;      
    if (ReducedTree->AK8JetsNum<1) continue;
    //    if (ReducedTree->AK8JetsNum < 1 || ReducedTree->AK8Jets_AK8isLooseJetId[0]==false) continue; //at least one jet

    //    if (ReducedTree->AK8Jets_PtCorr[0] < 150) continue; 
    //    if (ReducedTree->METPt < 50) continue; //50!

    //lepton Pt preselection
    //    if ( strcmp(leptonName.c_str(),"el")==0 && ReducedTree->ElectronsPt[0]<30) continue; 
    //    if ( strcmp(leptonName.c_str(),"mu")==0 && ReducedTree->MuonsPt[0]<30) continue; 
    
    //save event variables
    WWTree->run   = ReducedTree->RunNum;
    WWTree->event = ReducedTree->EvtNum;
    WWTree->lumi = ReducedTree->LumiBlockNum;
    WWTree->nPV  = ReducedTree->NVtx;

    /////////////////LEPTON
    WWTree->nLooseEle=0;
    WWTree->nLooseMu=0;

    float l_e=0.;
    float jet_e=0.;

    for (int i=0; i<ReducedTree->ElectronsNum; i++) 
      if (ReducedTree->Electrons_isHEEP[i]==true && ReducedTree->ElectronsPt[i]>35) WWTree->nLooseEle++;
    for (int i=0; i<ReducedTree->MuonsNum; i++) 
      if (ReducedTree->Muons_isHighPt[i]==true && (ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])<0.1 && ReducedTree->MuonsPt[i]>20 && fabs(ReducedTree->MuonsEta[i])<2.4) WWTree->nLooseMu++;

    if (strcmp(leptonName.c_str(),"el")==0) {
      float tempPt=0.;
      for (int i=0; i<ReducedTree->ElectronsNum; i++) {
	if (ReducedTree->Electrons_isHEEP[i]==false) continue;       
	if (ReducedTree->ElectronsPt[i]<=90) continue;       
	if (ReducedTree->ElectronsPt[i]<tempPt) continue;
	WWTree->l_pt  = ReducedTree->ElectronsPt[i];
	WWTree->l_eta = ReducedTree->ElectronsEta[i];
	WWTree->l_phi = ReducedTree->ElectronsPhi[i];	
	l_e = ReducedTree->ElectronsE[i];	
	tempPt = WWTree->l_pt;
      }
    }
    else if (strcmp(leptonName.c_str(),"mu")==0) {
      float tempPt=0.;
      for (int i=0; i<ReducedTree->MuonsNum; i++) {
	if (ReducedTree->Muons_isHighPt[i]==false) continue;
    if ((ReducedTree->Muons_trackIso[i]/ReducedTree->MuonsPt[i])>=0.1) continue;
	if (ReducedTree->MuonsPt[i]<40) continue;
	if (fabs(ReducedTree->MuonsEta[i])>=2.4) continue;
	if (ReducedTree->MuonsPt[i]<tempPt) continue;
	WWTree->l_pt  = ReducedTree->MuonsPt[i];
	WWTree->l_eta = ReducedTree->MuonsEta[i];
	WWTree->l_phi = ReducedTree->MuonsPhi[i];
	l_e = ReducedTree->MuonsE[i];
	tempPt = WWTree->l_pt;
      }
    }


    ///////////FAT JET
    float tempPt=0.;
    for (unsigned int i=0; i<ReducedTree->AK8JetsNum; i++)
      {
	bool isCleanedJet = true;

	//	if (ReducedTree->AK8Jets_PtCorr[i]<30 || ReducedTree->AK8JetsEta[i]>4.7)  continue;
	if (ReducedTree->AK8Jets_PtCorr[i]<=tempPt) continue; //save the jet with largest pt
	if (ReducedTree->AK8Jets_AK8isLooseJetId[i]==false) continue;

	//CLEANING FROM LEPTONS
	for (int j=0; j<ReducedTree->ElectronsNum; j++) {
	  if (ReducedTree->Electrons_isHEEP[j]==false) continue;       
	  if (ReducedTree->ElectronsPt[j]<=90) continue;       
	  if (deltaR(ReducedTree->ElectronsEta[j], ReducedTree->ElectronsPhi[j],
		     ReducedTree->AK8JetsEta[i],   ReducedTree->AK8JetsPhi[i]) <1.0)
	    isCleanedJet = false;
	}
	
	for (int j=0; j<ReducedTree->MuonsNum; j++) {
	  if (ReducedTree->Muons_isHighPt[j]==false) continue;       
	  if ((ReducedTree->Muons_trackIso[j]/ReducedTree->MuonsPt[j])>=0.1) continue;
	  if (ReducedTree->MuonsPt[j]<40) continue;
	  if (fabs(ReducedTree->MuonsEta[j])>=2.4) continue;
	  if (deltaR(ReducedTree->MuonsEta[j], ReducedTree->MuonsPhi[j],
		     ReducedTree->AK8JetsEta[i],   ReducedTree->AK8JetsPhi[i]) <1.0) {
	    isCleanedJet = false;
	  }
	}


	if (isCleanedJet==false) continue; //jet is overlapped with a lepton

	WWTree->jet_pt  = ReducedTree->AK8Jets_PtCorr[i];
	WWTree->jet_eta = ReducedTree->AK8JetsEta[i];
	WWTree->jet_phi = ReducedTree->AK8JetsPhi[i];
	jet_e = ReducedTree->AK8Jets_ECorr[i];
        WWTree->jet_mass_pruned   = ReducedTree->AK8Jets_prunedMass[i];
        WWTree->jet_mass_softdrop   = ReducedTree->AK8Jets_softDropMass[i];
        WWTree->jet_tau2tau1   = ReducedTree->AK8Jets_tau2[i]/ReducedTree->AK8Jets_tau1[i];
	tempPt = WWTree->jet_pt;
      }

    ////////AK4 jets
    WWTree->njets=0;
    WWTree->nbtag=0;

    for (unsigned int i=0; i<ReducedTree->JetsNum; i++) //loop on AK4 jet
      {
	bool isCleanedJet = true;

	if (ReducedTree->Jets_isLooseJetId[i]==false) continue;
	if (ReducedTree->Jets_PtCorr[i]<=30) continue;
	if (fabs(ReducedTree->JetsEta[i])>=2.4) continue;

	//CLEANING
	if (deltaR(WWTree->jet_eta, WWTree->jet_phi,
		       ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.8)
	  isCleanedJet = false;

	//CLEANING FROM LEPTONS
	for (int j=0; j<ReducedTree->ElectronsNum; j++) {
	  if (ReducedTree->Electrons_isHEEP[j]==false) continue;       
	  if (ReducedTree->ElectronsPt[j]<=90) continue;       
	  if (deltaR(ReducedTree->ElectronsEta[j], ReducedTree->ElectronsPhi[j],
		     ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.3)
	    isCleanedJet = false;
	}      
	for (int j=0; j<ReducedTree->MuonsNum; j++) {
	  if (ReducedTree->Muons_isHighPt[j]==false) continue;       
	  if ((ReducedTree->Muons_trackIso[j]/ReducedTree->MuonsPt[j])>=0.1) continue;
	  if (ReducedTree->MuonsPt[j]<40) continue;
	  if (fabs(ReducedTree->MuonsEta[j])>=2.4) continue;
	  if (deltaR(ReducedTree->MuonsEta[j], ReducedTree->MuonsPhi[j],
		     ReducedTree->JetsEta[i],   ReducedTree->JetsPhi[i]) <0.3)
	    isCleanedJet = false;
	}

	if (isCleanedJet==false) continue;

	WWTree->njets++;
	if (ReducedTree->Jets_bDiscriminatorICSV[i]>0.679)   WWTree->nbtag++;

	if (ReducedTree->Jets_PtCorr[i]<WWTree->jet3_pt) continue;

	if (ReducedTree->Jets_PtCorr[i]>=WWTree->jet3_pt) {
	  if (ReducedTree->Jets_PtCorr[i]>=WWTree->jet2_pt) {
	    WWTree->jet3_pt = WWTree->jet2_pt;
	    WWTree->jet3_eta = WWTree->jet2_eta;
	    WWTree->jet3_phi = WWTree->jet2_phi;
	    WWTree->jet3_btag = WWTree->jet2_btag;
	    WWTree->jet2_pt = ReducedTree->Jets_PtCorr[i];
	    WWTree->jet2_eta = ReducedTree->JetsEta[i];
	    WWTree->jet2_phi = ReducedTree->JetsPhi[i];
	    WWTree->jet2_btag = ReducedTree->Jets_bDiscriminatorICSV[i];
	  }
	  else {
	    WWTree->jet3_pt = ReducedTree->Jets_PtCorr[i];
	    WWTree->jet3_eta = ReducedTree->JetsEta[i];
	    WWTree->jet3_phi = ReducedTree->JetsPhi[i];
	    WWTree->jet3_btag = ReducedTree->Jets_bDiscriminatorICSV[i];
	  }
	}
      }

    WWTree->pfMET   = sqrt(ReducedTree->METPt*ReducedTree->METPt);
    WWTree->pfMETPhi = ReducedTree->METPhi;

    //////////////MET
    TLorentzVector *W = new TLorentzVector();
    TLorentzVector *LEP = new TLorentzVector();
    TLorentzVector *NU0  = new TLorentzVector();
    TLorentzVector *NU2  = new TLorentzVector();

    if (WWTree->l_pt>=0) {

    // Calculate Neutrino Pz using all the possible choices : 
    // type0 -> if real roots, pick the one nearest to the lepton Pz except when the Pz so chosen
    //               is greater than 300 GeV in which case pick the most central root.
    // type1 -> type = 1: if real roots, choose the one closest to the lepton Pz if complex roots, use only the real part.
    //          type = 2: if real roots, choose the most central solution. if complex roots, use only the real part. 
    //          type = 3: if real roots, pick the largest value of the cosine*

    float Wmass = 80.385;

    TLorentzVector W_mu, W_Met;

    W_mu.SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,l_e);
    W_Met.SetPxPyPzE(ReducedTree->METPt * TMath::Cos(ReducedTree->METPhi), ReducedTree->METPt * TMath::Sin(ReducedTree->METPhi), 0., sqrt(ReducedTree->METPt*ReducedTree->METPt));

    //    if(W_mu.Pt()<=0 || W_Met.Pt() <= 0 ){ std::cerr<<" Negative Lepton - Neutrino Pt "<<std::endl; continue ; }


    // type0 calculation of neutrino pZ
    METzCalculator NeutrinoPz_type0;
    NeutrinoPz_type0.SetMET(W_Met);
    NeutrinoPz_type0.SetLepton(W_mu);
    NeutrinoPz_type0.SetLeptonType(leptonName.c_str());

    double pz1_type0 = NeutrinoPz_type0.Calculate(); // Default one -> according to type0
    // double pz2_type0 = NeutrinoPz_type0.getOther();  // Default one

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
    // double pz2_type2 = NeutrinoPz_type2.getOther();   // Default one

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

    /*    if (NeutrinoPz_type2.IsComplex()) {// if this is a complix, change MET
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
    */
    //    W_mass_type2 = (W_mu+W_neutrino_type2).M();
    //    W_pz_type2 = (W_mu+W_neutrino_type2).Pz();
    //    W_nu1_pz_type2 = pz1_type2;
    //    W_nu2_pz_type2 = pz2_type2;

    /////////////////LEPTONIC W
    
    LEP->SetPtEtaPhiE(WWTree->l_pt,WWTree->l_eta,WWTree->l_phi,l_e);
    NU0->SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),pz1_type0,WWTree->pfMET);
    NU2->SetPxPyPzE(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi),ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi),pz1_type2,TMath::Sqrt(ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi)*ReducedTree->METPt*TMath::Cos(ReducedTree->METPhi)+ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi)*ReducedTree->METPt*TMath::Sin(ReducedTree->METPhi)+pz1_type2*pz1_type2));
    *W = *LEP + *NU2;

    WWTree->W_pt = W->Pt();
    WWTree->W_eta = W->Eta();
    WWTree->W_phi = W->Phi();
    //    WWTree->v_mt = TMath::Sqrt(2*LEP->Et()*NU0->Et()*(1-TMath::Cos(LEP->DeltaPhi(*NU0))));
    //    W_mt = W->Mt();
    }

    TLorentzVector *JET = new TLorentzVector();

    if (WWTree->jet_pt>=0) {
      JET->SetPtEtaPhiE(WWTree->jet_pt, WWTree->jet_eta, WWTree->jet_phi, jet_e);
    }

    //    WWTree->deltaR_lak8jet = JET->DeltaR(*LEP);
    //    WWTree->deltaphi_METak8jet = JET->DeltaPhi(*NU0);
    //    WWTree->deltaphi_Vak8jet = JET->DeltaPhi(*W);
    //    if (WWTree->deltaR_lak8jet>(TMath::Pi()/2.0) && WWTree->deltaphi_METak8jet>2.0 && WWTree->deltaphi_Vak8jet>2.0)
      //      WWTree->issignal=1;

    //FOUR-BODY INVARIANT MASS

    if (WWTree->jet_pt>=0 && WWTree->l_pt>=0) {
      WWTree->m_lvj = (*LEP + *NU2 + *JET).M();
    }
    //    WWTree->mass_lvj_type2 = (*LEP + *NU2 + *JET).M();

    //delete all the TLorentzVector before a new selection
    delete W;
    delete LEP;
    delete NU0;
    delete NU2;
    delete JET;
    
    //    if (WWTree->v_pt < 150) continue;
    //    if (WWTree->deltaR_lak8jet < (TMath::Pi()/2.0))   continue;


    
    /////////////////MC Infos
    if (isMC)
    /*
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
    */    
    //fill the tree
    outTree->Fill();
    }

  //--------close everything-------------
  ReducedTree->fChain->Delete();
  outTree->Write();
  outROOT->Close();

  return(0);
}
