#include "../interface/setOutputTreeSynch.h"

setOutputTreeSynch::setOutputTreeSynch(TTree* outTree){
  fTree = outTree;
  setBranches();
}

setOutputTreeSynch::~setOutputTreeSynch(){
  delete fTree;
}

void setOutputTreeSynch::initializeVariables()
{
  run=-99;
  event=-99;
  lumi=-99;
  nPV=-99;
  pfMET=-99;
  pfMETPhi=-99;
  nLooseEle=-99;
  nLooseMu=-99;
  l_pt=-99;
  l_eta=-99;
  l_phi=-99;
  jet_pt=-99;
  jet_eta=-99;
  jet_phi=-99;
  jet_mass_pruned=-99;
  jet_mass_softdrop=-99;
  jet_tau2tau1=-99;
  W_pt=-99;
  W_eta=-99;
  W_phi=-99;
  m_lvj=-99;
  njets=-99;
  nbtag=-99;
  jet2_pt=-99;
  jet2_eta=-99;
  jet2_phi=-99;
  jet2_btag=-99;
  jet3_pt=-99;
  jet3_eta=-99;
  jet3_phi=-99;
  jet3_btag=-99;
}

void setOutputTreeSynch::setBranches()
{
  fTree->Branch("run",&run,"run/I");
  fTree->Branch("event",&event,"event/I");
  fTree->Branch("lumi",&lumi,"lumi/I");
  fTree->Branch("nPV",&nPV,"nPV/I");
  fTree->Branch("pfMET",&pfMET,"pfMET/F");
  fTree->Branch("pfMETPhi",&pfMETPhi,"pfMETPhi/F");
  fTree->Branch("nLooseEle",&nLooseEle,"nLooseEle/I");
  fTree->Branch("nLooseMu",&nLooseMu,"nLooseMu/I");
  fTree->Branch("l_pt",&l_pt,"l_pt/F");
  fTree->Branch("l_eta",&l_eta,"l_eta/F");
  fTree->Branch("l_phi",&l_phi,"l_phi/F");
  fTree->Branch("jet_pt",&jet_pt,"jet_pt/F");
  fTree->Branch("jet_eta",&jet_eta,"jet_eta/F");
  fTree->Branch("jet_phi",&jet_phi,"jet_phi/F");
  fTree->Branch("jet_mass_pruned",&jet_mass_pruned,"jet_mass_pruned");
  fTree->Branch("jet_mass_softdrop",&jet_mass_softdrop,"jet_mass_softdrop");
  fTree->Branch("jet_tau2tau1",&jet_tau2tau1,"jet_tau2tau1");
  fTree->Branch("W_pt",&W_pt,"W_pt/F");
  fTree->Branch("W_eta",&W_eta,"W_eta/F");
  fTree->Branch("W_phi",&W_phi,"W_phi/F");
  fTree->Branch("m_lvj",&m_lvj,"m_lvj/F");
  fTree->Branch("njets",&njets,"njets/I");
  fTree->Branch("nbtag",&nbtag,"nbtag/I");
  fTree->Branch("jet2_pt",&jet2_pt,"jet2_pt/F");
  fTree->Branch("jet2_eta",&jet2_eta,"jet2_eta/F");
  fTree->Branch("jet2_phi",&jet2_phi,"jet2_phi/F");
  fTree->Branch("jet2_btag",&jet2_btag,"jet2_btag/F");
  fTree->Branch("jet3_pt",&jet3_pt,"jet3_pt/F");
  fTree->Branch("jet3_eta",&jet3_eta,"jet3_eta/F");
  fTree->Branch("jet3_phi",&jet3_phi,"jet3_phi/F");
  fTree->Branch("jet3_btag",&jet3_btag,"jet3_btag/F");
}


