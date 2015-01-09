#ifndef __initInputTree__
#define __initInputTree__

#include "TTree.h"
#include "TChain.h"

// Declaration of leaf types
extern unsigned int RunNum;
extern unsigned int EvtNum;
extern float MET;
extern unsigned int NJets;
extern unsigned int BTags;
extern unsigned int NVtx;
extern std::vector<float> *selectedIDIsoElectronsPt;
extern std::vector<float> *selectedIDIsoElectronsEta;
extern std::vector<float> *selectedIDIsoElectronsPhi;
extern std::vector<float> *selectedIDIsoElectronsE;
extern unsigned int selectedIDIsoElectronsNum;
extern std::vector<float> *selectedIDIsoMuonsPt;
extern std::vector<float> *selectedIDIsoMuonsEta;
extern std::vector<float> *selectedIDIsoMuonsPhi;
extern std::vector<float> *selectedIDIsoMuonsE;
extern unsigned int selectedIDIsoMuonsNum;
extern std::vector<float> *GenBosonPt;
extern std::vector<float> *GenBosonEta;
extern std::vector<float> *GenBosonPhi;
extern std::vector<float> *GenBosonE;
extern std::vector<int>   *GenBosonPDGId;
extern std::vector<float> *GenMuPt;
extern std::vector<float> *GenMuEta;
extern std::vector<float> *GenMuPhi;
extern std::vector<float> *GenMuE;
extern std::vector<float> *GenElecPt;
extern std::vector<float> *GenElecEta;
extern std::vector<float> *GenElecPhi;
extern std::vector<float> *GenElecE;
extern unsigned int MHTJetsNum;
extern std::vector<float> *MHTJetsPt;
extern std::vector<float> *MHTJetsEta;
extern std::vector<float> *MHTJetsPhi;
extern std::vector<float> *MHTJetsE;
extern std::vector<float> *MHTJets_bDiscriminator;

// List of branches
extern TBranch *b_RunNum; 
extern TBranch *b_EvtNum; 
extern TBranch *b_MET; 
extern TBranch *b_NJets;
extern TBranch *b_BTags;
extern TBranch *b_NVtx;
extern TBranch *b_selectedIDIsoElectronsPt;
extern TBranch *b_selectedIDIsoElectronsEta;
extern TBranch *b_selectedIDIsoElectronsPhi;
extern TBranch *b_selectedIDIsoElectronsE;
extern TBranch *b_selectedIDIsoElectronsNum;
extern TBranch *b_selectedIDIsoMuonsPt;
extern TBranch *b_selectedIDIsoMuonsEta;
extern TBranch *b_selectedIDIsoMuonsPhi;
extern TBranch *b_selectedIDIsoMuonsE;
extern TBranch *b_selectedIDIsoMuonsNum;
extern TBranch *b_GenBosonPt;
extern TBranch *b_GenBosonEta;
extern TBranch *b_GenBosonPhi;
extern TBranch *b_GenBosonE;
extern TBranch *b_GenBosonPDGId;
extern TBranch *b_GenMuPt;
extern TBranch *b_GenMuEta;
extern TBranch *b_GenMuPhi;
extern TBranch *b_GenMuE;
extern TBranch *b_GenElecPt;
extern TBranch *b_GenElecEta;
extern TBranch *b_GenElecPhi;
extern TBranch *b_GenElecE;
extern TBranch *b_MHTJetsNum;
extern TBranch *b_MHTJetsPt;
extern TBranch *b_MHTJetsEta;
extern TBranch *b_MHTJetsPhi;
extern TBranch *b_MHTJetsE;
extern TBranch *b_MHTJets_bDiscriminator;

void InitTree(TChain* nt);

#endif
