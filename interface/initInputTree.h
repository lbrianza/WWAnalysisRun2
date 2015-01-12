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
extern float selectedIDIsoElectronsPt[100];
extern float selectedIDIsoElectronsEta[100];
extern float selectedIDIsoElectronsPhi[100];
extern float selectedIDIsoElectronsE[100];
extern unsigned int selectedIDIsoElectronsNum;
extern float selectedIDIsoMuonsPt[100];
extern float selectedIDIsoMuonsEta[100];
extern float selectedIDIsoMuonsPhi[100];
extern float selectedIDIsoMuonsE[100];
extern unsigned int selectedIDIsoMuonsNum;
extern float GenBosonPt[100];
extern float GenBosonEta[100];
extern float GenBosonPhi[100];
extern float GenBosonE[100];
extern int GenBoson_GenBosonPDGId[100];
extern float GenMuPt[100];
extern float GenMuEta[100];
extern float GenMuPhi[100];
extern float GenMuE[100];
extern float GenElecPt[100];
extern float GenElecEta[100];
extern float GenElecPhi[100];
extern float GenElecE[100];
extern unsigned int MHTJetsNum;
extern float MHTJetsPt[100];
extern float MHTJetsEta[100];
extern float MHTJetsPhi[100];
extern float MHTJetsE[100];
extern float MHTJets_bDiscriminator[100];

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
extern TBranch *b_GenBoson_GenBosonPDGId;
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
