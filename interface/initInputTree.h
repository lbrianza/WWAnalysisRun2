#ifndef __initInputTree__
#define __initInputTree__

#include "TTree.h"
#include "TChain.h"

// Declaration of leaf types
extern unsigned int RunNum;
extern unsigned int EvtNum;
extern float METPt;
extern float METPhi;
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
extern int GenBosonNum;
extern float GenBosonPt[100];
extern float GenBosonEta[100];
extern float GenBosonPhi[100];
extern float GenBosonE[100];
extern int GenBoson_GenBosonPDGId[100];
extern int GenMuNum;
extern float GenMuPt[100];
extern float GenMuEta[100];
extern float GenMuPhi[100];
extern float GenMuE[100];
extern int GenElecNum;
extern float GenElecPt[100];
extern float GenElecEta[100];
extern float GenElecPhi[100];
extern float GenElecE[100];
extern int GenNuNum;
extern float GenNuPt[100];
extern float GenNuEta[100];
extern float GenNuPhi[100];
extern float GenNuE[100];
extern unsigned int JetsNum;
extern float JetsPt[100];
extern float JetsEta[100];
extern float JetsPhi[100];
extern float JetsE[100];
extern float Jets_bDiscriminator[100];
extern unsigned int AK8JetsNum;
extern float AK8JetsPt[100];
extern float AK8JetsEta[100];
extern float AK8JetsPhi[100];
extern float AK8JetsE[100];
extern float AK8Jets_prunedMass[100];
extern float AK8Jets_trimmedMass[100];
extern float AK8Jets_filteredMass[100];
extern float AK8Jets_tau1[100];
extern float AK8Jets_tau2[100];
extern float AK8Jets_tau3[100];

// List of branches
extern TBranch *b_RunNum; 
extern TBranch *b_EvtNum; 
extern TBranch *b_METPt; 
extern TBranch *b_METPhi; 
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
extern TBranch *b_GenBosonNum;
extern TBranch *b_GenBosonPt;
extern TBranch *b_GenBosonEta;
extern TBranch *b_GenBosonPhi;
extern TBranch *b_GenBosonE;
extern TBranch *b_GenBoson_GenBosonPDGId;
extern TBranch *b_GenMuNum;
extern TBranch *b_GenMuPt;
extern TBranch *b_GenMuEta;
extern TBranch *b_GenMuPhi;
extern TBranch *b_GenMuE;
extern TBranch *b_GenElecNum;
extern TBranch *b_GenElecPt;
extern TBranch *b_GenElecEta;
extern TBranch *b_GenElecPhi;
extern TBranch *b_GenElecE;
extern TBranch *b_GenNuNum;
extern TBranch *b_GenNuPt;
extern TBranch *b_GenNuEta;
extern TBranch *b_GenNuPhi;
extern TBranch *b_GenNuE;
extern TBranch *b_JetsNum;
extern TBranch *b_JetsPt;
extern TBranch *b_JetsEta;
extern TBranch *b_JetsPhi;
extern TBranch *b_JetsE;
extern TBranch *b_Jets_bDiscriminator;
extern TBranch *b_AK8JetsNum;
extern TBranch *b_AK8JetsPt;
extern TBranch *b_AK8JetsEta;
extern TBranch *b_AK8JetsPhi;
extern TBranch *b_AK8JetsE;
extern TBranch *b_AK8Jets_prunedMass;
extern TBranch *b_AK8Jets_trimmedMass;
extern TBranch *b_AK8Jets_filteredMass;
extern TBranch *b_AK8Jets_tau1;
extern TBranch *b_AK8Jets_tau2;
extern TBranch *b_AK8Jets_tau3;

void InitTree(TChain* nt);

#endif
