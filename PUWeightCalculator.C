{
 
TFile *temp = new TFile("PU.root", "RECREATE");

 float inix= 0;
 float finx= 60;
 float nbinx= 60.0; 
 
TFile * fD = new TFile("EOS/cms/store/caf/user/lbrianza/WWReducedTree_run2/TTbar_for_pu.root");
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

TH1F *hD= new TH1F("hD","h500",nbinx,inix,finx);
TH1F *hRatio= new TH1F("hRatio","hRatio500",nbinx,inix,finx);
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

cout<< "nDCount: " << endl;
double nDError;
cout<< "nDCount1: " << endl;
double nDCount = hD->IntegralAndError(1,nbinx-1,nDError);
cout<< "nDCount: " << nDCount<<"nDError: "<<nDError<<endl;

nDCount = hD->IntegralAndError(1,nbinx,nDError);
cout<< "nDCount: " << nDCount<<"nDError: "<<nDError<<endl;

TFile * fMC = new TFile("EOS/cms/store/caf/user/lbrianza/WWReducedTree_run2/TTbar_amcatnlo_50ns/ReducedSelection_1.root");
TTree * treeMC = (TTree *)fMC->Get("TreeMaker2/PreSelection");
treeMC->SetBranchAddress("NVtx", &nVtx);
treeMC->SetBranchAddress("METPt", &METPt);
treeMC->SetBranchAddress("ElectronsPt", &ElectronsPt);
treeMC->SetBranchAddress("MuonsPt", &MuonsPt);

Long64_t nentriesMC = treeMC->GetEntriesFast();
TH1F *h1= new TH1F("h1","h500",nbinx,inix,finx);
h1->Sumw2();

      treeMC->Draw("NVtx>>h1","(ElectronsPt[0]>90 || MuonsPt[0]>40) && METPt>20 && AK8JetsPt[0]>100");

/*
for(int i=0; i<= nentriesMC; i++)
{
    treeMC->GetEntry(i);
    if ((ElectronsPt[0]>90 || MuonsPt[0]>40) && METPt>20) {
    h1->Fill(nVtx);
    }
}
*/

double n1Error;
double n1Count = h1->IntegralAndError(1,nbinx-1,n1Error);
cout<<"n1Count: "<<n1Count<<"n1Error: "<<n1Error<<endl;

hRatio->Scale(n1Count/nDCount);
hRatio->Divide(h1);  
 

temp->cd();
hRatio->Write();
temp->Close();
delete temp;
}
