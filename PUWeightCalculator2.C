{
 
TFile *temp = new TFile("PU.root", "RECREATE");

 float inix= 0;
 float finx= 60;
 float nbinx= 60.0; 
 
TFile * fD = new TFile("/afs/cern.ch/user/l/lbrianza/work/public/WWAnalysis/pileupDATA.root");
 TH1F* hD = (TH1F*)fD->Get("pileup");

TH1F *hRatio= new TH1F("hRatio","hRatio500",nbinx,inix,finx);
// hD->Clone(hRatio);

 hD->Sumw2();
 hRatio->Sumw2();

cout<< "nDCount: " << endl;
double nDError;
cout<< "nDCount1: " << endl;
double nDCount = hD->IntegralAndError(1,nbinx-1,nDError);
cout<< "nDCount: " << nDCount<<"nDError: "<<nDError<<endl;

nDCount = hD->IntegralAndError(1,nbinx,nDError);
cout<< "nDCount: " << nDCount<<"nDError: "<<nDError<<endl;

TFile * fMC = new TFile("/afs/cern.ch/user/l/lbrianza/work/public/WWAnalysis/pileupMC.root");
 TH1F* h1 = (TH1F*)fMC->Get("pileup");

h1->Sumw2();

 double n1Error;
 double n1Count = h1->IntegralAndError(1,nbinx-1,n1Error);
 cout<<"n1Count: "<<n1Count<<"n1Error: "<<n1Error<<endl;

hD->Scale(n1Count/nDCount);
hD->Divide(h1);  
 

temp->cd();
hD->Write();
temp->Close();
delete temp;
}
