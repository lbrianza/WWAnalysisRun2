#include "../interface//PUWeight.h"
 #include "TFile.h"
 #include <iostream>
 
 //==============================================================================================
 // Get MC pile-up scenario from string representation
 PUWeight::Scenario PUWeight::toScenario(const std::string& str) {
   PUWeight::Scenario sc = Summer12S10;
   if( str == "PUS10" ) sc = Summer12S10;
   else {
     std::cerr << "\n\nERROR unknown scenario '" << str << "'" << std::endl;
     throw std::exception();
   }
 
   return sc;
 }
 
 //==============================================================================================
 // MC pile-up scenario to string representation
 std::string PUWeight::toString(const PUWeight::Scenario sc) {
   std::string str;
   if( sc == Summer12S10 ) str = "PUS10";
   else {
     std::cerr << "\n\nERROR unknown scenario '" << sc << "'" << std::endl;
     throw std::exception();
   }
 
   return str;
 }
 
 //==============================================================================================
 // Constructor. Initializes default behaviour to return PU weight of 1
 PUWeight::PUWeight()
   : isInit_(false), nPUMax_(0) {}
 
 //==============================================================================================
 // Initialise weights for a given MC pile-up scenario. Can only be
 // called once.
 void PUWeight::initPUWeights(const std::string& nameOfDataDistribution, const PUWeight::Scenario sc) {
 
   if( isInit_ ) {
     std::cerr << "\n\nERROR in PUWeight: weights already initialised" << std::endl;
     throw std::exception();
   }
 
   // Get data distribution from file
   TFile file(nameOfDataDistribution.c_str(), "READ");
   TH1* h = NULL;
   file.GetObject("pileup",h);
   if( h == NULL ) {
     std::cerr << "\n\nERROR in PUWeight: Histogram 'pileup' does not exist in file '" << nameOfDataDistribution << "'\n.";
     throw std::exception();
   }
   h->SetDirectory(0);
   file.Close();
 
   // Computing weights
   puWeigths_ = generateWeights(sc,h);
   nPUMax_ = puWeigths_.size();
 
   // Clean up
   delete h;
 
   isInit_ = true;
 }
 
 //==============================================================================================
 // Get weight factor dependent on number of added PU interactions
 double PUWeight::getPUWeight(const int nPU) const {
   double w = 1.;
   if( isInit_ ) {
     if( nPU >= nPUMax_ ) {
 //std::cerr << "WARNING: Number of PU vertices = " << nPU << " out of histogram binning." << std::endl;
 // In case nPU is out-of data-profile binning,
 // use weight from last bin
 w = puWeigths_.back();
     } else {
 w = puWeigths_.at(nPU);
     }
   }
 
   return w;
 }
 
 //==============================================================================================
 // Generate weights for given data PU distribution
 // Scenarios from: https://twiki.cern.ch/twiki/bin/view/CMS/Pileup_MC_Gen_Scenarios
 // Code adapted from: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
 std::vector<double> PUWeight::generateWeights(const PUWeight::Scenario sc, const TH1* data_npu_estimated) const {
 
   // Store probabilites for each pu bin
   unsigned int nPUMax = 0;
   double *npuProbs = 0;
 
   if( sc == Summer12S10 ) {
     nPUMax = 60;
     double npuSummer12_S10[60] = {
 2.560E-06,
 5.239E-06,
 1.420E-05,
 5.005E-05,
 1.001E-04,
 2.705E-04,
 1.999E-03,
 6.097E-03,
 1.046E-02,
 1.383E-02,
 1.685E-02,
 2.055E-02,
 2.572E-02,
 3.262E-02,
 4.121E-02,
 4.977E-02,
 5.539E-02,
 5.725E-02,
 5.607E-02,
 5.312E-02,
 5.008E-02,
 4.763E-02,
 4.558E-02,
 4.363E-02,
 4.159E-02,
 3.933E-02,
 3.681E-02,
 3.406E-02,
 3.116E-02,
 2.818E-02,
 2.519E-02,
 2.226E-02,
 1.946E-02,
 1.682E-02,
 1.437E-02,
 1.215E-02,
 1.016E-02,
 8.400E-03,
 6.873E-03,
 5.564E-03,
 4.457E-03,
 3.533E-03,
 2.772E-03,
 2.154E-03,
 1.656E-03,
 1.261E-03,
 9.513E-04,
 7.107E-04,
 5.259E-04,
 3.856E-04,
 2.801E-04,
 2.017E-04,
 1.439E-04,
 1.017E-04,
 7.126E-05,
 4.948E-05,
 3.405E-05,
 2.322E-05,
 1.570E-05,
 5.005E-06};
     npuProbs = npuSummer12_S10;
   }
 
   // Check that binning of data-profile matches MC scenario
   if( nPUMax != static_cast<unsigned int>(data_npu_estimated->GetNbinsX()) ) {
     std::cerr << "\n\nERROR number of bins (" << data_npu_estimated->GetNbinsX() << ") in data PU-profile does not match number of bins (" << nPUMax << ") in MC scenario " << toString(sc) << std::endl;
     throw std::exception();
   }
 
   std::vector<double> result(nPUMax,0.);
   double s = 0.;
   for(unsigned int npu = 0; npu < nPUMax; ++npu) {
     const double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));
     result[npu] = npu_estimated / npuProbs[npu];
     s += npu_estimated;
   }
   // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
   for(unsigned int npu = 0; npu < nPUMax; ++npu) {
     result[npu] /= s;
   }
 
   return result;
 }
