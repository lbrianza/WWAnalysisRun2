 #ifndef PUWeight_H
  #define PUWeight_H
 
  #include <string>
  #include <vector>
 
  #include "TH1.h"
  #include "TAxis.h"
  #include "TString.h"

  // Compute pile-up weights to match data distribtution
  class PUWeight {
  public:
    enum Scenario { Summer12S10 };
 
    static Scenario toScenario(const std::string& str);
    static std::string toString(const Scenario sc);
 
    // Constructor. Initializes default behaviour to return PU weight of 1
    PUWeight();
 
    // Initialise weights for a given MC pile-up scenario. Can only be
    // called once.
    void initPUWeights(const TString& nameOfDataDistribution, const PUWeight::Scenario sc) {
  	initPUWeights( std::string( nameOfDataDistribution.Data() ), sc );
    }
    void initPUWeights(const std::string& nameOfDataDistribution, const PUWeight::Scenario sc);
 
    // Get weight factor dependent on number of added PU interactions
    // nPU is the 'true' number of interactions as explained in
    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#2012_Pileup_JSON_Files
    double getPUWeight(const int nPU) const;
 
 
  private:
    bool isInit_;
    int nPUMax_;
    std::vector<double> puWeigths_;
 
    // Generate weights for given data PU distribution
    std::vector<double> generateWeights(const Scenario sc, const TH1* data_npu_estimated) const;
  };
 
  #endif //PUWeight_H
