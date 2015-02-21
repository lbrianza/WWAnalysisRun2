#ifndef analysisUtils_h
#define analysisUtils_h

#include<iostream>
#include "TMath.h"

double getDeltaPhi(double phi1, double phi2 );

double deltaPhi(const double& phi1, const double& phi2);

double deltaEta(const double& eta1, const double& eta2);

double deltaR(const double& eta1, const double& phi1,
              const double& eta2, const double& phi2);

#endif
