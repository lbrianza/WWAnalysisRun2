#include "../interface/analysisUtils.h"

double getDeltaPhi(double phi1, double phi2 ){
  const double PI = 3.14159265;
  double result = phi1 - phi2;

  if(result > PI) result = result - 2 * PI;
  if(result <= (-1 * PI)) result = result + 2 * PI;

  result = TMath::Abs(result);
  return result;
}

double deltaPhi(const double& phi1, const double& phi2)
{
  double deltaphi = fabs(phi1 - phi2);
  if (deltaphi > 6.283185308) deltaphi -= 6.283185308;
  if (deltaphi > 3.141592654) deltaphi = 6.283185308 - deltaphi;
  return deltaphi;
}
double deltaEta(const double& eta1, const double& eta2)
{
  double deltaeta = fabs(eta1 - eta2);
  return deltaeta;
}

double deltaR(const double& eta1, const double& phi1,
	      const double& eta2, const double& phi2)
{
  double deltaphi = deltaPhi(phi1, phi2);
  double deltaeta = deltaEta(eta1, eta2);
  double deltar = sqrt(deltaphi*deltaphi + deltaeta*deltaeta);
  return deltar;
}
