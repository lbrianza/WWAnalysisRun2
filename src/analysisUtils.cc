#include "../interface/analysisUtils.h"

double getDeltaPhi(double phi1, double phi2 ){
  const double PI = 3.14159265;
  double result = phi1 - phi2;

  if(result > PI) result = result - 2 * PI;
  if(result <= (-1 * PI)) result = result + 2 * PI;

  result = TMath::Abs(result);
  return result;
}
