#include "../interface/METzCalculator_Run2.h"
#include "TMath.h"

/// constructor
METzCalculator_Run2::METzCalculator_Run2() {
  isComplex_ = false;
  otherSol_ = 0.;
  Ctype_ = -999;
  leptonMass_ = 0.105658367;
  newPtneutrino1_ = -1;
  newPtneutrino2_ = -1;
}

/// destructor
METzCalculator_Run2::~METzCalculator_Run2() {
}

/// member functions
double METzCalculator_Run2::Calculate(int type) {
  double M_W = 80.4;
  double emu  = lepton_.E();
  double etamu = lepton_.Eta();
  double pxmu = lepton_.Px();
  double pymu = lepton_.Py();
  double pzmu = lepton_.Pz();
  double pxnu = MET_.Px();
  double pynu = MET_.Py();
  double pznu = 0.;
  otherSol_ = 0.;
  double az = pzmu/emu;
  double d = (M_W*M_W/2. + pxnu*pxmu + pynu*pymu)/emu;
  double B = 2*d*az/(az*az-1.);
  double C = (d*d - pxnu*pxnu - pynu*pynu) / (az*az-1) ;
  double tmproot = B*B - 4.0*C;

  if (tmproot<0) {
    isComplex_= true;
    // solve for W mass that has a real root

    double pnu = MET_.E();
    double g =  TMath::Sqrt( (1. - az*az) * pnu*pnu );
    double gg = g - (pxmu*pxnu + pxmu*pxnu)/emu ;
    double mmw = TMath::Sqrt(2.*emu*gg);
   
    double dd = (mmw*mmw/2. + pxnu*pxmu + pynu*pymu)/emu;
    double BB = 2*dd*az/(az*az-1);
    pznu = - BB/2.;
    Ctype_=1;

    otherSol_ = pznu; 

    if( mmw < 76 || mmw > 84 ) {
      Ctype_=2;
      g = TMath::Sqrt( 1. - az*az);
      double gg = g - (pxmu*pxnu + pxmu*pxnu)/(emu*pnu) ;
      double metrec = M_W * M_W / ( 2.*emu*gg);
      double mmetx = metrec * pxnu / pnu ;
      double mmety = metrec * pynu / pnu ;
      d = (M_W*M_W/2. + mmetx*pxmu + mmety*pymu)/emu;
      B = 2*d*az/(az*az-1);
      pznu = - B/2.;
    }
  }
  else {
    isComplex_ = false;
    Ctype_=0;
    double tmpsol1 = (-B + TMath::Sqrt(tmproot))/(2.0);
    double tmpsol2 = (-B - TMath::Sqrt(tmproot))/(2.0);
    if(pzmu>0.) {pznu = tmpsol1; otherSol_ = tmpsol2;}
      else {pznu = tmpsol2; otherSol_ = tmpsol1;}

//  errors with low pz and low yu correlations
    if( TMath::Abs(etamu)<0.4 || TMath::Abs(pznu) <25. ) {
      Ctype_=-1;
      if(TMath::Abs(etamu)<0.4) Ctype_=-2;
      if(TMath::Abs(etamu)<0.4&&TMath::Abs(pznu) <25.) Ctype_=-3;

      if ( TMath::Abs(tmpsol2-pzmu) < TMath::Abs(tmpsol1-pzmu)) { pznu = tmpsol2; otherSol_ = tmpsol1;}
      else { pznu = tmpsol1; otherSol_ = tmpsol2; }
      
    }
  }


  return pznu;
}
