#include <iostream>
#include <math.h>
#include "TMath.h"

#include "RooGammaOC.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"

using namespace RooFit;

 ClassImp(RooGammaOC)

 Int_t RooGammaOC::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* range) const
 {
   if (matchArgs(allVars,analVars,x)) return 1;
   return 0;
 }

 Double_t RooGammaOC::analyticalIntegral(Int_t code, const char* rangeName) const
 {
   assert(code==1) ;
   double epsilon = 1.e-5;
   double upper_scale = 1e5;
   return TMath::Gamma(epsilon) * ( TMath::Gamma(epsilon, upper_scale) - TMath::Gamma(epsilon, lambda * x0) );

 }
