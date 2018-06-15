#include "RooFit.h"
#include "Riostream.h"

#include "RooPowerLaw.h"
#include "RooRealVar.h"
#include "TMath.h"

ClassImp(RooPowerLaw)


RooPowerLaw::RooPowerLaw(const char * name, const char * title, RooAbsReal& _x,
			 RooAbsReal& _power) : 
  RooAbsPdf(name, title),
  x("x", "dependent", this, _x), 
  power("power", "power", this, _power) {
}

RooPowerLaw::RooPowerLaw(const RooPowerLaw& other, const char * name) :
  RooAbsPdf(other, name), 
  x("x", this, other.x), 
  power("power", this, other.power) {
}

Int_t RooPowerLaw::getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& anVars, 
					 const char * rangeName) const {
  if ((power == -1.0) && (x.min(rangeName) < 0.))
    return 0;
  if (matchArgs(allVars, anVars, x)) return 1 ;
  return 0;
}

Double_t RooPowerLaw::analyticalIntegral(Int_t code, 
					 const char * rangeName) const {
  switch(code) {
  case 1:
    Double_t ret(0.);
    if (power == -1.) {
      ret = TMath::Log(x.max(rangeName)) - TMath::Log(x.min(rangeName));
    } else {
      ret = (TMath::Power(x.max(rangeName), power + 1) - 
	     TMath::Power(x.min(rangeName), power + 1))/(power + 1);
    }
    return ret;
    break;
  }

  assert(0);
  return 0.;
}

Double_t RooPowerLaw::evaluate() const {
  return TMath::Power(x, power);
}

