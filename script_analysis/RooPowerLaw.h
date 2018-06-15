// -*- mode: C++ -*-
#ifndef ROO_POWERLAW_H
#define ROO_POWERLAW_H

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;
class RooAbsReal;

class RooPowerLaw : public RooAbsPdf {
public:
  RooPowerLaw() { } ;
  RooPowerLaw(const char * name, const char * title, RooAbsReal& _x,
	      RooAbsReal& _power) ;
  RooPowerLaw(const RooPowerLaw& other, const char * name = 0) ;

  virtual TObject * clone(const char * newname) const { return new RooPowerLaw(*this, newname); }

  inline virtual ~RooPowerLaw() { }

  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& anVars, 
  			      const char * rangeName = 0) const ;
  Double_t analyticalIntegral(Int_t code, const char * rangeName = 0) const ;

protected:
  RooRealProxy x;
  RooRealProxy power;

  Double_t evaluate() const;

private:
  ClassDef(RooPowerLaw, 1) // Power law PDF
};

#endif
