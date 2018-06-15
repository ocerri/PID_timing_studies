#ifndef ROOGAMMAOC
#define ROOGAMMAOC

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"

class RooGammaOC : public RooAbsPdf {
public:
  RooGammaOC(){}

  RooGammaOC(const char *name, const char *title,
              RooAbsReal& _x,
              RooAbsReal& _lambda,
              RooAbsReal& _x0
           ):
             RooAbsPdf(name,title),
             x("x","x",this,_x),
             lambda("lambda","lambda",this,_lambda),
             x0("x0","x0",this,_x0)
           {};

  RooGammaOC(const RooGammaOC& other, const char* name=0):
    RooAbsPdf(other,name),
    x("x","x",this,other.x),
    lambda("lambda","lambda",this,other.lambda),
    x0("x0","x0",this,other.x0)
    { } ;
  virtual TObject* clone(const char* newname) const { return new RooGammaOC(*this,newname); }
  inline virtual ~RooGammaOC() { }
  Int_t getAnalyticalIntegral(RooArgSet& allVars, RooArgSet& analVars, const char* rangeName=0) const ;
  Double_t analyticalIntegral(Int_t code, const char* rangeName=0) const;

protected:

  RooRealProxy x;
  RooRealProxy lambda;
  RooRealProxy x0;

  Double_t evaluate() const {
    return exp(-lambda * x)/x;
  };

private:

  ClassDef(RooGammaOC,1)
};

#endif
