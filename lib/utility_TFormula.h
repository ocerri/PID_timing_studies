#include<TString.h>

TString ReduceAngle_TFormula(TString angle)
{
  TString out = "(";
  out += angle;
  out += Form("+2*TMath::Pi()*(TMath::Sign(0.5,-(%s)/TMath::Pi()-1)+0.5)", angle.Data());
  out += Form("-2*TMath::Pi()*(TMath::Sign(0.5,(%s)/TMath::Pi()-1)+0.5)", angle.Data());
  out += ")";
  return out;
}
