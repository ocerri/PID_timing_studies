#include <Riostream.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include "TFile.h"
#include "TError.h"
#include "TMath.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TClonesArray.h"
#include "TH1.h"
#include "TF1.h"
#include "TString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TDatime.h"
#include "TLegend.h"
#include "TColor.h"
#include "TPaveText.h"
#include "TPaveStats.h"
#include "TGraphAsymmErrors.h"
#include "TMacro.h"
#include <vector>
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooFitResult.h"
#include "RooArgList.h"
#include "RooCBShape.h"
#include "RooAddPdf.h"
#include "RooPolynomial.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "RooPlot.h"
#include "TStyle.h"
#include "TChain.h"
#include "TTree.h"
#include "TSystem.h"
#include "TAxis.h"
#include "TMath.h"
#include "RooRandom.h"
#include "THStack.h"

using namespace RooFit;

//Double_t Luminosity = 500; //fb^-1 come articolo LEP3
Double_t Luminosity = 3500; //fb^-1 come FCC-ee/year

TString std_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>86.2&&M_ll<95.2";


Double_t find_br_upper_limit(Int_t seed,TH1 *h_bg, TH1 *h_sg , TH1 *h_WW, TH1 *h_ZZ, TH1 *h_HZ)
{
  RooRealVar x("x","Missing Mass [GeV]",40.,160.);
  RooDataHist data_bg("data_bg","datahist_bg",x,Import(*h_bg));
  RooDataHist data_WW("data_WW","datahist_WW",x,Import(*h_WW));
  RooDataHist data_ZZ("data_ZZ","datahist_ZZ",x,Import(*h_ZZ));
  RooDataHist data_HZ("data_HZ","datahist_HZ",x,Import(*h_HZ));

  // S e t u p   f i t   m o d e l
  // ---------------------
  RooHistPdf pdf_ZZ("pdf_ZZ", "PDF from histo (N=2)", x, data_ZZ, 2) ; // PDF from histo for ZZ bg
  RooHistPdf pdf_WW("pdf_WW", "PDF from histo (N=2)", x, data_WW, 2) ; // PDF from histo for WW bg
  RooRealVar n_ZZ("n_ZZ","n_ZZ",3700*Luminosity/500, 0, 10000*Luminosity/500);
  RooRealVar n_WW("n_WW","n_WW",6100*Luminosity/500, 0, 10000*Luminosity/500);
  RooAddPdf pdf_bg("pdf_bg", "Sum of backgrounds PDF", RooArgList(pdf_ZZ,pdf_WW), RooArgList(n_ZZ,n_WW)); // Sum of WW and Z peak
  RooHistPdf pdf_sg("pdf_sg", "PDF from histo (N=2)", x, data_HZ, 2) ; // PDF from histo for HZ bg

  RooRealVar n_sg("n_sg","signal ev",500*Luminosity/500,0,10000*Luminosity/500) ;
  RooRealVar n_bg("n_bg","background ev",10000*Luminosity/500,0,50000*Luminosity/500) ;
  RooAddPdf pdf_sum("pdf_sum", "Sum PDF", RooArgList(pdf_sg,pdf_bg), RooArgList(n_sg,n_bg));

  //Generate toydata set------------------------------------------------------------------------
  pdf_bg.fitTo(data_bg,Extended(kTRUE),SumW2Error(kTRUE),Save());
  RooRandom::randomGenerator()->SetSeed(111+seed);
  RooDataSet *toydata_bg = pdf_bg.generate(x,h_bg->Integral(),Extended(kTRUE));


  Double_t min_L = 0;
  vector<Double_t> v_nsg;
  vector<Double_t> v_min2NLL;
  Double_t nsg_max = 300.0;
  Double_t nsg_step = 5;

  for(Int_t i=0; i< nsg_max+nsg_step; i+=nsg_step){

    n_sg.setVal(i);
    v_nsg.push_back(i);
    n_sg.setConstant(kTRUE);

    RooFitResult *fr_tmp = pdf_sum.fitTo(*toydata_bg,SumW2Error(kTRUE),Extended(kTRUE),Save());

    v_min2NLL.push_back(2*fr_tmp->minNll());
    if(2*fr_tmp->minNll() < min_L){
      min_L=2*fr_tmp->minNll();
    }
  }

  TH1 *h_posterior = new TH1F("h_posterior","h_posterior", 1 + nsg_max/nsg_step, 0.0, nsg_max);

  for(unsigned int j=0; j<v_nsg.size() ; j++){
    v_min2NLL[j]-=min_L;
    h_posterior->SetBinContent(j+1, std::exp(-v_min2NLL[j]/2));
  }

  h_posterior->Scale(1/h_posterior->Integral());

  Double_t sum_part = 0;
  Int_t index_min = 0;
  Double_t alpha_limit = 0.05;
  while (sum_part < 1-alpha_limit ) {
    sum_part = h_posterior->Integral(0,index_min);
    index_min++;
  }

  index_min -=1;
  Double_t br_limit = v_nsg[index_min]/h_sg->Integral();

  // h_posterior->SetXTitle("Signal events");
  // TCanvas *c3 = new TCanvas("c3","c3",700,700);
  // c3->cd();
  // c3->SetGrid();
  // h_posterior->Draw();

  return 100*br_limit;

}
