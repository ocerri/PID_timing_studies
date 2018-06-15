#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
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

using namespace RooFit;


void fit_hist_WW()
{

  // Define the variables

  RooRealVar x("x","Missing Mass",40.,160.);

  // Read the histogram from file .root
  TFile *f = TFile::Open("_root/histo/plot_M_miss_paper2012_87_95Mll.root");
  //TFile f = new TFile("../_root/histo/plot_M_miss_paper2012_and_93_86.root","READ");
  TH1* h_bg =(TH1*)f->Get("h_background"); // background histogram
  TH1* h_sg =(TH1*)f->Get("h_signal"); // Signal histogram
  TH1* h_WW =(TH1*)f->Get("h_WW_miss"); // WW histogram


  TH1* h_sum = new TH1F("h_sum","hsum",60,40.,160.); // Create total histogram for signal with chosen BR plus background
  Double_t br = 0.07; // Higgs to invisible branching ratio
  h_sum->Add(h_sg,h_bg,br);


  RooDataHist data_bg("data_bg","datahist_bg",x,Import(*h_bg));
  RooDataHist data_WW("data_WW","datahist_WW",x,Import(*h_WW));
  RooDataHist data_sg("data_sg","datahist_sg",x,Import(*h_sg));
  RooDataHist data_sum("data_sum","datahist_sum",x,Import(*h_sum));


  // S e t u p   f i t   m o d e l
  // ---------------------


  RooRealVar cbmean1("cbmean1","cbmean1",91.,80,100) ;
  RooRealVar cbsigma1("cbsigma1","cbsigma1",2.,0,4) ;
  RooRealVar alpha1("alpha1","alpha1",-1.,-10,10) ;
  RooRealVar n1("n1","n1",2.,0,5) ;
  RooCBShape pdf_ZZ("pdf_ZZ", "Crystall Ball PDF", x, cbmean1, cbsigma1, alpha1, n1); // Crystall ball to fit Z peak

  RooHistPdf pdf_WW("pdf_WW", "PDF from histo (N=2)", x, data_WW, 2) ; // PDF from histo for WW bg


  RooRealVar n_ZZ("n_ZZ","n_ZZ",10000, 0, 100000);
  RooRealVar n_WW("n_WW","n_WW",10000, 0, 100000);
  RooAddPdf pdf_bg("pdf_bg", "Sum of backgrounds PDF", RooArgList(pdf_ZZ,pdf_WW), RooArgList(n_ZZ,n_WW)); // Sum of WW and Z peak


  RooRealVar cbmean2("cbmean2","cbmean2",124.,120,130) ;
  RooRealVar cbsigma2("cbsigma2","cbsigma2",2.,0,4) ;
  RooRealVar alpha2("alpha2","alpha2",-1.,-10,10) ;
  RooRealVar n2("n2","n2",2.,0,5) ;
  RooCBShape pdf_sg("pdf_sg", "Crystall Ball PDF", x, cbmean2, cbsigma2, alpha2, n2); // Crystall ball to fit signal


  RooRealVar n_bg("n_bg","n_bg",10000, 1000, 100000);
  RooRealVar n_sg("n_sg","n_sg",4000*br, 0, 10000);
  //RooRealVar sgfrac("sgfrac","fraction of signal",0.5,0.,1.) ;


  RooAddPdf pdf_sum("pdf_sum", "Sum PDF", RooArgList(pdf_sg,pdf_bg), RooArgList(n_sg,n_bg));



  // P l o t   a n d   f i t   a   R o o D a t a H i s t
  // ---------------------------------------------------

  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot* frame = x.frame(Title(Form("Background + signal (Br=%.2f) (pdf_sg) Poisson error bars",br))) ;
  RooPlot* frame2 = x.frame(Title(Form("Background + signal (Br=%.2f) (pdf_bg) with Poisson error bars",br))) ;
  RooPlot* frame3 = x.frame(Title("Background (pdf_bg) Poisson error bars")) ;
  data_bg.plotOn(frame3);
// data_sg.plotOn(frame2);
  data_sum.plotOn(frame);
  data_sum.plotOn(frame2);

  // P l o t   a n d   f i t   a   R o o D a t a H i s t   w i t h   i n t e r n a l   e r r o r s
  // ---------------------------------------------------------------------------------------------

  // If histogram has custom error (i.e. its contents is does not originate from a Poisson process
  // but e.g. is a sum of weighted events) you can data with symmetric 'sum-of-weights' error instead
  // (same error bars as shown by ROOT)

//  RooPlot* frame2 = x.frame(Title("Imported TH1 with internal errors")) ;
//  data_bg.plotOn(frame2,DataError(RooAbsData::SumW2));
//  data_sg.plotOn(frame2,DataError(RooAbsData::SumW2));

  // Please note that error bars shown (Poisson or SumW2) are for visualization only, the are NOT used
  // in a maximum likelihood fit
  //
  // A (binned) ML fit will ALWAYS assume the Poisson error interpretation of data (the mathematical definition
  // of likelihood does not take any external definition of errors). Data with non-unit weights can only be correctly
  // fitted with a chi^2 fit (see rf602_chi2fit.C)


  // Fit p.d.f to the data and plot

  pdf_bg.fitTo(data_bg,Extended(kTRUE));
  pdf_bg.plotOn(frame3,LineColor(kBlue));
//  pdf_bg.plotOn(frame2,LineColor(kBlue));

//  pdf_sg.fitTo(data_sg);
//  pdf_sg.plotOn(frame2,LineColor(kRed));
//  pdf_sg.plotOn(frame,LineColor(kRed));

  RooFitResult *fr_sum = pdf_sum.fitTo(data_sum,Extended(kTRUE),Save()); // Fitting sum pdf to background plus signal with fixed BR
  Double_t min_nll_sum = 2* fr_sum->minNll(); // Minimum of the Negative Log Likelihood fitting sum to bg+sg
  pdf_sum.plotOn(frame);

  RooFitResult *fr_bg = pdf_bg.fitTo(data_sum,Extended(kTRUE),Save()); // Fitting bg only pdf to background plus signal with fixed BR !!!!!
  Double_t min_nll_bg =  2* fr_bg->minNll(); // Minimum of the Negative Log Likelihood fitting bg only to bg+sg !!!!!
  pdf_bg.plotOn(frame2);

  Double_t significance = TMath::Sqrt(std::fabs(min_nll_bg - min_nll_sum));
  //sgfrac.Print();
  n_sg.Print();
  n_bg.Print();
  std::cout << "Significance =" << significance << std::endl;


  std::cout << "Total events = " << n_sg.getValV()+ n_bg.getValV() << std::endl;

  std::cout << std::endl << "Injected parameters : " << std::endl;
  std::cout << "HZ events = " << br*h_sg->Integral() << std::endl << "Total events = " << h_sum->Integral() << std::endl;




  // Print values of a0 and a1 (that now reflect fitted values and errors)
// a0.Print();
// a1.Print();
// a2.Print();
// a3.Print();
//
// cbmean.Print();
// cbsigma.Print();
// alpha.Print();
// n.Print();

  // ----------------------------- Likelihood

// RooAbsReal* nll = pdf_sum.createNLL(data_bg,NumCPU(4));
//
// RooPlot* frame1 = n_bg.frame(Bins(10),Range(0.01,0.95));
// nll->plotOn(frame1,ShiftToZero());





  // ----------------------------------

  // Draw all frames on a canvas
  TCanvas* c1 = new TCanvas("M_miss","Missing Mass",800,400) ;
  c1->Divide(3) ;
  c1->cd(1) ; gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
  c1->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
  c1->cd(3) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.6) ; frame3->Draw() ;


// TCanvas* c2 = new TCanvas("Likelihood","Likelihood",800,400) ;
// c->Divide(2) ;
// c->cd(1) ; gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
// c->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;


}
