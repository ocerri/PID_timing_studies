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


void fit_sigma_scan()
{

  // Define the variables

  RooRealVar x("x","Missing Mass",40.,160.);

  // Read the histogram from file .root
  TFile *f = TFile::Open("_root/histo/plot_M_miss_paper2012_and_93_86.root");
  //TFile f = new TFile("../_root/histo/plot_M_miss_paper2012_and_93_86.root","READ");
  TH1* h_bg =(TH1*)f->Get("h_background"); // background histogram
  h_bg->Sumw2();
  TH1* h_sg =(TH1*)f->Get("h_signal"); // Signal histogram
  h_sg->Sumw2();
  TH1* h_WW =(TH1*)f->Get("h_WW_miss"); // WW histogram
  h_WW->Sumw2()


  RooDataHist data_bg("data_bg","datahist_bg",x,Import(*h_bg));
  RooDataHist data_sg("data_sg","datahist_sg",x,Import(*h_sg));
  RooDataHist data_WW("data_WW","datahist_WW",x,Import(*h_WW));


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


  RooRealVar n_sg("n_sg","signal ev",3500,0,10000) ;
  RooRealVar n_bg("n_bg","background ev",10000,0,50000) ;
  RooAddPdf pdf_sum("pdf_sum", "Sum PDF", RooArgList(pdf_sg,pdf_bg), RooArgList(n_sg,n_bg));



  Double_t min_L = 0;
  vector<Double_t> v_br;
  vector<Double_t> v_sigma;

  TH1* h_sum = 0;
  Int_t jj=10;
  Double_t sigma = 10;

  do{
    Double_t br = jj/100.0;
    //Put some line to re-initialize fit parameter to starting values
    n_sg.setVal(3500*br);
    n_bg.setVal(10000);

    v_br.push_back(jj);

    h_sum = new TH1F(Form("h_sum%d", jj),Form("h_sum%d", jj),60,40.,160.); // Create total histogram for signal with chosen BR plus backgroun
    h_sum->Add(h_sg,h_bg,br);
    RooDataHist data_sum("data_sum","datahist_sum",x,Import(*h_sum));

    RooFitResult *fr_bg = pdf_bg.fitTo(data_sum,Extended(kTRUE),NumCPU(2),Save());
    RooFitResult *fr_sum = pdf_sum.fitTo(data_sum,Extended(kTRUE),NumCPU(2),Save());

    sigma = TMath::Sqrt(2*std::fabs(fr_bg->minNll() - fr_sum->minNll()));
    v_sigma.push_back(sigma);
    jj -= 1;
    std::cout << std::endl << std::endl << "Iteration with br=" << jj+1 << "/100 is over" << std::endl << std::endl;
  } while( jj > 0 && sigma > 2);


  TGraph *g_pll = new TGraph(v_br.size(),&v_br[0], &v_sigma[0]);
  g_pll->GetYaxis()->SetTitle("Signifiance [#sigma]");
  g_pll->GetXaxis()->SetTitle("BR h->inv [per cent]");
  g_pll->SetTitle("Reachable signifiance at FCC-ee for H->inv");
  g_pll->SetLineStyle(2);
  g_pll->SetLineColor(2);
  TCanvas* c_pll = new TCanvas("c_pll","c_pll",700,700);
  c_pll->cd();
  g_pll->Draw("AC*");
  c_pll->SetGrid();


}
