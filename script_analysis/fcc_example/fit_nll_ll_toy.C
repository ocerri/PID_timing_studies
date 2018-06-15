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


void fit_nll_ll_toy()
{

  // Define the variables

  RooRealVar x("x","Missing Mass",40.,160.);

  // Read the histogram from file .root
  TFile *f = TFile::Open("_root/histo/plot_M_miss_paper2012_and_93_86.root");
  TH1* h_bg =(TH1*)f->Get("h_background"); // background histogram
  TH1* h_sg =(TH1*)f->Get("h_signal"); // Signal histogram
  TH1* h_WW =(TH1*)f->Get("h_WW_miss"); // WW histogram
  TH1* h_ZZ =(TH1*)f->Get("h_ZZ_miss"); // ZZ histogram
  TH1* h_HZ =(TH1*)f->Get("h_HZ_miss"); // ZZ histogram

  TH1* h_sum = new TH1F("h_sum","hsum",60,40.,160.); // Create total histogram for signal with chosen BR plus background
  Double_t br = 0.2; // Higgs to invisible branching ratio
  h_sum->Add(h_sg,h_bg,br);


  RooDataHist data_bg("data_bg","datahist_bg",x,Import(*h_bg));
  RooDataHist data_sg("data_sg","datahist_sg",x,Import(*h_sg));
  RooDataHist data_sum("data_sum","datahist_sum",x,Import(*h_sum));
  RooDataHist data_WW("data_WW","datahist_WW",x,Import(*h_WW));
  RooDataHist data_ZZ("data_ZZ","datahist_ZZ",x,Import(*h_ZZ));
  RooDataHist data_HZ("data_HZ","datahist_HZ",x,Import(*h_HZ));


  // S e t u p   f i t   m o d e l
  // ---------------------


  RooHistPdf pdf_ZZ("pdf_ZZ", "PDF from histo (N=2)", x, data_ZZ, 2) ; // PDF from histo for ZZ bg
  RooHistPdf pdf_WW("pdf_WW", "PDF from histo (N=2)", x, data_WW, 2) ; // PDF from histo for WW bg
  RooRealVar n_ZZ("n_ZZ","n_ZZ",3700, 0, 10000);
  RooRealVar n_WW("n_WW","n_WW",6100, 0, 10000);
  RooAddPdf pdf_bg("pdf_bg", "Sum of backgrounds PDF", RooArgList(pdf_ZZ,pdf_WW), RooArgList(n_ZZ,n_WW)); // Sum of WW and Z peak

  RooHistPdf pdf_sg("pdf_sg", "PDF from histo (N=2)", x, data_HZ, 2) ; // PDF from histo for HZ bg

  RooRealVar n_sg("n_sg","signal ev",500,0,10000) ;
  RooRealVar n_bg("n_bg","background ev",10000,0,50000) ;
  RooAddPdf pdf_sum("pdf_sum", "Sum PDF", RooArgList(pdf_sg,pdf_bg), RooArgList(n_sg,n_bg));


  // ---------------------------------------------------

  pdf_bg.fitTo(data_bg,Extended(kTRUE),Save());
  RooDataSet *toydata_bg = pdf_bg.generate(x,h_bg->Integral(),Extended(kTRUE));

  pdf_sum.fitTo(data_sum,Extended(kTRUE),Save());
  RooDataSet *toydata_sum = pdf_sum.generate(x,h_sum->Integral(),Extended(kTRUE));

  // Make plot of binned dataset showing Poisson error bars (RooFit default)
  RooPlot* frame = x.frame(Title(Form("Background + signal (Br=%.2f) (pdf_sg + pdf_bg) Poisson error bars",br))) ;
  RooPlot* frame2 = x.frame(Title(Form("Background + signal (Br=%.2f) (pdf_bg) with Poisson error bars",br))) ;

  toydata_sum->plotOn(frame);
  toydata_sum->plotOn(frame2);



  // Fit p.d.f to the data and plot

  RooFitResult *fr_sum = pdf_sum.fitTo(*toydata_sum,Extended(kTRUE),Save()); // Fitting sum pdf to background plus signal with fixed BR
  Double_t min_nll_sum = 2* fr_sum->minNll(); // Minimum of the Negative Log Likelihood fitting sum to bg+sg
  pdf_sum.plotOn(frame);

  RooFitResult *fr_bg = pdf_bg.fitTo(*toydata_sum,Extended(kTRUE),Save()); // Fitting bg only pdf to background plus signal with fixed BR !!!!!
  Double_t min_nll_bg =  2* fr_bg->minNll(); // Minimum of the Negative Log Likelihood fitting bg only to bg+sg !!!!!
  pdf_bg.plotOn(frame2);

  // TCanvas* c1 = new TCanvas("M_miss","Missing Mass",800,400) ;
  // c1->Divide(2) ;
  // c1->cd(1) ; gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
  // c1->cd(2) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;

  Double_t significance = TMath::Sqrt(std::fabs(min_nll_bg - min_nll_sum));

  std::cout << "Total events = " << n_sg.getValV()+ n_bg.getValV() << std::endl;

  std::cout << std::endl << "Injected parameters : " << std::endl;
  std::cout << "HZ events = " << br*h_sg->Integral() << std::endl << "Total events = " << h_sum->Integral() << std::endl;



//Likelihood profile ----------------------------------------------------------------------------


  RooPlot* frame_scan = x.frame(Title("Background (pdf_bg + pdf_sg)")) ;
  toydata_bg->plotOn(frame_scan,MarkerStyle(7));

  Double_t min_L = 0;
  vector<Double_t> v_nsg;
  vector<Double_t> v_min2NLL;
  Double_t nsg_max = 200.0;
  Double_t nsg_step = 2;
  for(Int_t i=0; i< nsg_max+nsg_step; i+=nsg_step){
    //Put some line to re-initialize fit parameter to starting values

    n_sg.setVal(i);
    v_nsg.push_back(i);
    n_sg.setConstant(kTRUE);

    RooFitResult *fr_tmp = pdf_sum.fitTo(*toydata_bg,Extended(kTRUE),Save());
    if(i%50==0){
      pdf_sum.plotOn(frame_scan,LineWidth(1),LineStyle(2),LineColor(i/100 +1));
    }

    v_min2NLL.push_back(2*fr_tmp->minNll());
    if(2*fr_tmp->minNll() < min_L){
      min_L=2*fr_tmp->minNll();
    }
  }

  TH1 *h_posterior = new TH1F("h_posterior","h_posterior", 1 + nsg_max/nsg_step, 0.0, nsg_max);

  Int_t index_min = 0;
  Double_t dfrom4 = 10;
  for(Int_t j=0; j<v_nsg.size() ; j++){

    v_min2NLL[j]-=min_L;
    h_posterior->SetBinContent(j+1, std::exp(-v_min2NLL[j]/2));

    if(std::fabs(v_min2NLL[j]-4.0) < dfrom4){
      dfrom4 = std::fabs(v_min2NLL[j]-4);
      index_min = j;
    }
  }

  h_posterior->Scale(1/h_posterior->Integral());
  h_posterior->SetXTitle("Signal events");
  TCanvas *c3 = new TCanvas("c3","c3",700,700);
  c3->cd();
  c3->SetGrid();
  h_posterior->Draw();

  Double_t sum_part = 0;
  Int_t index_min2 = 0;
  Double_t alpha_limit = 0.05;
  while (sum_part < 1-alpha_limit ) {
    //sum_part+= h_posterior->GetBinContent(index_min2);
    sum_part = h_posterior->Integral(0,index_min2);
    index_min2++;
  }

  index_min =index_min2-1;
  Double_t br_limit = v_nsg[index_min]/h_sg->Integral();

  TGraph *g_pll = new TGraph(v_nsg.size(),&v_nsg[0], &v_min2NLL[0]);
  g_pll->GetYaxis()->SetTitle("-2ln(L)");
  g_pll->GetXaxis()->SetTitle("N_{signal}");
  g_pll->SetTitle("Likelihood scan");
  g_pll->SetLineStyle(2);
  g_pll->SetLineColor(2);
  TCanvas* c_pll = new TCanvas("c_pll","c_pll",700,700);
  c_pll->cd();
  g_pll->Draw("AC*");
  c_pll->SetGrid();
  TLegend* leg1 = new TLegend(0.3318928,0.8514908,0.6110631,0.9008028,Form("N_{limit} = %.0f", v_nsg[index_min]));
  TLegend* leg = new TLegend(0.0997151,0.8118012,0.3323837,0.9006211,Form("Br limit = %.2f %%", 100*br_limit));
  leg1->Draw();
  leg->Draw();

  TCanvas *c2 = new TCanvas("c2","c2", 700,700);
  c2->cd();
  c2->SetGrid();
  frame_scan->Draw();


  std::cout << "N_limit = " << v_nsg[index_min] << std::endl << "Br Limit = " << Form("%.2f",100*br_limit) << " %"<< std::endl;






}
