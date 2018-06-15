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
#define failure_limit 50



using namespace RooFit;

Double_t Luminosity = 500; //fb^-1 come articolo LEP3
//Double_t Luminosity = 3500; //fb^-1 come FCC-ee/year
//TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>76.2&&M_ll<106.2"; //original paper cuts
TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>86&&M_ll<93"; // modified cuts

Double_t find_br_upper_limit(Int_t seed,TH1 *h_bg, TH1 *h_sg , TH1 *h_WW, TH1 *h_ZZ, TH1 *h_HZ, TGraph **g_nll_profile, Int_t print_canvas =0)
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
  Double_t nsg_step = 10;

  TLegend *leg_scan = new TLegend(0.0997151,0.7118012,0.3323837,0.9006211);
  RooPlot* frame_scan = x.frame(Title("Background (pdf_bg + pdf_sg)")) ;
  toydata_bg->plotOn(frame_scan,MarkerStyle(7),RooFit::Name("toydata_bg"));
  leg_scan->AddEntry(frame_scan->findObject("toydata_bg"),"Toydata","lep");

  for(Int_t i=0; i< nsg_max+nsg_step; i+=nsg_step){

    n_sg.setVal(i);
    v_nsg.push_back(i);
    n_sg.setConstant(kTRUE);

    RooFitResult *fr_tmp = pdf_sum.fitTo(*toydata_bg,SumW2Error(kTRUE),Extended(kTRUE),Save());
    if(i%50==0){
      pdf_sum.plotOn(frame_scan,LineWidth(2),LineStyle(2),LineColor(i/100 +1),RooFit::Name(Form("pdf_sum%d",i)));
      leg_scan->AddEntry(frame_scan->findObject(Form("pdf_sum%d",i)),Form("N signal = %d",i),"l");
    }

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

  if(print_canvas==1){
    h_posterior->SetXTitle("Signal events");
    TCanvas *c3 = new TCanvas("c3","c3",700,700);
    c3->cd();
    c3->SetGrid();
    h_posterior->Draw();

    TCanvas *c2 = new TCanvas("c2","c2", 700,700);
    c2->cd();
    c2->SetGrid();
    frame_scan->Draw();
    leg_scan->Draw();

    TGraph *g_pll = new TGraph(v_nsg.size(),&v_nsg[0], &v_min2NLL[0]);
    g_pll->GetYaxis()->SetTitle("-2log L");
    g_pll->GetXaxis()->SetTitle("N_{signal}");
    g_pll->SetTitle("Likelihood profile");
    g_pll->SetLineStyle(2);
    g_pll->SetLineColor(2);
    TCanvas* c_pll = new TCanvas("c_pll","c_pll",700,700);
    c_pll->cd();
    gPad->SetLeftMargin(0.15);
    g_pll->GetYaxis()->SetTitleOffset(1.6) ;
    g_pll->Draw("AC*");
    c_pll->SetGrid();
    TLegend* leg = new TLegend(0.3318928,0.8514908,0.6810631,0.9008028,Form("N_{limit} = %.0f, Br limit = %.2f %%", v_nsg[index_min],100*br_limit));
    leg->Draw();
    *g_nll_profile = g_pll;
  }

  return 100*br_limit;

  // return 50.0; // Debug

}


Double_t find_significance_limit(Int_t seed,TH1 *h_bg, TH1 *h_sg, TH1 *h_WW, TH1 *h_ZZ, TH1 *h_HZ,
                            TGraphErrors **g_sig_br, Int_t print_canvas =0, int toy_count=5, int points=5)
{

  RooRealVar x("x","Missing Mass",40.,160.);
  RooDataHist data_bg("data_bg","datahist_bg",x,Import(*h_bg));
  RooDataHist data_sg("data_sg","datahist_sg",x,Import(*h_sg));
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


  //Generate toydata set----------------------------------------------------------------------


  Double_t v_br[toy_count][points];
  Double_t v_br_mean[points]; // Mean on toy models
  memset(v_br_mean, 0, points);
  Double_t v_br_rms[points]; // Zero: no errors on br
  memset(v_br_rms, 0, points);
  Double_t v_sigma[toy_count][points];
  Double_t v_sigma_mean[points];  // Significance for a given br is the mean on the toy models.
  memset(v_sigma_mean, 0, points);
  Double_t v_sigma_ms[points]; // Mean of squares from toy models
  memset(v_sigma_ms, 0, points);
  Double_t v_sigma_rms[points]; // Sigma root mean square evaluated from toy models.

  int i = 0;
  int failure_count = 0; // stop the cycle after failure_limit abort in fitting
  while ((i < toy_count)&&(failure_count<failure_limit)) // generate toy_count toy models, stop after failure_limit failures in fitting
  {
    pdf_bg.fitTo(data_bg,Extended(kTRUE),SumW2Error(kTRUE),Save());
    pdf_sg.fitTo(data_sg,SumW2Error(kTRUE),Save());
    RooRandom::randomGenerator()->SetSeed(111+seed);
    RooDataSet *toydata_bg = pdf_bg.generate(x,h_bg->Integral(),Extended(kTRUE));
    RooDataSet *toydata_sg = pdf_sg.generate(x,h_sg->Integral(),Extended(kTRUE));

    TH1* h_toy_bg =toydata_bg->createHistogram("h_toy_bg",x,Binning(60,40.,160.));
    TH1* h_toy_sg =toydata_sg->createHistogram("h_toy_sg",x,Binning(60,40.,160.));

    vector<Double_t> v_br_aux;
    vector<Double_t> v_sigma_aux;

    bool fit_status = true; // control fit convergence status
    int jj=0;

    while ((jj<points)&&(fit_status==true)) // find significance for various ("points") branching ratios beetween 0.01 and 0.1
    {

      Double_t br_min = 0.01;
      Double_t br_max = 0.1;
      Double_t br = br_min + (br_max-br_min)*jj/(points-1);

      TH1* h_toy_sum = new TH1F(Form("h_toy_sum%d", jj),Form("h_sum%d", jj),60,40.,160.); // Create total histogram for signal with chosen BR plus backgroun
      h_toy_sum->Add(h_toy_sg,h_toy_bg,br);
      RooDataHist toydata_sum("toydata_sum","toydata_sum",x,Import(*h_toy_sum));

      RooFitResult *fr_bg = pdf_bg.fitTo(toydata_sum,Extended(kTRUE),SumW2Error(kTRUE),NumCPU(4),Save());
      RooFitResult *fr_sum = pdf_sum.fitTo(toydata_sum,Extended(kTRUE),SumW2Error(kTRUE),NumCPU(4),Save());

      Double_t sigma = TMath::Sqrt(2*std::fabs(fr_bg->minNll() - fr_sum->minNll()));

      if ((sigma<0.001)||((jj>0)&&(sigma<v_sigma[i][jj-1]))) // break cycle if fit does not converge
      {
        fit_status = false;
        failure_count++;
      }

      if (fit_status == true) // Save br and significance for each point for this toy model
      {
        v_br[i][jj]=br*100;
        v_sigma[i][jj]=sigma;

        v_br_aux.push_back(br*100); // vector for intermediate printing
        v_sigma_aux.push_back(sigma);
      }

      jj++;
    }

    if (fit_status==true) // if fit did not converge generate a new toy model without increasing counter
    {
      // Plot graph for each toy model

      // TGraph *g_br = new TGraph(v_br_aux.size(),&v_br_aux[0], &v_sigma_aux[0]);
      // g_br->GetYaxis()->SetTitle("Signifiance [#sigma]");
      // g_br->GetXaxis()->SetTitle("BR h->inv [per cent]");
      // g_br->SetTitle(Form("Reachable signifiance at FCC-ee for H->inv (Toy model #%d)", i));
      // g_br->SetLineStyle(2);
      // g_br->SetLineColor(2);
      // TCanvas* c_br = new TCanvas("c_br","c_br",700,700);
      // c_br->cd();
      // g_br->Draw("AC*");
      // c_br->SetGrid();
      // c_br->SaveAs(Form("log/toy_model_%d.pdf",i));

      i++;
    }
  }

  if (failure_count<failure_limit)
  {
    int lim_5 = 0; // counter to determine the last j (br) with significance smaller than 5

    for (int j = 0; j < points; j++) // find mean on toy model of br and sigma for each point and lim_5
    {
      for (int i=0; i< toy_count; i++)
      {
        v_br_mean[j]= v_br_mean[j]+ v_br[i][j]/toy_count;
        v_sigma_mean[j]= v_sigma_mean[j]+ v_sigma[i][j]/toy_count;
      }
      if (v_sigma_mean[j]<5)
      {
        lim_5 = j;
      }
    }

    Double_t br_lim_5; // branching ratio limit at 5 sigma from linear interpolation
    br_lim_5 = (5 - v_sigma_mean[lim_5])*(v_br_mean[lim_5+1] - v_br_mean[lim_5])/(v_sigma_mean[lim_5+1] - v_sigma_mean[lim_5]) + v_br_mean[lim_5];


    for (int j = 0; j < points; j++) // find standard deviation on toy models of significance for each point
    {
      for (int i=0; i< toy_count; i++)
      {
        v_sigma_ms[j] = v_sigma_ms[j] + (v_sigma[i][j]-v_sigma_mean[j])*(v_sigma[i][j]-v_sigma_mean[j]);
      }
      v_sigma_rms[j] = TMath::Sqrt(v_sigma_ms[j]/toy_count);
    }

    if(print_canvas==1) // print
    {
      TGraphErrors *g_br_mean = new TGraphErrors(points,&v_br_mean[0], &v_sigma_mean[0],&v_br_rms[0],&v_sigma_rms[0]);
      g_br_mean->GetYaxis()->SetTitle("Significance [#sigma]");
      g_br_mean->GetXaxis()->SetTitle("Branching ratio H->inv [per cent]");
      g_br_mean->SetTitle("Significance from simulation, L = 500 fb^{-1}");
      g_br_mean->SetFillStyle(3003);
      g_br_mean->SetFillColor(9);
      g_br_mean->SetLineStyle(1);
      g_br_mean->SetLineWidth(2);
      g_br_mean->SetLineColor(9);
      g_br_mean->SetMarkerColor(1);
      g_br_mean->SetMarkerSize(0.8);
      g_br_mean->SetMarkerStyle(8);
      TCanvas* c_br_mean = new TCanvas("c_br","c_br",700,700);
      c_br_mean->cd();
      g_br_mean->Draw("AC3P");

      TLegend* leg = new TLegend(0.0997151,0.7118012,0.3323837,0.9006211,Form("BR 5 sigma = %.2f %%", br_lim_5));
      leg->Draw();

      c_br_mean->SetGrid();
      *g_sig_br = g_br_mean;
    }

    return br_lim_5;
  }

  else
  {
    std::cout << std::endl << std::endl << std:: endl << "Failure!";
    std::cout << std::endl << std::endl << std:: endl;

    Double_t br_lim_5 = 100.0;
    return br_lim_5;
  }

  // return 100.0; // Debug
}
