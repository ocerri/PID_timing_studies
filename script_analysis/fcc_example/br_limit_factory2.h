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
  Double_t nsg_max = 3000.0;
  Double_t nsg_step = 500;

//----------------------plot estetici------------------------------------------
  TLegend *leg_scan = new TLegend(0.7,0.6118012,0.9,0.9006211);
  // RooPlot* frame_scan = x.frame(Title("Background (pdf_{bg} + pdf_{sg})")) ;
  RooPlot* frame_scan = x.frame(Title(" ")) ;
  toydata_bg->plotOn(frame_scan,MarkerStyle(7),RooFit::Name("toydata_bg"));
  leg_scan->AddEntry(frame_scan->findObject("toydata_bg"),"Toy data","lep");
  //---------------------------------------------------------------------------
  Int_t colori[] = {1,2,3,4,6,46,38,8,9};
  Int_t icol = 0;
  for(Int_t i=0; i< nsg_max+nsg_step; i+=nsg_step){

    n_sg.setVal(i);
    v_nsg.push_back(i);
    n_sg.setConstant(kTRUE);

    RooFitResult *fr_tmp = pdf_sum.fitTo(*toydata_bg,SumW2Error(kTRUE),Extended(kTRUE),Save());

    //----------------------plot estetici------------------------------------------
    if(i>0 && i%50==0){
      // pdf_sum.plotOn(frame_scan,LineWidth(2),LineStyle(2),LineColor(i/100 +1),RooFit::Name(Form("pdf_sum%d",i)));
      pdf_sum.plotOn(frame_scan,LineWidth(2),LineStyle(2),LineColor(colori[icol]),RooFit::Name(Form("pdf_sum%d",i)));
      leg_scan->AddEntry(frame_scan->findObject(Form("pdf_sum%d",i)),Form("N signal = %d",i),"l");
      icol++;
    }
    //--------------------------------------------------------------------------

    v_min2NLL.push_back(2*fr_tmp->minNll());
    if(2*fr_tmp->minNll() < min_L){
      min_L=2*fr_tmp->minNll();
    }
  }

  TH1 *h_posterior = new TH1F("h_posterior","Posterior distribution, pdf_{post} = L #times pdf_{prior}", 1 + nsg_max/nsg_step, 0.0, nsg_max);

  for(unsigned int j=0; j<v_nsg.size() ; j++){
    v_min2NLL[j]-=min_L;
    h_posterior->SetBinContent(j+1, std::exp(-v_min2NLL[j]/2));
  }

  h_posterior->Scale(1/h_posterior->Integral());
  TF1 *gfit = new TF1("gfit","gaus",0,500);
  h_posterior->Fit(gfit);


  Double_t sum_part = 0;
  Double_t upper_limit = 50;
  // Double_t alpha_limit = 0.05;
  Double_t alpha_limit = 0.1;

  Double_t gf_norm = gfit->Integral(0,500);
  while (sum_part < 1-alpha_limit ) {
    upper_limit+=1;
    sum_part = gfit->Integral(0,upper_limit) / gf_norm;
  }

  Double_t br_limit = upper_limit/h_sg->Integral();

  h_posterior->SetXTitle("N_{signal}");
  h_posterior->SetStats(0);
  // TCanvas *c3 = new TCanvas("c3","c3",700,700);
  // c3->cd();
  // c3->SetGrid();
  // h_posterior->DrawCopy();
  // Int_t linea = 19;
  // TLegend *leg4 = new TLegend(0.3318928,0.8514908,0.6110631,0.9008028);
  // leg4->SetLineWidth(0);
  // leg4->SetLineColor(0);
  // leg4->SetFillStyle(0);
  // leg4->AddEntry(gfit,"Gaussian fit","l");
  // h_posterior->SetFillColor(4);
  // h_posterior->SetFillStyle(3003);
  // h_posterior->GetXaxis()->SetRange(0,linea);
  // h_posterior->DrawCopy("same");
  // leg4->AddEntry(h_posterior,Form("95%% CL zone"),"f");
  // leg4->Draw();

  TCanvas *c2 = new TCanvas("c2","c2", 700,700);
  c2->cd();
  c2->SetGrid();
  gPad->SetLeftMargin(0.15);
  frame_scan->GetYaxis()->SetTitleOffset(1.8);
  frame_scan->GetYaxis()->SetTitle("Events / 1.2 GeV");
  frame_scan->Draw();
  leg_scan->Draw();
  // c2->SaveAs("bg_fit_varying_Ns5.C");
  c2->SaveAs("bg_fit_varying_Ns5.pdf");
  // c2->SaveAs("bg_fit_varying_Ns5.root");
  c2->SaveAs("bg_fit_varying_Ns5.C");
  

  // TGraph *g_pll = new TGraph(v_nsg.size(),&v_nsg[0], &v_min2NLL[0]);
  // g_pll->GetYaxis()->SetTitle("-2log L");
  // g_pll->GetXaxis()->SetTitle("N_{signal}");
  // g_pll->SetTitle("Likelihood scan");
  // g_pll->SetLineStyle(2);
  // g_pll->SetLineColor(2);
  // TCanvas* c_pll = new TCanvas("c_pll","c_pll",700,700);
  // c_pll->cd();
  // gPad->SetLeftMargin(0.15);
  // g_pll->GetYaxis()->SetTitleOffset(1.6) ;
  // g_pll->Draw("AC*");
  // c_pll->SetGrid();
  // TLegend* leg5 = new TLegend(0.3318928,0.8514908,0.6810631,0.9008028,Form("N_{limit} = %.0f, Br limit = %.2f %%", upper_limit,100*br_limit));
  // leg5->Draw();

  // TFile *f_out1 = TFile::Open("_root/plots.root","RECREATE");
  // f_out1->cd();
  // c_pll->Write();
  // c2->Write();
  // c3->Write();

  //cout << upper_limit << "/" << h_sg->Integral() << endl;

  return 100*br_limit;

}
