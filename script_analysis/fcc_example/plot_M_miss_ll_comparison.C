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
#include "THStack.h"
#include "TLegendEntry.h"



//Create histos-----------------------------------------------------------------
Double_t Luminosity = 3500; //FCCee
TCanvas* c_M_miss = new TCanvas("c_M_miss", "c_M_miss",0,800,1250,900);

TLegend* leg1 = new TLegend(0.7407705,0.6983945,0.9911717,0.928899);

TH1 *h_sum_ILD = new TH1F("h_sum_ILD","Missing Mass in Z#rightarrow l^{+}l^{-} tagged events (M_{Z}+/- 4 GeV)",60,40.0,160.0 );
TH1 *h_sum_CMS = new TH1F("h_sum_CMS","Missing Mass in Z#rightarrow l^{+}l^{-} tagged events (M_{Z}+/- 4 GeV)",60,40.0,160.0 );
THStack *h_sum_miss = new THStack("h_sum_miss","Missing Mass in Z#rightarrow l^{+}l^{-} tagged events (M_{Z}+/- 4 GeV)");
THStack *h_sum_miss_ILD = new THStack("h_sum_miss_ILD","Missing Mass in Z#rightarrow l^{+}l^{-} tagged events (M_{Z}+/- 4 GeV)");
//TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>76.2&&M_ll<106.2";
TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>87&&M_ll<95"; // modified cuts


//Merge code -------------------------------------------------------------------
THStack* plot_M_miss_ll_ILD(TString cuts = "paper2012", const char* outtag="cuts"){

  if( cuts.Contains("paper2012") ){
    cuts.ReplaceAll("paper2012",paper2012_cuts);
  }

  const char* cut_in = cuts;


  //load data from trees--------------------------------------------------------
  TFile *f_WW = TFile::Open("_root/tree_WW_W2lep_2M_ILDdelphes_ll.root");
  TTree *tree_WW = (TTree*) f_WW->Get("tree_ll");
  Double_t scale_WW = 1.715e-09/2e6; //Cross section divided total events generated

  TFile *f_ZZ = TFile::Open("_root/tree_ZZ_Z2lep_2M_ILDdelphes_ll.root");
  TTree *tree_ZZ = (TTree*) f_ZZ->Get("tree_ll");
  Double_t scale_ZZ = 1.370e-10/2e6; //Cross section divided total events generated

  TFile *f_HZ = TFile::Open("_root/tree_HZ_H2inv_Z2lep_2M_ILDdelphes_ll.root");
  TTree *tree_HZ = (TTree*) f_HZ->Get("tree_ll");
  Double_t scale_HZ = 5.720e-11/2e6; //Cross section divided total events generated

  //Fill the histos-------------------------------------------------------------
  TH1* h_WW_miss = new TH1D("h_WW_miss","h_WW_miss",60,40.0,160.0);
  TH1* h_ZZ_miss = new TH1D("h_ZZ_miss","h_ZZ_miss",60,40.0,160.0);
  TH1* h_HZ_miss = new TH1D("h_HZ_miss","h_HZ_miss",60,40.0,160.0);

  Int_t n_surv_WW;
  n_surv_WW = tree_WW->Project("h_WW_miss", "M_miss", cut_in);

  Int_t n_surv_ZZ;
  n_surv_ZZ = tree_ZZ->Project("h_ZZ_miss", "M_miss", cut_in);

  Int_t n_surv_HZ;
  n_surv_HZ = tree_HZ->Project("h_HZ_miss", "M_miss", cut_in);


  //Scaling histos and plotting final M_miss plot------------------------------
  h_WW_miss->Scale(Luminosity*1e12*scale_WW);
  h_WW_miss->SetLineColor(6);
  h_WW_miss->SetLineWidth(3);
  h_WW_miss->SetFillColor(0);
  h_WW_miss->SetFillStyle(0);

  h_ZZ_miss->Scale(Luminosity*1e12*scale_ZZ);
  h_ZZ_miss->SetLineColor(4);
  h_ZZ_miss->SetLineWidth(3);
  h_ZZ_miss->SetFillColor(0);
  h_ZZ_miss->SetFillStyle(0);

  h_HZ_miss->Scale(Luminosity*1e12*scale_HZ);
  h_HZ_miss->SetLineColor(2);
  h_HZ_miss->SetLineWidth(3);
  h_HZ_miss->SetFillColor(0);
  h_HZ_miss->SetFillStyle(0);

  h_sum_ILD->Add(h_WW_miss,h_ZZ_miss);
  h_sum_ILD->Add(h_HZ_miss);

  h_sum_miss_ILD->Add(h_WW_miss);
  h_sum_miss_ILD->Add(h_ZZ_miss);
  h_sum_miss_ILD->Add(h_HZ_miss);

  h_ZZ_miss->SetStats(0);
  h_HZ_miss->SetStats(0);
  h_WW_miss->SetStats(0);

  leg1->AddEntry(h_HZ_miss,Form("Signal ILD"),"l");
  leg1->AddEntry(h_ZZ_miss,Form("ZZ ILD"),"l");
  leg1->AddEntry(h_WW_miss,Form("WW ILD"),"l");

  return h_sum_miss_ILD;
}



THStack* plot_M_miss_ll(TString cuts = "paper2012", const char* outtag="cuts"){

  if( cuts.Contains("paper2012") ){
    cuts.ReplaceAll("paper2012",paper2012_cuts);
  }

  const char* cut_in = cuts;


  //load data from trees--------------------------------------------------------
  TFile *f_WW = TFile::Open("_root/tree_WW_W2lep_1M_CMSdelphes_ll.root");
  TTree *tree_WW = (TTree*) f_WW->Get("tree_ll");
  Double_t scale_WW = 1.715e-09/1e6; //Cross section divided total events generated

  TFile *f_ZZ = TFile::Open("_root/tree_ZZ_Z2lep_1M_CMSdelphes_ll.root");
  TTree *tree_ZZ = (TTree*) f_ZZ->Get("tree_ll");
  Double_t scale_ZZ = 1.370e-10/1e6; //Cross section divided total events generated

  TFile *f_HZ = TFile::Open("_root/tree_HZ_H2inv_Z2lep_1M_CMSdelphes_ll.root");
  TTree *tree_HZ = (TTree*) f_HZ->Get("tree_ll");
  Double_t scale_HZ = 5.720e-11/1e6; //Cross section divided total events generated

  //Fill the histos-------------------------------------------------------------
  TH1* h_WW_miss = new TH1D("h_WW_miss","h_WW_miss",60,40.0,160.0);
  TH1* h_ZZ_miss = new TH1D("h_ZZ_miss","h_ZZ_miss",60,40.0,160.0);
  TH1* h_HZ_miss = new TH1D("h_HZ_miss","h_HZ_miss",60,40.0,160.0);


  tree_WW->Project("h_WW_miss", "M_miss", cut_in);

  tree_ZZ->Project("h_ZZ_miss", "M_miss", cut_in);

  tree_HZ->Project("h_HZ_miss", "M_miss", cut_in);


  //Scaling histos and plotting final M_miss plot------------------------------
  h_WW_miss->Scale(Luminosity*1e12*scale_WW);
  h_WW_miss->SetLineColor(40);
  h_WW_miss->SetLineWidth(3);
  h_WW_miss->SetFillColor(0);
  h_WW_miss->SetFillStyle(0);

  h_ZZ_miss->Scale(Luminosity*1e12*scale_ZZ);
  h_ZZ_miss->SetLineColor(38);
  h_ZZ_miss->SetLineWidth(3);
  h_ZZ_miss->SetFillColor(0);
  h_ZZ_miss->SetFillStyle(0);

  h_HZ_miss->Scale(Luminosity*1e12*scale_HZ);
  h_HZ_miss->SetLineColor(46);
  h_HZ_miss->SetLineWidth(3);
  h_HZ_miss->SetFillColor(0);
  h_HZ_miss->SetFillStyle(0);

  h_sum_CMS->Add(h_WW_miss,h_ZZ_miss);
  h_sum_CMS->Add(h_HZ_miss);

  h_sum_miss->Add(h_WW_miss);
  h_sum_miss->Add(h_ZZ_miss);
  h_sum_miss->Add(h_HZ_miss);

  leg1->AddEntry(h_HZ_miss,Form("Signal CMS"),"l");
  leg1->AddEntry(h_ZZ_miss,Form("ZZ CMS"),"l");
  leg1->AddEntry(h_WW_miss,Form("WW CMS"),"l");

  return h_sum_miss;
}



void plot_M_miss_ll_comparison(){

  THStack *hs_ILD = plot_M_miss_ll_ILD();
  THStack *hs_CMS = plot_M_miss_ll();

  TLegend *leg = new TLegend(0.1380417,0.787844,0.6845907,0.891055,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextSize(0.03555046);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  TLegendEntry *entry=leg->AddEntry("NULL","FCC-ee, 3.5ab^{-1}, #sqrt{s}=240GeV, BR_{H->inv}=100%","h");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry->SetTextFont(42);

  //Draw on Canvas--------------------------------------------------------------
  c_M_miss->cd();
  hs_ILD->Draw();
  c_M_miss->SetGridx();
  c_M_miss->SetGridy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.07);
  hs_ILD->GetHistogram()->GetYaxis()->SetTitleOffset(1.6);
  hs_ILD->GetHistogram()->GetXaxis()->SetTitle("Missing Mass [GeV]");
  hs_ILD->GetHistogram()->GetYaxis()->SetTitle("Events / 2 GeV");
  hs_CMS->Draw("SAME");

  leg1->Draw();
  leg->Draw();

  c_M_miss->Modified();

  TCanvas* c_miss_tot = new TCanvas("c_miss_tot", "c_miss_tot",0,800,1250,900);
  c_miss_tot->SetGrid();
  gPad->SetLeftMargin(0.12);
  h_sum_ILD->SetXTitle("Missing Mass [GeV]");
  h_sum_ILD->SetYTitle("Events / 2 GeV");
  h_sum_ILD->GetYaxis()->SetTitleOffset(1.6);
  h_sum_ILD->SetLineColor(2);
  h_sum_ILD->SetLineWidth(3);
  h_sum_ILD->SetStats(0);
  h_sum_ILD->Draw();

  h_sum_CMS->SetXTitle("Missing Mass [GeV]");
  h_sum_CMS->SetYTitle("Events / 2 GeV");
  h_sum_CMS->SetLineColor(4);
  h_sum_CMS->SetLineWidth(3);
  h_sum_CMS->SetStats(0);
  h_sum_CMS->Draw("SAME");

  TLegend* leg2 = new TLegend(0.7407705,0.6983945,0.9911717,0.928899);
  //leg->SetHeader("Legend");
  leg2->AddEntry(h_sum_ILD,Form("ILD"),"l");
  leg2->AddEntry(h_sum_CMS,Form("CMS"),"l");
  leg2->Draw();
  leg->Draw();

  cout << h_sum_ILD->Integral() << endl << h_sum_CMS->Integral() << endl ;

}
