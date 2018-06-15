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


//Create histos-----------------------------------------------------------------
//Double_t Luminosity = 500; //fb^-1 come articolo
Double_t Luminosity = 3500; //FCCee

THStack *h_sum_miss = new THStack("h_sum_miss","Missing Mass in Z#rightarrow l^{+}l^{-} tagged events (M_{Z}+/- 4 GeV)");
//TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>76.2&&M_ll<106.2";
TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>87&&M_ll<95"; // modified cuts

//Merge code -------------------------------------------------------------------
void plot_M_miss_ll(TString cuts = "paper2012", const char* outtag="cuts"){

  if( cuts.Contains("paper2012") ){
    cuts.ReplaceAll("paper2012",paper2012_cuts);
  }

  const char* cut_in = cuts;

  //Create file in to store data
  TMacro *infos = new TMacro("infos","infos");
  infos->AddLine("Applied cuts:");
  infos->AddLine(cut_in);

  //load data from trees--------------------------------------------------------
  TFile *f_WW = TFile::Open("_root/tree_WW_W2lep_2M_CMSdelphes_ll.root");
  TTree *tree_WW = (TTree*) f_WW->Get("tree_ll");
  Double_t scale_WW = 1.715e-09/1e6; //Cross section divided total events generated

  TFile *f_ZZ = TFile::Open("_root/tree_ZZ_Z2lep_2M_CMSdelphes_ll.root");
  TTree *tree_ZZ = (TTree*) f_ZZ->Get("tree_ll");
  Double_t scale_ZZ = 1.370e-10/1e6; //Cross section divided total events generated

  TFile *f_HZ = TFile::Open("_root/tree_HZ_H2inv_Z2lep_2M_CMSdelphes_ll.root");
  TTree *tree_HZ = (TTree*) f_HZ->Get("tree_ll");
  Double_t scale_HZ = 5.720e-11/1e6; //Cross section divided total events generated

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

  //Create Canvas and histos before and after cuts------------------------------
  TCanvas *c_WW_diff = new TCanvas("c_WW_diff","c_WW_diff",600,600);
  c_WW_diff->cd();
  c_WW_diff->SetGrid();
  Int_t n_tot_WW = tree_WW->Draw("M_miss", "1.");
  h_WW_miss->SetLineColor(2);
  h_WW_miss->DrawCopy("same");

  TCanvas *c_ZZ_diff = new TCanvas("c_ZZ_diff","c_ZZ_diff",600,600);
  c_ZZ_diff->cd();
  c_ZZ_diff->SetGrid();
  Int_t n_tot_ZZ = tree_ZZ->Draw("M_miss", "1.");
  h_ZZ_miss->SetLineColor(2);
  h_ZZ_miss->DrawCopy("same");

  TCanvas *c_HZ_diff = new TCanvas("c_HZ_diff","c_HZ_diff",600,600);
  c_HZ_diff->cd();
  c_HZ_diff->SetGrid();
  Int_t n_tot_HZ = tree_HZ->Draw("M_miss", "1.");
  h_HZ_miss->SetLineColor(2);
  h_HZ_miss->DrawCopy("same");

  //Computing efficency & printig output----------------------------------------
  cout << endl <<"Requested cut : " << cut_in << endl << endl;
  cout << "WW efficency = " << 100.0*(float)n_surv_WW/n_tot_WW << " %"<< endl;
  infos->AddLine(Form("WW efficency = %.3f %%",100.0*(float)n_surv_WW/n_tot_WW));
  cout << "ZZ efficency = " << 100.0*(float)n_surv_ZZ/n_tot_ZZ << " %"<< endl;
  infos->AddLine(Form("ZZ efficency = %.3f %%",100.0*(float)n_surv_ZZ/n_tot_ZZ));
  cout << "HZ efficency = " << 100.0*(float)n_surv_HZ/n_tot_HZ << " %"<< endl;
  infos->AddLine(Form("HZ efficency = %.3f %%",100.0*(float)n_surv_HZ/n_tot_HZ));

  //Scaling histos and plotting final M_miss plot------------------------------
  h_WW_miss->Scale(Luminosity*1e12*scale_WW);
  h_WW_miss->SetLineColor(6);
  h_WW_miss->SetLineWidth(3);
  h_WW_miss->SetFillColor(1);
  h_WW_miss->SetFillStyle(3004);

  h_ZZ_miss->Scale(Luminosity*1e12*scale_ZZ);
  h_ZZ_miss->SetLineColor(4);
  h_ZZ_miss->SetLineWidth(3);
  h_ZZ_miss->SetFillColor(1);
  h_ZZ_miss->SetFillStyle(3005);

  h_HZ_miss->Scale(Luminosity*1e12*scale_HZ);
  h_HZ_miss->SetLineColor(2);
  h_HZ_miss->SetLineWidth(3);
  h_HZ_miss->SetFillColor(10);

  h_sum_miss->Add(h_WW_miss);
  h_sum_miss->Add(h_ZZ_miss);
  h_sum_miss->Add(h_HZ_miss);

  //Draw on Canvas--------------------------------------------------------------
  TCanvas* c_M_miss = new TCanvas("c_M_miss", "c_M_miss",0,800,1250,900);
  c_M_miss->cd();
  h_sum_miss->Draw();
  c_M_miss->SetGridx();
  c_M_miss->SetGridy();
  gPad->SetLeftMargin(0.12);
  gPad->SetRightMargin(0.07);
  h_sum_miss->GetHistogram()->GetYaxis()->SetTitleOffset(1.6);
  h_sum_miss->GetHistogram()->GetXaxis()->SetTitle("Missing Mass [GeV]");
  h_sum_miss->GetHistogram()->GetYaxis()->SetTitle("Events / 2 GeV");

  TLegend* leg = new TLegend(0.1495246,0.6616972,0.3820225,0.8509174);
  //leg->SetHeader("Legend");
  leg->AddEntry(h_HZ_miss,Form("Signal"),"l");
  leg->AddEntry(h_ZZ_miss,Form("ZZ"),"l");
  leg->AddEntry(h_WW_miss,Form("WW"),"l");
  leg->AddEntry(h_ZZ_miss,Form("All backgrounds"),"f");
  leg->Draw();

  TLegend* leg1 = new TLegend(0.1495246,0.8486239,0.5294555,0.9002294,"FCC-ee, 3.5ab^{-1}, #sqrt{s}=240GeV, BR_{H->inv}=100%");
  leg1->Draw();

  // the following line is needed to avoid that the automatic redrawing of stats
  h_ZZ_miss->SetStats(0);
  h_HZ_miss->SetStats(0);
  h_WW_miss->SetStats(0);

  c_M_miss->Modified();

  //Saving all plots and histo to .root-----------------------------------------
  TH1 *h_signal = (TH1F*) h_HZ_miss->Clone("h_signal");
  h_signal->Sumw2();
  TH1 *h_background = new TH1F("h_background","h_background",60,40.0,160.0);
  h_background->Add(h_ZZ_miss,h_WW_miss);
  h_background->Sumw2();

  Double_t sig_ev = h_signal->Integral(41,43);
  Double_t bg_ev = h_background->Integral(41,43);
  Double_t br_limit = 100*TMath::Sqrt(bg_ev)/sig_ev;
  cout << "Limit on H2inv BR (raw) : " << br_limit << endl;
  infos->AddLine(Form("Limit on H2inv BR (raw) : %.2f", br_limit));


  char *out_name = Form("_root/plot_M_miss_%s.root",outtag);
  TFile *f_out = TFile::Open(out_name,"RECREATE");
  f_out->cd();
  c_M_miss->Write();
  c_WW_diff->Write();
  c_ZZ_diff->Write();
  c_HZ_diff->Write();
  h_ZZ_miss->Write();
  h_HZ_miss->Write();
  h_WW_miss->Write();
  h_background->Write();
  h_signal->Write();
  infos->Write();

  cout << endl << "Output file saved in : < " << out_name << " >" << endl << endl;



}
