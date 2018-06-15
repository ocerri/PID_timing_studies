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

TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10"; //NO M_ll cuts


//Create histos-----------------------------------------------------------------
Double_t Luminosity = 500; //fb^-1 come articolo
Double_t M_Z = 91.19;
Double_t M_H = 125.0;

//Merge code -------------------------------------------------------------------
void plot_process_miss(){

  //Decide cuts-----------------------------------------------------------------
  vector<Double_t> dmass = {90, 2, 5 ,10};
  vector<Double_t> dmiss = {100, 2, 5 ,10 ,15};

  //load data from trees--------------------------------------------------------
  // TFile *f = TFile::Open("_root/tree_WW_W2lep_1M_CMSdelphes_ll.root");
  // TTree *tree = (TTree*) f->Get("tree_ll");
  // Double_t scale = 1.715e-09/1e6; //Cross section divided total events generated
  // const char* ev_name = "WW";

  // TFile *f = TFile::Open("_root/tree_ZZ_Z2lep_1M_CMSdelphes_ll.root");
  // TTree *tree = (TTree*) f->Get("tree_ll");
  // Double_t scale = 1.370e-10/1e6; //Cross section divided total events generated
  // const char* ev_name = "ZZ";

  TFile *f = TFile::Open("_root/tree_HZ_H2inv_Z2lep_1M_CMSdelphes_ll.root");
  TTree *tree = (TTree*) f->Get("tree_ll");
  Double_t scale = 5.720e-11/1e6; //Cross section divided total events generated
  const char* ev_name = "HZ";

  //Fill the histos-------------------------------------------------------------
  TCanvas *c_miss = new TCanvas("c_WW_diff","c_WW_diff",0,45,1100,900);
  c_miss->cd();
  c_miss->SetGrid();
  c_miss->SetLeftMargin(0.12);
  Float_t ntot = (float)tree->Draw("M_miss", "1.");


  TLegend* leg = new TLegend(0.09900091,0.6534091,0.390554,0.9011364);

  for(unsigned short i=0;i<dmass.size()+1;++i){
    const char* title= Form("h_miss%d",i);
    TH1* h_mass = new TH1D(title,title,110,0.0,220.0);
    TString cuts;
    if(i!=dmass.size()){
      cuts = Form("M_ll>%.0f&&M_ll<%.0f", M_Z-dmass[i], M_Z+dmass[i]);
    }
    else{
      cuts = Form("M_ll>%.0f&&M_ll<%.0f&&paper2012", M_Z-4, M_Z+4);
    }
    if( cuts.Contains("paper2012") ){
      cuts.ReplaceAll("paper2012",paper2012_cuts);
    }
    const char* cut = cuts;

    Int_t n_surv = tree->Project(title, "M_miss", cut);
    h_mass->Scale(Luminosity*1e12*scale);
    h_mass->SetLineColor(1+i);
    if(i==4) h_mass->SetLineColor(6);
    h_mass->SetTitle(Form("Missing mass in %s events, L=%.0ffb^{-1}",ev_name, Luminosity));
    h_mass->GetXaxis()->SetTitle("Missing Mass [GeV]");
    h_mass->GetYaxis()->SetTitle("Events / 2 GeV");
    h_mass->GetYaxis()->SetTitleOffset(1.6);
    h_mass->SetStats(0);

    if(i==0){
      h_mass->DrawCopy();
    }
    else {
      h_mass->DrawCopy("same");
    }

    Double_t eff = 100* n_surv/ntot;

    if(i!=dmass.size()){
      leg->AddEntry(h_mass,Form("M_{ll} = M_{Z}#pm%.0f, eff=%.2f %%", dmass[i], eff),"l");
    }
    else{
      leg->AddEntry(h_mass,Form("M_{ll} = M_{Z}#pm 4 and paper2012, eff=%.2f %%", eff),"l");
    }

  }
  leg->Draw();


}
