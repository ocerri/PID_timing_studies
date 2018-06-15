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
Int_t binMcms = 10; //Multiplies the number of bins
Int_t binMild = 10; //Multiplies the number of bins


//TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>76.2&&M_ll<106.2";
TString paper2012_cuts = "dtheta_ll_lab>100&&Pt_ll>10&&Pl_ll<50&&acoplanarity_angle>10&&M_ll>87&&M_ll<95"; // modified cuts


void plot_M_miss_tails_HZ(TString cuts = "paper2012"){

  if( cuts.Contains("paper2012") ){
    cuts.ReplaceAll("paper2012",paper2012_cuts);
  }

  const char* cut_in = cuts;


  TFile *f_HZ_ILD = TFile::Open("_root/tree_HZ_H2inv_Z2lep_2M_ILDdelphes_ll.root");
  TTree *tree_HZ_ILD = (TTree*) f_HZ_ILD->Get("tree_ll");

  TFile *f_HZ_CMS = TFile::Open("_root/tree_HZ_H2inv_Z2lep_1M_CMSdelphes_ll.root");
  TTree *tree_HZ_CMS = (TTree*) f_HZ_CMS->Get("tree_ll");

  Double_t scale_HZ = 5.720e-11/1e6; //Cross section divided total events generated

  //Fill the histos-------------------------------------------------------------
  TH1 *h_ILD = new TH1F("h_ILD","Missing mass in HZ and Z#rightarrow l^{+}l^{-} tagged events (M_{Z}+/- 4 GeV)",60*binMild,40.0,160.0 );
  TH1 *h_CMS = new TH1F("h_CMS","Missing mass in HZ and Z#rightarrow l^{+}l^{-} tagged events (M_{Z}+/- 4 GeV)",60*binMcms,40.0,160.0 );

  tree_HZ_ILD->Project("h_ILD", "M_miss", cut_in);
  tree_HZ_CMS->Project("h_CMS", "M_miss", cut_in);

  //Scaling histos and plotting final M_miss plot------------------------------
  h_ILD->Scale(Luminosity*0.5*1e12*scale_HZ);
  h_CMS->Scale(Luminosity*1e12*scale_HZ);

  //Draw total distribution-----------------------------------------------------

  h_ILD->Scale(1/h_ILD->Integral());
  h_CMS->Scale(1/h_CMS->Integral());


  TCanvas* c_miss_tot = new TCanvas("c_miss_tot", "c_miss_tot",0,800,1250,900);
  c_miss_tot->SetGrid();
  gPad->SetLeftMargin(0.12);
  h_ILD->SetXTitle("Missing Mass [GeV]");
  h_ILD->SetYTitle(Form("Events / %.2f GeV",2/(float)binMild));
  h_ILD->GetYaxis()->SetTitleOffset(1.6);
  h_ILD->SetLineColor(2);
  h_ILD->SetLineWidth(3);
  h_ILD->SetStats(0);
  h_ILD->Draw();

  h_CMS->SetXTitle("Missing Mass [GeV]");
  h_CMS->SetYTitle(Form("Events / %.2f GeV",2/(float)binMcms));
  h_CMS->SetLineColor(4);
  h_CMS->SetLineWidth(3);
  h_CMS->SetStats(0);
  h_CMS->Draw("SAME");

  TLegend* leg2 = new TLegend(0.7407705,0.6983945,0.9911717,0.928899);
  //leg->SetHeader("Legend");
  leg2->AddEntry(h_ILD,Form("ILD"),"l");
  leg2->AddEntry(h_CMS,Form("CMS"),"l");
  leg2->Draw();

  TLegend *leg = new TLegend(0.1380417,0.787844,0.6845907,0.891055,NULL,"brNDC");
  leg->SetBorderSize(1);
  leg->SetTextSize(0.03555046);
  leg->SetLineColor(0);
  leg->SetLineStyle(1);
  leg->SetLineWidth(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  TLegendEntry *entry=leg->AddEntry("NULL","3.5ab^{-1}, #sqrt{s}=240GeV, BR_{H->inv}=100%","h");
  entry->SetLineColor(1);
  entry->SetLineStyle(1);
  entry->SetLineWidth(1);
  entry->SetMarkerColor(1);
  entry->SetMarkerStyle(21);
  entry->SetMarkerSize(1);
  entry->SetTextFont(42);
  leg->Draw();

  cout << h_ILD->Integral() << endl << h_CMS->Integral() << endl ;

}
