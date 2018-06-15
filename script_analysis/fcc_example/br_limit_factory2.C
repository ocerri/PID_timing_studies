#include "br_limit_factory2.h"
Int_t binM = 10; //Multiply the number of bin

void br_limit_factory2(unsigned int trials = 100, Int_t seed_st = 0){
  std::ofstream logfile;
  logfile.open(Form("log/br_limit_factory2_%d_%d.log",seed_st,seed_st+trials-1));

  const char* cut_in = std_cuts;
  logfile << "Cuts: " << cut_in << std::endl << std::flush;

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
  TH1* h_WW_miss = new TH1D("h_WW_miss","h_WW_miss",60*binM,40.0,160.0);
  TH1* h_ZZ_miss = new TH1D("h_ZZ_miss","h_ZZ_miss",60*binM,40.0,160.0);
  TH1* h_HZ_miss = new TH1D("h_HZ_miss","h_HZ_miss",60*binM,40.0,160.0);

  tree_WW->Project("h_WW_miss", "M_miss", cut_in);
  tree_ZZ->Project("h_ZZ_miss", "M_miss", cut_in);
  tree_HZ->Project("h_HZ_miss", "M_miss", cut_in);

  //Scale histos----------------------------------------------------------------
  h_WW_miss->Scale(Luminosity*1e12*scale_WW);
  h_ZZ_miss->Scale(Luminosity*1e12*scale_ZZ);
  h_HZ_miss->Scale(Luminosity*1e12*scale_HZ);

  //Create bg and sg histos-----------------------------------------------------
  TH1 *h_sg = (TH1F*) h_HZ_miss->Clone("h_sg");
  h_sg->Sumw2();
  TH1 *h_bg = new TH1F("h_bg","h_bg",60*binM,40.0,160.0);
  h_bg->Add(h_ZZ_miss,h_WW_miss);
  h_bg->Sumw2();

  logfile << "Data loaded" << std::endl << std::flush;

  //Create BR distrib histo-----------------------------------------------------
  TH1 *h_br = new TH1F("h_br","BR limit distribution, CMS card", 50, 0.0,5.0);
  h_br->SetXTitle("BR [per cent]");
  h_br->SetYTitle("Number of toymodels");

  //Histo filling loop----------------------------------------------------------
  for (size_t i = 0; i < trials; i++) {
    logfile << "Iteration " << i+1 << "/" << trials << std::endl << std::flush;

    Double_t br_tmp = find_br_upper_limit(seed_st+i, h_bg, h_sg, h_WW_miss, h_ZZ_miss, h_HZ_miss);
    h_br->Fill(br_tmp);

  }

  TFile *f_out = TFile::Open(Form("_root/br_limit_factory2_%d_%d_%.0f.root",seed_st,seed_st+trials-1, Luminosity),"RECREATE");
  f_out->cd();
  h_br->Write();

  logfile << "Histo saved to file" << std::endl << std::flush;


  TLegend *leg = new TLegend(0.3318928,0.8514908,0.6110631,0.9008028,Form("#int Ldt=%.1f fb^{-1}, #sqrt{s}=240GeV, std cuts",Luminosity));
  leg->SetLineWidth(0);
  leg->SetLineColor(0);
  leg->SetFillStyle(0);

  TCanvas *c_br = new TCanvas("c_br","c_br",900,1200);
  c_br->cd();
  c_br->SetGrid();
  h_br->Draw();
  leg->Draw();

  logfile << "Canvas printed" << std::endl << std::flush;

  f_out->cd();
  c_br->Write();

}
