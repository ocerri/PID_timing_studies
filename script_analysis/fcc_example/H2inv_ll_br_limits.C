#include "H2inv_ll_br_limits.h"

void H2inv_ll_br_limits(Double_t *br_limit_out, Double_t *br_sign5_out, TString cuts = "", Int_t save_output=0, const char* outtag="cuts"){

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

  Int_t n_surv_WW;
  n_surv_WW = tree_WW->Project("h_WW_miss", "M_miss", cut_in);

  Int_t n_surv_ZZ;
  n_surv_ZZ = tree_ZZ->Project("h_ZZ_miss", "M_miss", cut_in);

  Int_t n_surv_HZ;
  n_surv_HZ = tree_HZ->Project("h_HZ_miss", "M_miss", cut_in);

  TCanvas c_WW_diff("c_WW_diff","c_WW_diff",600,600);
  TCanvas c_ZZ_diff("c_ZZ_diff","c_ZZ_diff",600,600);
  TCanvas c_HZ_diff("c_HZ_diff","c_HZ_diff",600,600);
  TMacro *infos = new TMacro("infos","infos");
  if(save_output==1){

    c_WW_diff.cd();
    c_WW_diff.SetGrid();
    Int_t n_tot_WW = tree_WW->Draw("M_miss", "1.");
    h_WW_miss->SetLineColor(2);
    h_WW_miss->DrawCopy("same");

    c_ZZ_diff.cd();
    c_ZZ_diff.SetGrid();
    Int_t n_tot_ZZ = tree_ZZ->Draw("M_miss", "1.");
    h_ZZ_miss->SetLineColor(2);
    h_ZZ_miss->DrawCopy("same");

    c_HZ_diff.cd();
    c_HZ_diff.SetGrid();
    Int_t n_tot_HZ = tree_HZ->Draw("M_miss", "1.");
    h_HZ_miss->SetLineColor(2);
    h_HZ_miss->DrawCopy("same");

    infos->AddLine("Applied cuts:");
    infos->AddLine(cut_in);
    infos->AddLine(Form("WW efficency = %.3f %%",100.0*(float)n_surv_WW/n_tot_WW));
    infos->AddLine(Form("ZZ efficency = %.3f %%",100.0*(float)n_surv_ZZ/n_tot_ZZ));
    infos->AddLine(Form("HZ efficency = %.3f %%",100.0*(float)n_surv_HZ/n_tot_HZ));
  }

  //Scaling histos and building final M_miss plot------------------------------
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

  THStack *h_sum_miss = new THStack("h_sum_miss","Missing Mass in Z#rightarrow l^{+}l^{-} tagged events");
  TCanvas c_M_miss("c_M_miss", "c_M_miss",0,800,1161,900);
  TLegend* leg = new TLegend(0.0997151,0.7118012,0.3323837,0.9006211);
  TLegend* leg1 = new TLegend(0.3318928,0.8514908,0.6110631,0.9008028,Form("FCC-ee, %.1f fb^{-1}, #sqrt{s}=240GeV, BR_{H->inv}=100%",Luminosity));


  if (save_output==1){
    h_sum_miss->Add(h_WW_miss);
    h_sum_miss->Add(h_ZZ_miss);
    h_sum_miss->Add(h_HZ_miss);

    c_M_miss.cd();
    gPad->SetLeftMargin(0.15);
    h_sum_miss->Draw();
    c_M_miss.SetGrid();
    h_sum_miss->GetHistogram()->GetXaxis()->SetTitle("Missing Mass [GeV]");
    h_sum_miss->GetHistogram()->GetYaxis()->SetTitle("Events / 2 GeV");
    h_sum_miss->GetHistogram()->GetYaxis()->SetTitleOffset(1.6) ;

    leg->SetHeader("Legend");
    leg->AddEntry(h_HZ_miss,Form("Signal"),"l");
    leg->AddEntry(h_ZZ_miss,Form("ZZ"),"l");
    leg->AddEntry(h_WW_miss,Form("WW"),"l");
    leg->AddEntry(h_ZZ_miss,Form("All backgrounds"),"f");
    leg->Draw();

    leg1->Draw();

    // the following line is needed to avoid that the automatic redrawing of stats
    h_ZZ_miss->SetStats(0);
    h_HZ_miss->SetStats(0);
    h_WW_miss->SetStats(0);

    c_M_miss.Modified();

  }

  TH1 *h_sg = (TH1F*) h_HZ_miss->Clone("h_sg");
  h_sg->Sumw2();
  TH1 *h_bg = new TH1F("h_bg","h_bg",60,40.0,160.0);
  h_bg->Add(h_ZZ_miss,h_WW_miss);
  h_bg->Sumw2();


  TGraph **g_nll_temp = 0;
  TGraph *g_nll_tmp =0;
  g_nll_temp = &g_nll_tmp;
  TGraphErrors **g_sig_temp = 0;
  TGraphErrors *g_sig_tmp = 0;
  g_sig_temp = &g_sig_tmp;
  Double_t br_upper_limit = find_br_upper_limit(time(NULL), h_bg, h_sg, h_WW_miss, h_ZZ_miss, h_HZ_miss, g_nll_temp, save_output);

  int toy_count;
  int points;
  if (save_output==0){
    toy_count = 10;
    points = 10;
  }
  else {
    toy_count = 10;
    points = 10;
  }
  Double_t br_significance_limit_5 = find_significance_limit(time(NULL),h_bg, h_sg, h_WW_miss, h_ZZ_miss, h_HZ_miss, g_sig_temp, save_output,toy_count,points);

  std::cout << std::endl << "Br limit = " << br_upper_limit << " %"  << std::endl;
  std::cout << "Br for 5 sigma significance = " << br_significance_limit_5 << " %" << std::endl << std::endl;


  *br_limit_out=br_upper_limit;
  *br_sign5_out=br_significance_limit_5;


  if (save_output==1){
    Double_t sig_ev = h_sg->Integral(41,43);
    Double_t bg_ev = h_bg->Integral(41,43);
    Double_t br_limit = 100*TMath::Sqrt(bg_ev)/sig_ev;
    infos->AddLine(Form("Luminosity : %.2f",Luminosity));
    infos->AddLine("Extended maximum likelihood fit, with separate template pdf for signal, ZZ background and WW background");
    infos->AddLine(Form("Limit on H2inv BR (raw) : %.2f", br_limit));
    infos->AddLine(Form("Limit on H2inv BR (1 toy model): %.2f", br_upper_limit));
    infos->AddLine(Form("H2inv BR for 5 sigma significance (%d toy models): %.2f", toy_count, br_significance_limit_5));
    char *out_name = Form("_root/H2inv_ll_br_limits_%s.root",outtag);
    TFile *f_out = TFile::Open(out_name,"RECREATE");
    f_out->cd();
    c_M_miss.Write();
    c_WW_diff.Write();
    c_ZZ_diff.Write();
    c_HZ_diff.Write();
    h_ZZ_miss->Write();
    h_HZ_miss->Write();
    h_WW_miss->Write();
    h_bg->Write();
    h_sg->Write();
    g_nll_tmp->Write();
    g_sig_tmp->Write();
    infos->Write();
  }


}
