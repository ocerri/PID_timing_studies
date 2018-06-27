#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"

void PlotEfficiency() {
    auto f_HT = TFile::Open("_fig/HTtrigger/SigmaEff.root", "READ");
    auto g1_HT = (TGraphErrors*) f_HT->Get("g_e1");
    auto g2_HT = (TGraphErrors*) f_HT->Get("g_e2");
    auto f_TOF = TFile::Open("_fig/TOFtrigger/SigmaEff.root", "READ");
    auto g1_TOF = (TGraphErrors*) f_TOF->Get("g_e1");
    g1_TOF->SetLineStyle(9);
    g1_TOF->SetMarkerStyle(25);
    auto g2_TOF = (TGraphErrors*) f_TOF->Get("g_e2");
    g2_TOF->SetLineStyle(9);
    g2_TOF->SetMarkerStyle(25);

    auto leg = new TLegend(0.8, 0.2, 0.98, 0.48);
    leg->AddEntry(g1_HT, "Single particle category, HT trigger", "lep");
    leg->AddEntry(g1_TOF, "Single particle category, TOF trigger", "lep");
    leg->AddEntry(g2_HT, "Two particles category, HT trigger", "lep");
    leg->AddEntry(g2_TOF, "Two particles category, TOF trigger", "lep");

    auto c = new TCanvas("c", "c", 800, 600);
    c->SetLogy();
    c->SetGrid();
    g1_HT->Draw("APC");
    g1_HT->GetYaxis()->SetRangeUser(5e-3, 1);
    g2_HT->Draw("PC");
    g1_TOF->Draw("PC");
    g2_TOF->Draw("PC");
    leg->Draw();

    gPad->SetLeftMargin(.13);
  	gPad->SetRightMargin(.05);
    gPad->SetBottomMargin(.14);
  	gPad->SetTopMargin(.05);
}
