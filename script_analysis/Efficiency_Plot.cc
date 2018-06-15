#include <iostream>
#include <math.h>
#include <vector>

#include "TMath.h"
#include "TGraph.h"
#include "TCanvas.h"

using namespace std;

void Efficiency_Plot() {
  gStyle->SetTitleSize(.05,"X");//.055
  gStyle->SetTitleOffset(1.4,"X");//1.2,0.9
  gStyle->SetLabelSize(.05,"X");

  gStyle->SetTitleSize(.05,"Y");//.055
  gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(.05,"Y");

  gStyle->SetPadLeftMargin(.12);
  gStyle->SetPadRightMargin(.02);
  gStyle->SetPadBottomMargin(.15);
  gStyle->SetPadTopMargin(.02);

  vector<double> mass = {100, 200, 350, 500, 750, 1000, 1250, 1500, 2000};

  vector<double> eff2 = {1.67e-02, 7.45e-02, 1.39e-01, 1.67e-01, 1.81e-01, 1.85e-01, 1.87e-01, 1.92e-01, 1.94e-01};
  vector<double> eff1 = {9.31e-03, 5.21e-02, 1.52e-01, 2.42e-01, 3.45e-01, 4.02e-01, 4.33e-01, 4.57e-01, 4.81e-01};

  vector<double> eff_newT = {1.12e-01, 1.25e-01, 1.34e-01, 1.44e-01, 1.45e-01, 1.50e-01};
  double m_new[] = {100, 200, 350, 500, 750, 1000};

  auto gr1 = new TGraph(mass.size(), &mass[0], &eff1[0]);
  gr1->SetTitle("");
  gr1->SetMarkerStyle(4);
  gr1->SetLineStyle(7);
  gr1->SetLineWidth(2);
  gr1->SetMarkerColor(4);
  gr1->SetLineColor(4);


  auto gr2 = new TGraph(mass.size(), &mass[0], &eff2[0]);
  gr2->SetMarkerStyle(4);
  gr2->SetLineStyle(7);
  gr2->SetLineWidth(2);
  gr2->SetMarkerColor(2);
  gr2->SetLineColor(2);

  auto gr_new = new TGraph(eff_newT.size(), &m_new[0], &eff_newT[0]);
  gr_new->SetMarkerStyle(4);
  gr_new->SetLineStyle(7);
  gr_new->SetLineWidth(2);
  gr_new->SetMarkerColor(1);
  gr_new->SetLineColor(1);


  auto leg = new TLegend(0.6, 0.2, 0.98, 0.5);
  leg->AddEntry(gr1, "Single particle category", "lep");
  leg->AddEntry(gr2, "Two particles category", "lep");
  leg->AddEntry(gr_new, "TOF trigger", "lep");


  auto c = new TCanvas("c", "c", 800, 600);
  c->SetGrid();
  c->SetLogy();
  c->SetLogx();

  gr1->GetXaxis()->SetTitle("#tilde{t}_{1} Mass [GeV]");
  gr1->GetYaxis()->SetTitle("Selection efficiency");
  gr1->GetXaxis()->SetRangeUser(50, 3000);
  gr1->GetYaxis()->SetRangeUser(5e-3, 1);

  gr1->Draw("APC");
  gr2->Draw("PC");
  gr_new->Draw("PC");
  leg->Draw();
}
