#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"

#include "RooGausDoubleExp.cc"
#include "HypoTestInv.cc"

using namespace RooFit;

void Set_GausDoubleExpPdf_pars(RooWorkspace* w, double mass, TString label="", bool is_fix = true) {
  w->var("mu"+label)->setVal(mass*9.98e-01 + 1.67e-01);
  w->var("sigma"+label)->setVal(-2.89e-08*pow(mass,3) + 9.40e-05*mass*mass + -9.87e-03*mass + 4.15e+00);
  w->var("aL"+label)->setVal(2.16e-05*mass + 7.35e-01);
  w->var("aR"+label)->setVal(8.99e-07*mass + 7.69e-01);

  w->var("mu"+label)->setConstant(is_fix);
  w->var("sigma"+label)->setConstant(is_fix);
  w->var("aL"+label)->setConstant(is_fix);
  w->var("aR"+label)->setConstant(is_fix);
}

double GetEfficiency(double mass) {
  double val_m[] = {50., 100., 200., 350., 500., 750., 1000., 1500.};
  double val_eff[] = {0.00277778, 0.03997778, 0.19394, 0.42635556, 0.57726, 0.73302, 0.8118, 0.87264};

  auto gr = new TGraph(8, val_m, val_eff);

  return TMath::Max(0., (double)gr->Eval(mass, 0, "S"));
}

TCanvas* DrawExampleMassSpectrum(int nsig = 1500, double lumi = 1e6, double test_mass = 500) {
  int nbkg = (int) (3e7 * lumi/1.0e6);
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  double binning[] = {200, 35, 2000};

  RooWorkspace* w = new RooWorkspace("w");
  w->addClassDeclImportDir("/Users/olmo/cernbox/PID_timing_studies/script_analysis");
  w->importClassCode("RooGausDoubleExp",kTRUE);

  // Create the fitting model
  // w->factory(Form("Exponential:bkg_pdf(Mass[%f,%f], lambda[-0.025,-5,-0.001])",binning[1], binning[2]));
  w->factory(Form("Gamma:bkg_pdf(Mass[%f,%f], gamma[0.000, 0.000, 0.00], lambda[64.4899, 30, 80], mu_Gamma[0., 0., 0.])",binning[1], binning[2]));
  w->factory("GausDoubleExp:sig_pdf(Mass, mu[200., 100, 4000], sigma[30., 10, 500], aL[2, 0.1, 3], aR[1.5, 0.1, 3])");
  w->factory("SUM:model(nsig[0,10000]*sig_pdf, nbkg[3000,0,1000000000]*bkg_pdf)");  // for extended model

  // Settina names
  w->var("Mass")->SetTitle("Mass [GeV]");


  //=========== Generate the data
  Set_GausDoubleExpPdf_pars(w, test_mass, "");
  RooRandom::randomGenerator()->SetSeed(111);

  RooRealVar * m = w->var("Mass");  // the observable
  m->setBins(binning[0]);
  // set the desired value of signal and background events
  w->var("nsig")->setVal(nsig);
  w->var("nbkg")->setVal(nbkg);
  RooDataHist * data = w->pdf("model")->generateBinned(*m);
  data->SetName("data");
  w->import(*data);

  data->Print();

  // Fit the bkg to the exponential
  RooFitResult * r = w->pdf("model")->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  r->Print();

  RooPlot * frame = m->frame(Title(Form("Mass spectrum of QCD + RHad")));
  data->plotOn(frame, Name("Pseudo-data"));
  w->pdf("model")->plotOn(frame, Name("S+B model"));
  w->pdf("model")->plotOn(frame, RooFit::Components("bkg_pdf"), RooFit::LineStyle(kDashed), Name("Bkg component"));
  w->pdf("model")->plotOn(frame, RooFit::Components("sig_pdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), Name("Sig component") );
  frame->GetYaxis()->SetRangeUser(1, 6e7);
  frame->GetYaxis()->SetTitleSize(.05);
  frame->GetYaxis()->SetLabelSize(.05);
  frame->GetXaxis()->SetLimits(0.2*test_mass, 2*test_mass);
  frame->GetXaxis()->SetTitleSize(.05);
  frame->GetXaxis()->SetLabelSize(.05);



  // auto h_pull_bkg = frame->residHist(0,0,true);
  // auto frame2 = m->frame(Title(""));
  // frame2->addPlotable(h_pull_bkg, "P");

  TCanvas* c_out = new TCanvas("c_out", "c_out", 800, 600);
  // c_out->Divide(1,2);
  //
  // c_out->cd(1);
  frame->Draw();
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->SetLeftMargin(.13);
	gPad->SetRightMargin(.05);
  gPad->SetBottomMargin(.14);
	gPad->SetTopMargin(.07);

  // c_out->cd(2);
  // frame2->Draw();

  auto leg = new TLegend(0.7,0.7,0.95,0.93);
  leg->AddEntry(frame->findObject("Pseudo-data"), "Pseudo-data", "lep");
  leg->AddEntry(frame->findObject("S+B model"), "S+B model", "l");
  leg->AddEntry(frame->findObject("Bkg component"), "Bkg component", "l");
  leg->AddEntry(frame->findObject("Sig component"), "Sig component", "l");
  leg->Draw();

  auto note = new TLatex();
  double xsec = nsig / (GetEfficiency(test_mass) * lumi);
  note->DrawLatexNDC(0.18, 0.85, Form("#splitline{#sqrt{s} 14TeV, L = 1ab^{-1}}{M_{#tilde{t}_{1}} = %.0f GeV,  xsec = %.2e fb}", test_mass, xsec));

  return c_out;
}

RooWorkspace* ModelTOFAnalysis_RHad(int nsig = 0, int nbkg = 0 )
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  double binning[] = {200, 35, 2000};

  RooWorkspace* w = new RooWorkspace("w");
  w->addClassDeclImportDir("/Users/olmo/cernbox/PID_timing_studies/script_analysis");
  w->importClassCode("RooGausDoubleExp",kTRUE);

  // Create the fitting model
  // w->factory(Form("Exponential:bkg_pdf(Mass[%f,%f], lambda[-0.025,-5,-0.001])",binning[1], binning[2]));
  w->factory(Form("Gamma:bkg_pdf(Mass[%f,%f], gamma[1.e-6, 1.e-6, 1.e-6], lambda[64.4899, 30, 90], mu_Gamma[0., 0., 0.])",binning[1], binning[2]));
  w->factory("GausDoubleExp:sig_pdf(Mass, mu[200., 100, 4000], sigma[30., 10, 500], aL[2, 0.1, 3], aR[1.5, 0.1, 3])");
  w->factory("SUM:model(nsig[0,10000]*sig_pdf, nbkg[3000,0,1000000000]*bkg_pdf)");  // for extended model

  // Create the data distribution

  TFile* f_in = TFile::Open("/Users/olmo/cernbox/PID_timing_studies/_root/histos/RecoMassTree_bkg.root");
  TTree* tree = (TTree*) f_in->Get("T");
  auto bkg_MC = new RooDataSet("MC_bkg", "MC_bkg", tree, RooArgSet(*w->var("Mass")));
  // auto bkg_MC_pdf = new RooHistPdf("bkg_MC_pdf", "bkg_MC_pdf", *w->var("Mass"), *bkg_MC, 2);
  //
  // w->import(*bkg_MC_pdf);
  // w->factory("GausDoubleExp:sig_MC_pdf(Mass, mu_MC[200.], sigma_MC[30.], aL_MC[2], aR_MC[1.5])");
  // w->factory("SUM:data_pdf(nsig_MC[0,10000]*sig_MC_pdf, nbkg_MC[0,10000]*bkg_MC_pdf)");  // for extended model

  // Settina names
  w->var("Mass")->SetTitle("Mass [GeV]");


  //=========== Generate the data
  double test_mass = 300.;
  Set_GausDoubleExpPdf_pars(w, test_mass, "");
  // use fixed random numbers for reproducibility (use 0 for changing every time)
  RooRandom::randomGenerator()->SetSeed(111);

  RooRealVar * m = w->var("Mass");  // the observable
  m->setBins(binning[0]);

  // Fit the bkg to the exponential
  // RooFitResult * r_bkg = w->pdf("bkg_pdf")->fitTo(*bkg_MC, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  // r_bkg->Print();
  //
  // TCanvas* c_bkg = new TCanvas("c_bkg", "c_bkg", 800, 600);
  // RooPlot * plot_bkg = m->frame(Title("Mass spectrum of QCD "));
  // bkg_MC->plotOn(plot_bkg);
  // w->pdf("bkg_pdf")->plotOn(plot_bkg, RooFit::LineStyle(kDashed));
  // plot_bkg->Draw();
  // c_bkg->SetGrid();


  // set the desired value of signal and background events
  w->var("nsig")->setVal(0);
  w->var("nbkg")->setVal(nbkg);

  // RooDataSet * data = w->pdf("model")->generate(*m);
  RooDataHist * data = w->pdf("model")->generateBinned(*m);
  data->SetName("data");
  w->import(*data);

  data->Print();

  TCanvas* c_out = new TCanvas("c_out", "c_out", 800, 600);
  RooPlot * plot = m->frame(Title(Form("Mass spectrum of QCD + RHad (%.0f GeV)", test_mass)));
  data->plotOn(plot);
  plot->Draw();

  // RooFitResult * r = w->pdf("model")->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  // r->Print();

  w->pdf("model")->plotOn(plot);
  w->pdf("model")->plotOn(plot, RooFit::Components("bkg_pdf"), RooFit::LineStyle(kDashed) );
  w->pdf("model")->plotOn(plot, RooFit::Components("sig_pdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed) );

  w->pdf("model")->paramOn(plot,Layout(0.7,0.99,0.95));

  plot->Draw();

  c_out->SetLogy();
  c_out->SetGrid();


  RooStats::ModelConfig mc("ModelConfig",w);
  mc.SetPdf(*w->pdf("model"));
  mc.SetParametersOfInterest(*w->var("nsig"));
  mc.SetObservables(*w->var("Mass"));
  // define set of nuisance parameters
  w->defineSet("nuisParams","lambda,nbkg");

  mc.SetNuisanceParameters(*w->set("nuisParams"));

  // import model in the workspace
  w->import(mc);

  // write the workspace in the file
  TString fileName = "/Users/olmo/cernbox/PID_timing_studies/_root/histos/ModelRHad.root";
  w->writeToFile(fileName,true);
  cout << "model written to file " << fileName << endl;

  return w;
}

TCanvas* MakeBrazilPlot(vector<double> RHad_mass_list, map<int, vector<double>> exp_upper_limit, double lumi = 1e6) {
  // double lumi = 12.; //[pb^-1]
  // double lumi = 1e6; //[pb^-1]

  auto gr = new TGraph(RHad_mass_list.size());
  auto gr_1s = new TGraphAsymmErrors(RHad_mass_list.size());
  auto gr_2s = new TGraphAsymmErrors(RHad_mass_list.size());

  for(unsigned int i=0; i<RHad_mass_list.size(); i++) {
    double mass = RHad_mass_list[i];
    double eff = GetEfficiency(mass);

    double scale = 1./(eff*lumi);
    gr->SetPoint(i, mass, exp_upper_limit[0][i]*scale);

    gr_1s->SetPoint(i, mass, exp_upper_limit[0][i]*scale);
    double aux_l = (exp_upper_limit[0][i] - exp_upper_limit[-1][i])*scale;
    double aux_u = (exp_upper_limit[1][i] - exp_upper_limit[0][i])*scale;
    gr_1s->SetPointError(i, 0, 0, aux_l, aux_u);

    gr_2s->SetPoint(i, mass, exp_upper_limit[0][i]*scale);
    aux_l = (exp_upper_limit[0][i] - exp_upper_limit[-2][i])*scale;
    aux_u = (exp_upper_limit[2][i] - exp_upper_limit[0][i])*scale;
    gr_2s->SetPointError(i, 0, 0, aux_l, aux_u);

    cout << " =============== Mass " << mass << " ========================" << endl;
    cout << Form("NSig: %.1f %.1f %.1f %.1f %.1f", exp_upper_limit[-2][i], exp_upper_limit[-1][i], exp_upper_limit[0][i], exp_upper_limit[1][i], exp_upper_limit[2][i]) << endl;
    cout << "Eff: " << eff << endl;
    cout << "Scale: " << scale << endl;
  }

  auto c_mass = new TCanvas("c_mass", "c_mass", 800, 600);
  gr_2s->SetTitle("");
  gr_2s->SetFillColor(5);
  gr_2s->SetFillStyle(1001);
  gr_2s->Draw("A3");
  gr_2s->GetXaxis()->SetTitle("#tilde{t}_{1} mass [GeV]");
  gr_2s->GetYaxis()->SetTitle("Excluded xsec @95% CLs [pb]");
  gr_2s->GetYaxis()->SetRangeUser(1e-7, 1);
  gr_2s->GetYaxis()->SetTitleSize(.05);
	gr_2s->GetYaxis()->SetTitleOffset(1.25);
	gr_2s->GetYaxis()->SetLabelSize(.05);
  gr_2s->GetXaxis()->SetTitleSize(.05);
	gr_2s->GetXaxis()->SetTitleOffset(1.2);
	gr_2s->GetXaxis()->SetLabelSize(.05);

  gr_1s->SetFillColor(8);
  gr_1s->SetFillStyle(1001);
  gr_1s->Draw("3");

  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  gr->SetLineColor(1);
  gr->SetLineStyle(7);
  gr->SetLineWidth(2);
  gr->Draw("LP");



  // Add current CMS limits
  double masses[] = { 99.1816, 398.5258, 599.9218, 795.8401, 1195.8373, 1590.3867, 2001.2604, 2393.0807 };
  double xsec_up[] = { 0.0107, 4.68e-3, 1.4103e-3, 1.4103e-3, 1.5569e-3, 2.0691e-3, 2.9617e-3, 5.1667e-3 };
  auto gr_CMS = new TGraph(8, masses, xsec_up);
  gr_CMS->SetLineColor(4);
  gr_CMS->SetLineWidth(2);
  gr_CMS->SetMarkerColor(4);
  gr_CMS->SetMarkerStyle(20);
  gr_CMS->Draw("CP");

  for(uint i = 0; i < 8; i++) {
    xsec_up[i] *= sqrt(12./1000.);
  }
  auto gr_CMS_exp = new TGraph(8, masses, xsec_up);
  gr_CMS_exp->SetLineColor(2);
  gr_CMS_exp->SetLineWidth(2);
  gr_CMS_exp->SetMarkerColor(2);
  gr_CMS_exp->SetMarkerStyle(20);
  gr_CMS_exp->Draw("CP");

  auto leg = new TLegend(0.6, 0.68, 0.95, 0.95);
  leg->AddEntry(gr_CMS, "CMS EXO-16-36 (12 fb^{-1})", "lp");
  leg->AddEntry(gr_CMS_exp, "CMS expected @ 1 ab^{-1}", "lp");
  leg->AddEntry(gr, "Expected median", "lp");
  leg->AddEntry(gr_1s, "Expected #pm 1 #sigma", "F");
  leg->AddEntry(gr_2s, "Expected #pm 2 #sigma", "F");
  leg->Draw();

  auto note = new TLatex();
  note->DrawLatexNDC(0.2, 0.88, Form("#splitline{HL-LHC Simulation}{#sqrt{s} 14TeV, L = %.0f fb^{-1}}", lumi*1e-3));

  c_mass->SetGrid();
  c_mass->SetLogy();
  gPad->SetLeftMargin(.13);
	gPad->SetRightMargin(.05);
  gPad->SetBottomMargin(.14);
	gPad->SetTopMargin(.05);

  return c_mass;
}

void TOFAnalysis_RHad(){
  // double Luminosity = 1e6; //[pb^-1]
  double Luminosity = 12e3; //[pb^-1]
  // DrawExampleMassSpectrum();
  // return;

  int nbkg = 3.39e7 * Luminosity * 1e-6;

  cout << "Building the model" << endl;
  auto w = ModelTOFAnalysis_RHad(0, nbkg);

  cout << "Computing the exclusion limit\n";
  vector<double> RHad_mass_list;
  vector<int> sigma_eval = {-2, -1, 0, 1, 2};
  map<int, vector<double>> exp_upper_limit;
  for( auto s : sigma_eval ) {
    exp_upper_limit[s] = {};
  }


  for (double mass = 150; mass <= 1700; mass += 150) {
    Set_GausDoubleExpPdf_pars(w, mass);

    HypoTestInvTool calc(&optHTInv);
    int nsig_max = 15 + 4e3*exp(-mass/150.);
    cout << "=========== NSIG max = " << nsig_max << endl;
    HypoTestInverterResult * r = calc.RunInverter(w, "ModelConfig", "",
                                                  "data", 2, 3, true,
                                                  30, 0, nsig_max,
                                                  1000);

    calc.AnalyzeResult( r, 2, 3, true, 10, "/Users/olmo/cernbox/PID_timing_studies/_root/results/", Form("%.0f", mass));

    RHad_mass_list.push_back(mass);
    for( auto s : sigma_eval ) {
      exp_upper_limit[s].push_back(r->GetExpectedUpperLimit(s));
    }
  }

  auto c_mass = MakeBrazilPlot(RHad_mass_list, exp_upper_limit, Luminosity);
  c_mass->SaveAs("~/Desktop/ExclusionLimits.root");

}
