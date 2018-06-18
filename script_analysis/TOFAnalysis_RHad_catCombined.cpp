#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"

#include "RooGausDoubleExp.cc"
#include "RooPowerLaw.cc"
#include "HypoTestInv.cc"

using namespace RooFit;

double test_mass = 400;
double test_xsec = 100.; //[fb]

double Lumi_Nbkg = 1.23e+01; //[fb-1]
double Nevs_bkg_HT[] = {2.24e+05, 7.64e+03};
double Mstar_bkg_HT[] = {67.0, 22.9};
double Nevs_bkg_TOF[] = {4.47e+05, 2.00e+04};
double Mstar_bkg_TOF[] = {66.6, 24.9};

void Set_GausDoubleExpPdf_pars(RooWorkspace* w, double mass, bool TOFtrigger = false) {
  vector<double> M_st = {100.000, 200.000, 350.000, 500.000, 750.000, 1000.000, 1250.000, 1500.000, 2000.000};

  double aL_HT1[] = {0.579, 0.644, 0.780, 0.781, 0.894, 0.724, 0.698, 0.812, 0.850};
  double aR_HT1[] = {0.500, 0.588, 0.719, 0.669, 0.802, 0.644, 0.572, 0.655, 0.661};
  double mu_HT1[] = {99.953, 199.664, 349.487, 498.793, 748.687, 996.947, 1242.569, 1493.386, 1983.368};
  double sigma_HT1[] = {3.487, 6.181, 13.705, 21.230, 41.974, 54.583, 71.972, 108.938, 171.725};

  double aL_HT2[] = {1.029, 0.908, 0.934, 0.897, 0.885, 0.908, 0.992, 0.970, 1.076};
  double aR_HT2[] = {1.126, 0.980, 0.960, 0.966, 0.958, 1.036, 0.973, 0.954, 1.055};
  double mu_HT2[] = {100.046, 200.032, 349.727, 499.744, 748.320, 998.664, 1244.701, 1493.342, 1986.155};
  double sigma_HT2[] = {2.914, 4.730, 9.363, 15.778, 28.561, 46.456, 61.688, 82.566, 131.018};

  double aL_TOF1[] = {0.534, 0.510, 0.607, 0.625, 0.683, 0.582, 0.592, 0.672, 0.717};
  double aR_TOF1[] = {0.551, 0.556, 0.632, 0.600, 0.649, 0.527, 0.501, 0.567, 0.580};
  double mu_TOF1[] = {100.280, 200.393, 350.002, 498.981, 748.393, 995.516, 1242.172, 1492.932, 1982.999};
  double sigma_TOF1[] = {2.079, 3.517, 8.761, 14.787, 30.429, 41.554, 59.130, 90.762, 148.094};

  double aL_TOF2[] = {0.688, 0.806, 0.826, 0.839, 0.854, 0.886, 0.964, 0.941, 1.062};
  double aR_TOF2[] = {0.789, 0.916, 0.855, 0.905, 0.931, 1.010, 0.950, 0.925, 1.047};
  double mu_TOF2[] = {100.396, 200.284, 349.787, 499.844, 748.470, 998.663, 1244.760, 1493.372, 1986.291};
  double sigma_TOF2[] = {1.339, 3.325, 7.686, 14.235, 27.394, 45.285, 60.213, 80.647, 130.084};

  TGraph *g_mu, *g_sigma, *g_aL, *g_aR;
  if (!TOFtrigger) {
    g_mu = new TGraph(M_st.size(), &M_st[0], mu_HT1);
    g_sigma = new TGraph(M_st.size(), &M_st[0], sigma_HT1);
    g_aL = new TGraph(M_st.size(), &M_st[0], aL_HT1);
    g_aR = new TGraph(M_st.size(), &M_st[0], aR_HT1);
  }
  else {
    g_mu = new TGraph(M_st.size(), &M_st[0], mu_TOF1);
    g_sigma = new TGraph(M_st.size(), &M_st[0], sigma_TOF1);
    g_aL = new TGraph(M_st.size(), &M_st[0], aL_TOF1);
    g_aR = new TGraph(M_st.size(), &M_st[0], aR_TOF1);
  }

  w->var("mu_1")->setVal(g_mu->Eval(mass));
  w->var("sigma_1")->setVal(g_sigma->Eval(mass));
  w->var("aL_1")->setVal(g_aL->Eval(mass));
  w->var("aR_1")->setVal(g_aR->Eval(mass));

  w->var("mu_1")->setConstant(true);
  w->var("sigma_1")->setConstant(true);
  w->var("aL_1")->setConstant(true);
  w->var("aR_1")->setConstant(true);

  TGraph *g2_mu, *g2_sigma, *g2_aL, *g2_aR;
  if (!TOFtrigger) {
    g2_mu = new TGraph(M_st.size(), &M_st[0], mu_HT2);
    g2_sigma = new TGraph(M_st.size(), &M_st[0], sigma_HT2);
    g2_aL = new TGraph(M_st.size(), &M_st[0], aL_HT2);
    g2_aR = new TGraph(M_st.size(), &M_st[0], aR_HT2);
  }
  else {
    g2_mu = new TGraph(M_st.size(), &M_st[0], mu_TOF2);
    g2_sigma = new TGraph(M_st.size(), &M_st[0], sigma_TOF2);
    g2_aL = new TGraph(M_st.size(), &M_st[0], aL_TOF2);
    g2_aR = new TGraph(M_st.size(), &M_st[0], aR_TOF2);
  }

  w->var("mu_2")->setVal(g2_mu->Eval(mass));
  w->var("sigma_2")->setVal(g2_sigma->Eval(mass));
  w->var("aL_2")->setVal(g2_aL->Eval(mass));
  w->var("aR_2")->setVal(g2_aR->Eval(mass));

  w->var("mu_2")->setConstant(true);
  w->var("sigma_2")->setConstant(true);
  w->var("aL_2")->setConstant(true);
  w->var("aR_2")->setConstant(true);

  delete g_mu, g_sigma, g_aL, g_aR;
  delete g2_mu, g2_sigma, g2_aL, g2_aR;
}

void Set_effL(RooWorkspace* w, double lumi, double mass, bool TOFtrigger = false) {
    vector<double> M_st = {100.000, 200.000, 350.000, 500.000, 750.000, 1000.000, 1250.000, 1500.000, 2000.000};

    double eff_HT1[] = {0.011, 0.061, 0.168, 0.258, 0.359, 0.412, 0.440, 0.463, 0.483};
    double eff_HT2[] = {0.019, 0.081, 0.145, 0.170, 0.188, 0.193, 0.194, 0.186, 0.181};
    double eff_TOF1[] = {0.091, 0.242, 0.365, 0.413, 0.452, 0.468, 0.479, 0.490, 0.500};
    double eff_TOF2[] = {0.090, 0.148, 0.176, 0.185, 0.193, 0.195, 0.196, 0.187, 0.181};


    TGraph* gr_1 = nullptr;
    TGraph* gr_2 = nullptr;
    if (!TOFtrigger) {
      gr_1 = new TGraph(M_st.size(), &M_st[0], eff_HT1);
      gr_2 = new TGraph(M_st.size(), &M_st[0], eff_HT2);
    }
    else {
      gr_1 = new TGraph(M_st.size(), &M_st[0], eff_TOF1);
      gr_2 = new TGraph(M_st.size(), &M_st[0], eff_TOF2);
    }


    w->var("effL_1")->setVal(gr_1->Eval(mass) * lumi);
    w->var("effL_2")->setVal(gr_2->Eval(mass) * lumi);

    delete gr_1, gr_2;
}

TCanvas* DrawExampleMassSpectrum(RooWorkspace* w, double lumi, double xsec_injected = 0, TString label = "") {
  // w->Print();

  RooRealVar * m = w->var("Mass");  // the observable
  RooCategory * index = w->cat("index");
  RooAbsData * data = w->data("data");
  RooAbsPdf * pdf = w->pdf("jointSBModel");

  // Fit the bkg to the exponential
  RooFitResult * r = pdf->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  r->Print();

  RooPlot * plot1 = m->frame(Title("Single particle category (1)"));
  data->plotOn(plot1,RooFit::Cut("index==index::ch1"), Name("Pseudo-data"));
  pdf->plotOn(plot1, RooFit::Components("model_1"), RooFit::ProjWData(*data), RooFit::Slice(*index,"ch1"), Name("S+B model"));
  pdf->plotOn(plot1, RooFit::Components("bkg_pdf_1"), RooFit::ProjWData(*data), RooFit::Slice(*index,"ch1"), RooFit::LineStyle(kDashed), Name("Bkg component"));
  pdf->plotOn(plot1, RooFit::Components("sig_pdf_1"), RooFit::ProjWData(*data), RooFit::Slice(*index,"ch1"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), Name("Sig component") );
  // plot1->GetYaxis()->SetRangeUser(0.5, 6e5);
  plot1->GetXaxis()->SetLimits(50, 2.5*test_mass);
  plot1->GetYaxis()->SetTitleSize(.05);
  plot1->GetYaxis()->SetLabelSize(.05);
  plot1->GetXaxis()->SetTitleSize(.05);
  plot1->GetXaxis()->SetLabelSize(.05);


  RooPlot * plot2 = m->frame(Title("Two particles category (2)"));
  data->plotOn(plot2,RooFit::Cut("index==index::ch2"));
  pdf->plotOn(plot2, RooFit::Components("model_2"), RooFit::ProjWData(*data), RooFit::Slice(*index,"ch2"), Name("S+B model"));
  pdf->plotOn(plot2, RooFit::Components("bkg_pdf_2"), RooFit::ProjWData(*data), RooFit::Slice(*index,"ch2"), RooFit::LineStyle(kDashed), Name("Bkg component"));
  pdf->plotOn(plot2, RooFit::Components("sig_pdf_2"), RooFit::ProjWData(*data), RooFit::Slice(*index,"ch2"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed), Name("Sig component") );
  // plot2->GetYaxis()->SetRangeUser(0.5, 6e3);
  plot2->GetXaxis()->SetLimits(50, 2*test_mass);
  plot2->GetYaxis()->SetTitleSize(.05);
  plot2->GetYaxis()->SetLabelSize(.05);
  plot2->GetXaxis()->SetTitleSize(.05);
  plot2->GetXaxis()->SetLabelSize(.05);
  pdf->paramOn(plot2,Layout(0.4,0.97,0.93), Format("NEU",AutoPrecision(1)));

  // auto h_pull_bkg = frame->residHist(0,0,true);
  // auto frame2 = m->frame(Title(""));
  // frame2->addPlotable(h_pull_bkg, "P");

  TCanvas* c_out = new TCanvas("c_out_"+label, "c_out_"+label, 1200, 600);
  c_out->Divide(2,1);

  c_out->cd(1);
  plot1->Draw();
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->SetLeftMargin(.15);
	gPad->SetRightMargin(.03);
  gPad->SetBottomMargin(.125);
	gPad->SetTopMargin(.07);
  auto leg = new TLegend(0.6,0.6,0.97,0.93);
  leg->AddEntry(plot1->findObject("Pseudo-data"), "Pseudo-data", "lep");
  leg->AddEntry(plot1->findObject("S+B model"), "S+B model", "l");
  leg->AddEntry(plot1->findObject("Bkg component"), "Bkg component", "l");
  leg->AddEntry(plot1->findObject("Sig component"), "Sig component", "l");
  leg->Draw();
  auto note = new TLatex();
  note->SetTextSize(0.04);
  note->SetLineColor(1);
  note->SetLineWidth(1);
  note->DrawLatexNDC(0.18, 0.85, Form("#splitline{#sqrt{s} 14TeV, L = %.0ffb^{-1}}{M_{#tilde{t}_{1}} = %.0f GeV,  #sigma_{MC} = %.0f fb}", lumi, test_mass, xsec_injected));


  c_out->cd(2);
  plot2->Draw();
  gPad->SetGrid();
  gPad->SetLogy();
  gPad->SetLeftMargin(.15);
	gPad->SetRightMargin(.03);
  gPad->SetBottomMargin(.125);
	gPad->SetTopMargin(.07);

  c_out->Update();
  return c_out;
}

RooWorkspace* ModelTOFAnalysis_RHad(double lumi, double xsec_data = 0, bool TOFtrigger = false) //lumi [pb]
{
  double binning[] = {1150, 50, 3500};

  RooWorkspace* w = new RooWorkspace("w");
  w->addClassDeclImportDir("/Users/olmo/cernbox/PID_timing_studies/script_analysis");
  w->importClassCode("RooGausDoubleExp",kTRUE);
  w->importClassCode("RooPowerLaw",kTRUE);

  cout << "Creating the pdfs\n";
  // Create the model for category 1 particle
  w->factory(Form("Gamma:bkg_pdf_1(Mass[%f,%f], gamma[1.e-6, 1.e-6, 1.e-6], lambda_1[%f, 10, 90], mu_Gamma[0., 0., 0.])",binning[1], binning[2], TOFtrigger?Mstar_bkg_TOF[0]:Mstar_bkg_HT[0]));
  w->var("gamma")->setConstant(true);
  w->var("mu_Gamma")->setConstant(true);
  w->factory(Form("GausDoubleExp:sig_pdf_1(Mass, mu_1[200., 10, 4000], sigma_1[30., 1., 3000], aL_1[2, 0.1, 3], aR_1[1.5, 0.1, 3])"));
  w->factory("prod:nsig_1(xsec[1., 1e-5, 2e3], effL_1[12e2])");
  w->factory("SUM:model_1(nsig_1*sig_pdf_1, nbkg_1[1e6, 1e5, 5e7]*bkg_pdf_1)");
  w->factory("PROD:model_bkg_1(nbkg_1, bkg_pdf_1)");

  // Create model for 2 particle category
  double l_aux = TOFtrigger?Mstar_bkg_TOF[1]:Mstar_bkg_HT[1];
  l_aux = -1./l_aux;
  w->factory(Form("Exponential:bkg_pdf_2(Mass, lambda_2[%f, -5, -0.001])", l_aux));
  // w->factory(Form("Gamma:bkg_pdf_2(Mass, gamma_2[1.e-6, 1.e-6, 1.e-6], lambda_2[%f, 10, 50], mu_Gamma2[0., 0., 0.])", TOFtrigger?Mstar_bkg_TOF[1]:Mstar_bkg_HT[1]));
  // w->var("gamma_2")->setConstant(true);
  // w->var("mu_Gamma2")->setConstant(true);

  w->factory(Form("GausDoubleExp:sig_pdf_2(Mass, mu_2[200., 10, 4000], sigma_2[30., 1., 1000], aL_2[2, 0.1, 3], aR_2[1.5, 0.1, 3])"));
  w->factory("prod:nsig_2(xsec, effL_2[12e2])");
  w->factory("SUM:model_2(nsig_2*sig_pdf_2, nbkg_2[1e5, 1e3, 1e7]*bkg_pdf_2)");
  w->factory("PROD:model_bkg_2(nbkg_2, bkg_pdf_2)");


  // Create joint model
  w->var("Mass")->SetTitle("Mass");
  w->var("Mass")->setUnit("GeV");
  w->var("xsec")->SetTitle("#sigma");
  w->var("xsec")->setUnit("fb");
  w->var("nbkg_1")->SetTitle("N_{bkg}^{(1)}");
  w->var("nbkg_2")->SetTitle("N_{bkg}^{(2)}");
  w->factory("index[ch1, ch2]");
  w->factory("SIMUL:jointSBModel(index, ch1=model_1, ch2=model_2)");   //RooSimultaneous fit
  w->factory("SIMUL:jointBkgModel(index, ch1=model_bkg_1, ch2=model_bkg_2)");   //RooSimultaneous fit

  //=========== Generate the data =====================================
  cout << "Generating the data\n";
  Set_GausDoubleExpPdf_pars(w, test_mass, TOFtrigger);
  Set_effL(w, lumi, test_mass, TOFtrigger);

  // Set the correct amount of events for exclusion limits
  w->var("xsec")->setVal(xsec_data);
  if (!TOFtrigger) {
    w->var("nbkg_1")->setVal( (lumi/Lumi_Nbkg) * Nevs_bkg_HT[0] );
    w->var("nbkg_2")->setVal( (lumi/Lumi_Nbkg) * Nevs_bkg_HT[1] );
  }
  else {
    w->var("nbkg_1")->setVal( (lumi/Lumi_Nbkg) * Nevs_bkg_TOF [0] );
    w->var("nbkg_2")->setVal( (lumi/Lumi_Nbkg) * Nevs_bkg_TOF [1] );
  }

  RooRealVar * m = w->var("Mass");  // the observable
  m->setBins(binning[0]);
  RooCategory * index = w->cat("index");

  RooRandom::randomGenerator()->SetSeed(111);
  RooDataHist * data = w->pdf("jointSBModel")->generateBinned(RooArgSet(*m, *index));
  data->SetName("data");
  w->import(*data);

  //=========== Create mdoel config =====================================
  cout << "Creating the model config\n";
  RooStats::ModelConfig mc("ModelSBConfig",w);
  mc.SetPdf(*w->pdf("jointSBModel"));
  mc.SetParametersOfInterest(*w->var("xsec"));
  mc.SetObservables(*w->var("Mass"));
  // define set of nuisance parameters
  w->defineSet("nuisParams","lambda_1,nbkg_1,lambda_2,nbkg_2");
  mc.SetNuisanceParameters(*w->set("nuisParams"));
  // import model in the workspace
  w->import(mc);

  RooStats::ModelConfig mc_bkg("ModelBConfig",w);
  mc_bkg.SetPdf(*w->pdf("jointBkgModel"));
  mc_bkg.SetParametersOfInterest(*w->var("xsec"));
  mc_bkg.SetObservables(*w->var("Mass"));
  // define set of nuisance parameters
  // w->defineSet("nuisParams","lambda_1,nbkg_1,lambda_2,nbkg_2");
  mc_bkg.SetNuisanceParameters(*w->set("nuisParams"));
  // import model in the workspace
  w->import(mc_bkg);

  // write the workspace in the file
  TString fileName = "/Users/olmo/cernbox/PID_timing_studies/_root/histos/ModelRHad.root";
  w->writeToFile(fileName,true);
  cout << "model written to file " << fileName << endl;

  return w;
}

TCanvas* MakeBrazilPlot(vector<double> RHad_mass_list, map<int, vector<double>> exp_upper_limit_now, map<int, vector<double>> exp_upper_limit_HL, map<int, vector<double>> exp_upper_limit_HL_TOFtrigger) {

  auto gr = new TGraph(RHad_mass_list.size());
  auto gr_1s = new TGraphAsymmErrors(RHad_mass_list.size());
  auto gr_2s = new TGraphAsymmErrors(RHad_mass_list.size());

  for(unsigned int i=0; i<RHad_mass_list.size(); i++) {
    double mass = RHad_mass_list[i];

    gr->SetPoint(i, mass, exp_upper_limit_HL[0][i]);

    double scale = 1.;
    gr_1s->SetPoint(i, mass, exp_upper_limit_HL[0][i]);
    double aux_l = scale*(exp_upper_limit_HL[0][i] - exp_upper_limit_HL[-1][i]);
    double aux_u = scale*(exp_upper_limit_HL[1][i] - exp_upper_limit_HL[0][i]);
    gr_1s->SetPointError(i, 0, 0, aux_l, aux_u);

    gr_2s->SetPoint(i, mass, exp_upper_limit_HL[0][i]);
    aux_l = scale*(exp_upper_limit_HL[0][i] - exp_upper_limit_HL[-2][i]);
    aux_u = scale*(exp_upper_limit_HL[2][i] - exp_upper_limit_HL[0][i]);
    gr_2s->SetPointError(i, 0, 0, aux_l, aux_u);

    cout << " =============== Mass " << mass << " ========================" << endl;
    cout << Form("xsec: %.1e %.1e %.1e %.1e %.1e", exp_upper_limit_HL[-2][i], exp_upper_limit_HL[-1][i], exp_upper_limit_HL[0][i], exp_upper_limit_HL[1][i], exp_upper_limit_HL[2][i]) << endl;
  }

  auto c_mass = new TCanvas("c_mass", "c_mass", 800, 600);
  gr_2s->SetTitle("");
  gr_2s->SetFillColor(5);
  // gr_2s->SetFillColorAlpha(kYellow);//, 0.8);
  gr_2s->SetFillStyle(1001);
  gr_2s->Draw("A3");

  gr_1s->SetFillColor(8);
  // gr_1s->SetFillColorAlpha(kGreen);//, 0.7);
  gr_1s->SetFillStyle(1001);
  gr_1s->Draw("3");

  gr->SetMarkerStyle(20);
  gr->SetMarkerColor(1);
  gr->SetLineColor(1);
  gr->SetLineStyle(1);
  gr->SetLineWidth(2);
  gr->Draw("LP");

  // Limits with 12 fb-1
  auto gr_now = new TGraph(RHad_mass_list.size(), &RHad_mass_list[0], &(exp_upper_limit_now[0][0]));
  gr_now->SetMarkerStyle(20);
  gr_now->SetMarkerColor(2);
  gr_now->SetLineColor(2);
  gr_now->SetLineStyle(1);
  gr_now->SetLineWidth(2);
  gr_now->Draw("LP");

  //TOF trigger
  auto gr_toftrigger = new TGraph(exp_upper_limit_HL_TOFtrigger[0].size(), &RHad_mass_list[0], &(exp_upper_limit_HL_TOFtrigger[0][0]));
  gr_toftrigger->SetMarkerStyle(4);
  gr_toftrigger->SetMarkerColor(1);
  gr_toftrigger->SetLineColor(1);
  gr_toftrigger->SetLineStyle(7);
  gr_toftrigger->SetLineWidth(2);
  gr_toftrigger->Draw("LP");




  // Add current CMS limits
  double masses[] = { 99.1816, 398.5258, 599.9218, 795.8401, 1195.8373, 1590.3867, 2001.2604, 2393.0807 };
  double xsec_up[] = { 10.7, 4.68, 1.4103, 1.4103, 1.5569, 2.0691, 2.9617, 5.1667 }; //[fb]
  auto gr_CMS = new TGraph(8, masses, xsec_up);
  gr_CMS->SetLineColor(4);
  gr_CMS->SetLineWidth(2);
  gr_CMS->SetMarkerColor(4);
  gr_CMS->SetMarkerStyle(20);
  gr_CMS->Draw("CP");

  // for(uint i = 0; i < 8; i++) {
  //   xsec_up[i] *= sqrt(12./1000.);
  // }
  // auto gr_CMS_exp = new TGraph(8, masses, xsec_up);
  // gr_CMS_exp->SetLineColor(2);
  // gr_CMS_exp->SetLineWidth(2);
  // gr_CMS_exp->SetMarkerColor(2);
  // gr_CMS_exp->SetMarkerStyle(20);
  // gr_CMS_exp->Draw("CP");

  auto leg = new TLegend(0.63, 0.65, 0.95, 0.95);
  leg->AddEntry(gr_CMS, "CMS EXO-16-36 (12 fb^{-1})", "lp");
  // leg->AddEntry(gr_CMS_exp, "CMS expected @ 1 ab^{-1}", "lp");
  leg->AddEntry(gr_now, "Expected (12 fb^{-1})", "lp");
  leg->AddEntry(gr, "Expected (1 ab^{-1})", "lp");
  leg->AddEntry(gr_toftrigger, "Expected (1 ab^{-1}) - TOF trigger", "lp");
  leg->AddEntry(gr_1s, "Expected #pm 1 #sigma", "F");
  leg->AddEntry(gr_2s, "Expected #pm 2 #sigma", "F");
  leg->Draw();

  auto note = new TLatex();
  note->DrawLatexNDC(0.15, 0.88, "#splitline{Delphes CMS@HL-LHC Sim.}{#sqrt{s} = 14TeV, <PU> = 140}");

  c_mass->SetGrid();
  c_mass->SetLogy();
  c_mass->SetLogx();
  gr_2s->GetXaxis()->SetTitle("#tilde{t}_{1} mass [GeV]");
  gr_2s->GetYaxis()->SetTitle("Excluded xsec @95% CLs [fb]");
  gr_2s->GetYaxis()->SetRangeUser(5e-4, 1e4);
  gr_2s->GetYaxis()->SetTitleSize(.05);
  gr_2s->GetYaxis()->SetTitleOffset(1.25);
  gr_2s->GetYaxis()->SetLabelSize(.05);
  gr_2s->GetXaxis()->SetTitleSize(.05);
  gr_2s->GetXaxis()->SetTitleOffset(1.2);
  gr_2s->GetXaxis()->SetLabelSize(.05);
  gPad->SetLeftMargin(.13);
	gPad->SetRightMargin(.05);
  gPad->SetBottomMargin(.14);
	gPad->SetTopMargin(.05);
  c_mass->Update();
  return c_mass;
}

map<int, vector<double>> Compute_limit_band(vector<double> masses_scan, double Luminosity = 12, bool TOFtrigger = false){
  cout << "Building the model" << endl;
  auto w_sig = ModelTOFAnalysis_RHad(Luminosity, test_xsec);
  DrawExampleMassSpectrum(w_sig, Luminosity, test_xsec, Form("M%.0f_L%.0f", test_mass, Luminosity));

  auto w_bkg = ModelTOFAnalysis_RHad(Luminosity, 0, TOFtrigger);
  DrawExampleMassSpectrum(w_bkg, Luminosity, 0, Form("bkg_L%.0f", Luminosity));
  // map<int, vector<double>> a; return a;

  cout << "Computing the exclusion limit\n";
  vector<int> sigma_eval = {-2, -1, 0, 1, 2};
  map<int, vector<double>> exp_upper_limit;
  for( auto s : sigma_eval ) {
    exp_upper_limit[s] = {};
  }


  float points[4][9] = {
    {100, 150,  200,  300,  500,   600,  1000,   6000},
    {500,  90,   35,   10,    6,     4,     2,      2},
    { 60,   8,    2,  0.3, 0.04,  0.05, 0.022, 0.01},
    { 10,   5,  1.6,  0.3, 0.04,  0.05, 0.025, 0.01},
  };

  TGraph *gr_max_xsec;
  if (!TOFtrigger && Luminosity == (double)12) gr_max_xsec = new TGraph(9, points[0], points[1]);
  else if (!TOFtrigger && Luminosity == (double)1e3) gr_max_xsec = new TGraph(9, points[0], points[2]);
  else if (TOFtrigger && Luminosity == (double)1e3) gr_max_xsec = new TGraph(9, points[0], points[3]);
  else {
    cout << "Lumi " << Luminosity << " not expected\n";
    exit(0);
  }

  for (auto mass : masses_scan) {
    Set_GausDoubleExpPdf_pars(w_bkg, mass, TOFtrigger);
    Set_effL(w_bkg, Luminosity, mass, TOFtrigger);

    double max_xsec = gr_max_xsec->Eval(mass);

    HypoTestInvTool calc(&optHTInv);

    HypoTestInverterResult * r = calc.RunInverter(w_bkg, "ModelSBConfig", "ModelBConfig",
                                                  "data", 2, 3, true,
                                                  50, 0, max_xsec,
                                                  100);

    calc.AnalyzeResult( r, 2, 3, true, 50, "/Users/olmo/cernbox/PID_timing_studies/_root/results/", Form("_L%.0ffb_%s_M%.0f", Luminosity, TOFtrigger?"TOF":"HT", mass));

    for( auto s : sigma_eval ) {
      exp_upper_limit[s].push_back(r->GetExpectedUpperLimit(s));
    }
  }

  return exp_upper_limit;
}

void TOFAnalysis_RHad_catCombined(){
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  vector<double> masses_scan = {200};
  // vector<double> masses_scan = {100, 150, 200, 300, 400, 600, 800, 1000, 1200, 1600, 2000, 2500};
  // vector<double> masses_scan_low = {200};
  // vector<double> masses_scan_low = {100, 150, 200, 300, 400, 600, 800, 1000};

  auto exp_upper_limit_now = Compute_limit_band(masses_scan, 12, false);
  // auto exp_upper_limit_HL = Compute_limit_band(masses_scan, 1e3, false);

  // auto exp_upper_limit_HL_TOFtrigger = Compute_limit_band(masses_scan_low, 1e3, true);


  // auto c_mass = MakeBrazilPlot(masses_scan, exp_upper_limit_now, exp_upper_limit_HL, exp_upper_limit_HL_TOFtrigger);
  // c_mass->SaveAs("./_fig/ExclusionLimits.root");
  // c_mass->SaveAs("./_fig/ExclusionLimits.png");

}
