#include "TCanvas.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooDataSet.h"
#include "RooFitResult.h"
#include "RooPlot.h"
#include "RooRealVar.h"
#include "RooRandom.h"

#include "RooStats/ModelConfig.h"

using namespace RooFit;



void GausExpModel(int nsig = 200,    // number of signal events
                  int nbkg = 3100 )  // number of background events
{
  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING) ;

  double binning[] = {50, 0, 400};

  RooWorkspace w("w");
  w.factory(Form("Exponential:bkg_pdf(Mass[%f,%f], a[-0.023,-5,-0.001])",binning[1], binning[2]));
  w.factory("Gaussian:sig_pdf(Mass, mu[200.], sigma[30.])");

  w.factory("SUM:model(nsig[0,10000]*sig_pdf, nbkg[0,10000]*bkg_pdf)");  // for extended model

  RooAbsPdf * pdf = w.pdf("model");
  RooRealVar * m = w.var("Mass");  // the observable

  // set the desired value of signal and background events
  w.var("nsig")->setVal(nsig);
  w.var("nbkg")->setVal(nbkg);

  // generate the data

  // use fixed random numbers for reproducibility (use 0 for changing every time)
  RooRandom::randomGenerator()->SetSeed(111);

  // fix number of bins to 50 to plot or to generate data (default is 100 bins)
  m->setBins(binning[0]);

  RooDataSet * data = pdf->generate( *m);  // will generate accordint to total S+B events
  //RooDataSet * data = pdf->generate( *x, AllBinned());  // will generate accordint to total S+B events
  data->SetName("data");
  w.import(*data);

  data->Print();

  TCanvas* c_out = new TCanvas("c_out", "c_out", 800, 600);

  RooPlot * plot = m->frame(Title("Gaussian Signal over Exponential Background"));
  data->plotOn(plot);
  plot->Draw();

  RooFitResult * r = pdf->fitTo(*data, RooFit::Save(true), RooFit::Minimizer("Minuit2","Migrad"));
  r->Print();

  pdf->plotOn(plot);
  //draw the two separate pdf's
  pdf->plotOn(plot, RooFit::Components("bkg_pdf"), RooFit::LineStyle(kDashed) );
  pdf->plotOn(plot, RooFit::Components("sig_pdf"), RooFit::LineColor(kRed), RooFit::LineStyle(kDashed) );

  pdf->paramOn(plot,Layout(0.5,0.9,0.85));

  plot->Draw();

  c_out->SetLogy();
  c_out->SetGrid();


  RooStats::ModelConfig mc("ModelConfig",&w);
  mc.SetPdf(*pdf);
  mc.SetParametersOfInterest(*w.var("nsig"));
  mc.SetObservables(*w.var("Mass"));
  // define set of nuisance parameters
  w.defineSet("nuisParams","a,nbkg");

  mc.SetNuisanceParameters(*w.set("nuisParams"));

  // import model in the workspace
  w.import(mc);

  // write the workspace in the file
  TString fileName = "_root/histos/GausExpModel.root";
  w.writeToFile(fileName,true);
  cout << "model written to file " << fileName << endl;
}
