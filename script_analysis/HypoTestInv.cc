/// \file
/// \ingroup tutorial_roostats
/// \notebook
/// Standard tutorial macro for performing an inverted  hypothesis test for computing an interval
///
/// This macro will perform a scan of the p-values for computing the interval or limit
///
/// Usage:
///
/// ~~~{.cpp}
/// root>.L StandardHypoTestInv.C
/// root> StandardHypoTestInv("fileName","workspace name","S+B modelconfig name","B model name","data set name",calculator type, test statistic type, use CLS,
///                                number of points, xmin, xmax, number of toys, use number counting)
///
/// type = 0 Freq calculator
/// type = 1 Hybrid calculator
/// type = 2 Asymptotic calculator
/// type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)
///
/// testStatType = 0 LEP
///              = 1 Tevatron
///              = 2 Profile Likelihood two sided
///              = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
///              = 4 Profile Likelihood signed ( pll = -pll if mu < mu_hat)
///              = 5 Max Likelihood Estimate as test statistic
///              = 6 Number of observed event as test statistic
/// ~~~
///
/// \macro_image
/// \macro_output
/// \macro_code
///
/// \author Lorenzo Moneta

#include "TFile.h"
#include "RooWorkspace.h"
#include "RooAbsPdf.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooStats/ModelConfig.h"
#include "RooRandom.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TROOT.h"
#include "TSystem.h"

#include "RooStats/AsymptoticCalculator.h"
#include "RooStats/HybridCalculator.h"
#include "RooStats/FrequentistCalculator.h"
#include "RooStats/ToyMCSampler.h"
#include "RooStats/HypoTestPlot.h"

#include "RooStats/NumEventsTestStat.h"
#include "RooStats/ProfileLikelihoodTestStat.h"
#include "RooStats/SimpleLikelihoodRatioTestStat.h"
#include "RooStats/RatioOfProfiledLikelihoodsTestStat.h"
#include "RooStats/MaxLikelihoodEstimateTestStat.h"
#include "RooStats/NumEventsTestStat.h"

#include "RooStats/HypoTestInverter.h"
#include "RooStats/HypoTestInverterResult.h"
#include "RooStats/HypoTestInverterPlot.h"

using namespace RooFit;
using namespace RooStats;
using namespace std;

// structure defining the options
struct HypoTestInvOptions {

  bool plotHypoTestResult = true;          // plot test statistic result at each point
  bool writeResult = true;                 // write HypoTestInverterResult in a file
  TString resultFileName;                  // file with results (by default is built automatically using the workspace input file name)
  bool optimize = true;                    // optimize evaluation of test statistic
  bool useVectorStore = true;              // convert data to use new roofit data store
  bool generateBinned = true;             // generate binned data sets
  bool noSystematics = false;              // force all systematics to be off (i.e. set all nuisance parameters as constat
                                           // to their nominal values)
  double nToysRatio = 2;                   // ratio Ntoys S+b/ntoysB
  double maxPOI = -1;                      // max value used of POI (in case of auto scan)
  bool useProof = false;                   // use Proof Lite when using toys (for freq or hybrid)
  int nworkers = 0;                        // number of worker for ProofLite (default use all available cores)
  bool enableDetailedOutput = true;       // enable detailed output with all fit information for each toys (output will be written in result file)
  int initialFit = 1;                     // do a first  fit to the model (-1 : default, 0 skip fit, 1 do always fit)
  int randomSeed = 0;                     // random seed (if = -1: use default value, if = 0 always random )
                                           // NOTE: Proof uses automatically a random seed

  int nAsimovBins = 0;                     // number of bins in observables used for Asimov data sets (0 is the default and it is given by workspace, typically is 100)

  bool reuseAltToys = false;                // reuse same toys for alternate hypothesis (if set one gets more stable bands)
  double confLevel = 0.95;                // confidence level value



  std::string  minimizerType = "";                  // minimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()
  int   printLevel = 0;                    // print level for debugging PL test statistics and calculators

  bool useNLLOffset = false;               // use NLL offset when fitting (this increase stability of fits)

};
HypoTestInvOptions optHTInv;

class HypoTestInvTool{

   public:
      HypoTestInvTool(HypoTestInvOptions* opt = nullptr) {
        if (opt != nullptr) {
          PlotHypoTestResult = opt->plotHypoTestResult;
          WriteResult = opt->writeResult;
          Optimize = opt->optimize;
          UseVectorStore = opt->useVectorStore;
          GenerateBinned = opt->generateBinned;
          UseProof = opt->useProof;
          ReuseAltToys = opt->reuseAltToys;
          EnableDetOutput = opt->enableDetailedOutput;
          NWorkers = opt->nworkers;
          Verbose = opt->printLevel;
          InitialFit = opt->initialFit;
          RandomSeed = opt->randomSeed;
          NToysRatio = opt->nToysRatio;
          MaxPoi = opt->maxPOI;
          AsimovBins = opt->nAsimovBins;
          MinimizerType = opt->minimizerType;
          ResultFileName = opt->resultFileName;

          if (optHTInv.useNLLOffset) RooStats::UseNLLOffset(true);

          if (MinimizerType.size()==0) MinimizerType = ROOT::Math::MinimizerOptions::DefaultMinimizerType();
          else ROOT::Math::MinimizerOptions::SetDefaultMinimizer(MinimizerType.c_str());
        }
      };
      ~HypoTestInvTool(){};

      HypoTestInverterResult * RunInverter(RooWorkspace * w,
                  const char * modelSBName, const char * modelBName,
                  const char * dataName,
                  int type,  int testStatType,
                  bool useCLs,
                  int npoints, double poimin, double poimax, int ntoys,
                  bool useNumberCounting = false,
                  const char * nuisPriorName = 0);



      void AnalyzeResult( HypoTestInverterResult * r,
                     int calculatorType,
                     int testStatType,
                     bool useCLs,
                     int npoints,
                     const char * fileNameBase = 0,
                     std::string MassValue = "");

      bool        PlotHypoTestResult = false;
      bool        WriteResult = false;
      bool        Optimize = true;
      bool        UseVectorStore = true;
      bool        GenerateBinned = false;
      bool        UseProof = false;
      bool        ReuseAltToys = false;
      bool        EnableDetOutput = false;
      int         NWorkers = 4;
      int         Verbose = 0;
      int         InitialFit = -1;
      int         RandomSeed = -1;
      double      NToysRatio = 1;
      double      MaxPoi = -1;
      int         AsimovBins = 0;
      std::string MinimizerType = "";  // inimizer type (default is what is in ROOT::Math::MinimizerOptions::DefaultMinimizerType()
      TString     ResultFileName = "";
   };

void StandardHypoTestInv(const char * infile = 0,
                        const char * wsName = "w",
                        const char * modelSBName = "ModelConfig",
                        const char * modelBName = "",
                        const char * dataName = "data",
                        int calculatorType = 2,
                        int testStatType = 3,
                        bool useCLs = true ,
                        int npoints = 25,
                        double poimin = 0,
                        double poimax = 50,
                        int ntoys=1000,
                        bool useNumberCounting = false,
                        const char * nuisPriorName = 0)
{
  /*

    Other Parameter to pass in tutorial
    apart from standard for filename, ws, modelconfig and data

    type = 0 Freq calculator
    type = 1 Hybrid calculator
    type = 2 Asymptotic calculator
    type = 3 Asymptotic calculator using nominal Asimov data sets (not using fitted parameter values but nominal ones)

    testStatType = 0 LEP
    = 1 Tevatron
    = 2 Profile Likelihood
    = 3 Profile Likelihood one sided (i.e. = 0 if mu < mu_hat)
    = 4 Profiel Likelihood signed ( pll = -pll if mu < mu_hat)
    = 5 Max Likelihood Estimate as test statistic
    = 6 Number of observed event as test statistic

    useCLs          scan for CLs (otherwise for CLs+b)

    npoints:        number of points to scan , for autoscan set npoints = -1

    poimin,poimax:  min/max value to scan in case of fixed scans
    (if min >  max, try to find automatically)

    ntoys:         number of toys to use

    useNumberCounting:  set to true when using number counting events

    nuisPriorName:   name of prior for the nuisance. This is often expressed as constraint term in the global model
    It is needed only when using the HybridCalculator (type=1)
    If not given by default the prior pdf from ModelConfig is used.

    extra options are available as global parameters of the macro. They major ones are:

    plotHypoTestResult   plot result of tests at each point (TS distributions) (default is true)
    useProof             use Proof   (default is true)
    writeResult          write result of scan (default is true)
    rebuild              rebuild scan for expected limits (require extra toys) (default is false)
    generateBinned       generate binned data sets for toys (default is false) - be careful not to activate with
    a too large (>=3) number of observables
    nToyRatio            ratio of S+B/B toys (default is 2)


  */



   TString filename(infile);
   TFile *file = TFile::Open(filename);
   if(!file ){
      cout <<"StandardRooStatsDemoMacro: Input file " << filename << " not found" << endl;
      exit(0);
   }

   RooWorkspace * w = dynamic_cast<RooWorkspace*>( file->Get(wsName) );
   if ( w == NULL ) {
     std::cerr << "File " << filename << " does not contain a workspace named " << wsName << std::endl;
     exit(0);
   }

    HypoTestInvTool calc(&optHTInv);
    HypoTestInverterResult * r = calc.RunInverter(w, modelSBName, modelBName,
                         dataName, calculatorType, testStatType, useCLs,
                         npoints, poimin, poimax,
                         ntoys, useNumberCounting, nuisPriorName );


   calc.AnalyzeResult( r, calculatorType, testStatType, useCLs, npoints, infile );

   return;
}



void HypoTestInvTool::AnalyzeResult( HypoTestInverterResult * r,
                                          int calculatorType,
                                          int testStatType,
                                          bool useCLs,
                                          int npoints,
                                          const char * fileNameBase,
                                          std::string MassValue)
{

   // analyze result produced by the inverter, optionally save it in a file
   double lowerLimit = 0;
   double llError = 0;
   if (r->IsTwoSided()) {
      lowerLimit = r->LowerLimit();
      llError = r->LowerLimitEstimatedError();
   }

   double upperLimit = r->UpperLimit();
   double ulError = r->UpperLimitEstimatedError();

   if (lowerLimit < upperLimit*(1.- 1.E-4) && lowerLimit != 0)
      std::cout << "The computed lower limit is: " << lowerLimit << " +/- " << llError << std::endl;
   std::cout << "The computed upper limit is: " << upperLimit << " +/- " << ulError << std::endl;


   // compute expected limit
   std::cout << "Expected upper limits, using the B (alternate) model : " << std::endl;
   std::cout << " expected limit (median) " << r->GetExpectedUpperLimit(0) << std::endl;
   std::cout << " expected limit (-1 sig) " << r->GetExpectedUpperLimit(-1) << std::endl;
   std::cout << " expected limit (+1 sig) " << r->GetExpectedUpperLimit(1) << std::endl;
   std::cout << " expected limit (-2 sig) " << r->GetExpectedUpperLimit(-2) << std::endl;
   std::cout << " expected limit (+2 sig) " << r->GetExpectedUpperLimit(2) << std::endl;



   // plot the result ( p values vs scan points)
   std::string typeName = "";
   if (calculatorType == 0 )
      typeName = "Frequentist";
   if (calculatorType == 1 )
      typeName = "Hybrid";
   else if (calculatorType == 2 || calculatorType == 3) {
      typeName = "Asymptotic";
      PlotHypoTestResult = false;
   }

   const char * resultName = r->GetName();
   TString plotTitle = TString::Format("%s CL Scan for workspace %s",typeName.c_str(),resultName);
   HypoTestInverterPlot *plot = new HypoTestInverterPlot("HTI_Result_Plot",plotTitle,r);

   // plot in a new canvas with style
   TString c1Name = TString::Format("%s_Scan",typeName.c_str());
   TCanvas * c1 = new TCanvas(c1Name);
   c1->SetLogy(false);

   plot->Draw("CLb 2CL");  // plot all and Clb

   TLine * ln = new TLine();
   ln->SetLineStyle(7);
   ln->SetLineColor(46);
   double eval_sigma[] = {-2, -1, 0, 1, 2};
   for( auto s : eval_sigma) {
     double x = r->GetExpectedUpperLimit(s);
     ln->DrawLine(x, 1-optHTInv.confLevel, x, 0);
   }

   c1->SetLogy();

   const int nEntries = r->ArraySize();

   // plot test statistics distributions for the two hypothesis
   if (PlotHypoTestResult) {
      TCanvas * c2 = new TCanvas("c2");
      if (nEntries > 1) {
         int ny = TMath::CeilNint(TMath::Sqrt(nEntries));
         int nx = TMath::CeilNint(double(nEntries)/ny);
         c2->Divide( nx,ny);
      }
      for (int i=0; i<nEntries; i++) {
         if (nEntries > 1) c2->cd(i+1);
         SamplingDistPlot * pl = plot->MakeTestStatPlot(i);
         pl->SetLogYaxis(true);
         pl->Draw();
      }
   }
   gPad = c1;

   // write result in a file
   if (r != NULL && WriteResult) {

     // write to a file the results
     const char *  calcType = (calculatorType == 0) ? "Freq" : (calculatorType == 1) ? "Hybr" : "Asym";
     const char *  limitType = (useCLs) ? "CLs" : "Cls+b";
     const char * scanType = (npoints < 0) ? "auto" : "grid";
     if (ResultFileName.IsNull()) {
       ResultFileName = TString::Format("%s_%s_%s_ts%d",calcType,limitType,scanType,testStatType);
       //strip the / from the filename
       if (MassValue.size()>0) ResultFileName += Form("_M%sGeV", MassValue.c_str());

       ResultFileName += ".root";

       ResultFileName = fileNameBase + ResultFileName;
     }

     TString PlotFileName = ResultFileName.Data();
     PlotFileName.ReplaceAll(".root", "_pval_scan.png");
     c1->SaveAs(PlotFileName);

     // get (if existing) rebuilt UL distribution
     TString uldistFile = "RULDist.root";
     TObject * ulDist = 0;
     bool existULDist = !gSystem->AccessPathName(uldistFile);
     if (existULDist) {
       TFile * fileULDist = TFile::Open(uldistFile);
       if (fileULDist) ulDist= fileULDist->Get("RULDist");
     }


     TFile * fileOut = new TFile(ResultFileName,"RECREATE");
     r->Write();
     if (ulDist) ulDist->Write();
     Info("StandardHypoTestInv","HypoTestInverterResult has been written in the file %s",ResultFileName.Data());

     fileOut->Close();
   }
}



// internal routine to run the inverter
HypoTestInverterResult *HypoTestInvTool::RunInverter(RooWorkspace * w,
                                       const char * modelSBName, const char * modelBName,
                                       const char * dataName, int type,  int testStatType,
                                       bool useCLs, int npoints, double poimin, double poimax,
                                       int ntoys,
                                       bool useNumberCounting,
                                       const char * nuisPriorName )
{
   std::cout << "Running HypoTestInverter on the workspace " << w->GetName() << std::endl;
   w->Print();

   RooAbsData * data = w->data(dataName);
   if (!data) {
      std::cerr << "Data with name " << dataName << " not found." << std::endl;
      exit(0);
   }
   else { std::cout << "Using data set " << dataName << std::endl;}

   if (UseVectorStore) {
      RooAbsData::setDefaultStorageType(RooAbsData::Vector);
      data->convertToVectorStore() ;
   }


   // get models from WS
   // get the modelConfig out of the file
   ModelConfig* bModel = (ModelConfig*) w->obj(modelBName);
   ModelConfig* sbModel = (ModelConfig*) w->obj(modelSBName);

   //Sanity check
   {
     if (!sbModel) {
        Error("StandardHypoTestDemo","Not existing ModelConfig %s",modelSBName);
        return 0;
     }
     if (!sbModel->GetPdf()) {
        Error("StandardHypoTestDemo","Model %s has no pdf ",modelSBName);
        return 0;
     }
     if (!sbModel->GetParametersOfInterest()) {
        Error("StandardHypoTestDemo","Model %s has no poi ",modelSBName);
        return 0;
     }
     if (!sbModel->GetObservables()) {
        Error("StandardHypoTestInv","Model %s has no observables ",modelSBName);
        return 0;
     }
     if (!sbModel->GetSnapshot() ) {
        Info("StandardHypoTestInv","Model %s has no snapshot  - make one using model poi",modelSBName);
        sbModel->SetSnapshot( *sbModel->GetParametersOfInterest() );
     }
   }

   // In case of no systematics remove nuisance parameters from model
   if (optHTInv.noSystematics) {
      const RooArgSet * nuisPar = sbModel->GetNuisanceParameters();
      if (nuisPar && nuisPar->getSize() > 0) {
         std::cout << "StandardHypoTestInv" << "  -  Switch off all systematics by setting them constant to their initial values" << std::endl;
         RooStats::SetAllConstant(*nuisPar);
      }
      if (bModel) {
         const RooArgSet * bnuisPar = bModel->GetNuisanceParameters();
         if (bnuisPar)
            RooStats::SetAllConstant(*bnuisPar);
      }
   }

   if (!bModel || bModel == sbModel) {
      Info("StandardHypoTestInv","The background model %s does not exist",modelBName);
      Info("StandardHypoTestInv","Copy it from ModelConfig %s and set POI to zero",modelSBName);
      bModel = (ModelConfig*) sbModel->Clone();
      bModel->SetName(TString(modelSBName)+TString("_with_poi_0"));
      RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
      if (!var) return 0;
      double oldval = var->getVal();
      var->setVal(0);
      bModel->SetSnapshot( RooArgSet(*var)  );
      var->setVal(oldval);
   }
   else {
      if (!bModel->GetSnapshot() ) {
         Info("StandardHypoTestInv","Model %s has no snapshot  - make one using model poi and 0 values ",modelBName);
         RooRealVar * var = dynamic_cast<RooRealVar*>(bModel->GetParametersOfInterest()->first());
         if (var) {
            double oldval = var->getVal();
            var->setVal(0);
            bModel->SetSnapshot( RooArgSet(*var)  );
            var->setVal(oldval);
         }
         else {
            Error("StandardHypoTestInv","Model %s has no valid poi",modelBName);
            return 0;
         }
      }
   }

   // check model  has global observables when there are nuisance pdf
   // for the hybrid case the globals are not needed
   if (type != 1 ) {
      bool hasNuisParam = (sbModel->GetNuisanceParameters() && sbModel->GetNuisanceParameters()->getSize() > 0);
      bool hasGlobalObs = (sbModel->GetGlobalObservables() && sbModel->GetGlobalObservables()->getSize() > 0);
      if (hasNuisParam && !hasGlobalObs ) {
         // try to see if model has nuisance parameters first
         RooAbsPdf * constrPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisanceConstraintPdf_sbmodel");
         if (constrPdf) {
            Warning("StandardHypoTestInv","Model %s has nuisance parameters but no global observables associated",sbModel->GetName());
            Warning("StandardHypoTestInv","\tThe effect of the nuisance parameters will not be treated correctly ");
         }
      }
   }

   // save all initial parameters of the model including the global observables
   RooArgSet initialParameters;
   RooArgSet * allParams = sbModel->GetPdf()->getParameters(*data);
   allParams->snapshot(initialParameters);
   delete allParams;

   // run first a data fit

   const RooArgSet * poiSet = sbModel->GetParametersOfInterest();
   RooRealVar *poi = (RooRealVar*)poiSet->first();

   std::cout << "StandardHypoTestInv : POI initial value:   " << poi->GetName() << " = " << poi->getVal()   << std::endl;

   // fit the data first (need to use constraint )
   TStopwatch tw;

   bool doFit = InitialFit;
   if (testStatType == 0 && InitialFit == -1) doFit = false;  // case of LEP test statistic
   if (type == 3  && InitialFit == -1) doFit = false;         // case of Asymptoticcalculator with nominal Asimov
   double poihat = 0;

   if (doFit)  {

      // do the fit : By doing a fit the POI snapshot (for S+B)  is set to the fit value
      // and the nuisance parameters nominal values will be set to the fit value.
      // This is relevant when using LEP test statistics

      RooArgSet constrainParams;
      if (sbModel->GetNuisanceParameters() ) constrainParams.add(*sbModel->GetNuisanceParameters());
      RooStats::RemoveConstantParameters(&constrainParams);
      tw.Start();
      RooFitResult * fitres = sbModel->GetPdf()->fitTo(*data,InitialHesse(false), Hesse(false),
                                                       Minimizer(MinimizerType.c_str(),"Migrad"), Strategy(0), PrintLevel(Verbose), Constrain(constrainParams), Save(true), Offset(RooStats::IsNLLOffset()) );
      if (fitres->status() != 0) {
         Warning("StandardHypoTestInv","Fit to the model failed - try with strategy 1 and perform first an Hesse computation");
         fitres = sbModel->GetPdf()->fitTo(*data,InitialHesse(true), Hesse(false),Minimizer(MinimizerType.c_str(),"Migrad"), Strategy(1), PrintLevel(Verbose+1), Constrain(constrainParams),
                                           Save(true), Offset(RooStats::IsNLLOffset()) );
      }
      if (fitres->status() != 0)
         Warning("StandardHypoTestInv"," Fit still failed - continue anyway.....");


      poihat  = poi->getVal();
      std::cout << "StandardHypoTestInv - Best Fit value : " << poi->GetName() << " = "
                << poihat << " +/- " << poi->getError() << std::endl;
      std::cout << "Time for fitting : "; tw.Print();

      //save best fit value in the poi snapshot
      sbModel->SetSnapshot(*sbModel->GetParametersOfInterest());
      std::cout << "StandardHypoTestInvo: snapshot of S+B Model " << sbModel->GetName()
                << " is set to the best fit value" << std::endl;

   }

   // build test statistics and hypotest calculators for running the inverter
   SimpleLikelihoodRatioTestStat slrts(*sbModel->GetPdf(),*bModel->GetPdf());

   // null parameters must includes snapshot of poi plus the nuisance values
   RooArgSet nullParams(*sbModel->GetSnapshot());
   if (sbModel->GetNuisanceParameters()) nullParams.add(*sbModel->GetNuisanceParameters());
   if (sbModel->GetSnapshot()) slrts.SetNullParameters(nullParams);
   RooArgSet altParams(*bModel->GetSnapshot());
   if (bModel->GetNuisanceParameters()) altParams.add(*bModel->GetNuisanceParameters());
   if (bModel->GetSnapshot()) slrts.SetAltParameters(altParams);
   if (EnableDetOutput) slrts.EnableDetailedOutput();

   // ratio of profile likelihood - need to pass snapshot for the alt
   RatioOfProfiledLikelihoodsTestStat
      ropl(*sbModel->GetPdf(), *bModel->GetPdf(), bModel->GetSnapshot());
   ropl.SetSubtractMLE(false);
   if (testStatType == 11) ropl.SetSubtractMLE(true);
   ropl.SetPrintLevel(Verbose);
   ropl.SetMinimizer(MinimizerType.c_str());
   if (EnableDetOutput) ropl.EnableDetailedOutput();

   ProfileLikelihoodTestStat profll(*sbModel->GetPdf());
   if (testStatType == 3) profll.SetOneSided(true);
   if (testStatType == 4) profll.SetSigned(true);
   profll.SetMinimizer(MinimizerType.c_str());
   profll.SetPrintLevel(Verbose);
   if (EnableDetOutput) profll.EnableDetailedOutput();

   profll.SetReuseNLL(Optimize);
   slrts.SetReuseNLL(Optimize);
   ropl.SetReuseNLL(Optimize);

   if (Optimize) {
      profll.SetStrategy(0);
      ropl.SetStrategy(0);
      ROOT::Math::MinimizerOptions::SetDefaultStrategy(0);
   }

   if (MaxPoi > 0) poi->setMax(MaxPoi);  // increase limit

   MaxLikelihoodEstimateTestStat maxll(*sbModel->GetPdf(),*poi);
   NumEventsTestStat nevtts;

   AsymptoticCalculator::SetPrintLevel(Verbose);

   // create the HypoTest calculator class
   HypoTestCalculatorGeneric *  hc = 0;
   if (type == 0) hc = new FrequentistCalculator(*data, *bModel, *sbModel);
   else if (type == 1) hc = new HybridCalculator(*data, *bModel, *sbModel);
   else if (type == 2 ) hc = new AsymptoticCalculator(*data, *bModel, *sbModel, false );
   else if (type == 3 ) hc = new AsymptoticCalculator(*data, *bModel, *sbModel, true );  // for using Asimov data generated with nominal values
   else {
      Error("StandardHypoTestInv","Invalid - calculator type = %d supported values are only :\n\t\t\t 0 (Frequentist) , 1 (Hybrid) , 2 (Asymptotic) ",type);
      return 0;
   }

   // set the test statistic
   TestStatistic * testStat = 0;
   if (testStatType == 0) testStat = &slrts;
   if (testStatType == 1 || testStatType == 11) testStat = &ropl;
   if (testStatType == 2 || testStatType == 3 || testStatType == 4) testStat = &profll;
   if (testStatType == 5) testStat = &maxll;
   if (testStatType == 6) testStat = &nevtts;

   if (testStat == 0) {
      Error("StandardHypoTestInv","Invalid - test statistic type = %d supported values are only :\n\t\t\t 0 (SLR) , 1 (Tevatron) , 2 (PLR), 3 (PLR1), 4(MLE)",testStatType);
      return 0;
   }


   ToyMCSampler *toymcs = (ToyMCSampler*)hc->GetTestStatSampler();
   if (toymcs && (type == 0 || type == 1) ) {
      // look if pdf is number counting or extended
      if (sbModel->GetPdf()->canBeExtended() ) {
         if (useNumberCounting)   Warning("StandardHypoTestInv","Pdf is extended: but number counting flag is set: ignore it ");
      }
      else {
         // for not extended pdf
         if (!useNumberCounting  )  {
            int nEvents = data->numEntries();
            Info("StandardHypoTestInv","Pdf is not extended: number of events to generate taken  from observed data set is %d",nEvents);
            toymcs->SetNEventsPerToy(nEvents);
         }
         else {
            Info("StandardHypoTestInv","using a number counting pdf");
            toymcs->SetNEventsPerToy(1);
         }
      }

      toymcs->SetTestStatistic(testStat);

      if (data->isWeighted() && !GenerateBinned) {
         Info("StandardHypoTestInv","Data set is weighted, nentries = %d and sum of weights = %8.1f but toy generation is unbinned - it would be faster to set mGenerateBinned to true\n",data->numEntries(), data->sumEntries());
      }
      toymcs->SetGenerateBinned(GenerateBinned);

      toymcs->SetUseMultiGen(Optimize);

      if (GenerateBinned &&  sbModel->GetObservables()->getSize() > 2) {
         Warning("StandardHypoTestInv","generate binned is activated but the number of observable is %d. Too much memory could be needed for allocating all the bins",sbModel->GetObservables()->getSize() );
      }

      // set the random seed if needed
      if (RandomSeed >= 0) RooRandom::randomGenerator()->SetSeed(RandomSeed);

   }

   // specify if need to re-use same toys
   if (ReuseAltToys) {
      hc->UseSameAltToys();
   }

   if (type == 1) {
      HybridCalculator *hhc = dynamic_cast<HybridCalculator*> (hc);
      assert(hhc);

      hhc->SetToys(ntoys,ntoys/NToysRatio); // can use less ntoys for b hypothesis

      // remove global observables from ModelConfig (this is probably not needed anymore in 5.32)
      bModel->SetGlobalObservables(RooArgSet() );
      sbModel->SetGlobalObservables(RooArgSet() );


      // check for nuisance prior pdf in case of nuisance parameters
      if (bModel->GetNuisanceParameters() || sbModel->GetNuisanceParameters() ) {

         // fix for using multigen (does not work in this case)
         toymcs->SetUseMultiGen(false);
         ToyMCSampler::SetAlwaysUseMultiGen(false);

         RooAbsPdf * nuisPdf = 0;
         if (nuisPriorName) nuisPdf = w->pdf(nuisPriorName);
         // use prior defined first in bModel (then in SbModel)
         if (!nuisPdf)  {
            Info("StandardHypoTestInv","No nuisance pdf given for the HybridCalculator - try to deduce  pdf from the model");
            if (bModel->GetPdf() && bModel->GetObservables() )
               nuisPdf = RooStats::MakeNuisancePdf(*bModel,"nuisancePdf_bmodel");
            else
               nuisPdf = RooStats::MakeNuisancePdf(*sbModel,"nuisancePdf_sbmodel");
         }
         if (!nuisPdf ) {
            if (bModel->GetPriorPdf())  {
               nuisPdf = bModel->GetPriorPdf();
               Info("StandardHypoTestInv","No nuisance pdf given - try to use %s that is defined as a prior pdf in the B model",nuisPdf->GetName());
            }
            else {
               Error("StandardHypoTestInv","Cannot run Hybrid calculator because no prior on the nuisance parameter is specified or can be derived");
               return 0;
            }
         }
         assert(nuisPdf);
         Info("StandardHypoTestInv","Using as nuisance Pdf ... " );
         nuisPdf->Print();

         const RooArgSet * nuisParams = (bModel->GetNuisanceParameters() ) ? bModel->GetNuisanceParameters() : sbModel->GetNuisanceParameters();
         RooArgSet * np = nuisPdf->getObservables(*nuisParams);
         if (np->getSize() == 0) {
            Warning("StandardHypoTestInv","Prior nuisance does not depend on nuisance parameters. They will be smeared in their full range");
         }
         delete np;

         hhc->ForcePriorNuisanceAlt(*nuisPdf);
         hhc->ForcePriorNuisanceNull(*nuisPdf);


      }
   }
   else if (type == 2 || type == 3) {
      if (testStatType == 3) ((AsymptoticCalculator*) hc)->SetOneSided(true);
      if (testStatType != 2 && testStatType != 3)
         Warning("StandardHypoTestInv","Only the PL test statistic can be used with AsymptoticCalculator - use by default a two-sided PL");
   }
   else if (type == 0 ) {
      ((FrequentistCalculator*) hc)->SetToys(ntoys,ntoys/NToysRatio);
      // store also the fit information for each poi point used by calculator based on toys
      if (EnableDetOutput) ((FrequentistCalculator*) hc)->StoreFitInfo(true);
   }
   else if (type == 1 ) {
      ((HybridCalculator*) hc)->SetToys(ntoys,ntoys/NToysRatio);
   }

   // Get the result
   RooMsgService::instance().getStream(1).removeTopic(RooFit::NumIntegration);



   HypoTestInverter calc(*hc);
   calc.SetConfidenceLevel(optHTInv.confLevel);


   calc.UseCLs(useCLs);
   calc.SetVerbose(true);

   // can speed up using proof-lite
   if (UseProof) {
      ProofConfig pc(*w, NWorkers, "", kFALSE);
      toymcs->SetProofConfig(&pc);    // enable proof
   }


   if (npoints > 0) {
      if (poimin > poimax) {
         // if no min/max given scan between MLE and +4 sigma
         poimin = int(poihat);
         poimax = int(poihat +  4 * poi->getError());
      }
      std::cout << "Doing a fixed scan  in interval : " << poimin << " , " << poimax << std::endl;
      calc.SetFixedScan(npoints,poimin,poimax);
   }
   else {
      //poi->setMax(10*int( (poihat+ 10 *poi->getError() )/10 ) );
      std::cout << "Doing an  automatic scan  in interval : " << poi->getMin() << " , " << poi->getMax() << std::endl;
   }

   tw.Start();
   HypoTestInverterResult * r = calc.GetInterval();
   std::cout << "Time to perform limit scan \n";
   tw.Print();

   return r;
}



void ReadResult(const char * infile, const char * resultName="", bool useCLs=true) {
  TString filename(infile);
  TFile *file = TFile::Open(filename);
  if(!file ){
     cout <<"StandardRooStatsDemoMacro: Input file " << filename << " not found" << endl;
     exit(0);
  }

  HypoTestInverterResult * r = dynamic_cast<HypoTestInverterResult*>( file->Get(resultName) );
  if ( r == NULL ) {
    std::cerr << "File " << filename << " does not contain a HypoTestInverterResult named " << resultName << std::endl;
    exit(0);
  }

  HypoTestInvTool calc(&optHTInv);
  calc.WriteResult = false;
  cout << "Calculator type?\n";
  int calculatorType = 0;
  cin >> calculatorType;
  cout << "Test stat type?\n";
  int testStatType = 0;
  cin >> testStatType;

  calc.AnalyzeResult( r, calculatorType, testStatType, useCLs, 0, "");
}


#ifdef USE_AS_MAIN
int main() {
    StandardHypoTestInv();
}
#endif
