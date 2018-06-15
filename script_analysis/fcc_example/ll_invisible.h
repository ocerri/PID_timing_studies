#include <cstdio>
#include <stdio.h>
#include <stdlib.h>
#include "TFile.h"
#include <TCanvas.h>
#include <iostream>
#include <fstream>  // necessario per scrivere su file
#include <string>
#include <TRandom3.h>
#include <TMinuit.h>
#include "TH1.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TStyle.h"
#include "TMath.h"
#include "TMatrixDSym.h"
#include "TFitResult.h"
#include <vector>
#include "TArc.h"
#include "TFitter.h"
#include "TVirtualFitter.h"
#include "TLegend.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include <TRandom3.h>

using namespace std;



//------------------------------------------------------------------------------
Double_t M_Z = 91.19;
Double_t dM_Z = 5.; //Same as Patrick paper
Double_t M_Hp = 135.0;
Double_t M_Hm = 115.0;
TLorentzVector P_in = TLorentzVector(0.,0.,0.,240.0);

typedef struct {
  Int_t is_mumu;
  Int_t is_ee;
  Int_t is_eeph;
  Int_t is_mumuph;

  Double_t M_ll;
  Double_t M_miss;
  Double_t Pt_ll;
  Double_t Pl_ll;

  Double_t dtheta_ll_lab; //angle between leptons in the lab frame
  Double_t helicity_angle;
  Double_t acoplanarity_angle;
  Double_t sin2dthetahalf;

  Double_t pmaxpmin; //product of muons momentum
  Double_t pmax;
  Double_t pmin;
} Event_ll_t;

Int_t Jet_lepton_tagging(TLorentzVector jet, TLorentzVector lep){
  Double_t angle = lep.Angle(jet.Vect());
  Double_t dPRatio = 2*std::fabs(jet.P()-lep.P())/(jet.P()+lep.P());

  if (angle < 0.3 && dPRatio < 0.2) return 1;
  return 0;
}

void show_progress_bar(int entry, int numberOfEntries){
	Int_t step = numberOfEntries>200 ?  numberOfEntries/ 200 : 2;
	if (entry==0 || entry%step==0 || entry==numberOfEntries-1){
      if (entry==0){
    	cout << "[--------------------]  0%";
      }
      else if((entry+1)/(float)numberOfEntries == 1){
    	cout << '\r' << "[####################]  100%"<< endl;
      }
      else if((entry+1)/(float)numberOfEntries >= 0.95){
    	cout << '\r' << "[###################-]  95%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.90){
    	cout << '\r' << "[##################--]  90%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.85){
    	cout << '\r' << "[#################---]  85%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.80){
    	cout << '\r' << "[################----]  80%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.75){
    	cout << '\r' << "[###############-----]  75%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.70){
    	cout << '\r' << "[##############------]  70%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.65){
    	cout << '\r' << "[#############-------]  65%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.60){
    	cout << '\r' << "[############--------]  60%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.55){
    	cout << '\r' << "[###########---------]  55%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.50){
    	cout << '\r' << "[##########----------]  50%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.45){
    	cout << '\r' << "[#########-----------]  45%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.40){
    	cout << '\r' << "[########------------]  40%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.35){
    	cout << '\r' << "[#######-------------]  35%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.30){
    	cout << '\r' << "[######--------------]  30%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.25){
    	cout << '\r' << "[#####---------------]  25%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.20){
    	cout << '\r' << "[####----------------]  20%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.15){
    	cout << '\r' << "[###-----------------]  15%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.10){
    	cout << '\r' << "[##------------------]  10%";
      }
      else if((entry+1)/(float)numberOfEntries >= 0.05){
    	cout << '\r' << "[#-------------------]  5%";
      }
      else if((entry+1)/(float)numberOfEntries < 0.05){
    	cout << '\r' << "[--------------------]  0%";
      }
    }
    cout << flush;
}
