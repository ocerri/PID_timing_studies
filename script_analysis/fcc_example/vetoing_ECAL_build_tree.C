#include "mumu_invisible.h"


/*
root -l Analysis/mumu_invisible_analysis.C'(process)'
*/

//------------------------------------------------------------------------------


//pythia xsec must be in millibarn
void vetoing_ECAL_build_tree(int process=1, TString namein = "", TString nameout = "")
{

// Setting input and output files and cross sections ------------------------------

  if (process==1){
    const char *inputFile = "_root/HZ_H2inv_Z2lep_1M_CMSdelphes.root";
    const char *outputFile = "_root/tree_HZ_H2inv_Z2lep_1M_CMSdelphes_analysis_vetoing.root";
  }
  else{
    const char *inputFile = namein;
    const char *outputFile = nameout;
  }


//Setting style-----------------------------------------------

  //labstyle();

  cout << "Input file = " << inputFile << endl;


  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  Double_t mass, miss, angle, energy, mass_lp;
  Int_t n_mu;

  TTree t("tree_veto","Tree with HZ events data");
  t.Branch("mass", &mass,"mass/D");
  t.Branch("miss", &miss,"miss/D");
  t.Branch("angle", &angle,"angle/D");
  t.Branch("energy", &energy,"energy/D");
  t.Branch("n_mu", &n_mu, "n_mu/I");
  t.Branch("mass_lp", &mass_lp,"mass_lp/D");



  // Loop over all events--------------------------------------------------------------------------------------------------------
  cout << "-----------Starting loop over events--------------" << endl;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    //Progress bar generator
    show_progress_bar(entry,numberOfEntries);

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // If event contains 2 muons
    if(branchMuon->GetEntries() == 2 && branchPhoton->GetEntries()==1)
    {
      Muon *muon1 = (Muon *) branchMuon->At(0);
      Muon *muon2 = (Muon *) branchMuon->At(1);
      Photon *gamma1 = (Photon *) branchPhoton->At(0);

      TLorentzVector v1 = muon1->P4();
      TLorentzVector v2 = muon2->P4();
      TLorentzVector g1 = gamma1->P4();

      // Plot their invariant mass
      mass_lp = (v1 + v2 + g1).M();
      mass = (v1+v2).M();
      miss = (P_in - v1 - v2).M();
      energy = g1.E();
      angle = min(v1.Angle(g1.Vect()) , v2.Angle(g1.Vect()) );
      n_mu = 2;
      t.Fill();
    }
    else if(branchElectron->GetEntries() == 2 && branchPhoton->GetEntries()==1)
    {
      Electron *electron1 = (Electron *) branchElectron->At(0);
      Electron *electron2 = (Electron *) branchElectron->At(1);
      Photon *gamma1 = (Photon *) branchPhoton->At(0);

      TLorentzVector v1 = electron1->P4();
      TLorentzVector v2 = electron2->P4();
      TLorentzVector g1 = gamma1->P4();

      mass_lp = (v1 + v2 + g1).M();
      mass = (v1+v2).M();
      miss = (P_in - v1 - v2).M();
      energy = g1.E();
      angle = min( v1.Angle(g1.Vect()) , v2.Angle(g1.Vect()) );
      n_mu = 0;
      t.Fill();
    }
    else if(branchMuon->GetEntries() == 1 && branchElectron->GetEntries()==1 && branchPhoton->GetEntries()==0)
    {
      Muon *muon1 = (Muon *) branchMuon->At(0);
      Electron *electron1 = (Electron *) branchElectron->At(0);

      TLorentzVector v1 = muon1->P4();
      TLorentzVector v2 = electron1->P4();

      mass_lp = 0;
      mass = (v1+v2).M();
      miss = (P_in - v1 - v2).M();
      energy = -1;
      angle = -1;
      n_mu = 1;
      t.Fill();
    }
  }
  //End loop----------------------------------------------------------------------------------------------------------------------

  cout << "-----------Loop on event terminated---------------" << endl << endl << endl;
  cout << "Total event : " << numberOfEntries/1000 << " k" << endl;



  //Saving on root file
  TFile *file = TFile::Open(outputFile, "RECREATE");
  file->cd();
  t.Write();
  file->Close();

  cout << "Output .root file wrote in :  < " << outputFile << " >" << endl<< endl;


}
