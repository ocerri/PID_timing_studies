#include "mumu_invisible.h"


/*
root -l Analysis/mumu_invisible_analysis.C'(process)'
*/

//------------------------------------------------------------------------------


//pythia xsec must be in millibarn
void HCAL_ECAL_vetoing_analysis(int process=1, TString namein = "", Double_t xsec_in =1, TString nameout = "")
{

// Setting input and output files and cross sections ------------------------------

  if (process==1){
    const char *inputFile = "_root/HZ_H2inv_Z2lep_1M_CMSdelphes.root";
    Double_t pythiaXsec = 5.720e-11;
    const char *outputFile = "_root/HZ_H2inv_Z2lep_1M_CMSdelphes_analysis_vetoing.root";
  }
  else{
    const char *inputFile = namein;
    Double_t pythiaXsec = xsec_in;
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
  Long64_t numberOfEntries = treeReader->GetEntries()/10;
  Double_t scale_factor = Luminosity*1e12*pythiaXsec/numberOfEntries;

  // Get pointers to branches used in this analysis
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");

  // Book histograms
  {
    TH1 *h_mass_mu = new TH1F("h_mass_mu","h_mass_mu", 120, 0.0, 240.0);
    h_mass_mu->SetXTitle("M_{#mu#mu}");
    h_mass_mu->SetYTitle("Events / 2GeV");

    TH1 *h_miss_mu = new TH1F("h_miss_mu","h_miss_mu", 120, 0.0, 240.0);
    h_miss_mu->SetXTitle("Missing mass [GeVß]");
    h_miss_mu->SetYTitle("Events / 2GeV");

	  TH1 *h_photon_energy_mu = new TH1F("h_photon_energy_mu", "photon_energy_mu", 120, 0.0, 240.0);
	  h_photon_energy_mu->SetXTitle("E_{photon} [GeV]");
	  h_photon_energy_mu->SetYTitle("Events / 2GeV");
	  //h_photon_energy->Sumw2();

    TH1 *h_photon_angle_mu = new TH1F("h_photon_angle_mu", "photon_angle_mu", 90, 0.0, 180.0);
	  h_photon_angle_mu->SetXTitle("E_{photon} [GeV]");
	  h_photon_angle_mu->SetYTitle("Events / 2GeV");

    TH1 *h_mass_ele = new TH1F("h_mass_ele","h_mass_ele", 120, 0.0, 240.0);
    h_mass_ele->SetXTitle("M_{ee}");
    h_mass_ele->SetYTitle("Events / 2GeV");

    TH1 *h_miss_ele = new TH1F("h_miss_ele","h_miss_ele", 120, 0.0, 240.0);
    h_miss_ele->SetXTitle("Missing mass [GeV]");
    h_miss_ele->SetYTitle("Events / 2GeV");

    TH1 *h_photon_energy_ele = new TH1F("h_photon_energy_ele", "photon_energy_ele", 120, 0.0, 240.0);
	  h_photon_energy_ele->SetXTitle("E_{photon} [GeV]");
	  h_photon_energy_ele->SetYTitle("Events / 2GeV");
	  //h_photon_energy->Sumw2();

    TH1 *h_photon_angle_ele = new TH1F("h_photon_angle_ele", "photon_angle_ele", 90, 0.0, 180.0);
	  h_photon_angle_ele->SetXTitle("E_{photon} [GeV]");
	  h_photon_angle_ele->SetYTitle("Events / 2GeV");

    TH1 *h_mass_emu = new TH1F("h_mass_emu","h_mass_emu", 120,0.0,240.0);
    h_mass_emu->SetXTitle("Invariant mass e#mu pair [GeV}]");
    h_mass_emu->SetYTitle("Events / 2GeV");

    TH1 *h_miss_emu = new TH1F("h_M_emu","h_M_emu", 120,0.0,240.0);
    h_miss_emu->SetXTitle("Missing mass e#mu pair ev [GeV}]");
    h_miss_emu->SetYTitle("Events / 2GeV");

	}



  // Loop over all events--------------------------------------------------------------------------------------------------------
  cout << "-----------Starting loop over events--------------" << endl;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {

    //Progress bar generator
    show_progress_bar(entry,numberOfEntries);

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    Muon *muon1, *muon2;

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
      Double_t Mass = (v1+v2).M();
      h_mass_mu->Fill(Mass);

      Double_t M_miss = (P_in - v1 - v2).M();
      h_miss_mu->Fill(M_miss);

      h_photon_energy_mu->Fill(g1.E());

      Double_t angle = min(v1.Angle(g1.Vect()) , v2.Angle(g1.Vect()) );
      h_photon_angle_mu->Fill(angle);
    }
    else if(branchElectron->GetEntries() == 2 && branchPhoton->GetEntries()==1)
    {
      Electron *electron1 = (Electron *) branchElectron->At(0);
      Electron *electron2 = (Electron *) branchElectron->At(1);
      Photon *gamma1 = (Photon *) branchPhoton->At(0);

      TLorentzVector v1 = electron1->P4();
      TLorentzVector v2 = electron2->P4();
      TLorentzVector g1 = gamma1->P4();

      // Plot their invariant mass
      Double_t Mass = (v1+v2).M();
      h_mass_ele->Fill(Mass);

      Double_t M_miss = (P_in - v1 - v2).M();
      h_miss_ele->Fill(M_miss);

      h_photon_energy_ele->Fill(g1.E());

      Double_t angle = min(v1.Angle(g1.Vect()) , v2.Angle(g1.Vect()) );
      h_photon_angle_ele->Fill(angle);
    }
    else if(branchMuon->GetEntries() == 1 && branchElectron->GetEntries()==1 && branchPhoton->GetEntries()==0)
    {
      Muon *muon1 = (Muon *) branchMuon->At(0);
      Electron *electron1 = (Electron *) branchElectron->At(0);

      TLorentzVector v1 = muon1->P4();
      TLorentzVector v2 = electron1->P4();

      // Plot their invariant mass
      Double_t Mass = (v1+v2).M();
      h_mass_emu->Fill(Mass);

      Double_t M_miss = (P_in - v1 - v2).M();
      h_miss_emu->Fill(M_miss);
    }
  }
  //End loop----------------------------------------------------------------------------------------------------------------------

  cout << "-----------Loop on event terminated---------------" << endl << endl << endl;
  cout << "Total event : " << numberOfEntries/1000 << " k" << endl;


  /*
  TCanvas *c_2pp = new TCanvas("2pp","pp",600,600);
  c_2pp->cd();
  h2_pp->Draw("LEGO2Z 0");

  TCanvas *c_sin2vspp = new TCanvas("sin2vspp","sin2vspp",600,600);
  c_sin2vspp->cd();
  h2_sin2vspp->Draw("LEGO2Z 0");

  //TCanvas *c_missing_mass = new TCanvas("Missing_Mass","Missing Mass",600,600);
  //c_missing_mass->cd();
  //h_missing_mass->Draw();

  TCanvas *c_missing_mass_Ztag = new TCanvas("Missing_Mass_Ztag","Missing Mass in Z tagged events",600,600);
  c_missing_mass_Ztag->cd();
  h_missing_mass_Ztag->Draw();

  TCanvas *c_muon_angle = new TCanvas("Muons_angle","Anlge between muons",600,600);
  c_muon_angle->cd();
  h_angle_muons->Draw();

  TCanvas *c_pvsp = new TCanvas("pvsp","pvsp",600,600);
  c_pvsp->cd();
  g_pmaxvspmin->SetLineColor(0);
  g_pmaxvspmin->Draw("ap");

  TCanvas *c_pvssin = new TCanvas("pvssin","pvssin",600,600);
  c_pvssin->cd();
  g_sinvspp->SetLineColor(0);
  g_sinvspp->Draw("ap");

  TCanvas *c_elicity_angle = new TCanvas("elicity_angle","Elicity angle",600,600);
  c_elicity_angle->cd();
  h_elicity_angle->Draw();
  */


  //Saving on root file
  TFile *file = TFile::Open(outputFile, "RECREATE");
  file->cd();
  h_mass_mu->Write();
  h_miss_mu->Write();
  h_photon_angle_mu->Write();
  h_photon_energy_mu->Write();
  h_mass_ele->Write();
  h_miss_ele->Write();
  h_photon_angle_ele->Write();
  h_photon_energy_ele->Write();
  h_mass_emu->Write();
  h_miss_emu->Write();
  file->Close();

  cout << "Output .root file wrote in :  < " << outputFile << " >" << endl<< endl;


}
