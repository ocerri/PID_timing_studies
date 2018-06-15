#include "ll_invisible.h"

//Macro to build tree from delphes output----------------------------------------
//this is the principal macro to analyze ll events
// beam_en_spread espresso in MeV
void ll_invisible_build_tree(Double_t beam_en_spread=288, const char* inputFile="_root/HZ_H2inv_Z2lep_2M_ILDdelphes.root")
{
  TRandom3 randGen(0);

  TString out_name = inputFile;
  out_name.Insert(6,"tree_");
  // if(beam_en_spread >= 0)
  // {
  out_name.Insert(out_name.Length()-5,Form("BES%.0fMeV_ll",beam_en_spread));
  // }
  // else
  // {
    // out_name.Insert(out_name.Length()-5,"_ll_noBES"); //no beam energy spread
  // }
  const char* outputFile = out_name;
  cout << "Output file: " << outputFile << endl;

  cout << "Applied cuts:" << endl;
  cout << "- Only 2 electrons or only 2 muons" << endl;
  cout << "- Only 2 electrons or only 2 muons and only 1 photon" << endl;
  vector<Int_t> event_ll;

  gSystem->Load("libDelphes");

  // Create chain of root trees and create object of class ExRootTreeReader
  TChain chain("Delphes");
  chain.Add(inputFile);
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

  // Get pointers to branches used in this analysis
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchJet = treeReader->UseBranch("Jet");


  // Loop over all events--------------------------------------------------------------------------------------------------------
  cout << "---------------Starting loop over all events----------------" << endl << endl;
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    //Progress bar generator
    show_progress_bar(entry,numberOfEntries);

    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    Int_t Nmu = branchMuon->GetEntries();
    Int_t Ne = branchElectron->GetEntries();
    Int_t Npho = branchPhoton->GetEntries();
    Int_t Njet = branchJet->GetEntries();

    Bool_t ev2mu = (Nmu == 2 && Njet==0 && Ne == 0 && Npho < 2);
    Bool_t ev2e = (Nmu == 0 && Njet==0 && Ne == 2 && Npho < 2);

    if( ev2mu || ev2e ){
      event_ll.push_back(entry);}
  }
  //End loop----------------------------------------------------------------------------------------------------------------------

  cout << endl << "-----------------------Loop terminated----------------------" << endl << endl << endl;
  cout << "Total event : " << numberOfEntries/1000 << " k" << endl;
  cout << "Total event with 2 leptons detected : " << event_ll.size() << endl;
  cout << "Trigger efficency : " << 100*event_ll.size()/(float)numberOfEntries << " %" << endl;

  TTree t("tree_ll","Tree with ll events data");
  Event_ll_t ev;

  t.Branch("is_mumu", &ev.is_mumu, "is_mumu/I");
  t.Branch("is_ee", &ev.is_ee, "is_ee/I");
  t.Branch("is_mumuph", &ev.is_mumuph, "is_mumuph/I");
  t.Branch("is_eeph", &ev.is_eeph, "is_eeph/I");

  t.Branch("M_ll", &ev.M_ll,"M_ll/D");
  t.Branch("M_miss", &ev.M_miss,"M_miss/D");
  t.Branch("Pt_ll", &ev.Pt_ll,"Pt_ll/D");
  t.Branch("Pl_ll", &ev.Pl_ll,"Pl_ll/D");

  t.Branch("dtheta_ll_lab", &ev.dtheta_ll_lab,"dtheta_ll_lab/D");
  t.Branch("helicity_angle", &ev.helicity_angle,"helicity_angle/D");
  t.Branch("acoplanarity_angle", &ev.acoplanarity_angle,"acoplanarity_angle/D");
  t.Branch("sin2dthetahalf", &ev.sin2dthetahalf,"sin2dthetahalf/D");

  t.Branch("pmaxpmin", &ev.pmaxpmin,"pmaxpmin/D");
  t.Branch("pmax", &ev.pmax,"pmax/D");
  t.Branch("pmin", &ev.pmin,"pmin/D");



  cout << endl << "------Loop over HZ candidate events and tree building-------"<< endl << endl;
  //Start LOOP over HZ candidate event only
  //storing loop variable
  Int_t Nmu; //Number of muons
  Int_t Ne; //Number of electrons
  Int_t Npho; //Number of photond
  TLorentzVector v1; //positive lepton
  TLorentzVector v2;//neative lepton
  TLorentzVector g1;//photon if any
  TLorentzVector P_vis; //Total visible 4-vector
  for(Int_t i=0; i<event_ll.size(); ++i){

    show_progress_bar(i, event_ll.size());

    treeReader->ReadEntry(event_ll[i]);
    Nmu = branchMuon->GetEntries();
    Ne = branchElectron->GetEntries();
    Npho = branchPhoton->GetEntries();

    ev.is_mumu = (Nmu == 2 && Ne == 0 && Npho == 0);
    ev.is_ee = (Nmu == 0 && Ne == 2 && Npho == 0);
    ev.is_mumuph = (Nmu == 2 && Ne == 0 && Npho == 1);
    ev.is_eeph = (Nmu == 0 && Ne == 2 && Npho == 1);

    //setting local variable
    if(ev.is_mumu){
      Muon *muon1 = (Muon *) branchMuon->At(0);
      Muon *muon2 = (Muon *) branchMuon->At(1);
      if(muon1->Charge < 0) {
        Muon *muon_tmp= muon1;
        muon1=muon2;
        muon2=muon_tmp;
      }
      v1 = muon1->P4();
      v2 = muon2->P4();
      P_vis = v1 +v2;
    }
    else if(ev.is_ee){
      Electron *el1 = (Electron *) branchElectron->At(0);
      Electron *el2 = (Electron *) branchElectron->At(1);
      if(el1->Charge < 0) {
        Electron *el_tmp= el1;
        el1=el2;
        el2=el_tmp;
      }
      v1 = el1->P4();
      v2 = el2->P4();
      P_vis = v1 +v2;
    }
    else if(ev.is_mumuph){
      Muon *muon1 = (Muon *) branchMuon->At(0);
      Muon *muon2 = (Muon *) branchMuon->At(1);
      Photon *ph1 = (Photon *) branchPhoton->At(0);
      if(muon1->Charge < 0) {
        Muon *muon_tmp= muon1;
        muon1=muon2;
        muon2=muon_tmp;
      }
      v1 = muon1->P4();
      v2 = muon2->P4();
      g1 = ph1->P4();
      P_vis = v1 + v2 + g1;
    }
    else if(ev.is_eeph){
      Electron *el1 = (Electron *) branchElectron->At(0);
      Electron *el2 = (Electron *) branchElectron->At(1);
      Photon *ph1 = (Photon *) branchPhoton->At(0);
      if(el1->Charge < 0) {
        Electron *el_tmp= el1;
        el1=el2;
        el2=el_tmp;
      }
      v1 = el1->P4();
      v2 = el2->P4();
      g1 = ph1->P4();
      P_vis = v1 + v2 + g1;
    }

    //Create momentum variable

    if(beam_en_spread > 0)
    {
      //Double_t sqrt_s_smeared = 240.0;
      Double_t sqrt_s_smeared = 240.0 + randGen.Gaus(0,beam_en_spread/1e3); //1.2/1000 gaussian smearing (288MeV)
    }
    else {
      Double_t sqrt_s_smeared = 240.0; //non smeared
    }
    P_in.SetE(sqrt_s_smeared);

    ev.M_ll = P_vis.M();
    ev.M_miss = (P_vis - P_in).M();
    ev.Pt_ll = P_vis.Pt();
    ev.Pl_ll = P_vis.Pz();

    ev.pmax = TMath::Max(v1.P() , v2.P());
    ev.pmin = TMath::Min(v1.P() , v2.P());
    ev.pmaxpmin = v1.P()*v2.P();


    //Create angular Variable
    Double_t angle = v1.Angle(v2.Vect());
    ev.dtheta_ll_lab = angle *180/TMath::Pi();
    ev.sin2dthetahalf = sin(angle/2)*sin(angle/2);

    TVector3 vp = v1.Vect().Cross( v2.Vect() );
    ev.acoplanarity_angle = 180*vp.Theta() / 3.14159;

    v1.Boost(-P_vis.BoostVector());
    ev.helicity_angle = v1.Angle(P_vis.Vect()) *180/3.14159;

    //Fill the tree
    t.Fill();
  }
  cout << endl << "-----------------------Loop terminated----------------------" << endl;

  TFile tree_file(outputFile,"recreate");
  t.Write();
  //tree_file.Close();
  cout << "Output tree.root file wrote in :  < " << outputFile << " >" << endl<< endl;
}
