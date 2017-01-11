#include <iostream>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "ObjectID.C"
#include "TLorentzVector.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z

int getEtaBin(float abseta){
  int bin=-1;
  if(abseta>2.0) bin=6;
  else if(abseta>1.556) bin=5;
  else if(abseta>1.442) bin=4;
  else if(abseta>0.8) bin=3;
  else if(abseta>0.4) bin=2;
  else if(abseta>=0.0) bin=1;

  return bin;
}


//Feed this method a given eta and  your weight array, and this will return the corresponding weight
double EtaWeight(double* weights,double eta,double pt){
  int fbin = getEtaBin( abs(eta) );
  if (pt > 100.0){ //high pt 
    fbin+=6;
  }
  //std::cout<<"eta-pt bin "<<fbin<<std::endl;
  return weights[fbin-1];
};

void Ex_1p3(){
  /*add in charge misID rate you measured earlier, weights should conform to the binning for eta from getEtaBin:
    0.0-0.4 |  0.4-0.8 | 0.8-1.442 | 1.442-1.556 | 1.556-2.0 | 2.0 - infinity
    The six eta region weights for low pt (<100 GeV) should come first, followed by the high pt (>100 GeV) weights 
   */
  //uncomment when ready to use

  double weights[12] = {0.00026, 0.00021, 0.00106, 0.0, 0.00806, 0.00970,
			0.00240, 0.00202, 0.01165, 0.0, 0.02621, 0.03616};

  //load data file - This file is a sum of all data files with opposite-signed leptons
  //TFile* f = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  TFile* f = new TFile("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/ChargeMisID_MuCBTightMiniIsoTight_ElMVATightRC_2016B-H.root");
  TTree* t = (TTree*)f->Get("tEvts_ssdl");

  //output files
  TFile* fbg = new TFile("bg_chargeMisID.root","RECREATE");
  TTree* T = new TTree("T","test");

  int nEntries = t->GetEntries();

  int Nmetcut=0, NJetscut=0, NJetpt1cut=0, NJetpt2cut=0, NHTcut=0, NZcut=0, NMasscut=0, NlepPt1cut=0, NlepPt2cut=0;

  float lepPts1, lepEtas1, lepPhis1, lepEs1;
  float lepPts2, lepEtas2, lepPhis2, lepEs2;
  float jetPts1, jetEtas1, jetPhis1, jetEs1;
  float jetPts2, jetEtas2, jetPhis2, jetEs2;
  int lepFlavor1, lepFlavor2;
  int lepCharge1, lepCharge2;
  int lepTight1, lepTight2;
  int lepLoose1, lepLoose2;
  float assocMass, dilepMass, HT;
  float met = 0;

  //Set branch addresses
  t->SetBranchAddress("Lep1Pt", &lepPts1);
  t->SetBranchAddress("Lep1Eta", &lepEtas1);
  t->SetBranchAddress("Lep1Phi", &lepPhis1);
  t->SetBranchAddress("Lep1Energy", &lepEs1);
  t->SetBranchAddress("Lep1Charge", &lepCharge1);
  t->SetBranchAddress("Lep1Flavor", &lepFlavor1);
  //t->SetBranchAddress("Lep1Tight", &lepTight1);
  //t->SetBranchAddress("Lep1Loose", &lepLoose1);
  t->SetBranchAddress("Lep2Pt", &lepPts2);
  t->SetBranchAddress("Lep2Eta", &lepEtas2);
  t->SetBranchAddress("Lep2Phi", &lepPhis2);
  t->SetBranchAddress("Lep2Energy", &lepEs2);
  t->SetBranchAddress("Lep2Charge", &lepCharge2);
  t->SetBranchAddress("Lep2Flavor", &lepFlavor2);
  //t->SetBranchAddress("Lep2Tight", &lepTight2);
  //t->SetBranchAddress("Lep2Loose", &lepLoose2);
  t->SetBranchAddress("cleanAK4Jet1Pt", &jetPts1);
  t->SetBranchAddress("cleanAK4Jet1Eta", &jetEtas1);
  t->SetBranchAddress("cleanAK4Jet1Phi", &jetPhis1);
  t->SetBranchAddress("cleanAK4Jet1Energy", &jetEs1);
  t->SetBranchAddress("cleanAK4Jet2Pt", &jetPts2);
  t->SetBranchAddress("cleanAK4Jet2Eta", &jetEtas2);
  t->SetBranchAddress("cleanAK4Jet2Phi", &jetPhis2);
  t->SetBranchAddress("cleanAK4Jet2Energy", &jetEs2);
  t->SetBranchAddress("MET",&met);
  t->SetBranchAddress("AssocMass",&assocMass); //any two electrons that whose invariant mass is closest to Z mass 
  t->SetBranchAddress("DilepMass",&dilepMass); //invariant mass of two opposite sign leptons
  t->SetBranchAddress("cleanAK4HT",&HT); //invariant mass of two opposite sign leptons

  //Histograms
  //TH1F* ssEtaHist = new TH1F("ssEtaHist","#eta",30,-3,3);
  //TH1F* osEtaHist = new TH1F("osEtaHist","#eta",30,-3,3);

  //TH1F* ssPtHist = new TH1F("ssPtHist","p_{T}",100,0.,200.);
  //TH1F* osPtHist = new TH1F("osPtHist","p_{T}",100,0.,200.);

  TH1F* ssHTHist = new TH1F("ssHTHist","Weight H_{T} of charge misID events",40,0.,2000);
  TH1F* osHTHist = new TH1F("osHTHist","H_{T} of opposite sign dilepton events",40,0.,2000);


  //Loop over the tree and look for an electron pair that makes a Z
  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;

   //require more than one jet
    if(!jetPts2) continue;
    NJetscut +=1;

    if( HT < 900) continue;
    NHTcut+=1;

    if (assocMass > M_Z - dM && assocMass < M_Z + dM) continue;
    if (dilepMass > M_Z - dM && dilepMass < M_Z + dM) continue;
    NZcut+=1;
    
    if(lepPts1<40) continue;
    NlepPt1cut+=1;

    if(lepPts2<30) continue;
    NlepPt2cut+=1;

    //cut on dilepton mass
    if (dilepMass < 20) continue;
    NMasscut+=1;

    //check channel
    bool ee = false;
    bool emu = false;
    bool mumu = false;

    if( (lepFlavor1 == 0) && (lepFlavor2 == 0) ) ee = true;
    if( (lepFlavor1 == 0) && (lepFlavor2 == 1) ) emu = true;
    if( (lepFlavor1 == 1) && (lepFlavor2 == 0) ) emu = true;
    if( (lepFlavor1 == 1) && (lepFlavor2 == 1) ) mumu = true;

    //skip any di-mu events since we assume they don't contribute to ss via charge misID
    if(mumu) continue;

    /*now we want to fill our predicted HT distribution with the os events weighted by the appropriate weights
      I've written a simple function that gives the correct element of the weight array if you pass it the eta of the electron.
      Usage is simply EtaWeight(weightsarray, eta,pt). Take advantage of it to make life easier :)
     */

    float ee_weight;
    float emu_weight;

    // ADD CODE HERE TO CALCULATE PROBABLITY FOR EACH EVENT


    //fill ss HT histogram
    if(ee){
      ssHTHist->Fill(HT,ee_weight);
    }
    if(emu){
      ssHTHist->Fill(HT,emu_weight);
    }

    //fill os hist for reference
    osHTHist->Fill(HT);

    //end event loop
  }

  printf("Num of events passing Njet cut         : %d\n", NJetscut);
  printf("Num of events passing lep1pt cut      : %d\n", NlepPt1cut);
  printf("Num of events passing lep2pt cut      : %d\n", NlepPt2cut);
  printf("Num of events passing HT cut          : %d\n", NHTcut);
  printf("Num of events passing Z cut           : %d\n", NZcut);
  printf("Num of events passing dileptonMass cut: %d\n", NMasscut);

  TCanvas c1;
  osHTHist->Draw();
  c1.Print("HT_oppositeSignEvents.pdf");
  //Note: sample does not have any detected same sign lepton. Instead, these are the predicted number of same sign events that have been misID as opposite sign events
  TCanvas c2; 
  ssHTHist->Draw();
  c2.Print("HT_sameSignEvents.pdf");
  //finally let's save our predictions to a root file for ease of use later
  T->Branch("osHTHist","TH1F",&osHTHist,32000,0);
  T->Branch("ssHTHist","TH1F",&ssHTHist,32000,0);
  T->Fill();
  T->Print();
  fbg->Write();
  fbg->Close();
}


