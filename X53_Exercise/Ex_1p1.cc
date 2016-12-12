#include <iostream>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "ObjectID.C"
#include "TLorentzVector.h"
#include "TChain.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z

void Ex_1p1(){
  /*
    The first part of this script is just a lot of logistics: initializing variables, getting information from the ttree, etc. To complete the exercise you don't need to understand all of it,
    but please feel free to use it as a reference for how to get information from a ttree and loop over all the entries.
  */

  //load file and tree
  //TChain* t = new TChain("ljmet");
  TChain* t = new TChain("ChargeMisID");
  t->Add("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/ChargeMisID_Data_Run2016_Electrons_MVATightRC.root");
  

  int nEntries = t->GetEntries();
  //kinematic variables

  float elPts1, elEtas1, elPhis1, elEs1;
  float elPts2, elEtas2, elPhis2, elEs2;
  int elCharge1, elCharge2;
  int lep1etabin, lep2etabin;

  //Set branch addresses
  t->SetBranchAddress("Lep1Pt", &elPts1);
  t->SetBranchAddress("Lep1Eta", &elEtas1);
  t->SetBranchAddress("Lep1Phi", &elPhis1);
  t->SetBranchAddress("Lep1E", &elEs1);
  t->SetBranchAddress("Lep1Charge", &elCharge1);
  t->SetBranchAddress("Lep2Pt", &elPts2);
  t->SetBranchAddress("Lep2Eta", &elEtas2);
  t->SetBranchAddress("Lep2Phi", &elPhis2);
  t->SetBranchAddress("Lep2E", &elEs2);
  t->SetBranchAddress("Lep2Charge", &elCharge2);

  //Histograms
  TH1F* allmass = new TH1F("allmass","DiElectron Invariant Mass All Events",100,0.,200.);
  TH1F* osmass = new TH1F("osmass","DiElecton Invariant Mass for Opposite-Sign Events",100,0.,200.);
  TH1F* ssmass = new TH1F("ssmass","DiElecton Invariant Mass for Same-Sign Events",100,0.,200.);

  TH1F* ssEtaHist = new TH1F("ssEtaHist","#eta",15,-3,3);
  TH1F* totEtaHist = new TH1F("totEtaHist","#eta",15,-3,3);

  TH1F* ssPtHist = new TH1F("ssPtHist","p_{T}",100,0.,200.);
  TH1F* totPtHist = new TH1F("totPtHist","p_{T}",100,0.,200.);

  //Loop over the tree and look for an electron pair that makes a Z
  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;

    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(elPts1, elEtas1, elPhis1, M_EL);
    v2.SetPtEtaPhiM(elPts2, elEtas2, elPhis2, M_EL);
    double mass = (v1+v2).M();

    allmass->Fill(mass);

    totEtaHist->Fill(v1.Eta());
    totEtaHist->Fill(v2.Eta());
	  
    totPtHist->Fill(v1.Pt());
    totPtHist->Fill(v2.Pt());

    if (elCharge1 == elCharge2 ){
      ssmass->Fill(mass);	    
      ssEtaHist->Fill(v1.Eta());
      ssEtaHist->Fill(v2.Eta());
      ssPtHist->Fill(v1.Pt());
      ssPtHist->Fill(v2.Pt());
    }
    else{
      osmass->Fill(mass);
    }
  }

  //Draw the histograms
  TCanvas c1;
  allmass->Draw();

  c1.Print("DiElectronInvariantMass_all.pdf");

  /* ADD CODE TO PLOT MASS DISTRIBUTION FOR OPPOSITE SIGN AND SAME SIGN. 
     MOST OF IT IS THERE, YOU JUST NEED TO NORMALIZE OSMASS AND SSMASS BEFORE DRAWING THEM
  */ 

  TCanvas c2;
  osmass->Draw();

  c2.Print("DiElectronInvariantMass_os.pdf");

  TCanvas c3;
  ssmass->Draw();
  
  c3.Print("DiElectronInvariantMass_ss.pdf");  


  //ADD CODE TO PRINT OUT CHARGE MISID RATE


  std::cout<<"fake rate: "<< ssmass->Integral() /  allmass->Integral()<<std::endl;

  /* ADD CODE HERE TO MAKE PLOTS OF THE CHARGE MISID RATE VS PT, ETA. SOME OF IT IS THERE SO THAT WE ALL HAVE THE SAME OUTPUT FILES
     HINT1: MAKE USE OF TGRAPHASYMMERRORS, 
     HINT2: SET YOUR VS ETA BINS ACCORDING TO THE getEtaBin METHOD PROVIDED AT THE BOTTOM OF THIS FILE
   */
  //tcanvas for eta plot
  TCanvas c4;

  //add code!


  c4.Print("chargeMisID_vEta.pdf");

  //same as above but for pt
  TCanvas c5;



  c5.Print("chargeMisID_vPt.pdf");

}

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
