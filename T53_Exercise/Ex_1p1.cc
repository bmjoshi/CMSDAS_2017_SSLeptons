#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"

void Ex_1p1(){

  //load file and tree
  TFile* f = new TFile("ljmet_tree_DY.root");
  TTree* t = f->Get("ljmet");

  //Draw the mass histogram for all events
  TH1F* allmass = new TH1F("allmass","DiElectron Invariant Mass All Events",100,0.,200.);
  t->Draw("diElMass_DileptonCalc>>allmass");


  ///TODO: add parts about IDs, but first need to modify dilepton calc to do this easily

  //make strings for cuts
  string oscuts = "elCharge1_DileptonCalc != elCharge2_DileptonCalc";
  string sscuts = "elCharge1_DileptonCalc == elCharge2_DileptonCalc";

  //histogram for os events
  TH1F* osmass = new TH1F("osmass","DiElecton Invariant Mass for Opposite-Sign Events",100,0.,200.);
  t->Draw("diElMass_DileptonCalc>>osmass",oscuts.c_str());

  //histogram for ss events
  TH1F* ssmass = new TH1F("ssmass","DiElecton Invariant Mass for Same-Sign Events",100,0.,200.);
  t->Draw("diElMass_DileptonCalc>>ssmass",sscuts.c_str());

  TCanvas c1;
  allmass->Draw();

  c1.Print("DiElectronInvariantMass_all.pdf");

  TCanvas c2;
  osmass->Draw();

  c2.Print("DiElectronInvariantMass_os.pdf");

  TCanvas c3;
  ssmass->Draw();
  
  c3.Print("DiElectronInvariantMass_ss.pdf");
}
