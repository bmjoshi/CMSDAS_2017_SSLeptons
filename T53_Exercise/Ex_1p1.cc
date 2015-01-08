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

  //Now draw fake rate vs eta, pt

  //define cuts
  string masscuts = "(76 < diElMass_DileptonCalc) && (diElMass_DileptonCalc < 111)";
  string SSmasscuts = "(76 < diElMass_DileptonCalc) && (diElMass_DileptonCalc < 111) && (elCharge1_DileptonCalc == elCharge2_DileptonCalc)";

  //first get histograms
  TH1F* ssEtaHist = new TH1F("ssEtaHist","#eta",30,-3,3);
  TH1F* totEtaHist = new TH1F("totEtaHist","#eta",30,-3,3);

  t->Draw("elEta_DileptonCalc>>ssEtaHist",SSmasscuts.c_str());
  t->Draw("elEta_DileptonCalc>>totEtaHist",masscuts.c_str());

  //make new TGraph to get ratio
  TGraphAsymmErrors* etaGraph = new TGraphAsymmErrors(ssEtaHist,totEtaHist);

  TCanvas c5;
  etaGraph->Draw("apl");
  c5.Print("chargeMisID_vEta.pdf");

  //same as above but for pt
  //first get histograms
  TH1F* ssPtHist = new TH1F("ssPtHist","p_{T}",100,0.,200.);
  TH1F* totPtHist = new TH1F("totPtHist","p_{T}",100,0.,200.);

  t->Draw("elPt_DileptonCalc>>ssPtHist",SSmasscuts.c_str());
  t->Draw("elPt_DileptonCalc>>totPtHist",masscuts.c_str());

  //make new TGraph to get ratio
  TGraphAsymmErrors* ptGraph = new TGraphAsymmErrors(ssPtHist,totPtHist);

  TCanvas c5;
  ptGraph->Draw("apl");
  c5.Print("chargeMisID_vPt.pdf");

}
