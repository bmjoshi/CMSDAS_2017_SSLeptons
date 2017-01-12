#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TStyle.h"
#include "ObjectID.C"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z



void Ex_2p2p5(){

  gStyle->SetOptStat(kFALSE);

  TFile* f = new TFile("/uscms_data/d3/clint/public/FakeRate_TTbar.root");
  TH1F* el_all = (TH1F*) f->Get("elNumHist_all");
  TH1F* el_lpt = (TH1F*) f->Get("elNumHist_lpt");
  TH1F* el_hpt = (TH1F*) f->Get("elNumHist_hpt");

  TH1F* mu_all = (TH1F*) f->Get("muNumHist_all");
  TH1F* mu_lpt = (TH1F*) f->Get("muNumHist_lpt");
  TH1F* mu_hpt = (TH1F*) f->Get("muNumHist_hpt");

  TH1F* elden_all = (TH1F*) f->Get("elDenHist_all");
  TH1F* elden_lpt = (TH1F*) f->Get("elDenHist_lpt");
  TH1F* elden_hpt = (TH1F*) f->Get("elDenHist_hpt");

  TH1F* muden_all = (TH1F*) f->Get("muDenHist_all");
  TH1F* muden_lpt = (TH1F*) f->Get("muDenHist_lpt");
  TH1F* muden_hpt = (TH1F*) f->Get("muDenHist_hpt");


  el_all->Divide(elden_all);
  el_lpt->Divide(elden_lpt);
  el_hpt->Divide(elden_hpt);

  mu_all->Divide(muden_all);
  mu_lpt->Divide(muden_lpt);
  mu_hpt->Divide(muden_hpt);

  TCanvas* c = new TCanvas();
  el_all->Draw("pe");
  el_all->GetYaxis()->SetRangeUser(0,1.1);
  c->Print("Electron_FakeRate_TTbar_All.pdf");

  el_lpt->Draw("pe");
  el_lpt->GetYaxis()->SetRangeUser(0,1.1);
  c->Print("Electron_FakeRate_TTbar_pt25-35.pdf");

  el_hpt->Draw("pe");
  el_hpt->GetYaxis()->SetRangeUser(0,1.1);
  c->Print("Electron_FakeRate_TTbar_pt35-Inf.pdf");

  mu_all->Draw("pe");
  mu_all->GetYaxis()->SetRangeUser(0,1.1);
  c->Print("Muon_FakeRate_TTbar_All.pdf");

  mu_lpt->Draw("pe");
  mu_lpt->GetYaxis()->SetRangeUser(0,1.1);
  c->Print("Muon_FakeRate_TTbar_pt25-35.pdf");

  mu_hpt->Draw("pe");
  mu_hpt->GetYaxis()->SetRangeUser(0,1.1);
  c->Print("Muon_FakeRate_TTbar_pt35-Inf.pdf");


}
