#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "ObjectID.C"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z


void Ex_2p2(){

  /*
    The first part of this script is just a lot of logistics: initializing variables, getting information from the ttree, etc. To complete the exercise you don't need to understand all of it,
    but please feel free to use it as a reference for how to get information from a ttree and loop over all the entries.
   */

  //load the files
  TFile* fEl = new TFile("/uscms_data/d3/clint/public/FakeRate_Data_Electrons_MVATightRC_modifiedLepPt.root");
  TTree* tEl = (TTree*)fEl->Get("FakeRate");
  TFile* fMu = new TFile("/uscms_data/d3/clint/public/FakeRate_Data_Muons_CBTightMiniIso_modifiedLepPt.root");
  TTree* tMu = (TTree*)fMu->Get("FakeRate");

  //initialize need histograms
  TH1F* numPtHist_e = new TH1F("numPtHist_e","Lepton p_{T} - TIGHT ID",80,0.,40.);
  TH1F* denPtHist_e = new TH1F("denPtHist_e","Lepton p_{T} - LOOSE ID",80,0.,40.);
  TH1F* numEtaHist_e = new TH1F("numEtaHist_e","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist_e = new TH1F("denEtaHist_e","Lepton #eta - LOOSE ID",30,-3.,3.);

  TH1F* numPtHist_m = new TH1F("numPtHist_m","Lepton p_{T} - TIGHT ID",80,0.,40.);
  TH1F* denPtHist_m = new TH1F("denPtHist_m","Lepton p_{T} - LOOSE ID",80,0.,40.);
  TH1F* numEtaHist_m = new TH1F("numEtaHist_m","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist_m = new TH1F("denEtaHist_m","Lepton #eta - LOOSE ID",30,-3.,3.);

  
  // * ADD CODE HERE TO GET TOTAL PROMPT RATE *
  
  int nEntriesMu = tMu->GetEntries();
  //set branch addresses
  float muPt,muEta,muPhi,muEnergy,muMiniIso;
  int muIsTight;
  tMu->SetBranchAddress("LepPt",&muPt);
  tMu->SetBranchAddress("LepEta",&muEta);
  tMu->SetBranchAddress("LepPhi",&muPhi);
  tMu->SetBranchAddress("LepE",&muEnergy);
  tMu->SetBranchAddress("LepPt",&muPt);
  tMu->SetBranchAddress("LepIsTight",&muIsTight);
  tMu->SetBranchAddress("LepPt",&muPt);
  tMu->SetBranchAddress("LepMiniIso",&muMiniIso);

  float elPt,elEta,elPhi,elEnergy,elMiniIso;
  int elIsTight;
  tEl->SetBranchAddress("LepPt",&elPt);
  tEl->SetBranchAddress("LepEta",&elEta);
  tEl->SetBranchAddress("LepPhi",&elPhi);
  tEl->SetBranchAddress("LepE",&elEnergy);
  tEl->SetBranchAddress("LepPt",&elPt);
  tEl->SetBranchAddress("LepIsTight",&elIsTight);
  tEl->SetBranchAddress("LepPt",&elPt);
  tEl->SetBranchAddress("LepMiniIso",&elMiniIso);

  for(int imu=0; imu<tMu->GetEntries();imu++){
    tMu->GetEntry(imu);
    if(imu % 100000 ==0) std::cout<<"Completed "<<imu<<" out of "<<tMu->GetEntries()<<" muon events"<<std::endl;

    //fill denominator histogram
    denPtHist_m->Fill(muPt);
    denEtaHist_m->Fill(muEta);
    //check tight
    if(muIsTight){

      // *Currently miniIso requirement is doing nothing* - PLAY WITH THIS NUMBER (0-0.4) TO SEE IT'S EFFECTS ON THE PROMPT RATE
      if(muMiniIso<0.1){
	numEtaHist_m->Fill(muEta);
	numPtHist_m->Fill(muPt);
       }
    }
  }//end loop on muons

  for(int iel=0; iel<tEl->GetEntries();iel++){
    tEl->GetEntry(iel);
    if(iel % 100000 ==0) std::cout<<"Completed "<<iel<<" out of "<<tEl->GetEntries()<<" electron events"<<std::endl;

    //fill denominator histogram
    denPtHist_e->Fill(elPt);
    denEtaHist_e->Fill(elEta);
    //check tight
    if(elIsTight){

      // *Currently miniIso requirement is doing nothing* - PLAY WITH THIS NUMBER (0-0.4) TO SEE IT'S EFFECTS ON THE PROMPT RATE
      if(elMiniIso<0.1){
	numEtaHist_e->Fill(elEta);
	numPtHist_e->Fill(elPt);
       }
    }
  }//end loop on electrons


  
  //masshist->Draw();  
  TGraphAsymmErrors* ptgraph_e = new TGraphAsymmErrors(numPtHist_e,denPtHist_e);
  TGraphAsymmErrors* etagraph_e = new TGraphAsymmErrors(numEtaHist_e,denEtaHist_e);
  
  TCanvas c1;
  ptgraph_e->Draw("apl");
  c1.Print("FakeRate_v_pT_e.pdf");
  
  TCanvas c2;
  etagraph_e->Draw("apl");
  c2.Print("FakeRate_v_Eta_e.pdf");
  
  TGraphAsymmErrors* ptgraph_m = new TGraphAsymmErrors(numPtHist_m,denPtHist_m);
  TGraphAsymmErrors* etagraph_m = new TGraphAsymmErrors(numEtaHist_m,denEtaHist_m);
  
  TCanvas c3;
  ptgraph_m->Draw("apl");
  c3.Print("FakeRate_v_pT_m.pdf");
  
  TCanvas c4;
  etagraph_m->Draw("apl");
  c4.Print("FakeRate_v_Eta_m.pdf");

  //* ADD CODE TO PRINT OUT OVERALL FAKE RATES *

}
