#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "ObjectID.C"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z


void Ex_2p1(){

  /*
    The first part of this script is just a lot of logistics: initializing variables, getting information from the ttree, etc. To complete the exercise you don't need to understand all of it,
    but please feel free to use it as a reference for how to get information from a ttree and loop over all the entries.
   */

  //load the files
  TFile* fEl = new TFile("/uscms_data/d3/clint/public/PromptRate_Data_All_Electrons_MVATightRC_SortByPhi.root");
  TTree* tEl = (TTree*)fEl->Get("FakeRate");
  TFile* fMu = new TFile("/uscms_data/d3/clint/public/PromptRate_Data_All_Muons_CBTightMiniIso_SortByPhi.root");
  TTree* tMu = (TTree*)fMu->Get("FakeRate");
  //  TTreeReader myReader("ljmet", f);

  //initialize need histograms
  //TH1F* masshist = new TH1F("masshist", "Dielectron Invariant Mass",100,0.,200.);
  TH1F* numPtHist_ee = new TH1F("numPtHist_ee","Lepton p_{T} - TIGHT ID",200,0.,400.);
  TH1F* denPtHist_ee = new TH1F("denPtHist_ee","Lepton p_{T} - LOOSE ID",200,0.,400.);
  TH1F* numEtaHist_ee = new TH1F("numEtaHist_ee","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist_ee = new TH1F("denEtaHist_ee","Lepton #eta - LOOSE ID",30,-3.,3.);

  TH1F* numPtHist_mumu = new TH1F("numPtHist_mumu","Lepton p_{T} - TIGHT ID",200,0.,400.);
  TH1F* denPtHist_mumu = new TH1F("denPtHist_mumu","Lepton p_{T} - LOOSE ID",200,0.,400.);
  TH1F* numEtaHist_mumu = new TH1F("numEtaHist_mumu","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist_mumu = new TH1F("denEtaHist_mumu","Lepton #eta - LOOSE ID",30,-3.,3.);

  // * ADD CODE  HERE TO PLOT PROMPT RATE VS. MINISOLATION REQUIREMENT *


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
    denPtHist_mumu->Fill(muPt);
    denEtaHist_mumu->Fill(muEta);
    //check tight
    if(muIsTight){

      // *Currently miniIso requirement is doing nothing* - PLAY WITH THIS NUMBER (0-0.4) TO SEE IT'S EFFECTS ON THE PROMPT RATE
      if(muMiniIso<999){
	numEtaHist_mumu->Fill(muEta);
	numPtHist_mumu->Fill(muPt);
       }
    }
  }//end loop on muons

  for(int iel=0; iel<tEl->GetEntries();iel++){
    tEl->GetEntry(iel);
    if(iel % 100000 ==0) std::cout<<"Completed "<<iel<<" out of "<<tEl->GetEntries()<<" electron events"<<std::endl;

    //fill denominator histogram
    denPtHist_ee->Fill(elPt);
    denEtaHist_ee->Fill(elEta);
    //check tight
    if(elIsTight){

      // *Currently miniIso requirement is doing nothing* - PLAY WITH THIS NUMBER (0-0.4) TO SEE IT'S EFFECTS ON THE PROMPT RATE
      if(elMiniIso<999){
	numEtaHist_ee->Fill(elEta);
	numPtHist_ee->Fill(elPt);
       }
    }
  }//end loop on electrons


  TGraphAsymmErrors* ptgraph_ee = new TGraphAsymmErrors(numPtHist_ee,denPtHist_ee);
  TGraphAsymmErrors* etagraph_ee = new TGraphAsymmErrors(numEtaHist_ee,denEtaHist_ee);  
  
  TCanvas c1;
  ptgraph_ee->GetYaxis()->SetRangeUser(0,1);
  ptgraph_ee->Draw("apl");
  // * FEEL FREE TO MAKE THE PLOTS PRETTIER (BETTER?) BY ADDING AXIS LABELS :) HINT: YOU CAN USE TGraphAsymmErrors::SetTitle() TO SET ALL THREE (GRAPH TITLE, X-AXIS LABEL, Y-AXIS LABEL) *
  c1.Print("PromptRate_v_pT_ee.pdf");

  TCanvas c2;
  etagraph_ee->GetYaxis()->SetRangeUser(0,1);
  etagraph_ee->Draw("apl");
  c2.Print("PromptRate_v_Eta_ee.pdf");

  TGraphAsymmErrors* ptgraph_mumu = new TGraphAsymmErrors(numPtHist_mumu,denPtHist_mumu);
  TGraphAsymmErrors* etagraph_mumu = new TGraphAsymmErrors(numEtaHist_mumu,denEtaHist_mumu);
  
  TCanvas c3;
  ptgraph_mumu->GetYaxis()->SetRangeUser(0,1);
  ptgraph_mumu->Draw("apl");
  c3.Print("PromptRate_v_pT_mumu.pdf");

  TCanvas c4;
  etagraph_mumu->GetYaxis()->SetRangeUser(0,1);
  etagraph_mumu->Draw("apl");
  c4.Print("PromptRate_v_Eta_mumu.pdf");

  //ADD YOUR CODE TO GIVE THE PROMPT RATE FOR ELECTRONS AND THE PROMPT RATE FOR MUONS



}
