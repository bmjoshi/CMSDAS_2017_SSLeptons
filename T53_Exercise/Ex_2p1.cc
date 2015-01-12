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
  TFile* fDY = new TFile("/uscms_data/d3/clint/public/ljmet_tree_DY.root");
  TTree* tDY = (TTree*)fDY->Get("ljmet");
  //  TTreeReader myReader("ljmet", f);

  //initialize need histograms
  //TH1F* masshist = new TH1F("masshist", "Dielectron Invariant Mass",100,0.,200.);
  TH1F* numPtHist = new TH1F("numPtHist","Lepton p_{T} - TIGHT ID",200,0.,400.);
  TH1F* denPtHist = new TH1F("denPtHist","Lepton p_{T} - LOOSE ID",200,0.,400.);
  TH1F* numEtaHist = new TH1F("numEtaHist","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist = new TH1F("denEtaHist","Lepton #eta - LOOSE ID",30,-3.,3.);

  
  int nEntries = tDY->GetEntries();
  vector<double> *diElMass = 0;
  //kinematic variables
  vector<double> *elPts = 0;
  vector<double> *elEtas =0;
  vector<double>* elPhis =0;
  vector<int> * elCharge = 0;
  //variables for tracking cuts
  vector<double> *elDeta =0;
  vector<double> *elDphi =0;
  //variables for primary vtx cuts
  vector<double> *elDZs = 0;
  vector<double> *elD0s = 0;
  // H over E
  vector<double> *elHoverEs = 0;
  //missing hits
  vector<int> *elMHits = 0;
  //ooemoop
  vector<double> *elOoEmooPs = 0;
  //charged isolation
  vector<double> *elRelIsos = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetas = 0;
  vector<int>* elChargeConsistency = 0;

  //set branch addresses
  tDY->SetBranchAddress("diElMass_DileptonCalc", &diElMass);
  tDY->SetBranchAddress("elPt_DileptonCalc", &elPts);
  tDY->SetBranchAddress("elEta_DileptonCalc", &elEtas);
  tDY->SetBranchAddress("elPhi_DileptonCalc", &elPhis);
  tDY->SetBranchAddress("elDeta_DileptonCalc", &elDeta);
  tDY->SetBranchAddress("elDphi_DileptonCalc", &elDphi);
  tDY->SetBranchAddress("elDZ_DileptonCalc", &elDZs);
  tDY->SetBranchAddress("elD0_DileptonCalc", &elD0s);
  tDY->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs);
  tDY->SetBranchAddress("elMHits_DileptonCalc",&elMHits);
  tDY->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs);
  tDY->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos);
  tDY->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas);
  tDY->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistency);
  tDY->SetBranchAddress("elCharge_DileptonCalc",&elCharge);

  for(int ient = 0; ient < nEntries; ient++){
    tDY->GetEntry(ient);
    bool masscut = false;
    //checks per event
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;
    //Put electrons back together into coherent objects
    vector <Electron*> vEl;
    for (unsigned int uiEl = 0; uiEl < elPts->size(); uiEl++){
      Electron* el = new Electron;

      el->pt                = elPts->at(uiEl);
      el->eta               = elEtas->at(uiEl);
      el->phi               = elPhis->at(uiEl);
      el->charge            = elCharge->at(uiEl);
      el->dEta              = elDeta->at(uiEl);
      el->dPhi              = elDphi->at(uiEl);
      el->dZ                = elDZs->at(uiEl);
      el->d0                = elD0s->at(uiEl);
      el->hOverE            = elHoverEs->at(uiEl);
      el->mHits             = elMHits->at(uiEl);
      el->ooEmooP           = elOoEmooPs->at(uiEl);
      el->relIso            = elRelIsos->at(uiEl);
      el->sigmaIetaIeta     = elSigmaIetaIetas->at(uiEl);
      el->chargeConsistency = elChargeConsistency->at(uiEl);
      
      vEl.push_back(el);
    }

    //vector for lepton pair
    vector<Electron*> vElPair;
   
    //checks per lepton:
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply loose selection to the electron
      if (!vEl.at(ui)->loose()) continue;
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->loose()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  masscut = true;
	  vElPair.push_back(vEl.at(ui));
	  vElPair.push_back(vEl.at(uj));
	  
	}	
      }//End loop over second lepton
      if (masscut) break;
    }//End loop over first lepton

    //skip events not passing mass cut
    if(!masscut) continue;
    
    //define id info
    bool elTight1 = vElPair.at(0)->tight();
    bool elLoose1 = vElPair.at(0)->loose();
    bool elTight2 = vElPair.at(1)->tight();
    bool elLoose2 = vElPair.at(1)->loose();
    
    //check to make sure there is at least one tight lepton
    if( elTight1 || elTight2){
      //then fill denominator with other lepton if it's a loose one
      if(elTight1){
	if(elLoose2){
	  denPtHist->Fill(vElPair.at(1)->pt);
	  denEtaHist->Fill(vElPair.at(1)->eta);
	  
	  if(elTight2){
	    numPtHist->Fill(vElPair.at(1)->pt);
	    numEtaHist->Fill(vElPair.at(1)->eta);
	  }
	}
	
      }
      /*      else{
	if(elLoose1){
	  denPtHist->Fill(vElPair.at(0)->pt);
	  denEtaHist->Fill(vElPair.at(0)->eta);
	  if(elTight1){
	    numPtHist->Fill(vElPair.at(0)->pt);
	    numEtaHist->Fill(vElPair.at(0)->eta);
	  }
	}
	
	}*/
    }
    
    
  }//finish loop over entries
  


  //masshist->Draw();
  
  TGraphAsymmErrors* ptgraph = new TGraphAsymmErrors(numPtHist,denPtHist);
  TGraphAsymmErrors* etagraph = new TGraphAsymmErrors(numEtaHist,denEtaHist);
  
  TCanvas c1;
  ptgraph->Draw("apl");
  c1.Print("PromptRate_v_pT.pdf");

  TCanvas c2;
  etagraph->Draw("apl");
  c2.Print("PromptRate_v_Eta.pdf");



}
