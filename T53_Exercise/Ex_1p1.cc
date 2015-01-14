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
  TChain* t = new TChain("ljmet");
  t->Add("/uscms_data/d3/clint/public/ljmet_tree_DY.root");
  //t->Add("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  

  int nEntries = t->GetEntries();
  //kinematic variables
  vector<double> *elPts = 0;
  vector<double> *elEtas =0;
  vector<double> *elPhis =0;

  vector<int>* elCharge = 0;

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

  //Set branch addresses
  t->SetBranchAddress("elPt_DileptonCalc", &elPts);
  t->SetBranchAddress("elEta_DileptonCalc", &elEtas);
  t->SetBranchAddress("elPhi_DileptonCalc", &elPhis);
  t->SetBranchAddress("elCharge_DileptonCalc", &elCharge);
  t->SetBranchAddress("elDeta_DileptonCalc", &elDeta);
  t->SetBranchAddress("elDphi_DileptonCalc", &elDphi);
  t->SetBranchAddress("elDZ_DileptonCalc", &elDZs);
  t->SetBranchAddress("elD0_DileptonCalc", &elD0s);
  t->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs);
  t->SetBranchAddress("elMHits_DileptonCalc",&elMHits);
  t->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs);
  t->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos);
  t->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas);
  t->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistency);

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

    bool foundPair = false;
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      // ADD CODE HERE - THE LINE BELOW DETERMINES THE QUALITY OF THE ELECTRON
      if (!vEl.at(ui)->tight()) continue;
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	  allmass->Fill(mass);

	  totEtaHist->Fill(v1.Eta());
	  totEtaHist->Fill(v2.Eta());

	  totPtHist->Fill(v1.Pt());
	  totPtHist->Fill(v2.Pt());
	  
	  if (vEl.at(ui)->charge == vEl.at(uj)->charge){
	    ssmass->Fill(mass);
	    
	    ssEtaHist->Fill(v1.Eta());
	    ssEtaHist->Fill(v2.Eta());
	    
	    ssPtHist->Fill(v1.Pt());
	    ssPtHist->Fill(v2.Pt());
	  }
	  else osmass->Fill(mass);
	}	
      }//End loop over second lepton
      if (foundPair) break;
    }//End loop over first lepton

  }//End loop over tree  

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
     HINT: MAKE USE OF TGRAPHASYMMERRORS
   */
  //tcanvas for eta plot
  TCanvas c4;

  //add code!


  c4.Print("chargeMisID_vEta.pdf");

  //same as above but for pt
  TCanvas c5;



  c5.Print("chargeMisID_vPt.pdf");

}
