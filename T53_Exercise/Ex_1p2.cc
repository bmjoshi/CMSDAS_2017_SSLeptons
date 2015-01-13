#include <iostream>
#include "TH1.h"
#include "TTree.h"
#include "TFile.h"
#include "ObjectID.C"
#include "TLorentzVector.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z

void Ex_1p2(){
  /*
    THIS SCRIPT WILL LOOP OVER ALL EVENTS IN THE DIFFERENT SAMPLES AND PRINT OUT HOW MANY
    PASS THE MASS WINDOW CUT, AND ALSO HOW MANY PASS THE MASS WINDOW CUT AND ARE SAME SIGN
   */


  //load TFiles
  TFile* fDY = new TFile("/uscms_data/d3/clint/public/ljmet_tree_DY.root");
  TFile* fWZ = new TFile("/uscms_data/d3/clint/public/ljmet_tree_WZ.root");
  TFile* fWJets = new TFile("/uscms_data/d3/clint/public/ljmet_tree_WJets.root");
  TFile* fTT = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  TFile* fTTZ = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TTZ.root");


  //load the TTrees
  TTree* tDY = (TTree*)fDY->Get("ljmet");
  TTree* tWZ = (TTree*)fWZ->Get("ljmet");
  TTree* tWJets = (TTree*)fWJets->Get("ljmet");
  TTree* tTT = (TTree*)fTT->Get("ljmet");
  TTree* tTTZ = (TTree*)fTTZ->Get("ljmet");

  //variables for DY
  //kinematic variables
  vector<double> *elPtsDY = 0;
  vector<double> *elEtasDY =0;
  vector<double> *elPhisDY =0;

  vector<int>* elChargeDY = 0;

  //variables for tracking cuts
  vector<double> *elDetaDY =0;
  vector<double> *elDphiDY =0;
  //variables for primary vtx cuts
  vector<double> *elDZsDY = 0;
  vector<double> *elD0sDY = 0;
  // H over E
  vector<double> *elHoverEsDY = 0;
  //missing hits
  vector<int> *elMHitsDY = 0;
  //ooemoop
  vector<double> *elOoEmooPsDY = 0;
  //charged isolation
  vector<double> *elRelIsosDY = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetasDY = 0;
  vector<int>* elChargeConsistencyDY = 0;

  //Set branch addresses
  tDY->SetBranchAddress("elPt_DileptonCalc", &elPtsDY);
  tDY->SetBranchAddress("elEta_DileptonCalc", &elEtasDY);
  tDY->SetBranchAddress("elPhi_DileptonCalc", &elPhisDY);
  tDY->SetBranchAddress("elCharge_DileptonCalc", &elChargeDY);
  tDY->SetBranchAddress("elDeta_DileptonCalc", &elDetaDY);
  tDY->SetBranchAddress("elDphi_DileptonCalc", &elDphiDY);
  tDY->SetBranchAddress("elDZ_DileptonCalc", &elDZsDY);
  tDY->SetBranchAddress("elD0_DileptonCalc", &elD0sDY);
  tDY->SetBranchAddress("elHoE_DileptonCalc",&elHoverEsDY);
  tDY->SetBranchAddress("elMHits_DileptonCalc",&elMHitsDY);
  tDY->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPsDY);
  tDY->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsosDY);
  tDY->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetasDY);
  tDY->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistencyDY);

  //variables for TT
  //kinematic variables
  vector<double> *elPtsTT = 0;
  vector<double> *elEtasTT =0;
  vector<double> *elPhisTT =0;

  vector<int>* elChargeTT = 0;

  //variables for tracking cuts
  vector<double> *elDetaTT =0;
  vector<double> *elDphiTT =0;
  //variables for primary vtx cuts
  vector<double> *elDZsTT = 0;
  vector<double> *elD0sTT = 0;
  // H over E
  vector<double> *elHoverEsTT = 0;
  //missing hits
  vector<int> *elMHitsTT = 0;
  //ooemoop
  vector<double> *elOoEmooPsTT = 0;
  //charged isolation
  vector<double> *elRelIsosTT = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetasTT = 0;
  vector<int>* elChargeConsistencyTT = 0;

  //Set branch addresses
  tTT->SetBranchAddress("elPt_DileptonCalc", &elPtsTT);
  tTT->SetBranchAddress("elEta_DileptonCalc", &elEtasTT);
  tTT->SetBranchAddress("elPhi_DileptonCalc", &elPhisTT);
  tTT->SetBranchAddress("elCharge_DileptonCalc", &elChargeTT);
  tTT->SetBranchAddress("elDeta_DileptonCalc", &elDetaTT);
  tTT->SetBranchAddress("elDphi_DileptonCalc", &elDphiTT);
  tTT->SetBranchAddress("elDZ_DileptonCalc", &elDZsTT);
  tTT->SetBranchAddress("elD0_DileptonCalc", &elD0sTT);
  tTT->SetBranchAddress("elHoE_DileptonCalc",&elHoverEsTT);
  tTT->SetBranchAddress("elMHits_DileptonCalc",&elMHitsTT);
  tTT->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPsTT);
  tTT->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsosTT);
  tTT->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetasTT);
  tTT->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistencyTT);

  //variables for TTZ
  //kinematic variables
  vector<double> *elPtsTTZ = 0;
  vector<double> *elEtasTTZ =0;
  vector<double> *elPhisTTZ =0;

  vector<int>* elChargeTTZ = 0;

  //variables for tracking cuts
  vector<double> *elDetaTTZ =0;
  vector<double> *elDphiTTZ =0;
  //variables for primary vtx cuts
  vector<double> *elDZsTTZ = 0;
  vector<double> *elD0sTTZ = 0;
  // H over E
  vector<double> *elHoverEsTTZ = 0;
  //missing hits
  vector<int> *elMHitsTTZ = 0;
  //ooemoop
  vector<double> *elOoEmooPsTTZ = 0;
  //charged isolation
  vector<double> *elRelIsosTTZ = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetasTTZ = 0;
  vector<int>* elChargeConsistencyTTZ = 0;

  //Set branch addresses
  tTTZ->SetBranchAddress("elPt_DileptonCalc", &elPtsTTZ);
  tTTZ->SetBranchAddress("elEta_DileptonCalc", &elEtasTTZ);
  tTTZ->SetBranchAddress("elPhi_DileptonCalc", &elPhisTTZ);
  tTTZ->SetBranchAddress("elCharge_DileptonCalc", &elChargeTTZ);
  tTTZ->SetBranchAddress("elDeta_DileptonCalc", &elDetaTTZ);
  tTTZ->SetBranchAddress("elDphi_DileptonCalc", &elDphiTTZ);
  tTTZ->SetBranchAddress("elDZ_DileptonCalc", &elDZsTTZ);
  tTTZ->SetBranchAddress("elD0_DileptonCalc", &elD0sTTZ);
  tTTZ->SetBranchAddress("elHoE_DileptonCalc",&elHoverEsTTZ);
  tTTZ->SetBranchAddress("elMHits_DileptonCalc",&elMHitsTTZ);
  tTTZ->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPsTTZ);
  tTTZ->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsosTTZ);
  tTTZ->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetasTTZ);
  tTTZ->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistencyTTZ);

  //variables for WZ
  //kinematic variables
  vector<double> *elPtsWZ = 0;
  vector<double> *elEtasWZ =0;
  vector<double> *elPhisWZ =0;

  vector<int>* elChargeWZ = 0;

  //variables for tracking cuts
  vector<double> *elDetaWZ =0;
  vector<double> *elDphiWZ =0;
  //variables for primary vtx cuts
  vector<double> *elDZsWZ = 0;
  vector<double> *elD0sWZ = 0;
  // H over E
  vector<double> *elHoverEsWZ = 0;
  //missing hits
  vector<int> *elMHitsWZ = 0;
  //ooemoop
  vector<double> *elOoEmooPsWZ = 0;
  //charged isolation
  vector<double> *elRelIsosWZ = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetasWZ = 0;
  vector<int>* elChargeConsistencyWZ = 0;

  //Set branch addresses
  tWZ->SetBranchAddress("elPt_DileptonCalc", &elPtsWZ);
  tWZ->SetBranchAddress("elEta_DileptonCalc", &elEtasWZ);
  tWZ->SetBranchAddress("elPhi_DileptonCalc", &elPhisWZ);
  tWZ->SetBranchAddress("elCharge_DileptonCalc", &elChargeWZ);
  tWZ->SetBranchAddress("elDeta_DileptonCalc", &elDetaWZ);
  tWZ->SetBranchAddress("elDphi_DileptonCalc", &elDphiWZ);
  tWZ->SetBranchAddress("elDZ_DileptonCalc", &elDZsWZ);
  tWZ->SetBranchAddress("elD0_DileptonCalc", &elD0sWZ);
  tWZ->SetBranchAddress("elHoE_DileptonCalc",&elHoverEsWZ);
  tWZ->SetBranchAddress("elMHits_DileptonCalc",&elMHitsWZ);
  tWZ->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPsWZ);
  tWZ->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsosWZ);
  tWZ->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetasWZ);
  tWZ->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistencyWZ);

    //variables for WJets
  //kinematic variables
  vector<double> *elPtsWJets = 0;
  vector<double> *elEtasWJets =0;
  vector<double> *elPhisWJets =0;

  vector<int>* elChargeWJets = 0;

  //variables for tracking cuts
  vector<double> *elDetaWJets =0;
  vector<double> *elDphiWJets =0;
  //variables for primary vtx cuts
  vector<double> *elDZsWJets = 0;
  vector<double> *elD0sWJets = 0;
  // H over E
  vector<double> *elHoverEsWJets = 0;
  //missing hits
  vector<int> *elMHitsWJets = 0;
  //ooemoop
  vector<double> *elOoEmooPsWJets = 0;
  //charged isolation
  vector<double> *elRelIsosWJets = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetasWJets = 0;
  vector<int>* elChargeConsistencyWJets = 0;

  //Set branch addresses
  tWJets->SetBranchAddress("elPt_DileptonCalc", &elPtsWJets);
  tWJets->SetBranchAddress("elEta_DileptonCalc", &elEtasWJets);
  tWJets->SetBranchAddress("elPhi_DileptonCalc", &elPhisWJets);
  tWJets->SetBranchAddress("elCharge_DileptonCalc", &elChargeWJets);
  tWJets->SetBranchAddress("elDeta_DileptonCalc", &elDetaWJets);
  tWJets->SetBranchAddress("elDphi_DileptonCalc", &elDphiWJets);
  tWJets->SetBranchAddress("elDZ_DileptonCalc", &elDZsWJets);
  tWJets->SetBranchAddress("elD0_DileptonCalc", &elD0sWJets);
  tWJets->SetBranchAddress("elHoE_DileptonCalc",&elHoverEsWJets);
  tWJets->SetBranchAddress("elMHits_DileptonCalc",&elMHitsWJets);
  tWJets->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPsWJets);
  tWJets->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsosWJets);
  tWJets->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetasWJets);
  tWJets->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistencyWJets);

  //define our taget luminosity - we are picking 5 fb^{-1} to correspond to roughly half-way through the 2015 run
  float targetlumi = 5.0;

  //To simplify things a little I've placed the number of events ran over for each sample below:
  float nRunDY = 1366703;
  float nRunWZ = 237484;
  float nRunWJets = 3828404;
  float nRunTT = 2206600;
  float nRunTTZ = 249275;
  
  /* Now, to figure out how to scale things correctly we need to remember the following equation:

     N_events = lumi * cross section

     Hence, if we want to normalize to 5 fb^{-1} for each sample we want

     N_events / cross section = 5.0

     NOTE: This equation is not currently true because we are forcing the events to be normalized to 5 inverse femtobarns.
     To make this equation true we will introduce a weight W_WZ for instance, which will make our equation true. What we really want then is to weight N_events as below:

     N_normalized = Weight_Sample * N_RunSample

     Where N_normalized is the normalized number of events we would see for that sample if we had collected 5^{-1] fb of data.

  */

  //To start, look up the cross sections for the various samples in MCFM and fill them out below:

  float xsecDY = 1.0;
  float xsecWZ = 1.0;
  float xsecWJets = 1.0;
  float xsecTT = 1.0;
  float xsecTTZ = 1.0;

  //Here some math is done for you :)

  float weightDY = (targetlumi*xsecDY) / (nRunDY);
  float weightWZ = (targetlumi*xsecWZ) / (nRunWZ);
  float weightWJets = (targetlumi*xsecWJets) / (nRunWJets);
  float weightTT = (targetlumi*xsecTT) / (nRunTT);
  float weightTTZ = (targetlumi*xsecTTZ) / (nRunTTZ);

  //Now that we have the weights we can find out how many events passed out selection:

  int Nos_DY =0;
  int Nos_WZ =0;
  int Nos_WJets =0;
  int Nos_TT =0;
  int Nos_TTZ=0;

  int Nss_DY =0;
  int Nss_WZ =0;
  int Nss_WJets =0;
  int Nss_TT =0;
  int Nss_TTZ=0;


  //loop over DY
  int nEntDY = tDY->GetEntries();
  for(int ient = 0; ient <nEntDY; ient++){
    tDY->GetEntry(ient);
    vector <Electron*> vEl;
    for (unsigned int uiEl = 0; uiEl < elPtsDY->size(); uiEl++){
      Electron* el = new Electron;
      
      el->pt                = elPtsDY->at(uiEl);
      el->eta               = elEtasDY->at(uiEl);
      el->phi               = elPhisDY->at(uiEl);
      el->charge            = elChargeDY->at(uiEl);
      el->dEta              = elDetaDY->at(uiEl);
      el->dPhi              = elDphiDY->at(uiEl);
      el->dZ                = elDZsDY->at(uiEl);
      el->d0                = elD0sDY->at(uiEl);
      el->hOverE            = elHoverEsDY->at(uiEl);
      el->mHits             = elMHitsDY->at(uiEl);
      el->ooEmooP           = elOoEmooPsDY->at(uiEl);
      el->relIso            = elRelIsosDY->at(uiEl);
      el->sigmaIetaIeta     = elSigmaIetaIetasDY->at(uiEl);
      el->chargeConsistency = elChargeConsistencyDY->at(uiEl);

      vEl.push_back(el);
    }

    bool foundPair = false;
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply tight selection to the electron
      if (!vEl.at(ui)->tight()) continue;
      
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	  Nos_DY+=1;

	  if (vEl.at(ui)->charge == vEl.at(uj)->charge){
	    Nss_DY+=1;
	  }

	}	
      }//End loop over second lepton
      if (foundPair) break;
    }//End loop over first lepton

  }//End loop over tree      

  //loop over TT
  int nEntTT = tTT->GetEntries();
  for( int ient = 0; ient <nEntTT; ient++){
    tTT->GetEntry(ient);
    vector <Electron*> vEl;
    for (unsigned int uiEl = 0; uiEl < elPtsTT->size(); uiEl++){
      Electron* el = new Electron;
      
      el->pt                = elPtsTT->at(uiEl);
      el->eta               = elEtasTT->at(uiEl);
      el->phi               = elPhisTT->at(uiEl);
      el->charge            = elChargeTT->at(uiEl);
      el->dEta              = elDetaTT->at(uiEl);
      el->dPhi              = elDphiTT->at(uiEl);
      el->dZ                = elDZsTT->at(uiEl);
      el->d0                = elD0sTT->at(uiEl);
      el->hOverE            = elHoverEsTT->at(uiEl);
      el->mHits             = elMHitsTT->at(uiEl);
      el->ooEmooP           = elOoEmooPsTT->at(uiEl);
      el->relIso            = elRelIsosTT->at(uiEl);
      el->sigmaIetaIeta     = elSigmaIetaIetasTT->at(uiEl);
      el->chargeConsistency = elChargeConsistencyTT->at(uiEl);

      vEl.push_back(el);
    }

    bool foundPair = false;
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply tight selection to the electron
      if (!vEl.at(ui)->tight()) continue;
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	  Nos_TT+=1;
	  if (vEl.at(ui)->charge == vEl.at(uj)->charge){
	    Nss_TT+=1;
	  }

	}	
      }//End loop over second lepton
      if (foundPair) break;
    }//End loop over first lepton

  }//End loop over tree      

  //loop over TTZ
  int nEntTTZ = tTTZ->GetEntries();
  for(int ient = 0; ient <nEntTTZ; ient++){
    tTTZ->GetEntry(ient);
    vector <Electron*> vEl;
    for (unsigned int uiEl = 0; uiEl < elPtsTTZ->size(); uiEl++){
      Electron* el = new Electron;
      
      el->pt                = elPtsTTZ->at(uiEl);
      el->eta               = elEtasTTZ->at(uiEl);
      el->phi               = elPhisTTZ->at(uiEl);
      el->charge            = elChargeTTZ->at(uiEl);
      el->dEta              = elDetaTTZ->at(uiEl);
      el->dPhi              = elDphiTTZ->at(uiEl);
      el->dZ                = elDZsTTZ->at(uiEl);
      el->d0                = elD0sTTZ->at(uiEl);
      el->hOverE            = elHoverEsTTZ->at(uiEl);
      el->mHits             = elMHitsTTZ->at(uiEl);
      el->ooEmooP           = elOoEmooPsTTZ->at(uiEl);
      el->relIso            = elRelIsosTTZ->at(uiEl);
      el->sigmaIetaIeta     = elSigmaIetaIetasTTZ->at(uiEl);
      el->chargeConsistency = elChargeConsistencyTTZ->at(uiEl);

      vEl.push_back(el);
    }

    bool foundPair = false;
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply tight selection to the electron
      if (!vEl.at(ui)->tight()) continue;
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	  Nos_TTZ+=1;
	  if (vEl.at(ui)->charge == vEl.at(uj)->charge){
	    Nss_TTZ+=1;
	  }

	}	
      }//End loop over second lepton
      if (foundPair) break;
    }//End loop over first lepton

  }//End loop over tree      
  
  //loop over WZ
  int nEntWZ = tWZ->GetEntries();
  for(int ient = 0; ient <nEntWZ; ient++){
    tWZ->GetEntry(ient);
    vector <Electron*> vEl;
    for (unsigned int uiEl = 0; uiEl < elPtsWZ->size(); uiEl++){
      Electron* el = new Electron;
      
      el->pt                = elPtsWZ->at(uiEl);
      el->eta               = elEtasWZ->at(uiEl);
      el->phi               = elPhisWZ->at(uiEl);
      el->charge            = elChargeWZ->at(uiEl);
      el->dEta              = elDetaWZ->at(uiEl);
      el->dPhi              = elDphiWZ->at(uiEl);
      el->dZ                = elDZsWZ->at(uiEl);
      el->d0                = elD0sWZ->at(uiEl);
      el->hOverE            = elHoverEsWZ->at(uiEl);
      el->mHits             = elMHitsWZ->at(uiEl);
      el->ooEmooP           = elOoEmooPsWZ->at(uiEl);
      el->relIso            = elRelIsosWZ->at(uiEl);
      el->sigmaIetaIeta     = elSigmaIetaIetasWZ->at(uiEl);
      el->chargeConsistency = elChargeConsistencyWZ->at(uiEl);

      vEl.push_back(el);
    }

    bool foundPair = false;
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply tight selection to the electron
      if (!vEl.at(ui)->tight()) continue;
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	  Nos_WZ+=1;
	  if (vEl.at(ui)->charge == vEl.at(uj)->charge){
	    Nss_WZ+=1;
	  }

	}	
      }//End loop over second lepton
      if (foundPair) break;
    }//End loop over first lepton

  }//End loop over tree      

  //loop over WJets
  int nEntWJets = tWJets->GetEntries();
  for( int ient = 0; ient <nEntWJets; ient++){
    tWJets->GetEntry(ient);
    vector <Electron*> vEl;
    for (unsigned int uiEl = 0; uiEl < elPtsWJets->size(); uiEl++){
      Electron* el = new Electron;
      
      el->pt                = elPtsWJets->at(uiEl);
      el->eta               = elEtasWJets->at(uiEl);
      el->phi               = elPhisWJets->at(uiEl);
      el->charge            = elChargeWJets->at(uiEl);
      el->dEta              = elDetaWJets->at(uiEl);
      el->dPhi              = elDphiWJets->at(uiEl);
      el->dZ                = elDZsWJets->at(uiEl);
      el->d0                = elD0sWJets->at(uiEl);
      el->hOverE            = elHoverEsWJets->at(uiEl);
      el->mHits             = elMHitsWJets->at(uiEl);
      el->ooEmooP           = elOoEmooPsWJets->at(uiEl);
      el->relIso            = elRelIsosWJets->at(uiEl);
      el->sigmaIetaIeta     = elSigmaIetaIetasWJets->at(uiEl);
      el->chargeConsistency = elChargeConsistencyWJets->at(uiEl);

      vEl.push_back(el);
    }

    bool foundPair = false;
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply tight selection to the electron
      if (!vEl.at(ui)->tight()) continue;
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	  Nos_WJets+=1;
	  if (vEl.at(ui)->charge == vEl.at(uj)->charge){
	    Nss_WJets+=1;
	  }

	}	
      }//End loop over second lepton
      if (foundPair) break;
    }//End loop over first lepton

  }//End loop over tree      



  //now weight them and print out the values
  float nNormDY = Nos_DY * weightDY;
  float nNormWZ = Nos_WZ * weightWZ;
  float nNormWJets = Nos_WJets * weightWJets;
  float nNormTT = Nos_TT * weightTT;
  float nNormTTZ = Nos_TTZ * weightTTZ;

  //now weight them and print out the values
  float nSSNormDY = Nss_DY * weightDY;
  float nSSNormWZ = Nss_WZ * weightWZ;
  float nSSNormWJets = Nss_WJets * weightWJets;
  float nSSNormTT = Nss_TT * weightTT;
  float nSSNormTTZ = Nss_TTZ * weightTTZ;
  
  std::cout<<"Number of events passing mass window cut from DY: "<<nNormDY<<std::endl;
  std::cout<<"Number of events passing mass window cut from WZ: "<<nNormWZ<<std::endl;
  std::cout<<"Number of events passing mass window cut from WJets: "<<nNormWJets<<std::endl;
  std::cout<<"Number of events passing mass window cut from TT: "<<nNormTT<<std::endl;
  std::cout<<"Number of events passing mass window cut from TTZ: "<<nNormTTZ<<std::endl;

  std::cout<<"Number of same-sign events passing mass window cut from DY: "<<nSSNormDY<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from WZ: "<<nSSNormWZ<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from WJets: "<<nSSNormWJets<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from TT: "<<nSSNormTT<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from TTZ: "<<nSSNormTTZ<<std::endl;

  /*ADD CODE TO PRINT OUT THE TOTAL PURITY YOU WILL HAVE TO DEFINE YOUR OWN VARIABLES AND PRINT THEM OUT. 
   YOU PROBABLY WANT ONE FOR TOTAL EVENTS AND ONE FOR NON DY EVENTS*/

}
