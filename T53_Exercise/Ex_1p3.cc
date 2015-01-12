#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"

void Ex_1p3(){

  //load 'data' file - in this case ttbar mc
  TFile* f = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT.root");
  TTree* t = (TTree*)f->Get("ljmet");

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
  //jets
  vector<double>* jetPts = 0;
  //met
  double met = 0;

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
  t->SetBranchAddress("corr_met_DileptonCalc",&met);
  t->SetBranchAddress("AK5JetPt_DileptonCalc",&jetPts);

  //Histograms
  TH1F* ssEtaHist = new TH1F("ssEtaHist","#eta",30,-3,3);
  TH1F* osEtaHist = new TH1F("osEtaHist","#eta",30,-3,3);

  TH1F* ssPtHist = new TH1F("ssPtHist","p_{T}",100,0.,200.);
  TH1F* osPtHist = new TH1F("osPtHist","p_{T}",100,0.,200.);

  TH1F* ssHTHist = new TH1F("ssHTHist","Weight H_{T} of charge misID events",500,0.,1500);
  TH1F* osHTHist = new TH1F("osHTHist","H_{T} of opposite sign dilepton events",500,0.,1500);


  //Loop over the tree and look for an electron pair that makes a Z
  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;

    //first apply cuts that don't need a loop
        //check met req
    if( met < 100) continue;
    //Nmetcut +=1;


   //require more than one jet
    if(jetPts->size() < 2) continue;
    //NJetscut +=1;

    //check for high pt jet
    if(jetPts->at(0) < 150) continue;
    //NJetpt1cut+=1;

    //check for second jet
    if(jetPts->at(1) < 50) continue;
    //    NJetpt2cut+=1;

    //check HT req
    float HT = 0;
    for(int ijet = 0; ijet<jetPts->size(); ijet++){
      HT+= fabs(jetPts->at(ijet));
    }

    if( HT < 800) continue;

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

    //vector for tight electrons
    vector<Electonr*> vTightEl;

    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply tight selection to the electron
      if (!vEl.at(ui)->tight()) continue;
      //save tight leptons
      vTightEl.push_back(vEl.at(ui));
      for(unsigned int uj = ui + 1; uj < vEl.size(); uj++){
	if (!vEl.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vEl.at(ui)->pt, vEl.at(ui)->eta, vEl.at(ui)->phi, M_EL);
	v2.SetPtEtaPhiM(vEl.at(uj)->pt, vEl.at(uj)->eta, vEl.at(uj)->phi, M_EL);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	}
      }
      if(foundPair) break; //no need to continue if a pair has been found

    }

    if(foundPair) continue; //skip events with a Zboson

    //require only two tight leptons
    if(vTightEl.size()!=2) continue;
    //require opposite charges
    if(vTightEl.at(0)->charge==vTightEl.at(1)->charge()) continue;

    //now fill histograms

    //now we can fill histograms because events have passed

    //end event loop
  }

}
