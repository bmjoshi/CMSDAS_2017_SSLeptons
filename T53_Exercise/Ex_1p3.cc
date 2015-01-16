#include <iostream>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "ObjectID.C"
#include "TLorentzVector.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z

double EtaWeight(double* weights,double eta){
  float fbin = fabs(-3.0-eta)/0.4;
  std::cout<<"eta diff "<<fbin<<std::endl;
  int bin = (int) fbin;
  return weights[bin];
};

void Ex_1p3(){


  /*add in charge misID rate you measured earlier, weights should conform to the following binning for eta:
    (-3.0 to -2.6, -2.6 to -2.2, -2.2 to -1.8, -1.8 to -1.4, -1.4 to -1.0, -1.0 to -0.6, -0.6 to -0.2, -0.2 to 0.2, 0.2 to 0.6, 0.6 to 1.0, 1.0 to 1.4, 1.4 to 1.8, 1.8 to 2.2, 2.2 to 2.6, 2.6 to 3.0)
   */
  //uncomment when ready to use
  double weights[15] = {0,0.00310209,0.00319882,0.002216,0.00107367,0.000806831, 0.000135953,0.000669693,0.000472845,0.000678656,0.000774127,0.00261718,0.00264455,0.00472953,0};



  //load 'data' file - in this case ttbar mc
  TFile* f = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  TTree* t = (TTree*)f->Get("ljmet");

  //output files
  TFile* fbg = new TFile("bg_chargeMisID.root","RECREATE");
  TTree* T = new TTree("T","test");

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
  //jets
  vector<double>* jetPts = 0;
  //met
  double met = 0;
  //muon pt;
  vector<double> *muPts = 0;
  vector<double> *muEtas = 0;
  vector<double> *muPhis = 0;
  vector<int> * muIsTight = 0;
  vector<int> * muIsLoose = 0;
  vector<int> * muCharge =0;

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
  t->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistency);
  t->SetBranchAddress("muPt_DileptonCalc",&muPts);
  t->SetBranchAddress("muEta_DileptonCalc",&muEtas);
  t->SetBranchAddress("muPhi_DileptonCalc",&muPhis);
  t->SetBranchAddress("muIsLoose_DileptonCalc",&muIsLoose);
  t->SetBranchAddress("muIsTight_DileptonCalc",&muIsTight);
  t->SetBranchAddress("muCharge_DileptonCalc",&muCharge);

  //Histograms
  //TH1F* ssEtaHist = new TH1F("ssEtaHist","#eta",30,-3,3);
  //TH1F* osEtaHist = new TH1F("osEtaHist","#eta",30,-3,3);

  //TH1F* ssPtHist = new TH1F("ssPtHist","p_{T}",100,0.,200.);
  //TH1F* osPtHist = new TH1F("osPtHist","p_{T}",100,0.,200.);

  TH1F* ssHTHist = new TH1F("ssHTHist","Weight H_{T} of charge misID events",40,0.,2000);
  TH1F* osHTHist = new TH1F("osHTHist","H_{T} of opposite sign dilepton events",40,0.,2000);


  //Loop over the tree and look for an electron pair that makes a Z
  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;

    //first apply cuts that don't need a loop
        //check met req
    if( met < 130) continue;
    //Nmetcut +=1;


   //require more than one jet
    if(jetPts->size() < 2) continue;
    //NJetscut +=1;

    //check for high pt jet
    if(jetPts->at(0) < 240) continue;
    //NJetpt1cut+=1;

    //check for second jet
    if(jetPts->at(1) < 140) continue;
    //    NJetpt2cut+=1;

    //check HT req
    float HT = 0;
    for(unsigned int ijet = 0; ijet<jetPts->size(); ijet++){
      HT+= fabs(jetPts->at(ijet));
    }

    if( HT < 700) continue;

    /*
    //now make ST cut;
    float ST = HT;
    for(unsigned int uiEl = 0; uiEl < elPts->size(); uiEl++){
      ST+=elPts->at(uiEl);
    }

    ST+=met;

    if(ST<900) continue; //ST cut
    */

    //Put electrons back together into coherent objects
    vector<Electron*> vEl;
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
    vector<Electron*> vTightEl;
    //vector for tight leptons
    vector<Lepton*> vTightLep;
    for(unsigned int ui = 0; ui < vEl.size(); ui++){
      //Apply tight selection to the electron
      if (!vEl.at(ui)->tight()) continue;
      //save tight leptons
      vTightEl.push_back(vEl.at(ui));
      Lepton* lep = vEl.at(ui);
      lep->isMu = false;
      lep->isEl = true;
      lep->Tight = true;
      vTightLep.push_back(lep);
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

    //make muon vector
    vector<Muon*> vTightMu;
    for (unsigned int uiMu = 0; uiMu <muPts->size(); uiMu++){
      Muon* mu = new Muon;

      mu->pt      = muPts->at(uiMu);
      mu->eta     = muEtas->at(uiMu);
      mu->phi     = muPhis->at(uiMu);
      mu->isLoose = muIsLoose->at(uiMu);
      mu->isTight = muIsTight->at(uiMu);
      mu->charge  = muCharge->at(uiMu);
      if(mu->isTight){
	//only save tight muons
	vTightMu.push_back(mu);
	//make tight lepton
	Lepton* lep = mu;
	lep->isEl    = false;
	lep->isMu    = true;
	lep->Tight   = mu->tight();
	lep->Loose   = mu->loose();
	vTightLep.push_back(lep);
      }
    }

    bool foundMuPair=false;
    //veto events with z boson
    for(unsigned int uimu = 0; uimu < vTightMu.size(); uimu++){
      for(unsigned int ujmu = uimu+1; ujmu < vTightMu.size(); ujmu++){

	TLorentzVector v1, v2;
	v1.SetPtEtaPhiM(vTightMu.at(uimu)->pt, vTightMu.at(uimu)->eta, vTightMu.at(uimu)->phi, M_MU);
	v2.SetPtEtaPhiM(vTightMu.at(ujmu)->pt, vTightMu.at(ujmu)->eta, vTightMu .at(ujmu)->phi, M_MU);
	
	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundMuPair = true;
	}
      }
      if(foundMuPair) break; //no need to continue if a pair has been found	
    }

    //skip if zboson
    if(foundMuPair) continue;

    //now check for only two tight leptons
    if(vTightLep.size()!=2) continue;
    //check OPPOSITE sign
    if(vTightLep.at(0)->charge==vTightLep.at(1)->charge) continue;

        
    //cut harder on lepton pt
    float temp_pt = 0;
    for(unsigned int i = 0; i < vTightLep.size();i++){
      if(vTightLep.at(i)->pt > temp_pt) temp_pt = vTightLep.at(i)->pt;
    }

    //check to make sure hardest lepton is above 100 GeV
    if(temp_pt<140) continue;

    float sec_pt = 0;
    unsigned int lep2 = 0;
    for(unsigned int i = 0 ; i < vTightLep.size(); i++){
      if((vTightLep.at(i)->pt > sec_pt) && (vTightLep.at(i)->pt!=temp_pt)) {
	lep2=i;
	sec_pt = vTightLep.at(i)->pt;
      }
    }

    if(sec_pt < 100) continue;

    TLorentzVector v1, v2;
    if(vTightLep.at(0)->isEl){
      v1.SetPtEtaPhiM(vTightLep.at(0)->pt, vTightLep.at(0)->eta, vTightLep.at(0)->phi, M_EL);
    }
    else{
      v1.SetPtEtaPhiM(vTightLep.at(0)->pt, vTightLep.at(0)->eta, vTightLep.at(0)->phi, M_MU);
    }
    if(vTightLep.at(1)->isEl){
      v2.SetPtEtaPhiM(vTightLep.at(1)->pt, vTightLep.at(1)->eta, vTightLep.at(1)->phi, M_EL);
    }
    else{
      v2.SetPtEtaPhiM(vTightLep.at(1)->pt, vTightLep.at(1)->eta, vTightLep.at(1)->phi, M_MU);
    }

    
    //cut on dilepton mass
    if( (v1+v2).M() <250) continue;



    //check channel
    bool ee = false;
    bool emu = false;
    bool mumu = false;

    if(vTightLep.at(0)->isEl && vTightLep.at(1)->isEl) ee = true;
    if(vTightLep.at(0)->isEl && vTightLep.at(1)->isMu) emu = true;
    if(vTightLep.at(0)->isMu && vTightLep.at(1)->isEl) emu = true;
    if(vTightLep.at(0)->isMu && vTightLep.at(1)->isMu) mumu = true;

    //skip any di-mu events since we assume they don't contribute to ss via charge misID
    if(mumu) continue;

    /*now we want to fill our predicted HT distribution with the os events weighted by the appropriate weights
      I've written a simple function that gives the correct element of the weight array if you pass it the eta of the electron.
      Usage is simply EtaWeight(weightsarray, eta). Take advantage or it to make life easier :)
     */

    // ADD CODE HERE TO CALCULATE PROBABLITY FOR EACH EVENT

    float ee_weight = EtaWeight(weights, vTightLep.at(0)->eta) + EtaWeight(weights, vTightLep.at(1)->eta) - EtaWeight(weights, vTightLep.at(0)->eta)*EtaWeight(weights, vTightLep.at(1)->eta);;// fill in
    float emu_weight;// = EtaWeight(weights,vTightLep.at(;// fill in
    if(vTightLep.at(0)->isEl)
      emu_weight = EtaWeight(weights, vTightLep.at(0)->eta); //misId rate for mu is zero                                                                                                                          
    else
      emu_weight = EtaWeight(weights, vTightLep.at(1)->eta);


    //fill HT histogram
    if(ee){
      ssHTHist->Fill(HT,ee_weight);
    }
    if(emu){
      ssHTHist->Fill(HT,emu_weight);
    }

    //fill os hist for reference
    osHTHist->Fill(HT);

    //end event loop
  }

  TCanvas c1;
  osHTHist->Draw();
  c1.Print("HT_oppositeSignEvents.pdf");

  TCanvas c2;
  ssHTHist->Draw();
  c2.Print("HT_sameSignEvents.pdf");

  //finally let's save our predictions to a root file for ease of use later
  T->Branch("osHTHist","TH1F",&osHTHist,32000,0);
  T->Branch("ssHTHist","TH1F",&ssHTHist,32000,0);
  T->Fill();
  T->Print();
  fbg->Write();
  fbg->Close();


}


