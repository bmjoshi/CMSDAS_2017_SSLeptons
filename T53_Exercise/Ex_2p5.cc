#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "ObjectID.C"
#include <vector>
#include "TLorentzVector.h"
#include "TChain.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z


void Ex_2p5(){


  //make tchain
  TChain* tsig = new TChain("ljmet");
  
  //add files to tchain
  tsig->Add("/uscms_data/d3/clint/public/ljmet_tree_T53_right_M-1000.root");
  tsig->Add("/uscms_data/d3/clint/public/ljmet_tree_T53_left_M-1000.root");
  //output file
  TFile* fsig = new TFile("sig_HT.root","RECREATE");




  // LOAD THINGS FOR DATA
  //histograms
  TH1F* HTHist = new TH1F("HTHist","H_{T} for events with two tight leptons",40,0.,2000);


  //count events for each cut
  int Nmasscut=0;
  int Nmetcut=0;
  int NHTcut=0;
  int Ndilepcut=0;
  int NSSdilepcut=0;
  int NJetscut=0;
  int NJetpt2cut=0;
  int NJetpt1cut=0;

 

  //initialize variables
  int nEntries = tsig->GetEntries();
  //  int nEntriesmc   = tmc->GetEntries();

  //electron charges
  vector<int> *elCharge = 0;
  //invariant mass
  //vector<double> *diElMass = 0;
  //kinematic variables
  vector<double> *elPts = 0;
  vector<double> *elPhis =0;
  vector<double> *elEtas =0;
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


  //set branch addresses
  //tsig->SetBranchAddress("diElMass_DileptonCalc", &diElMass);
  tsig->SetBranchAddress("elPt_DileptonCalc", &elPts);
  tsig->SetBranchAddress("elEta_DileptonCalc", &elEtas);
  tsig->SetBranchAddress("elPhi_DileptonCalc", &elPhis);
  tsig->SetBranchAddress("elDeta_DileptonCalc", &elDeta);
  tsig->SetBranchAddress("elDphi_DileptonCalc", &elDphi);
  tsig->SetBranchAddress("elDZ_DileptonCalc", &elDZs);
  tsig->SetBranchAddress("elD0_DileptonCalc", &elD0s);
  tsig->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs);
  tsig->SetBranchAddress("elMHits_DileptonCalc",&elMHits);
  tsig->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs);
  tsig->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos);
  tsig->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas);
  tsig->SetBranchAddress("corr_met_DileptonCalc",&met);
  tsig->SetBranchAddress("AK5JetPt_DileptonCalc",&jetPts);
  tsig->SetBranchAddress("elCharge_DileptonCalc",&elCharge);
  tsig->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistency);
  tsig->SetBranchAddress("muPt_DileptonCalc",&muPts);
  tsig->SetBranchAddress("muEta_DileptonCalc",&muEtas);
  tsig->SetBranchAddress("muPhi_DileptonCalc",&muPhis);
  tsig->SetBranchAddress("muIsLoose_DileptonCalc",&muIsLoose);
  tsig->SetBranchAddress("muIsTight_DileptonCalc",&muIsTight);
  tsig->SetBranchAddress("muCharge_DileptonCalc",&muCharge);
  

  //event loop
  for(int ient = 0; ient < nEntries; ient++){
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;
    tsig->GetEntry(ient);
    /*
      Since the event cuts not having to do with leptons do not require loops, while those for leptons do, let's apply all the non-loop needing
      cuts first in order to speed things up.
     */
    //check met req
    if( met < 100) continue;
    Nmetcut +=1;


   //require more than one jet
    if(jetPts->size() < 2) continue;
    NJetscut +=1;

    //check for high pt jet
    if(jetPts->at(0) < 150) continue;
    NJetpt1cut+=1;

    //check for second jet
    if(jetPts->at(1) < 50) continue;
    NJetpt2cut+=1;

    //check HT req
    float HT = 0;
    for(unsigned int ijet = 0; ijet<jetPts->size(); ijet++){
      HT+= fabs(jetPts->at(ijet));
    }

    if( HT < 800) continue;
    NHTcut+=1;

    //Put electrons back together into coherent objects and add to lepton vector

    vector <Lepton*> vLep;
    //vector <Electron*> vEl;
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
      Lepton* lep = el;
      lep->isEl = true;
      lep->isMu = false;
      lep->Tight = el->tight();
      lep->Loose = el->loose();
      vLep.push_back(lep);
    }

    //add muons to lepton vector
    for (unsigned int uiMu = 0; uiMu <muPts->size(); uiMu++){
      Muon* mu = new Muon;

      mu->pt      = muPts->at(uiMu);
      mu->eta     = muEtas->at(uiMu);
      mu->phi     = muPhis->at(uiMu);
      mu->isLoose = muIsLoose->at(uiMu);
      mu->isTight = muIsTight->at(uiMu);
      mu->charge  = muCharge->at(uiMu);

      Lepton* lep = mu;
      lep->isEl    = false;
      lep->isMu    = true;
      lep->Tight   = mu->tight();
      lep->Loose   = mu->loose();
      vLep.push_back(lep);
    }

    //loop over all leptons and make sure none come from Z
    bool foundPair = false; 
    for(unsigned int ui = 0; ui<vLep.size();ui++){
      //skip loose leptons
      if (!vLep.at(ui)->Tight) continue;

      for(unsigned int uj = ui + 1; uj < vLep.size(); uj++){
	//skip loose leptons
	if (!vLep.at(uj)->Tight) continue;

	TLorentzVector v1, v2;
	if(vLep.at(ui)->isEl){
	  v1.SetPtEtaPhiM(vLep.at(ui)->pt, vLep.at(ui)->eta, vLep.at(ui)->phi, M_EL);
	}
	else{
	  v1.SetPtEtaPhiM(vLep.at(ui)->pt, vLep.at(ui)->eta, vLep.at(ui)->phi, M_MU);
	}
	if(vLep.at(uj)->isEl){
	  v2.SetPtEtaPhiM(vLep.at(uj)->pt, vLep.at(uj)->eta, vLep.at(uj)->phi, M_EL);
	}
	else{
	  v2.SetPtEtaPhiM(vLep.at(uj)->pt, vLep.at(uj)->eta, vLep.at(uj)->phi, M_MU);
	}

	double mass = (v1+v2).M();
	if (mass > M_Z - dM && mass < M_Z + dM){
	  foundPair = true;
	}

      }//End loop over second lepton
      if (foundPair) break;
    }//End loop over first lepton

    //now skip event if there is a pair which comes from z boson
    if(foundPair) continue;

    //now check for same sign:
    bool samesign = false;   
    //vector to hold the leptons that are same sign
    vector<Lepton*> vSSLep;

    for(unsigned int uiLep = 0; uiLep<vLep.size(); uiLep++){
      //make sure first lepton is tight
      if(!vLep.at(uiLep)->Tight) continue;
      //get charge
      int charge1 = vLep.at(uiLep)->charge;
      for(unsigned int ujLep = uiLep+1; ujLep < vLep.size(); ujLep++){
	//make sure second lepton is tight
	if(!vLep.at(ujLep)->Tight) continue;
	if( charge1 == vLep.at(ujLep)->charge) samesign = true;
	if(samesign){
	  //	  std::cout<<"found samesign event"<<std::endl;
	  vSSLep.push_back(vLep.at(uiLep));
	  vSSLep.push_back(vLep.at(ujLep));
	}
      }

      if(samesign) break;
    }
    
    //skip event without same sign leptons
    if(!samesign) continue;
    
    //now only events this far in loop have passed all cuts so plot HT
    HTHist->Fill(HT);

    
    //end event loop
  }


  //  std::cout<<"Number of events passing dilepton cut: "<<Ndilepcut<<std::endl;
  //std::cout<<"Number of events passing SS dilepton cut: "<<NSSdilepcut<<std::endl;
  std::cout<<"Number of eventspassing mass cut: "<<Nmasscut<<std::endl;
  std::cout<<"Number of eventspassing met cut: "<<Nmetcut<<std::endl;
  std::cout<<"Number of eventspassing NJets cut: "<<NJetscut<<std::endl;
  std::cout<<"Number of eventspassing Jetpt1 cut: "<<NJetpt1cut<<std::endl;
  std::cout<<"Number of eventspassing Jetpt2 cut: "<<NJetpt2cut<<std::endl;
  std::cout<<"Number of events passing HT cut: "<<NHTcut<<std::endl;


  fsig->Write();
  fsig->Close();


    
}
