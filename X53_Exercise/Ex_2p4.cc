#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "ObjectID.C"
#include <vector>
#include "TLorentzVector.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z


void Ex_2p4(){


  ///PLACE TO ADD IN NUMBERS FOR PROMPT AND FAKE RATE
  float f_e  = 0.36;
  float p_e  = 0.8;
  float f_mu = 0.97;
  float p_mu = 1.0;

  //variables we will use later based on that:
  float eps_e  = (f_e) / (1 - f_e);
  float eta_e  = (1 - p_e) / (p_e);
  float eps_mu = (f_mu) / (1 - f_mu);
  float eta_mu = (1 - p_mu) / (p_mu);

  //total number of events
  float Npf_emu=0;
  float Nfp_emu=0;
  float Nff_emu=0;
  float Nfp_ee=0;
  float Nff_ee=0;
  float Nfp_mumu=0;
  float Nff_mumu=0;

  //load the 'data' and mc
  TFile* fdata = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  //TFile* fmc   = new TFile("ljmet_tree_TT2.root");

  TTree* tdata = (TTree*)fdata->Get("ljmet");
  //  TTree* tmc   = fmc->Get("ljmet");


  // LOAD THINGS FOR DATA
  //histograms
  TH1F* ttHThist = new TH1F("ttHThist","H_{T} for events with two tight leptons",200,1000.,2000);
  TH1F* tlHThist = new TH1F("tlHThist","H_{T} for events with one tight and one loose leptons",200,1000.,2000);
  TH1F* llHThist = new TH1F("llHThist","H_{T} for events with two loose leptons",200,1000.,2000);

  //count of events with tight/loose leptons
  int Ntt_ee = 0;
  int Ntt_emu = 0;
  int Ntt_mumu = 0;
  //emu channel, left/right subscripts for t/l correspond to e/mu
  int Ntl_emu = 0;
  int Nlt_emu = 0;
  int Nll_emu = 0;
  //mu mu channel
  int Ntl_mumu=0;
  int Nll_mumu=0;
  //ee channel
  int Ntl_ee=0;
  int Nll_ee=0;

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
  int nEntries = tdata->GetEntries();
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
  //tdata->SetBranchAddress("diElMass_DileptonCalc", &diElMass);
  tdata->SetBranchAddress("elPt_DileptonCalc", &elPts);
  tdata->SetBranchAddress("elEta_DileptonCalc", &elEtas);
  tdata->SetBranchAddress("elPhi_DileptonCalc", &elPhis);
  tdata->SetBranchAddress("elDeta_DileptonCalc", &elDeta);
  tdata->SetBranchAddress("elDphi_DileptonCalc", &elDphi);
  tdata->SetBranchAddress("elDZ_DileptonCalc", &elDZs);
  tdata->SetBranchAddress("elD0_DileptonCalc", &elD0s);
  tdata->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs);
  tdata->SetBranchAddress("elMHits_DileptonCalc",&elMHits);
  tdata->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs);
  tdata->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos);
  tdata->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas);
  tdata->SetBranchAddress("corr_met_DileptonCalc",&met);
  tdata->SetBranchAddress("AK5JetPt_DileptonCalc",&jetPts);
  tdata->SetBranchAddress("elCharge_DileptonCalc",&elCharge);
  tdata->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistency);
  tdata->SetBranchAddress("muPt_DileptonCalc",&muPts);
  tdata->SetBranchAddress("muEta_DileptonCalc",&muEtas);
  tdata->SetBranchAddress("muPhi_DileptonCalc",&muPhis);
  tdata->SetBranchAddress("muIsLoose_DileptonCalc",&muIsLoose);
  tdata->SetBranchAddress("muIsTight_DileptonCalc",&muIsTight);
  tdata->SetBranchAddress("muCharge_DileptonCalc",&muCharge);
  

  //event loop
  for(int ient = 0; ient < nEntries; ient++){
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;
    tdata->GetEntry(ient);
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
    for(int ijet = 0; ijet<jetPts->size(); ijet++){
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

      for(unsigned int uj = ui + 1; uj < vLep.size(); uj++){
	//	if (!VLep.at(uj)->tight()) continue;

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
      if(! ( vLep.at(uiLep)->Tight || vLep.at(uiLep)->Loose)) continue;
      //get charge
      int charge1 = vLep.at(uiLep)->charge;
      for(unsigned int ujLep = uiLep+1; ujLep < vLep.size(); ujLep++){
	//make sure second lepton is at least loose
	if(! ( vLep.at(ujLep)->Tight || vLep.at(ujLep)->Loose)) continue;
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
    


    //need to separate out different channels:
    bool ee   = false;
    bool emu  = false;
    bool mumu = false;

    if( vSSLep.at(0)->isEl && vSSLep.at(1)->isEl) ee   = true;
    if( vSSLep.at(0)->isMu && vSSLep.at(1)->isEl) emu  = true;
    if( vSSLep.at(0)->isEl && vSSLep.at(1)->isMu) emu  = true; //two ways to have electron/muon
    if( vSSLep.at(0)->isMu && vSSLep.at(1)->isMu) mumu = true;

    //now get lepton id information
    if(emu){

      //check flavor of first lepton
      if(vSSLep.at(0)->isEl){
	//check if first lepton is tight
	if(vSSLep.at(0)->Tight){
	  //check if second  is
	  if(vSSLep.at(1)->Tight) Ntt_emu +=1;
	  //if second isn't tight, it must be loose (requirement for leptons to be put in sslep vector)
	  else{
	    Ntl_emu+=1; // we know the electron is tight and muon is loose
	  }
	}
	//now, we have first lepton is loose electron
	else{
	  //check for second lepton being tight
	  if(vSSLep.at(1)->Tight) Nlt_emu+=1;
	  //else we have two loose leptons
	  else Nll_emu+=1;
	}
      }
      //now handle case if first lepton is muon
      else {
	//check for first lepton to be tight
	if(vSSLep.at(0)->Tight){
	  //check for second lepton to be tight
	  if(vSSLep.at(1)->Tight) Ntt_emu +=1; //note that this case is not overcounting with Ntt_emu above because we are preserving order of leptons
	  else Nlt_emu+=1; //still required muon to be tight, so lt event
	}
	//now handle when first lepton is loose muon
	else{
	  //check if electron is tight
	  if(vSSLep.at(1)->Tight) Ntl_emu +=1;
	  else Nll_emu+=1; //again, not overcounting with above because we preserve order of leptons
	}
      }
    } //end emu channel
      

    //now ee channel
    if(ee){
      //check if first lepton is tight
      if(vSSLep.at(0)->Tight){
	//check second lepton
	if(vSSLep.at(1)->Tight) Ntt_ee +=1; //tight-tight event
	else Ntl_ee+=1; //tight-loose event
      }
      //now if first lepton is loose
      else {
	//check second lepton
	if(vSSLep.at(1)->Tight) Ntl_ee+=1; //loose-tight event
	else Nll_ee +=1; //loose-loose event
      }
    }
    //now mumu channel
    if(mumu){
      //check if first lepton is tight
      if(vSSLep.at(0)->Tight){
	//check second lepton
	if(vSSLep.at(1)->Tight) Ntt_mumu +=1; //tight-tight event
	else Ntl_mumu+=1; //tight-loose event
      }
      //now if first lepton is loose
      else {
	//check second lepton
	if(vSSLep.at(1)->Tight) Ntl_mumu+=1; //loose-tight event
	else Nll_mumu +=1; //loose-loose event
      }
    }
    
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



  //NOW CONVERT THESE INTO HOW MANY PASSING EVENTS ARE FROM BACKGROUND WITHOUT A PROMPT LEPTON
  

  //emu channel
  Npf_emu = ( eps_mu / ( (1-eps_e*eta_e) * ( 1 - eps_mu *eta_mu))) * ( (-eps_e * Nll_emu) + (Ntl_emu) + (eps_e * eta_mu * Nlt_emu) - (eta_mu*Ntt_emu));

  Nfp_emu = ( eps_e / ( (1-eps_e*eta_e) * ( 1 - eps_mu *eta_mu))) * ( (-eps_mu * Nll_emu) + (eta_e * eps_mu * Ntl_emu) + (Nlt_emu) - (eta_e*Ntt_emu) );
  
  Nff_emu = ( (eps_mu * eps_e) / ( (1 - eps_e*eta_e)*(1 - eps_mu*eta_mu))) * ( (Nll_emu) - (eta_e*Ntl_emu) - (eta_mu*Nlt_emu) + (eta_e*eta_mu*Ntt_emu));

  //ee channel

  Nfp_ee = ( (eps_e) / ( pow( 1 - (eps_e*eta_e),2))) * ( (-2*eps_e*Nll_ee) + (fabs(1+eps_e*eta_e)*Ntl_ee) - (2*eta_e*Ntt_ee));
 
  Nff_ee = (pow( (eps_e) / (1- (eps_e*eta_e)),2)) * ( (Nll_ee) - (eta_e*Ntl_ee) + (eta_e*eta_e*Ntt_ee));

  //mumu channel

  Nfp_mumu = ( (eps_mu) / ( pow( 1 - (eps_mu*eta_mu),2))) * ( (-2*eps_mu*Nll_mumu) + (fabs(1+eps_mu*eta_mu)*Ntl_mumu) - (2*eta_mu*Ntt_mumu));
 
  Nff_mumu = (pow( (eps_mu) / (1- (eps_mu*eta_mu)),2)) * ( (Nll_mumu) - (eta_mu*Ntl_mumu) + (eta_mu*eta_mu*Ntt_mumu));



  std::cout<<"Nfp_emu_pass: "<<Nfp_emu<<std::endl;
  std::cout<<"Nff_emu_pass: "<<Nff_emu<<std::endl;


  //then, N_fp plus N_ff should equal the number of events we have with two tight leptons so let's check:
  float Nfake_emu = Nfp_emu + Nff_emu + Npf_emu;
  float Nfake_ee = Nfp_ee + Nff_ee;
  float Nfake_mumu = Nfp_mumu + Nff_ee;

  std::cout<<"predicted number of non-prompt events in emu channel: "<<Nfake_emu<<std::endl;
  std::cout<<"observed number of non-prompt events in emu channel: "<<Ntt_emu<<std::endl;

  std::cout<<"predicted number of non-prompt events in ee channel: "<<Nfake_ee<<std::endl;
  std::cout<<"observed number of non-prompt events in ee channel: "<<Ntt_ee<<std::endl;

  std::cout<<"predicted number of non-prompt events in mumu channel: "<<Nfake_mumu<<std::endl;
  std::cout<<"observed number of non-prompt events in mumu channel: "<<Ntt_mumu<<std::endl;



    
}
