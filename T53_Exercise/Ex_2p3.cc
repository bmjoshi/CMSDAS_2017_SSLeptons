#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "ObjectID.C"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z


void Ex_2p3(){


  ///PLACE TO ADD IN NUMBERS FOR PROMPT AND FAKE RATE
  float f_e  = 0.1;
  float p_e  = 1.0;
  float f_mu = 0.1;
  float p_mu = 1.0;

  //variables we will use later based on that:
  float eps_e  = (f_e) / (1 - f_e);
  float eta_e  = (1 - p_e) / (p_e);
  float eps_mu = (f_mu) / (1 - f_mu);
  float eta_mu = (1 - p_mu) / (p_mu);

  //total numver of events
  float Npf_data;
  float Nfp_data;
  float Nff_data;
  float Nfp_mc;
  float Nff_mc;

  //load the 'data' and mc
  TFile* fdata = new TFile("ljmet_tree_TT1_temp.root");
  //TFile* fmc   = new TFile("ljmet_tree_TT2.root");

  TTree* tdata = fdata->Get("ljmet");
  //  TTree* tmc   = fmc->Get("ljmet");


  // LOAD THINGS FOR DATA
  //histograms
  TH1F* ttHThist_data = new TH1F("ttHThist_data","H_{T} for events with two tight leptons",200,1000.,2000);
  TH1F* tlHThist_data = new TH1F("tlHThist_data","H_{T} for events with one tight and one loose leptons",200,1000.,2000);
  TH1F* llHThist_data = new TH1F("llHThist_data","H_{T} for events with two loose leptons",200,1000.,2000);

  //count of events with tight/loose leptons
  int Ntt_data = 0;
  int Ntl_data = 0;
  int Nll_data = 0;
  //count events for each cut
  int Nmasscut_data=0;
  int Nmetcut_data=0;
  int NHTcut_data=0;
  int Ndilepcut_data=0;
  int NSSdilepcut_data=0;
  int NJetscut_data=0;
  int NJetpt2cut_data=0;
  int NJetpt1cut_data=0;

 

  //initialize variables
  int nEntriesdata = tdata->GetEntries();
  //  int nEntriesmc   = tmc->GetEntries();

  //electron charges
  vector<int> *elCharge_data = 0;
  //invariant mass
  //vector<double> *diElMass_data = 0;
  //kinematic variables
  vector<double> *elpts_data = 0;
  vector<double> *elEtas_data =0;
  //variables for tracking cuts
  vector<double> *elDeta_data =0;
  vector<double> *elDphi_data =0;
  //variables for primary vtx cuts
  vector<double> *elDZs_data = 0;
  vector<double> *elD0s_data = 0;
  // H over E
  vector<double> *elHoverEs_data = 0;
  //missing hits
  vector<int> *elMHits_data = 0;
  //ooemoop
  vector<double> *elOoEmooPs_data = 0;
  //charged isolation
  vector<double> *elRelIsos_data = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetas_data = 0;
  //jets
  vector<double>* jetPts_data = 0;
  //met
  double met_data = 0;

  //set branch addresses
  //tdata->SetBranchAddress("diElMass_DileptonCalc", &diElMass_data);
  tdata->SetBranchAddress("elPt_DileptonCalc", &elpts_data);
  tdata->SetBranchAddress("elEta_DileptonCalc", &elEtas_data);
  tdata->SetBranchAddress("elDeta_DileptonCalc", &elDeta_data);
  tdata->SetBranchAddress("elDphi_DileptonCalc", &elDphi_data);
  tdata->SetBranchAddress("elDZ_DileptonCalc", &elDZs_data);
  tdata->SetBranchAddress("elD0_DileptonCalc", &elD0s_data);
  tdata->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs_data);
  tdata->SetBranchAddress("elMHits_DileptonCalc",&elMHits_data);
  tdata->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs_data);
  tdata->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos_data);
  tdata->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas_data);
  tdata->SetBranchAddress("corr_met_DileptonCalc",&met_data);
  tdata->SetBranchAddress("AK5JetPt_DileptonCalc",&jetPts_data);
  tdata->SetBranchAddress("elCharge_DileptonCalc",&elCharge_data);
  // tdata->SetBranchAddress("elCharge2_DileptonCalc",&elCharge2_data);
  

  //LOAD THINGS FOR MC

  vector<double> *diElMass_mc = 0;
  //kinematic variables
  vector<double> *elpts_mc = 0;
  vector<double> *elEtas_mc =0;
  //variables for tracking cuts
  vector<double> *elDeta_mc =0;
  vector<double> *elDphi_mc =0;
  //variables for primary vtx cuts
  vector<double> *elDZs_mc = 0;
  vector<double> *elD0s_mc = 0;
  // H over E
  vector<double> *elHoverEs_mc = 0;
  //missing hits
  vector<int> *elMHits_mc = 0;
  //ooemoop
  vector<double> *elOoEmooPs_mc = 0;
  //charged isolation
  vector<double> *elRelIsos_mc = 0;
  //sigmaIetaIeta
  vector<double>* elSigmaIetaIetas_mc = 0;

  //set branch addresses
  //tmc->SetBranchAddress("diElMass_DileptonCalc", &diElMass_mc);
  //tmc->SetBranchAddress("elPt_DileptonCalc", &elpts_mc);
  //tmc->SetBranchAddress("elEta_DileptonCalc", &elEtas_mc);
  //tmc->SetBranchAddress("elDeta_DileptonCalc", &elDeta_mc);
  //tmc->SetBranchAddress("elDphi_DileptonCalc", &elDphi_mc);
  //tmc->SetBranchAddress("elDZ_DileptonCalc", &elDZs_mc);
  //tmc->SetBranchAddress("elD0_DileptonCalc", &elD0s_mc);
  //tmc->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs_mc);
  //tmc->SetBranchAddress("elMHits_DileptonCalc",&elMHits_mc);
  //tmc->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs_mc);
  //tmc->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos_mc);
  //tmc->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas_mc);


  /*
    DATA PROCEDURE: WE JUST DO THE ANALYSIS BUT GET NUMBERS FOR HOW MANY LOOSE
    AND HOW MANY TIGHT LEPTONS THERE ARE
   */

  //event loop - plots to produce: HT distribution for all cuts except 
  for(int ient = 0; ient < nEntriesdata; ient++){

    tdata->GetEntry(ient);
    /*
      Since the event cuts not having to do with leptons do not require loops, while those for leptons do, let's apply all the non-loop needing
      cuts first in order to speed things up.
     */
    //check met req
    if( met_data < 100) continue;
    Nmetcut_data +=1;

    //check lep1pt req
    if(elpts_data->at(0) < 80) continue;

    //check subleading lep pt req
    if(elpts_data->at(1) < 30) continue;

   //require more than one jet
    if(jetPts_data->size() < 2) continue;
    NJetscut_data +=1;

    //check for high pt jet
    if(jetPts_data->at(0) < 150) continue;
    NJetpt1cut_data+=1;

    //check for second jet
    if(jetPts_data->at(1) < 50) continue;
    NJetpt2cut_data+=1;

    //check HT req
    float HT = 0;
    for(int ijet = 0; ijet<jetPts_data->size(); ijet++){
      HT+= fabs(jetPts_data->at(ijet));
    }

    if( HT < 800) continue;
    NHTcut_data+=1;

    //Put electrons back together into coherent objects and add to lepton vector
    vector <Lepton*> vLep;
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
      el->isEl              = true;
      vLep.push_back(el);
    }

    //add muons to lepton vector
    for (unsigned int uiMu = 0; uiMu <muPts->size(); uiMu++){
      Muon* mu = new Muon;

      mu->pt      = muPts->at(uiMu);
      mu->eta     = muEtas->at(iuMu);
      mu->phi     = muPhis->at(uiMu);
      mu->isLoose = muLooseIDs->at(uiMu);
      mu->isTight = muTightIDs->at(uiMu);
      mu->isMu    = true;

      vLep.push_back(mu);
    }

    for(unsigned int ui = 0; ui<VLep.size();ui++){
      //loop over all electrons and make sure none come from Z
      bool foundPair = false; 
      for(unsigned int uj = ui + 1; uj < VLep.size(); uj++){
	//	if (!VLep.at(uj)->tight()) continue;

	TLorentzVector v1, v2;
	if(VLep.at(ui)->isEl){
	  v1.SetPtEtaPhiM(VLep.at(ui)->pt, VLep.at(ui)->eta, VLep.at(ui)->phi, M_EL);
	}
	else{
	  v1.SetPtEtaPhiM(VLep.at(ui)->pt, VLep.at(ui)->eta, VLep.at(ui)->phi, M_MU);
	}
	if(VLep.at(uj)->isEl){
	  v2.SetPtEtaPhiM(VLep.at(uj)->pt, VLep.at(uj)->eta, VLep.at(uj)->phi, M_EL);
	}
	else{
	  v2.SetPtEtaPhiM(VLep.at(uj)->pt, VLep.at(uj)->eta, VLep.at(uj)->phi, M_MU);
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

    for(unsigned int uiLep = 0; uiLep<VLep.size(); uiLep++){
      //make sure lepton at least passes loose requirement
      if(! ( VLep.at(uiLep)->tight() || vLep.at(uiLep)->loose())) continue;
      //get charge
      int charge1 = VLep.at(uiLep)->charge;
      for(unsigned int ujLep = uiLep+1; ujLep < VLep.size(); ujLep++){
	//make sure second lepton is at least loose
	if(! ( VLep.at(ujLep)->tight() || vLep.at(ujLep)->loose())) continue;
	if( charge1 == VLep.at(ujLep)->charge) samesign = true;
      }
      if(samesign){
	vSSLep.push_back(VLep.at(uiLep));
	vSSLep.push_back(VLep.at(ujLep));
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
    if( vSSLep.at(0)->isMu && vSSLep.at(1)->isMu) mumu = true;

    //now get lepton id information
    if(emu){
      if(vSSLep.at(0)->tight()){
	if(vSSLep.at(1)->tight()) Ntt_data +=1;
	else{
	  if(vSSLep.at(0).isEl){
	    if(vSSLep.at(1).isEl){
	      //???
	    }
	    else if(vSSLep.at(1).isMu){
	      Ntl_data+=1;
	    }
	  }
	}
      }
    }

    /*


  //  std::cout<<"Number of events passing dilepton cut: "<<Ndilepcut_data<<std::endl;
  //std::cout<<"Number of events passing SS dilepton cut: "<<NSSdilepcut_data<<std::endl;
  std::cout<<"Number of eventspassing mass cut: "<<Nmasscut_data<<std::endl;
  std::cout<<"Number of eventspassing met cut: "<<Nmetcut_data<<std::endl;
  std::cout<<"Number of eventspassing NJets cut: "<<NJetscut_data<<std::endl;
  std::cout<<"Number of eventspassing Jetpt1 cut: "<<NJetpt1cut_data<<std::endl;
  std::cout<<"Number of eventspassing Jetpt2 cut: "<<NJetpt2cut_data<<std::endl;
  std::cout<<"Number of events passing HT cut: "<<NHTcut_data<<std::endl;


  std::cout<<"Number of tight-tight events: "<<Ntt_data<<std::endl;
  std::cout<<"Number of tight-loose events: "<<Ntl_data<<std::endl;
  std::cout<<"Number of loose-loose events: "<<Nll_data<<std::endl;



  //NOW CONVERT THESE INTO HOW MANY PASSING EVENTS ARE FROM BACKGROUND WITHOUT A PROMPT LEPTON
  
  Npf_data = ( eps_mu / ( (1-eps_e*eta_e) * ( 1 - eps_mu *eta_mu))) * ( (-eps_e * Nll_data) + (Ntl_data) + (eps_e * eta_mu * Nlt_data) - (eta_mu*Ntt_data));

  Nfp_data = ( eps_e / ( (1-eps_e*eta_e) * ( 1 - eps_mu *eta_mu))) * ( (-eps_mu * Nll_data) + (eta_e * eps_mu * Ntl_data) + (Nlt_data) - (eta_e*Ntt_data) );
  
  Nff_data = ( (eps_mu * eps_e) / ( (1 - eps_e*eta_e)*(1 - eps_mu*eta_mu))) * ( (Nll_data) - (eta_e*Ntldata) - (eta_mu*Nlt_data) + (eta_e*eta_mu*Ntt_data));

  std::cout<<"Nfp_pass: "<<Nfp_data<<std::endl;
  std::cout<<"Nff_pass: "<<Nff_data<<std::endl;


  //then, N_fp plus N_ff should equal the number of events we have with two tight leptons so let's check:
  float N_fake = Nfp_data + Nff_data;
  std::cout<<"predicted number of non-prompt events: "<<N_fake<<std::endl;
  std::cout<<"observed number of non-prompt events: "<<Ntt_data<<std::endl;
  /* 
     MONTE CARLO PART. HERE WE WANT TO APPLY ALL THE ANALYSIS CUTS AND SIMPLY SEE HOW MANY EVENTS ACTUALLY DO PASS OUR SELECTION
   */

  }
}
