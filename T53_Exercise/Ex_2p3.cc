
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


void Ex_2p3(){


  ///PLACE TO ADD IN NUMBERS FOR PROMPT AND FAKE RATE
  float f_e  = 0.419847;
  float p_e  = 0.83468;
  float f_mu = 0.078;
  float p_mu = 0.828028;

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

  //make tchain
  TChain* tdata = new TChain("ljmet");
  
  //add files to tchain
  //tdata->Add("/uscms_data/d3/clint/public/ljmet_tree_DY.root");
  //tdata->Add("/uscms_data/d3/clint/public/ljmet_tree_WZ.root");
  //tdata->Add("/uscms_data/d3/clint/public/ljmet_tree_WJets.root");
  tdata->Add("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  //tdata->Add("/uscms_data/d3/clint/public/ljmet_tree_TTZ.root");

  //output file
  TFile* fbg = new TFile("bg_nonPrompt.root","RECREATE");



  // LOAD THINGS FOR DATA
  //histograms
  TH1F* HTtt_emu = new TH1F("HTtt_emu","H_{T} for events with two tight leptons",40,0,2000);
  TH1F* HTlt_emu = new TH1F("HTlt_emu","H_{T} for events with one tight and one loose leptons",40,0,2000);
  TH1F* HTtl_emu = new TH1F("HTtl_emu","H_{T} for events with one tight and one loose",40,0,2000);
  TH1F* HTll_emu = new TH1F("HTll_emu","H_{T} for events with two loose leptons",40,0,2000);

  TH1F* HTtt_mumu = new TH1F("HTtt_mumu","H_{T} for events with two tight leptons",40,0,2000);
  TH1F* HTtl_mumu = new TH1F("HTtl_mumu","H_{T} for events with one tight and one loose leptons",40,0,2000);
  TH1F* HTll_mumu = new TH1F("HTll_mumu","H_{T} for events with two loose leptons",40,0,2000);

  TH1F* HTtt_ee = new TH1F("HTtt_ee","H_{T} for events with two tight leptons",40,0,2000);
  TH1F* HTtl_ee = new TH1F("HTtl_ee","H_{T} for events with one tight and one loose leptons",40,0,2000);
  TH1F* HTll_ee = new TH1F("HTll_ee","H_{T} for events with two loose leptons",40,0,2000);

  TH1F* leadlep_pt = new TH1F("leadlep_pt","leading lepton pt",100,0,1000);
  TH1F* sublep_pt = new TH1F("sublep_pt","subleading lepton pt",100,0,1000);

  TH1F* leadjet_pt = new TH1F("leadjet_pt","leading jet pt",100,0,1000);
  TH1F* subjet_pt = new TH1F("subjet_pt","sub-leading jet pt",100,0,1000);

  TH1F* met_h = new TH1F("met_h","MET",100,0,1000);
  TH1F* HT_h = new TH1F("HT_h","HT",50,0,2000);

  TH1F* dilep_mass = new TH1F("dilep_mass","DilepMasss",120,0,1200);
  TH1F* njets = new TH1F("njets","Number of Jets",10,0,10);
  TH1F* TriMass = new TH1F("TriMass","Mass of leptons and Closest Jet",120,0,1200);

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
  vector<double>* jetEtas = 0;
  vector<double>* jetPhis = 0;
  //met
  double met = 0;
  //muon pt;
  vector<double> *muPts = 0;
  vector<double> *muEtas = 0;
  vector<double> *muPhis = 0;
  vector<int> * muIsTight = 0;
  vector<int> * muIsLoose = 0;
  vector<int> * muCharge =0;

  vector<int>    *global = 0;
  vector<double> *chi2 = 0;
  vector<int>    *nValMuHits = 0;
  vector<int>    *nMatchedStations = 0;
  vector<double> *dxy = 0;
  vector<double> *dz = 0;
  vector<int>    *nValPixelHits = 0;
  vector<int>    *nTrackerLayers = 0;
  vector<double> *relIso = 0;


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
  tdata->SetBranchAddress("AK5JetPhi_DileptonCalc",&jetPhis);
  tdata->SetBranchAddress("AK5JetEta_DileptonCalc",&jetEtas);
  tdata->SetBranchAddress("elCharge_DileptonCalc",&elCharge);
  tdata->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistency);
  tdata->SetBranchAddress("muPt_DileptonCalc",&muPts);
  tdata->SetBranchAddress("muEta_DileptonCalc",&muEtas);
  tdata->SetBranchAddress("muPhi_DileptonCalc",&muPhis);
  tdata->SetBranchAddress("muIsLoose_DileptonCalc",&muIsLoose);
  tdata->SetBranchAddress("muIsTight_DileptonCalc",&muIsTight);
  tdata->SetBranchAddress("muCharge_DileptonCalc",&muCharge);

  tdata->SetBranchAddress("muGlobal_DileptonCalc",&global);
  tdata->SetBranchAddress("muChi2_DileptonCalc",&chi2);
  tdata->SetBranchAddress("muNValMuHits_DileptonCalc",&nValMuHits);
  tdata->SetBranchAddress("muNMatchedStations_DileptonCalc",&nMatchedStations);
  tdata->SetBranchAddress("muDxy_DileptonCalc",&dxy);
  tdata->SetBranchAddress("muDz_DileptonCalc",&dz);
  tdata->SetBranchAddress("muNValPixelHits_DileptonCalc",&nValPixelHits);
  tdata->SetBranchAddress("muNTrackerLayers_DileptonCalc",&nTrackerLayers);
  tdata->SetBranchAddress("muRelIso_DileptonCalc",&relIso);
  

  //event loop
  for(int ient = 0; ient < nEntries; ient++){
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;
    tdata->GetEntry(ient);
    /*
      Since the event cuts not having to do with leptons do not require loops, while those for leptons do, let's apply all the non-loop needing
      cuts first in order to speed things up.
     */
    //check met req
    if( met < 130) continue;
    Nmetcut +=1;


   //require more than one jet
    if(jetPts->size() < 2) continue;
    NJetscut +=1;

    //check for high pt jet
    if(jetPts->at(0) < 240) continue;
    NJetpt1cut+=1;

    //check for second jet
    if(jetPts->at(1) < 140) continue;
    NJetpt2cut+=1;

    //check HT req
    float HT = 0;
    for(unsigned int ijet = 0; ijet<jetPts->size(); ijet++){
      HT+= fabs(jetPts->at(ijet));
    }

    if( HT < 700) continue;
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

      mu->global           = global->at(uiMu);
      mu->chi2             = chi2->at(uiMu);
      mu->nValMuHits       = nValMuHits->at(uiMu);
      mu->nMatchedStations = nMatchedStations->at(uiMu);
      mu->dxy              = dxy->at(uiMu);
      mu->dz               = dz->at(uiMu);
      mu->nValPixelHits    = nValPixelHits->at(uiMu);
      mu->nTrackerLayers   = nTrackerLayers->at(uiMu);
      mu->relIso           = relIso->at(uiMu);


      Lepton* lep = mu;
      lep->isEl    = false;
      lep->isMu    = true;
      lep->Tight   = mu->cutBasedTight();
      lep->Loose   = mu->cutBasedLoose();
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
    
    //cut harder on lepton pt
    float temp_pt = 0;
    for(unsigned int i = 0; i < vSSLep.size();i++){
      if(vSSLep.at(i)->pt > temp_pt) temp_pt = vSSLep.at(i)->pt;
    }



    //check to make sure hardest lepton is above 100 GeV
    if(temp_pt<140) continue;

    //check second leading lepton

    float sec_pt = 0;
    unsigned int lep2 = 0;
    for(unsigned int i = 0 ; i < vSSLep.size(); i++){
      if((vSSLep.at(i)->pt > sec_pt) && (vSSLep.at(i)->pt!=temp_pt)) {
	lep2=i;
	sec_pt = vSSLep.at(i)->pt;
      }
    }

    if(sec_pt < 100) continue;

    leadlep_pt->Fill(temp_pt);
    sublep_pt->Fill(sec_pt);
    leadjet_pt->Fill(jetPts->at(0));
    subjet_pt->Fill(jetPts->at(1));

    met_h->Fill(met);
    HT_h->Fill(HT);

    int njet=0;
    for(unsigned int ijet = 0; ijet<jetPts->size(); ijet++){
      if(jetPts->at(ijet)>40) njet+=1;
    }

    njets->Fill(njet);
    TLorentzVector v1, v2;
    if(vSSLep.at(0)->isEl){
      v1.SetPtEtaPhiM(vSSLep.at(0)->pt, vSSLep.at(0)->eta, vSSLep.at(0)->phi, M_EL);
    }
    else{
      v1.SetPtEtaPhiM(vSSLep.at(0)->pt, vSSLep.at(0)->eta, vSSLep.at(0)->phi, M_MU);
    }
    if(vSSLep.at(1)->isEl){
      v2.SetPtEtaPhiM(vSSLep.at(1)->pt, vSSLep.at(1)->eta, vSSLep.at(1)->phi, M_EL);
    }
    else{
      v2.SetPtEtaPhiM(vSSLep.at(1)->pt, vSSLep.at(1)->eta, vSSLep.at(1)->phi, M_MU);
    }
    dilep_mass->Fill( (v1+v2).M());
    
    //cut on dilepton mass
    if( (v1+v2).M() <250) continue;
   
    /*
    //check mass with two electrons and jet
    float dR = 1000;
    unsigned int closejet = -1;
    for(unsigned int i = 0; i <jetPts->size();i++){
      if(jetPts->at(i)<100) continue;
      float dr= pow( pow(jetEtas->at(i)-vSSLep.at(lep2)->eta,2) + pow( jetPhis->at(i)-vSSLep.at(lep2)->phi,2),0.5);
      if(dr<dR) closejet =i;
    }

    //now make lorentz vectors
    TLorentzVector vjet;
    vjet.SetPtEtaPhiM(jetPts->at(closejet),jetEtas->at(closejet),jetPhis->at(closejet),0);

    //mass of three
    if(dR<1){
      TriMass->Fill( (v1 + v2 + vjet).M() );
    }
*/
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
	  if(vSSLep.at(1)->Tight){
	    Ntt_emu +=1;
	    HTtt_emu->Fill(HT);
	  }
	  //if second isn't tight, it must be loose (requirement for leptons to be put in sslep vector)
	  else{
	    Ntl_emu+=1; // we know the electron is tight and muon is loose
	    HTtl_emu->Fill(HT);
	  }
	}
	//now, we have first lepton is loose electron
	else{
	  //check for second lepton being tight
	  if(vSSLep.at(1)->Tight){
	    Nlt_emu+=1;
	    HTlt_emu->Fill(HT);
	  }
	    //else we have two loose leptons
	  else {
	    Nll_emu+=1;
	    HTll_emu->Fill(HT);   
	  }
	}
      }
      //now handle case if first lepton is muon
      else {
	//check for first lepton to be tight
	if(vSSLep.at(0)->Tight){
	  //check for second lepton to be tight
	  if(vSSLep.at(1)->Tight){ 
	    Ntt_emu +=1; //note that this case is not overcounting with Ntt_emu above because we are preserving order of leptons
	    HTtt_emu->Fill(HT);
	  }
	  else {
	    Nlt_emu+=1; //still required muon to be tight, so lt event
	    HTlt_emu->Fill(HT);
	  }
	}
	//now handle when first lepton is loose muon
	else{
	  //check if electron is tight
	  if(vSSLep.at(1)->Tight) {
	    Ntl_emu +=1;
	    HTtl_emu->Fill(HT);
	  }
	  else {
	    Nll_emu+=1; //again, not overcounting with above because we preserve order of leptons
	    HTll_emu->Fill(HT);
	  }
	}
      }
    } //end emu channel
      

    //now ee channel
    if(ee){
      //check if first lepton is tight
      if(vSSLep.at(0)->Tight){
	//check second lepton
	if(vSSLep.at(1)->Tight){
	  Ntt_ee +=1; //tight-tight event
	  HTtt_ee->Fill(HT);
	}
	else {
	  Ntl_ee+=1; //tight-loose event
	  HTtl_ee->Fill(HT);
	}
      }
      //now if first lepton is loose
      else {
	//check second lepton
	if(vSSLep.at(1)->Tight) {
	  Ntl_ee+=1; //loose-tight event
	  HTtl_ee->Fill(HT);
	}
	else{
	  Nll_ee +=1; //loose-loose event
	  HTll_ee->Fill(HT);
	}
      }
    }

    //now mu mu channel
    if(mumu){
      //check if first lepton is tight
      if(vSSLep.at(0)->Tight){
	//check second lepton
	if(vSSLep.at(1)->Tight){
	  Ntt_mumu +=1; //tight-tight event
	  HTtt_mumu->Fill(HT);
	}
	else {
	  Ntl_mumu+=1; //tight-loose event
	  HTtl_mumu->Fill(HT);
	}
      }
      //now if first lepton is loose
      else {
	//check second lepton
	if(vSSLep.at(1)->Tight) {
	  Ntl_mumu+=1; //loose-tight event
	  HTtl_mumu->Fill(HT);
	}
	else{
	  Nll_mumu +=1; //loose-loose event
	  HTll_mumu->Fill(HT);
	}
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

  std::cout<<"nll_mumu "<<Nll_mumu<<" ntl_mumu "<<Ntl_mumu<<" Ntt_mumu :"<<Ntt_mumu<<std::endl;

  //CONVERT TO NUMBER OF EVENTS PREDICTED FROM EACH BACKGROUND SOURCE

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



  //scale histograms by contribution rate;
  TH1F* HT_emu = new TH1F("HT_emu","HT electron-muon",40,0.,2000);
  for(int i = 1; i <= HT_emu->GetNbinsX(); i++){
    //loop over all bins in tl, tt, lt histograms and fill ht bin with them * weights from formulas
    float pf_ht = ( eps_mu / ( (1-eps_e*eta_e) * ( 1 - eps_mu *eta_mu))) * ( (-eps_e * HTll_emu->GetBinContent(i)) + (HTlt_emu->GetBinContent(i)) + (eps_e * eta_mu * HTlt_emu->GetBinContent(i)) - (eta_mu*HTtt_emu->GetBinContent(i))); 
    float fp_ht = ( eps_e / ( (1-eps_e*eta_e) * ( 1 - eps_mu *eta_mu))) * ( (-eps_mu * HTll_emu->GetBinContent(i)) + (eta_e * eps_mu * HTlt_emu->GetBinContent(i)) + (HTlt_emu->GetBinContent(i)) - (eta_e*HTtt_emu->GetBinContent(i)) );
    float ff_ht =( (eps_mu * eps_e) / ( (1 - eps_e*eta_e)*(1 - eps_mu*eta_mu))) * ( (HTll_emu->GetBinContent(i)) - (eta_e*HTlt_emu->GetBinContent(i)) - (eta_mu*HTlt_emu->GetBinContent(i)) + (eta_e*eta_mu*HTtt_emu->GetBinContent(i)));

    float weighted_ht = pf_ht + fp_ht + ff_ht;
    HT_emu->SetBinContent(i,weighted_ht);
  }

  TH1F* HT_ee = new TH1F("HT_ee","HT di-electron",40,0.,2000);
  for(int i = 0; i<=HTtt_ee->GetNbinsX(); i++){
    float fp_ht =  ( (eps_e) / ( pow( 1 - (eps_e*eta_e),2))) * ( (-2*eps_e*HTll_ee->GetBinContent(i)) + (fabs(1+eps_e*eta_e)*HTtl_ee->GetBinContent(i)) - (2*eta_e*HTtt_ee->GetBinContent(i)));
    float ff_ht = (pow( (eps_e) / (1- (eps_e*eta_e)),2)) * ( (HTll_ee->GetBinContent(i)) - (eta_e*HTtl_ee->GetBinContent(i)) + (eta_e*eta_e*HTtt_ee->GetBinContent(i)));
    float weighted_ht = ff_ht+fp_ht;
    HT_ee->SetBinContent(i,weighted_ht);
  }

  TH1F* HT_mumu = new TH1F("HT_mumu","HT di-electron",40,0.,2000);
  for(int i = 0; i<=HTtt_mumu->GetNbinsX(); i++){
    float fp_ht =  ( (eps_mu) / ( pow( 1 - (eps_mu*eta_mu),2))) * ( (-2*eps_mu*HTll_mumu->GetBinContent(i)) + (fabs(1+eps_mu*eta_mu)*HTtl_mumu->GetBinContent(i)) - (2*eta_mu*HTtt_mumu->GetBinContent(i)));
    float ff_ht = (pow( (eps_mu) / (1- (eps_mu*eta_mu)),2)) * ( (HTll_mumu->GetBinContent(i)) - (eta_mu*HTtl_mumu->GetBinContent(i)) + (eta_mu*eta_mu*HTtt_mumu->GetBinContent(i)));
    float weighted_ht = ff_ht+fp_ht;
    HT_mumu->SetBinContent(i,weighted_ht);
  }

  //now make total HT plot
  TH1F* HT_nonPrompt = new TH1F("HT_nonPrompt","HT from prompt-fake and fake-fake events",40,0.,2000);
  HT_nonPrompt->Add(HT_emu);
  HT_nonPrompt->Add(HT_ee);
  HT_nonPrompt->Add(HT_mumu);

  //then, N_fp plus N_ff should equal the number of events we have with two tight leptons so let's check:
  float Nfake_emu = Nfp_emu + Nff_emu + Npf_emu;
  float Nfake_ee = Nfp_ee + Nff_ee;
  float Nfake_mumu = Nfp_mumu + Nff_mumu;



  std::cout<<"predicted number of non-prompt events in emu channel: "<<Nfake_emu<<std::endl;


  std::cout<<"predicted number of non-prompt events in ee channel: "<<Nfake_ee<<std::endl;


  std::cout<<"predicted number of non-prompt events in mumu channel: "<<Nfake_mumu<<std::endl;


  fbg->Write();
  fbg->Close();

    
}
