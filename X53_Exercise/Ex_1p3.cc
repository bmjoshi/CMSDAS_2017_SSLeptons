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

int getEtaBin(float abseta){
  int bin=-1;
  if(abseta>2.0) bin=6;
  else if(abseta>1.556) bin=5;
  else if(abseta>1.442) bin=4;
  else if(abseta>0.8) bin=3;
  else if(abseta>0.4) bin=2;
  else if(abseta>=0.0) bin=1;

  return bin;
}


//Feed this method a given eta and  your weight array, and this will return the corresponding weight
double EtaWeight(double* weights,double eta,double pt){
  int fbin = getEtaBin( abs(eta) );
  std::cout<<"eta diff "<<fbin<<std::endl;
  if (pt > 100.0){ //high pt 
    fbin+6;
  }
  return weights[fbin-1];
};

void Ex_1p3(){
  /*add in charge misID rate you measured earlier, weights should conform to the binning for eta from getEtaBin:
    0.0-0.4, 0.4-0.8, 0.8-1.442, 1.442-1.556, 1.556-2.0, 2.0 - infinity
    The six eta region weights for low pt should come first, followed by the high pt weights 
   */
  //uncomment when ready to use

  double weights[12] = {0.00026, 0.00021, 0.00106, 0.0, 0.00806, 0.00970,
			0.00240, 0.00202, 0.01165, 0.0, 0.02621, 0.03616};

  //load data file - This file is a sum of all data files with opposite-signed leptons
  //TFile* f = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  TFile* f = new TFile("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/ChargeMisID_MuCBTightMiniIsoTight_ElMVATightRC_2016B-H.root");
  //TTree* t = (TTree*)f->Get("ChargeMisID");
  TTree* t = (TTree*)f->Get("tEvts_ssdl");

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
  //muon pt;
  vector<double> *muPts = 0;
  vector<double> *muEtas = 0;
  vector<double> *muPhis = 0;
  vector<int> * muIsTight = 0;
  vector<int> * muIsLoose = 0;
  vector<int> * muCharge =0;

  float lepPts1, lepEtas1, lepPhis1, lepEs1;
  float lepPts2, lepEtas2, lepPhis2, lepEs2;
  int lepFlavor1, lepFlavor2;
  int lepCharge1, lepCharge2;
  int lepTight1, lepTight2;
  int lepLoose1, lepLoose2;
  float assocMass, dilepMass
  //jets
  vector<double>* jetPts = 0;
  //met
  double met = 0;

  //Set branch addresses
  t->SetBranchAddress("Lep1Pt", &lepPts1);
  t->SetBranchAddress("Lep1Eta", &lepEtas1);
  t->SetBranchAddress("Lep1Phi", &lepPhis1);
  t->SetBranchAddress("Lep1E", &lepEs1);
  t->SetBranchAddress("Lep1Charge", &lepCharge1);
  t->SetBranchAddress("Lep1Flavor", &lepFlavor1);
  //t->SetBranchAddress("Lep1Tight", &lepTight1);
  //t->SetBranchAddress("Lep1Loose", &lepLoose1);
  t->SetBranchAddress("Lep2Pt", &lepPts2);
  t->SetBranchAddress("Lep2Eta", &lepEtas2);
  t->SetBranchAddress("Lep2Phi", &lepPhis2);
  t->SetBranchAddress("Lep2E", &lepEs2);
  t->SetBranchAddress("Lep2Charge", &lepCharge2);
  t->SetBranchAddress("Lep2Flavor", &lepFlavor2);
  //t->SetBranchAddress("Lep2Tight", &lepTight2);
  //t->SetBranchAddress("Lep2Loose", &lepLoose2);
  t->SetBranchAddress("MET",&met);
  t->SetBranchAddress("AssocMass",&assocMass); //any two electrons that whose invariant mass is closest to Z mass 
  t->SetBranchAddress("DilepMass",&dilepMass); //invariant mass of two opposite sign leptons

  /*
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
  t->SetBranchAddress("AK5JetPt_DileptonCalc",&jetPts);
  t->SetBranchAddress("elChargeConsistent_DileptonCalc",&elChargeConsistency);
  t->SetBranchAddress("muIsLoose_DileptonCalc",&muIsLoose);
  t->SetBranchAddress("muIsTight_DileptonCalc",&muIsTight);
  t->SetBranchAddress("muCharge_DileptonCalc",&muCharge);
  */
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

    if (assocMass > M_Z - dM && assocMass < M_Z + dM) continue;
    if (dilepMass > M_Z - dM && dilepMass < M_Z + dM) continue;
    

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
    /*
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
    */

    /*
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

*/

    //now check for only two tight leptons
    if(vTightLep.size()!=2) continue; //FIXME: delete - Lep1Tight/Lep2Tight sample bug - not filled 
    //check OPPOSITE sign
    //if(vTightLep.at(0)->charge==vTightLep.at(1)->charge) continue; //sample is all opposite sign - Lep1Charge/Lep2Charge sample bug - not filled

    /*        
    //cut harder on lepton pt
    float temp_pt = 0;
    for(unsigned int i = 0; i < vTightLep.size();i++){
      if(vTightLep.at(i)->pt > temp_pt) temp_pt = vTightLep.at(i)->pt;
    }
    */
    
    //check to make sure hardest lepton is above 100 GeV
    //if(temp_pt<140) continue;
    if(lepPt1<40) continue;

    /*
    float sec_pt = 0;
    unsigned int lep2 = 0;
    for(unsigned int i = 0 ; i < vTightLep.size(); i++){
      if((vTightLep.at(i)->pt > sec_pt) && (vTightLep.at(i)->pt!=temp_pt)) {
	lep2=i;
	sec_pt = vTightLep.at(i)->pt;
      }
    }
    */
    //if(sec_pt < 100) continue;
    if(lepPt2<30) continue;

    //FIXME: dzou - not sure what this SetPtEtaPhiM does - need to read method
    /*
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
*/
    //cut on dilepton mass
    //if( (v1+v2).M() <250) continue;
    if (dilepMass < 250) continue;

    //check channel
    bool ee = false;
    bool emu = false;
    bool mumu = false;

    if( (lepFlavor1 == 0) && (lepFlavor2 == 0) ) ee = true;
    if( (lepFlavor1 == 0) && (lepFlavor2 == 1) ) emu = true;
    if( (lepFlavor1 == 1) && (lepFlavor2 == 0) ) emu = true;
    if( (lepFlavor1 == 1) && (lepFlavor2 == 1) ) mumu = true;
    /*
    if(vTightLep.at(0)->isEl && vTightLep.at(1)->isEl) ee = true;
    if(vTightLep.at(0)->isEl && vTightLep.at(1)->isMu) emu = true;
    if(vTightLep.at(0)->isMu && vTightLep.at(1)->isEl) emu = true;
    if(vTightLep.at(0)->isMu && vTightLep.at(1)->isMu) mumu = true;
    */
    //skip any di-mu events since we assume they don't contribute to ss via charge misID
    if(mumu) continue;

    /*now we want to fill our predicted HT distribution with the os events weighted by the appropriate weights
      I've written a simple function that gives the correct element of the weight array if you pass it the eta of the electron.
      Usage is simply EtaWeight(weightsarray, eta). Take advantage or it to make life easier :)
     */

    // ADD CODE HERE TO CALCULATE PROBABLITY FOR EACH EVENT

    float ee_weight = EtaWeight(weights, lepEta1, lepPt1) + EtaWeight(weights, lepEta2, lepPt2) - EtaWeight(weights, lepEta1, lepPt1)*EtaWeight(weights, lepEta2, lepPt2);
    float emu_weight;
    if(lepFlavor1 == 0)
      emu_weight = EtaWeight(weights, lepEta1, lepPt1); //misId rate for mu is zero
    else
      emu_weight = EtaWeight(weights, lepEta2, lepPt2);

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
  /* //sample does not have Same sign events
  TCanvas c2; 
  ssHTHist->Draw();
  c2.Print("HT_sameSignEvents.pdf");
  */
  //finally let's save our predictions to a root file for ease of use later
  T->Branch("osHTHist","TH1F",&osHTHist,32000,0);
  //T->Branch("ssHTHist","TH1F",&ssHTHist,32000,0);
  T->Fill();
  T->Print();
  fbg->Write();
  fbg->Close();


}


