#include <iostream>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "ObjectID.C"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z


void Ex_2p2(){

  /*
    The first part of this script is just a lot of logistics: initializing variables, getting information from the ttree, etc. To complete the exercise you don't need to understand all of it,
    but please feel free to use it as a reference for how to get information from a ttree and loop over all the entries.
   */

  //load the files
  TFile* f = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  // TFile* f = new TFile("/uscms_data/d2/avetisya/CMSDAS/2015/v3/CMSSW_7_3_0/src/LJMet/Com/python/ljmet_tree_Status23and33.root");
  TTree* t = (TTree*)f->Get("ljmet");

  //initialize need histograms
  TH1F* numPtHist_e = new TH1F("numPtHist_e","Lepton p_{T} - TIGHT ID",80,0.,40.);
  TH1F* denPtHist_e = new TH1F("denPtHist_e","Lepton p_{T} - LOOSE ID",80,0.,40.);
  TH1F* numEtaHist_e = new TH1F("numEtaHist_e","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist_e = new TH1F("denEtaHist_e","Lepton #eta - LOOSE ID",30,-3.,3.);

  TH1F* numPtHist_m = new TH1F("numPtHist_m","Lepton p_{T} - TIGHT ID",80,0.,40.);
  TH1F* denPtHist_m = new TH1F("denPtHist_m","Lepton p_{T} - LOOSE ID",80,0.,40.);
  TH1F* numEtaHist_m = new TH1F("numEtaHist_m","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist_m = new TH1F("denEtaHist_m","Lepton #eta - LOOSE ID",30,-3.,3.);

  
  int nEntries = t->GetEntries();
  // int nEntries = 20;

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

  //Muon varaibles;
  vector<double> *muPts = 0;
  vector<double> *muEtas = 0;
  vector<double> *muPhis = 0;
  vector<int> * muIsTight = 0;
  vector<int> * muIsLoose = 0;

  vector<int>    *global = 0;
  vector<double> *chi2 = 0;
  vector<int>    *nValMuHits = 0;
  vector<int>    *nMatchedStations = 0;
  vector<double> *dxy = 0;
  vector<double> *dz = 0;
  vector<int>    *nValPixelHits = 0;
  vector<int>    *nTrackerLayers = 0;
  vector<double> *relIso = 0;

  //Jets
  vector<double> *jetPts = 0;
  vector<double> *jetEtas = 0;
  vector<double> *jetPhis = 0;
  vector<double> *jetEnergies = 0;


  //set branch addresses
  //Electrons
  t->SetBranchAddress("elPt_DileptonCalc", &elPts);
  t->SetBranchAddress("elEta_DileptonCalc", &elEtas);
  t->SetBranchAddress("elPhi_DileptonCalc", &elPhis);
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
  t->SetBranchAddress("elCharge_DileptonCalc",&elCharge);

  //Muons
  t->SetBranchAddress("muPt_DileptonCalc",&muPts);
  t->SetBranchAddress("muEta_DileptonCalc",&muEtas);
  t->SetBranchAddress("muPhi_DileptonCalc",&muPhis);
  t->SetBranchAddress("muIsLoose_DileptonCalc",&muIsLoose);
  t->SetBranchAddress("muIsTight_DileptonCalc",&muIsTight);

  t->SetBranchAddress("muGlobal_DileptonCalc",&global);
  t->SetBranchAddress("muChi2_DileptonCalc",&chi2);
  t->SetBranchAddress("muNValMuHits_DileptonCalc",&nValMuHits);
  t->SetBranchAddress("muNMatchedStations_DileptonCalc",&nMatchedStations);
  t->SetBranchAddress("muDxy_DileptonCalc",&dxy);
  t->SetBranchAddress("muDz_DileptonCalc",&dz);
  t->SetBranchAddress("muNValPixelHits_DileptonCalc",&nValPixelHits);
  t->SetBranchAddress("muNTrackerLayers_DileptonCalc",&nTrackerLayers);
  t->SetBranchAddress("muRelIso_DileptonCalc",&relIso);


  //Jets
  t->SetBranchAddress("AK5JetPt_DileptonCalc", &jetPts);
  t->SetBranchAddress("AK5JetEta_DileptonCalc", &jetEtas);
  t->SetBranchAddress("AK5JetPhi_DileptonCalc", &jetPhis);
  t->SetBranchAddress("AK5JetEnergy_DileptonCalc", &jetEnergies);

  //MET cut to avoid leptons from W bosons
  double met     = 0;
  double met_phi = 0;
  t->SetBranchAddress("met_DileptonCalc", &met);
  t->SetBranchAddress("met_phi_DileptonCalc", &met_phi);

  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    //checks per event
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;

    //MET cut to avoid W's
    if (met > 25) continue;
    
    /*
      ELECTRONS
     */

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

      el->isEl    = true;
      el->isMu    = false;
      
      if (!el->loose()) continue;
      el->Loose = el->loose();
      el->Tight = el->tight();
      
      vEl.push_back(el);

      
    }


    /*
      MUONS
     */
    vector<Muon*> vMu;
    for (unsigned int uiMu = 0; uiMu <muPts->size(); uiMu++){
      Muon* mu = new Muon;

      mu->pt      = muPts->at(uiMu);
      mu->eta     = muEtas->at(uiMu);
      mu->phi     = muPhis->at(uiMu);
      mu->isLoose = muIsLoose->at(uiMu);
      mu->isTight = muIsTight->at(uiMu);
      mu->isEl    = false;
      mu->isMu    = true;

      mu->global           = global->at(uiMu);
      mu->chi2             = chi2->at(uiMu);
      mu->nValMuHits       = nValMuHits->at(uiMu);
      mu->nMatchedStations = nMatchedStations->at(uiMu);
      mu->dxy              = dxy->at(uiMu);
      mu->dz               = dz->at(uiMu);
      mu->nValPixelHits    = nValPixelHits->at(uiMu);
      mu->nTrackerLayers   = nTrackerLayers->at(uiMu);
      mu->relIso           = relIso->at(uiMu);

      // mu->Loose = mu->loose();
      // mu->Tight = mu->tight();

      mu->Loose = mu->cutBasedLoose();
      mu->Tight = mu->cutBasedTight();
      
      if (!mu->Loose) continue;
      
      vMu.push_back(mu);
    }

    //Reject events with more than 1 loose lepton to avoid tops and Zs
    if (vMu.size() + vEl.size() != 1) continue;
    
    Lepton *lepton = 0;
    if (vEl.size() > 0) lepton = vEl.at(0);
    if (vMu.size() > 0) lepton = vMu.at(0);
    
    //Reject high pt leptons
    if (lepton->pt > 35) continue;

    //Reject leptons with high transverse mass to avoid W's
    double dphi = fabs(lepton->phi - met_phi);
    if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
    double mt2 = 2 * lepton->pt * met * (1 - cos(dphi));
    if (mt2 > 25) continue;
 
    //Away-jet and lepton-jet Z veto
    bool foundAwayJet   = false;
    bool leptonJetZVeto = false;

    TLorentzVector v1, v2;
    if (lepton->isMu) v1.SetPtEtaPhiM(lepton->pt, lepton->eta, lepton->phi, M_MU);
    else              v1.SetPtEtaPhiM(lepton->pt, lepton->eta, lepton->phi, M_EL);

    for (unsigned int uj = 0; uj < jetPts->size(); uj++){
      v2.SetPtEtaPhiE(jetPts->at(uj), jetEtas->at(uj), jetPhis->at(uj), jetEnergies->at(uj));
      
      if ((v1+v2).M() < M_Z + dM && (v1+v2).M() > M_Z - dM){
	leptonJetZVeto = true;
	break;
      }
      if (v2.Pt() > 40 && v1.DeltaR(v2) > 1) foundAwayJet = true;
    }
    
    if (leptonJetZVeto || !foundAwayJet) continue;

    if (lepton->isEl){
      denPtHist_e->Fill(v1.Pt());
      denEtaHist_e->Fill(v1.Eta());
      
      if (lepton->Tight){
	numPtHist_e->Fill(v1.Pt());
	numEtaHist_e->Fill(v1.Eta());
      }
    }

    if (lepton->isMu){
      denPtHist_m->Fill(v1.Pt());
      denEtaHist_m->Fill(v1.Eta());
      
      if (lepton->Tight){
	numPtHist_m->Fill(v1.Pt());
	numEtaHist_m->Fill(v1.Eta());
      }
    }
    
  }
  
  //masshist->Draw();  
  TGraphAsymmErrors* ptgraph_e = new TGraphAsymmErrors(numPtHist_e,denPtHist_e);
  TGraphAsymmErrors* etagraph_e = new TGraphAsymmErrors(numEtaHist_e,denEtaHist_e);
  
  TCanvas c1;
  ptgraph_e->Draw("apl");
  c1.Print("PromptRate_v_pT_e.pdf");
  
  TCanvas c2;
  etagraph_e->Draw("apl");
  c2.Print("PromptRate_v_Eta_e.pdf");
  
  TGraphAsymmErrors* ptgraph_m = new TGraphAsymmErrors(numPtHist_m,denPtHist_m);
  TGraphAsymmErrors* etagraph_m = new TGraphAsymmErrors(numEtaHist_m,denEtaHist_m);
  
  TCanvas c3;
  ptgraph_m->Draw("apl");
  c3.Print("PromptRate_v_pT_m.pdf");
  
  TCanvas c4;
  etagraph_m->Draw("apl");
  c4.Print("PromptRate_v_Eta_m.pdf");
  
  int nLooseMu = denEtaHist_m->GetEntries();
  int nTightMu = numEtaHist_m->GetEntries();
  double muFR  = (double)nTightMu;
  if (nLooseMu != 0) muFR /= nLooseMu;
  
  int nLooseEl = denEtaHist_e->GetEntries();
  int nTightEl = numEtaHist_e->GetEntries();
  double elFR  = (double)nTightEl;
  if (nLooseEl != 0) elFR /= nLooseEl;
  
  std::cout<<"Loose Electrons: "<<nLooseEl<<std::endl;
  std::cout<<"Tight Electrons: "<<nTightEl<<std::endl;
  std::cout<<"Electron fake rate: "<<elFR<<std::endl<<std::endl;
  
  std::cout<<"Loose Muons: "<<nLooseMu<<std::endl;
  std::cout<<"Tight Muons: "<<nTightMu<<std::endl;
  std::cout<<"Muon fake rate: "<<muFR<<std::endl<<std::endl;
}
