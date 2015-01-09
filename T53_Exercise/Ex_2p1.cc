#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>



void Ex_2p1(){

  /*
    The first part of this script is just a lot of logistics: initializing variables, getting information from the ttree, etc. To complete the exercise you don't need to understand all of it,
    but please feel free to use it as a reference for how to get information from a ttree and loop over all the entries.
   */

  //load the files
  TFile* fDY = new TFile("ljmet_tree_DY_v2.root");
  TTree* tDY = fDY->Get("ljmet");
  //  TTreeReader myReader("ljmet", f);

  //initialize need histograms
  TH1F* masshist = new TH1F("masshist", "Dielectron Invariant Mass",100,0.,200.);
  TH1F* numPtHist = new TH1F("numPtHist","Lepton p_{T} - TIGHT ID",200,0.,400.);
  TH1F* denPtHist = new TH1F("denPtHist","Lepton p_{T} - LOOSE ID",200,0.,400.);
  TH1F* numEtaHist = new TH1F("numEtaHist","Lepton #eta - TIGHT ID",30,-3.,3.);
  TH1F* denEtaHist = new TH1F("denEtaHist","Lepton #eta - LOOSE ID",30,-3.,3.);

  
  int nEntries = tDY->GetEntries();
  vector<double> *diElMass = 0;
  //kinematic variables
  vector<double> *elpts = 0;
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

  //set branch addresses
  tDY->SetBranchAddress("diElMass_DileptonCalc", &diElMass);
  tDY->SetBranchAddress("elPt_DileptonCalc", &elpts);
  tDY->SetBranchAddress("elEta_DileptonCalc", &elEtas);
  tDY->SetBranchAddress("elDeta_DileptonCalc", &elDeta);
  tDY->SetBranchAddress("elDphi_DileptonCalc", &elDphi);
  tDY->SetBranchAddress("elDZ_DileptonCalc", &elDZs);
  tDY->SetBranchAddress("elD0_DileptonCalc", &elD0s);
  tDY->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs);
  tDY->SetBranchAddress("elMHits_DileptonCalc",&elMHits);
  tDY->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs);
  tDY->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos);
  tDY->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas);

  //setup booleans for passing the various cuts and initialize to false:
  bool masscut = false; // 81 < M_ee <101
  bool elptcut1 = false; // >20
  bool elptcut2 = false; // >20
  //tracking cuts
  bool elDEtaLoosecut1 = false; // fabs() <0.181;
  bool elDEtaLoosecut2 = false; // fabs() <0.181;
  bool elDEtaTightcut1 = false; // fabs() < 0.0091;
  bool elDEtaTightcut2 = false; //fabs() < 0.0091;
  bool elDPhiLoosecut1 = false; // fabs() <0.181;
  bool elDPhiLoosecut2 = false; // fabs() <0.181;
  bool elDPhiTightcut1 = false; // fabs() < 0.0091;
  bool elDPhiTightcut2 = false; //fabs() < 0.0091;
  //shower cuts
  bool elSigmaIetaIetaLoosecut1 = false;
  bool elSigmaIetaIetaLoosecut2 = false;
  bool elSigmaIetaIetaTightcut1 = false;
  bool elSigmaIetaIetaTightcut2 = false;
  //H/E cuts
  bool elHELoosecut1 = false;
  bool elHELoosecut2 = false;
  bool elHETightcut1 = false;
  bool elHETightcut2 = false;
  //primary vertex cuts
  bool elD0Loosecut1 = false;
  bool elD0Loosecut2 = false;
  bool elD0Tightcut1 = false;
  bool elD0Tightcut2 = false;
  bool elDZLoosecut1 = false;
  bool elDZLoosecut2 = false;
  bool elDZTightcut1 = false;
  bool elDZTightcut2 = false;
  //ooEooP cuts
  bool elooEmooPLoosecut1 = false;
  bool elooEmooPLoosecut2 = false;
  bool elooEmooPTightcut1 = false;
  bool elooEmooPTightcut2 = false;
  //isolation cuts
  bool elPFIsoLoosecut1 = false;
  bool elPFIsoLoosecut2 = false;
  bool elPFIsoTightcut1 = false;
  bool elPFIsoTightcut2 = false;
  //conversion rejection cuts
  bool elVtxFitProbLoosecut1 = false;
  bool elVtxFitProbLoosecut2 = false;
  bool elVtxFitProbTightcut1 = false;
  bool elVtxFitProbTightcut2 = false;
  bool elMissingHitsLoosecut1 = false;
  bool elMissingHitsLoosecut2 = false;
  bool elMissingHitsTightcut1 = false;
  bool elMissingHitsTightcut2 = false;

  //booleans for final ids
  bool elLoose1 = false; //and of all loose cuts
  bool elTight1 = false; //and of all tight cuts
  bool elLoose2 = false; //and of all loose cuts
  bool elTight2 = false; //and of all tight cuts

  for(int ient = 0; ient < nEntries; ient++){
    tDY->GetEntry(ient);
    // if(diElMass->size()>0) std::cout<<" mass: "<<diElMass->at(0)<<std::endl;
    //checks per event

    if( diElMass->size()>0 && ( diElMass->at(0) > 81) && ( diElMass->at(0) <101) ){ 
      //std::cout<<"diElMass is: "<<diElMass->at(0)<<std::endl;
      //      masshist->Fill(diElMass->at(0));
      masscut=true;
    }


    //checks per lepton:

    //std::cout<<" size of electron pt vector: "<<elpts->size()<<std::endl;
    if(elpts->size()>1){
      if(elpts->at(0)) elptcut1 = true;
      if(elpts->at(1)) elptcut2 = true;

      //need to separate between barrel and endcap
      //first lepton
      if(elEtas->at(0) <= 1.479){
	//check delta Eta between track and supercluster
	if(fabs(elDeta->at(0)) < 0.0181) elDEtaLoosecut1 = true;
	if(fabs(elDeta->at(0)) < 0.0091) elDEtaTightcut1 = true;
	//check delta phi between track and supercluster
	if(fabs(elDphi->at(0)) < 0.0936) elDPhiLoosecut1 = true;
	if(fabs(elDphi->at(0)) < 0.031) elDPhiTightcut1 = true;
	//check sigmaIetaIeta
	if(elSigmaIetaIetas->at(0) < 0.0123) elSigmaIetaIetaLoosecut1 = true;
	if(elSigmaIetaIetas->at(0) < 0.0106) elSigmaIetaIetaTightcut1 = true;
	//check H over E
	if(elHoverEs->at(0) < 0.141) elHELoosecut1 = true;
	if(elHoverEs->at(0) < 0.0532) elHETightcut1 = true;
	//check vertex cuts
	if(fabs(elD0s->at(0)) < 0.0166) elD0Loosecut1 = true;
	if(fabs(elD0s->at(0)) < 0.0126) elD0Tightcut1 = true;
	if(fabs(elDZs->at(0)) < 0.54342) elDZLoosecut1 = true;
	if(fabs(elDZs->at(0)) < 0.0116) elDZTightcut1 = true;
	//ooEmooP cuts
	if(fabs(elOoEmooPs->at(0)) < 0.1353) elooEmooPLoosecut1 = true;
	if(fabs(elOoEmooPs->at(0)) < 0.0609) elooEmooPTightcut1 = true;
	//pfIsoCuts
	if(elRelIsos->at(0) < 0.24) elPFIsoLoosecut1 = true;
	if(elRelIsos->at(0) < 0.1649) elPFIsoTightcut1 = true;
	//missing hits cuts
	if(elMHits->at(0) <= 1) elMissingHitsLoosecut1 = true;
	if(elMHits->at(0) <= 1) elMissingHitsTightcut1 = true;
      }
      else{
	//check delta Eta between track and supercluster
	if(fabs(elDeta->at(0)) < 0.0124) elDEtaLoosecut1 = true;
	if(fabs(elDeta->at(0)) < 0.0106) elDEtaTightcut1 = true;
	//check delta phi between track and supercluster
	if(fabs(elDphi->at(0)) < 0.0642) elDPhiLoosecut1 = true;
	if(fabs(elDphi->at(0)) < 0.0359) elDPhiTightcut1 = true;
	//check sigmaIetaIeta
	if(elSigmaIetaIetas->at(0) < 0.035) elSigmaIetaIetaLoosecut1 = true;
	if(elSigmaIetaIetas->at(0) < 0.0305) elSigmaIetaIetaTightcut1 = true;
	//check H over E
	if(elHoverEs->at(0) < 0.1115) elHELoosecut1 = true;
	if(elHoverEs->at(0) < 0.0835) elHETightcut1 = true;
	//check vertex cuts
	if(fabs(elD0s->at(0)) < 0.098) elD0Loosecut1 = true;
	if(fabs(elD0s->at(0)) < 0.0163) elD0Tightcut1 = true;
	if(fabs(elDZs->at(0)) < 0.9187) elDZLoosecut1 = true;
	if(fabs(elDZs->at(0)) < 0.5999) elDZTightcut1 = true;
	//ooEmooP cuts
	if(fabs(elOoEmooPs->at(0)) < 0.1443) elooEmooPLoosecut1 = true;
	if(fabs(elOoEmooPs->at(0)) < 0.1126) elooEmooPTightcut1 = true;
	//pfIsoCuts
	if(elRelIsos->at(0) < 0.3529) elPFIsoLoosecut1 = true;
	if(elRelIsos->at(0) < 0.2075) elPFIsoTightcut1 = true;
	//missing hits cuts
	if(elMHits->at(0) <= 1) elMissingHitsLoosecut1 = true;
	if(elMHits->at(0) <= 1) elMissingHitsTightcut1 = true;
      }

      //now for second lepton
      if(elEtas->at(1)<=1.479){
	//check delta Eta between track and supercluster	
	if(fabs(elDeta->at(1)) < 0.0181) elDEtaLoosecut2 = true;
	if(fabs(elDeta->at(1)) < 0.0091) elDEtaTightcut2 = true;
	//check delta phi between track and supercluster
	if(fabs(elDphi->at(1)) < 0.0936) elDPhiLoosecut2 = true;
	if(fabs(elDphi->at(1)) < 0.031) elDPhiTightcut2 = true;
	//check sigmaIetaIeta
	if(elSigmaIetaIetas->at(1) < 0.0123) elSigmaIetaIetaLoosecut1 = true;
	if(elSigmaIetaIetas->at(1) < 0.0106) elSigmaIetaIetaTightcut1 = true;
	//check H over E
	if(elHoverEs->at(1) < 0.141) elHELoosecut2 = true;
	if(elHoverEs->at(1) < 0.0532) elHETightcut2 = true;
	//check vertex cuts
	if(fabs(elD0s->at(1)) < 0.0166) elD0Loosecut2 = true;
	if(fabs(elD0s->at(1)) < 0.0126) elD0Tightcut2 = true;
	if(fabs(elDZs->at(1)) < 0.54342) elDZLoosecut2 = true;
	if(fabs(elDZs->at(1)) < 0.0116) elDZTightcut2 = true;
	//ooEmooP cuts
	if(fabs(elOoEmooPs->at(1)) < 0.1353) elooEmooPLoosecut2 = true;
	if(fabs(elOoEmooPs->at(1)) < 0.0609) elooEmooPTightcut2 = true;
	//pfIsoCuts
	if(elRelIsos->at(1) < 0.24) elPFIsoLoosecut2 = true;
	if(elRelIsos->at(1) < 0.1649) elPFIsoTightcut2 = true;
	//missing hits cuts
	if(elMHits->at(1) <= 1) elMissingHitsLoosecut2 = true;
	if(elMHits->at(1) <= 1) elMissingHitsTightcut2 = true;
      }

      else{
	//check delta Eta between track and supercluster	
	if(fabs(elDeta->at(1)) < 0.0124) elDEtaLoosecut2 = true;
	if(fabs(elDeta->at(1)) < 0.0106) elDEtaTightcut2 = true;
	//check delta phi between track and supercluster
	if(fabs(elDphi->at(1)) < 0.0642) elDPhiLoosecut2 = true;
	if(fabs(elDphi->at(1)) < 0.0359) elDPhiTightcut2 = true;
	//check sigmaIetaIeta
	if(elSigmaIetaIetas->at(1) < 0.035) elSigmaIetaIetaLoosecut2 = true;
	if(elSigmaIetaIetas->at(1) < 0.0305) elSigmaIetaIetaTightcut2 = true;
	//check H over E
	if(elHoverEs->at(1) < 0.1115) elHELoosecut2 = true;
	if(elHoverEs->at(1) < 0.0835) elHETightcut2 = true;
	//check vertex cuts
	if(fabs(elD0s->at(1)) < 0.098) elD0Loosecut2 = true;
	if(fabs(elD0s->at(1)) < 0.0163) elD0Tightcut2 = true;
	if(fabs(elDZs->at(1)) < 0.9187) elDZLoosecut2 = true;
	if(fabs(elDZs->at(1)) < 0.5999) elDZTightcut2 = true;
	//ooEmooP cuts
	if(fabs(elOoEmooPs->at(1)) < 0.1443) elooEmooPLoosecut2 = true;
	if(fabs(elOoEmooPs->at(1)) < 0.1126) elooEmooPTightcut2 = true;
	//pfIsoCuts
	if(elRelIsos->at(1) < 0.3529) elPFIsoLoosecut2 = true;
	if(elRelIsos->at(1) < 0.2075) elPFIsoTightcut2 = true;
	//missing hits cuts
	if(elMHits->at(1) <= 1) elMissingHitsLoosecut2 = true;
	if(elMHits->at(1) <= 1) elMissingHitsTightcut2 = true;
      }

    

      //now add them all up
      elTight1 = (elDPhiTightcut1 * elDEtaTightcut1 * elSigmaIetaIetaTightcut1 * elHETightcut1 * elD0Tightcut1 * elDZTightcut1 * elooEmooPTightcut1 * elPFIsoTightcut1 * elMissingHitsTightcut1);
      elTight2 = (elDPhiTightcut2 * elDEtaTightcut2 * elSigmaIetaIetaTightcut2 * elHETightcut2 * elD0Tightcut2 * elDZTightcut2 * elooEmooPTightcut2 * elPFIsoTightcut2 * elMissingHitsTightcut2);
      elLoose1 = (elDPhiLoosecut1 * elDEtaLoosecut1 * elSigmaIetaIetaLoosecut1 * elHELoosecut1 * elD0Loosecut1 * elDZLoosecut1 * elooEmooPLoosecut1 * elPFIsoLoosecut1 * elMissingHitsLoosecut1);
      elLoose2 = (elDPhiLoosecut2 * elDEtaLoosecut2 * elSigmaIetaIetaLoosecut2 * elHELoosecut2 * elD0Loosecut2 * elDZLoosecut2 * elooEmooPLoosecut2 * elPFIsoLoosecut2 * elMissingHitsLoosecut2);
      
      //check to make sure there is at least one tight lepton, but that the event passes the z mass window cut
      if( (elTight1 || elTight2) && masscut ){
	//then fill denominator with other lepton if it's a loose one
	if(elTight1){
	  if(elLoose2){
	    denPtHist->Fill(elpts->at(1));
	    denEtaHist->Fill(elEtas->at(1));
	    
	    if(elTight2){
	      numPtHist->Fill(elpts->at(1));
	      numEtaHist->Fill(elEtas->at(1));
	    }
	  }
	  
	}
	else{
	  if(elLoose1){
	    denPtHist->Fill(elpts->at(0));
	    denEtaHist->Fill(elEtas->at(0));
	    if(elTight1){
	      numPtHist->Fill(elpts->at(0));
	      numEtaHist->Fill(elEtas->at(0));
	    }
	  }
	  
	}
      }
      //finish requirement of two leptons
    }
    //finish loop over entries
  }

  //masshist->Draw();
  
  TGraphAsymmErrors* ptgraph = new TGraphAsymmErrors(numPtHist,denPtHist);
  TGraphAsymmErrors* etagraph = new TGraphAsymmErrors(numEtaHist,denEtaHist);
  
  TCanvas c1;
  ptgraph->Draw("apl");
  c1.Print("PromptRate_v_pT.pdf");

  TCanvas c2;
  etagraph->Draw("apl");
  c2.Print("PromptRate_v_Eta.pdf");



}
