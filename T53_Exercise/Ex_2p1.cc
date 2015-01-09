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
  TH1F* numEtaHist = new TH1F("numEtaHist","Lepton #eta - TIGHT ID",200,0.,400.);
  TH1F* denEtaHist = new TH1F("denEtaHist","Lepton #eta - LOOSE ID",200,0.,400.);

  
  int nEntries = tDY->GetEntries();
  vector<double> *diElMass = 0;
  //kinematic variables
  vector<double> *elpts = 0;
  vector<double> *elEtas =0;
  //variables for tracking cuts
  vector<double> *eleDeta =0;
  vector<double> *elDphi =0;
  //variables for primary vtx cuts
  vector<double> *elDZs = 0;
  vector<double> *elD0s = 0;
  // H over E
  vector<double> *elHoverEs = 0;
  //missing hits
  vector<double> *elMHits = 0;
  //ooemoop
  vector<double> *elOoEmooPs = 0;
  //charged isolation
  vector<double> *elRelIsos = 0;
  //sigmaIetaIeta
  vector<double>* elSigIetaIetas = 0;

  //set branch addresses
  tDY->SetBranchAddress("diElMass_DileptonCalc", &diElMass);
  tDY->SetBranchAddress("elPt_DileptonCalc", &elpts);
  tDY->SetBranchAddress("elEta_DileptonCalc", &elEtas);
  tDY->SetBranchAddress("elDeta_DileptonCalc", &elDeta);
  tDY->SetBranchAddress("elDphi_DileptonCalc", &elDphi);
  tDY->SetBranchAddress("elDZ_DileptonCalc", &elDZs);
  tDY->SetBranchAddress("elD0_DileptonCalc", &elD0s);
  tDY->SetBranchAddress("elHoE_DileptonCalc",&elD0s);
  tDY->SetBranchAddress("elMHits_DileptonCalc",&elMHits);
  tDY->SetBranchAddress("elOooemoop_DileptonCalc",&elOoEmooPs);
  tDY->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos);
  tDY->SetBranchAddress("elSihih_DileptonCalc",&elSigIetIetas);

  //setup booleans for passing the various cuts and initialize to false:
  bool masscut = false; // 81 < M_ee <101
  bool elptcut1 = false; // >20
  bool eptcut2 = false; // >20
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
  bool elHETighcut2 = false;
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
  bool elooEooPLoosecut1 = false;
  bool elooEooPLoosecut2 = false;
  bool elooEooPTightcut1 = false;
  bool elooEooPTightcut2 = false;
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
      }

    }

    //check to make sure there is at least one tight lepton
    if( elTight1 && elTight2){
      //then fill denominator with other lepton if it's a loose one
      if(elTight1){
	if(elLoose2){
	  denPtHist->Fill(elpts->at(1));
	  denEtaHist->Fill(elEtas->at(1));

	  if(elTight2){
	    numPtHist->Fill(elpts->at(1));
	    numEtaHist0>Fill(elEtas->at(1));
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
    //finish loop over entries
  }

  //masshist->Draw();
  //while (myReader.Next()){
  //std::cout<<"dEta is: "<<



}
