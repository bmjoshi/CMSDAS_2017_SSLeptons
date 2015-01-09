#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"

void Ex_2p3(){

  //load the 'data' and mc
  TFile* fdata = new TFile("ljmet_tree_TT1.root");
  TFile* fmc   = new TFile("ljmet_tree_TT2.root");

  TTree* tdata = fdata->Get("ljmet");
  TTree* tmc   = fmc->Get("ljmet");


  // LOAD THINGS FOR DATA

  //initialize variables
  int nEntriesdata = tdata->GetEntries();
  int nEntriesmc   = tmc->GetEntries();

  vector<double> *diElMass_data = 0;
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

  //set branch addresses
  tdata->SetBranchAddress("diElMass_DileptonCalc", &diElMass_data);
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


  //bools for data ids
  //tracking cuts
  bool elDEtaLoosecut1data = false; // fabs() <0.181;
  bool elDEtaLoosecut2data = false; // fabs() <0.181;
  bool elDEtaTightcut1data = false; // fabs() < 0.0091;
  bool elDEtaTightcut2data = false; //fabs() < 0.0091;
  bool elDPhiLoosecut1data = false; // fabs() <0.181;
  bool elDPhiLoosecut2data = false; // fabs() <0.181;
  bool elDPhiTightcut1data = false; // fabs() < 0.0091;
  bool elDPhiTightcut2data = false; //fabs() < 0.0091;
  //shower cuts
  bool elSigmaIetaIetaLoosecut1data = false;
  bool elSigmaIetaIetaLoosecut2data = false;
  bool elSigmaIetaIetaTightcut1data = false;
  bool elSigmaIetaIetaTightcut2data = false;
  //H/E cuts
  bool elHELoosecut1data = false;
  bool elHELoosecut2data = false;
  bool elHETightcut1data = false;
  bool elHETightcut2data = false;
  //primary vertex cuts
  bool elD0Loosecut1data = false;
  bool elD0Loosecut2data = false;
  bool elD0Tightcut1data = false;
  bool elD0Tightcut2data = false;
  bool elDZLoosecut1data = false;
  bool elDZLoosecut2data = false;
  bool elDZTightcut1data = false;
  bool elDZTightcut2data = false;
  //ooEooP cuts
  bool elooEmooPLoosecut1data = false;
  bool elooEmooPLoosecut2data = false;
  bool elooEmooPTightcut1data = false;
  bool elooEmooPTightcut2data = false;
  //isolation cuts
  bool elPFIsoLoosecut1data = false;
  bool elPFIsoLoosecut2data = false;
  bool elPFIsoTightcut1data = false;
  bool elPFIsoTightcut2data = false;
  //conversion rejection cuts
  bool elVtxFitProbLoosecut1data = false;
  bool elVtxFitProbLoosecut2data = false;
  bool elVtxFitProbTightcut1data = false;
  bool elVtxFitProbTightcut2data = false;
  bool elMissingHitsLoosecut1data = false;
  bool elMissingHitsLoosecut2data = false;
  bool elMissingHitsTightcut1data = false;
  bool elMissingHitsTightcut2data = false;

  //booleans for final ids
  bool elLoose1 = false; //and of all loose cuts
  bool elTight1 = false; //and of all tight cuts
  bool elLoose2 = false; //and of all loose cuts
  bool elTight2 = false; //and of all tight cuts

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
  tmc->SetBranchAddress("diElMass_DileptonCalc", &diElMass_mc);
  tmc->SetBranchAddress("elPt_DileptonCalc", &elpts_mc);
  tmc->SetBranchAddress("elEta_DileptonCalc", &elEtas_mc);
  tmc->SetBranchAddress("elDeta_DileptonCalc", &elDeta_mc);
  tmc->SetBranchAddress("elDphi_DileptonCalc", &elDphi_mc);
  tmc->SetBranchAddress("elDZ_DileptonCalc", &elDZs_mc);
  tmc->SetBranchAddress("elD0_DileptonCalc", &elD0s_mc);
  tmc->SetBranchAddress("elHoE_DileptonCalc",&elHoverEs_mc);
  tmc->SetBranchAddress("elMHits_DileptonCalc",&elMHits_mc);
  tmc->SetBranchAddress("elOoemoop_DileptonCalc",&elOoEmooPs_mc);
  tmc->SetBranchAddress("elRelIso_DileptonCalc",&elRelIsos_mc);
  tmc->SetBranchAddress("elSihih_DileptonCalc",&elSigmaIetaIetas_mc);


  /*
    DATA PROCEDURE: WE JUST DO THE ANALYSIS BUT GET NUMBERS FOR HOW MANY LOOSE
    AND HOW MANY TIGHT LEPTONS THERE ARE
   */

  //event loop - plots to produce: HT distribution for all cuts except 
  for(int iel = 0; iel < nEntriesdata; iel++){


    //check dilepton cut
    if(elPts->size() < 2) continue;
    //z mass veto
    if( (diElMass->at(0) > 76) || (diElMass->at(0) <116) ) continue;
    //check HT req
    if( HT_data->at(0) < 1400) continue;
    //check met req
    if( met_data->at(0) < 100) continue;
    //check lep1pt req
    if(elPts_data->at(0) < 80) continue;
    //check subleading lep pt req
    if(elPts_data->at(1) < 30) continue;
    //require more than one jet
    if(jetPts_data->size() < 2) continue;
    //check for high pt jet
    if(jetPts_data->at(0) < 150) continue;
    //check for second jet
    if(jetPts_data->at(1) < 50) continue;
    
    //check samesign lepton req
    if(elCharge1 != elCharge2) continue;
    
    //now that we have made the objects pass all the analysis cuts other than ID we can see how many events there are with 0, 1, and 2 tight leptons
    
    //loop for id
    //need to separate between barrel and endcap
    //first lepton
    if(elEtas->at(0) <= 1.479){
      //check delta Eta between track and supercluster
      if(fabs(elDeta->at(0)) < 0.0181) elDEtaLoosecut1data = true;
      if(fabs(elDeta->at(0)) < 0.0091) elDEtaTightcut1data = true;
      //check delta phi between track and supercluster
      if(fabs(elDphi->at(0)) < 0.0936) elDPhiLoosecut1data = true;
      if(fabs(elDphi->at(0)) < 0.031) elDPhiTightcut1data = true;
      //check sigmaIetaIeta
      if(elSigmaIetaIetas->at(0) < 0.0123) elSigmaIetaIetaLoosecut1data = true;
      if(elSigmaIetaIetas->at(0) < 0.0106) elSigmaIetaIetaTightcut1data = true;
      //check H over E
      if(elHoverEs->at(0) < 0.141) elHELoosecut1data = true;
      if(elHoverEs->at(0) < 0.0532) elHETightcut1data = true;
      //check vertex cuts
      if(fabs(elD0s->at(0)) < 0.0166) elD0Loosecut1data = true;
      if(fabs(elD0s->at(0)) < 0.0126) elD0Tightcut1data = true;
      if(fabs(elDZs->at(0)) < 0.54342) elDZLoosecut1data = true;
      if(fabs(elDZs->at(0)) < 0.0116) elDZTightcut1data = true;
      //ooEmooP cuts
      if(fabs(elOoEmooPs->at(0)) < 0.1353) elooEmooPLoosecut1data = true;
      if(fabs(elOoEmooPs->at(0)) < 0.0609) elooEmooPTightcut1data = true;
      //pfIsoCuts
      if(elRelIsos->at(0) < 0.24) elPFIsoLoosecut1data = true;
      if(elRelIsos->at(0) < 0.1649) elPFIsoTightcut1data = true;
      //missing hits cuts
      if(elMHits->at(0) <= 1) elMissingHitsLoosecut1data = true;
      if(elMHits->at(0) <= 1) elMissingHitsTightcut1data = true;
    }
    else{
      //check delta Eta between track and supercluster
      if(fabs(elDeta->at(0)) < 0.0124) elDEtaLoosecut1data = true;
      if(fabs(elDeta->at(0)) < 0.0106) elDEtaTightcut1data = true;
      //check delta phi between track and supercluster
      if(fabs(elDphi->at(0)) < 0.0642) elDPhiLoosecut1data = true;
      if(fabs(elDphi->at(0)) < 0.0359) elDPhiTightcut1data = true;
      //check sigmaIetaIeta
      if(elSigmaIetaIetas->at(0) < 0.035) elSigmaIetaIetaLoosecut1data = true;
      if(elSigmaIetaIetas->at(0) < 0.0305) elSigmaIetaIetaTightcut1data = true;
      //check H over E
      if(elHoverEs->at(0) < 0.1115) elHELoosecut1data = true;
      if(elHoverEs->at(0) < 0.0835) elHETightcut1data = true;
      //check vertex cuts
      if(fabs(elD0s->at(0)) < 0.098) elD0Loosecut1data = true;
      if(fabs(elD0s->at(0)) < 0.0163) elD0Tightcut1data = true;
      if(fabs(elDZs->at(0)) < 0.9187) elDZLoosecut1data = true;
      if(fabs(elDZs->at(0)) < 0.5999) elDZTightcut1data = true;
      //ooEmooP cuts
      if(fabs(elOoEmooPs->at(0)) < 0.1443) elooEmooPLoosecut1data = true;
      if(fabs(elOoEmooPs->at(0)) < 0.1126) elooEmooPTightcut1data = true;
      //pfIsoCuts
      if(elRelIsos->at(0) < 0.3529) elPFIsoLoosecut1data = true;
      if(elRelIsos->at(0) < 0.2075) elPFIsoTightcut1data = true;
      //missing hits cuts
      if(elMHits->at(0) <= 1) elMissingHitsLoosecut1data = true;
      if(elMHits->at(0) <= 1) elMissingHitsTightcut1data = true;
    }
    
    //now for second lepton
    if(elEtas->at(1)<=1.479){
      //check delta Eta between track and supercluster	
      if(fabs(elDeta->at(1)) < 0.0181) elDEtaLoosecut2data = true;
      if(fabs(elDeta->at(1)) < 0.0091) elDEtaTightcut2data = true;
      //check delta phi between track and supercluster
      if(fabs(elDphi->at(1)) < 0.0936) elDPhiLoosecut2data = true;
      if(fabs(elDphi->at(1)) < 0.031) elDPhiTightcut2data = true;
      //check sigmaIetaIeta
      if(elSigmaIetaIetas->at(1) < 0.0123) elSigmaIetaIetaLoosecut1data = true;
      if(elSigmaIetaIetas->at(1) < 0.0106) elSigmaIetaIetaTightcut1data = true;
      //check H over E
      if(elHoverEs->at(1) < 0.141) elHELoosecut2data = true;
      if(elHoverEs->at(1) < 0.0532) elHETightcut2data = true;
      //check vertex cuts
      if(fabs(elD0s->at(1)) < 0.0166) elD0Loosecut2data = true;
      if(fabs(elD0s->at(1)) < 0.0126) elD0Tightcut2data = true;
      if(fabs(elDZs->at(1)) < 0.54342) elDZLoosecut2data = true;
      if(fabs(elDZs->at(1)) < 0.0116) elDZTightcut2data = true;
      //ooEmooP cuts
      if(fabs(elOoEmooPs->at(1)) < 0.1353) elooEmooPLoosecut2data = true;
      if(fabs(elOoEmooPs->at(1)) < 0.0609) elooEmooPTightcut2data = true;
      //pfIsoCuts
      if(elRelIsos->at(1) < 0.24) elPFIsoLoosecut2data = true;
      if(elRelIsos->at(1) < 0.1649) elPFIsoTightcut2data = true;
      //missing hits cuts
      if(elMHits->at(1) <= 1) elMissingHitsLoosecut2data = true;
      if(elMHits->at(1) <= 1) elMissingHitsTightcut2data = true;
    }
    
    else{
      //check delta Eta between track and supercluster	
      if(fabs(elDeta->at(1)) < 0.0124) elDEtaLoosecut2data = true;
      if(fabs(elDeta->at(1)) < 0.0106) elDEtaTightcut2data = true;
      //check delta phi between track and supercluster
      if(fabs(elDphi->at(1)) < 0.0642) elDPhiLoosecut2data = true;
      if(fabs(elDphi->at(1)) < 0.0359) elDPhiTightcut2data = true;
      //check sigmaIetaIeta
      if(elSigmaIetaIetas->at(1) < 0.035) elSigmaIetaIetaLoosecut2data = true;
      if(elSigmaIetaIetas->at(1) < 0.0305) elSigmaIetaIetaTightcut2data = true;
      //check H over E
      if(elHoverEs->at(1) < 0.1115) elHELoosecut2data = true;
      if(elHoverEs->at(1) < 0.0835) elHETightcut2data = true;
      //check vertex cuts
      if(fabs(elD0s->at(1)) < 0.098) elD0Loosecut2data = true;
      if(fabs(elD0s->at(1)) < 0.0163) elD0Tightcut2data = true;
      if(fabs(elDZs->at(1)) < 0.9187) elDZLoosecut2data = true;
      if(fabs(elDZs->at(1)) < 0.5999) elDZTightcut2data = true;
      //ooEmooP cuts
      if(fabs(elOoEmooPs->at(1)) < 0.1443) elooEmooPLoosecut2data = true;
      if(fabs(elOoEmooPs->at(1)) < 0.1126) elooEmooPTightcut2data = true;
      //pfIsoCuts
      if(elRelIsos->at(1) < 0.3529) elPFIsoLoosecut2data = true;
      if(elRelIsos->at(1) < 0.2075) elPFIsoTightcut2data = true;
      //missing hits cuts
      if(elMHits->at(1) <= 1) elMissingHitsLoosecut2data = true;
      if(elMHits->at(1) <= 1) elMissingHitsTightcut2data = true;
    }
    
    //now add them all up
    elTight1 = (elDPhiTightcut1data * elDEtaTightcut1data * elSigmaIetaIetaTightcut1data * elHETightcut1data * elD0Tightcut1data * elDZTightcut1data * elooEmooPTightcut1data * elPFIsoTightcut1data * elMissingHitsTightcut1data);
    elTight2 = (elDPhiTightcut2data * elDEtaTightcut2data * elSigmaIetaIetaTightcut2data * elHETightcut2data * elD0Tightcut2data * elDZTightcut2data * elooEmooPTightcut2data * elPFIsoTightcut2data * elMissingHitsTightcut2data);
    elLoose1 = (elDPhiLoosecut1data * elDEtaLoosecut1data * elSigmaIetaIetaLoosecut1data * elHELoosecut1data * elD0Loosecut1data * elDZLoosecut1data * elooEmooPLoosecut1data * elPFIsoLoosecut1data * elMissingHitsLoosecut1data);
    elLoose2 = (elDPhiLoosecut2data * elDEtaLoosecut2data * elSigmaIetaIetaLoosecut2data * elHELoosecut2data * elD0Loosecut2data * elDZLoosecut2data * elooEmooPLoosecut2data * elPFIsoLoosecut2data * elMissingHitsLoosecut2data);
  }


}
