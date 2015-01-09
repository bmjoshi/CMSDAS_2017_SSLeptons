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
  


}
