#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"



void Ex_2p1(){

  TFile* fDY = new TFile("ljmet_tree_DY_v2.root");
  TTree* tDY = fDY->Get("ljmet");
  //  TTreeReader myReader("ljmet", f);
  
  TH1F* masshist = new TH1F("masshist", "Dielectron Invariant Mass",100,0.,200.);
  int nEntries = tDY->GetEntries();
  vector<double> *diElMass_DileptonCalc = 0;
  vector<double> *elpts = 0;
  vector<double> *elDEta =0;

  tDY->SetBranchAddress("diElMass_DileptonCalc", &diElMass_DileptonCalc);
  tDY->SetBranchAddress("elPt_DileptonCalc", &elpts);
  tDY->SetBranchAddress("elDeta_DileptonCalc", &elDeta);

  //setup booleans for passing the various cuts and initialize to false:
  bool masscut = false; // 81 < M_ee <101
  bool elptcut1 = false; // >20
  bool eptcut2 = false; // >20
  bool elDEtaLoosecut1 = false; // fabs() <0.181;
  bool elDEtaLoosecut2 = false; // fabs() <0.181;
  bool elDEtaTightcut1 = false; // fabs() < 0.0091;
  bool elDEtaTightcut2 = false; //fabs() < 0.0091;

  for(int ient = 0; ient < nEntries; ient++){
    tDY->GetEntry(ient);
    // if(diElMass_DileptonCalc->size()>0) std::cout<<" mass: "<<diElMass_DileptonCalc->at(0)<<std::endl;
    //checks per event

    if( diElMass_DileptonCalc->size()>0 && ( diElMass_DileptonCalc->at(0) > 81) && ( diElMass_DileptonCalc->at(0) <101) ){ 
      //std::cout<<"diElMass is: "<<diElMass_DileptonCalc->at(0)<<std::endl;
      //      masshist->Fill(diElMass_DileptonCalc->at(0));
      masscut=true;
    }


    //checks per lepton:

    //std::cout<<" size of electron pt vector: "<<elpts->size()<<std::endl;
    if(elpts->size()>1){
      if(elpts->at(0)) elptcut1 = true;
      if(elpts->at(1)) elptcut2 = true;
      if(fabs(elDeta->at(0)) < 0.0181) elDEtaLoosecut1 = true;
      if(fabs(elDeta->at(1)) < 0.0181) elDEtaLoosecut2 = true;
      if(fabs(elDeta->at(0)) < 0.0091) elDEtaTightcut1 = true;
      if(fabs(elDeta->at(1)) < 0.0091) elDEtaTightcut2 = true;
      

    }

  }
  masshist->Draw();
  //while (myReader.Next()){
  //std::cout<<"dEta is: "<<



}
