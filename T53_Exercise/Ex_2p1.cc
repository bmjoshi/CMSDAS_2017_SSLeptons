#include <iostream>
#include "TFile.h"
#include "TTree.h"
#include "TH1.h"
#include <vector>
//#include "TTreeReader.h"
//#include "TTreeReaderValue.h"



void Ex_2p1(){

  TFile* fDY = new TFile("ljmet_tree_DY.root");
  TTree* tDY = fDY->Get("ljmet");
  //  TTreeReader myReader("ljmet", f);
  
  TH1F* masshist = new TH1F("masshist", "Dielectron Invariant Mass",100,0.,200.);
  int nEntries = tDY->GetEntries();
  vector<double> *diElMass_DileptonCalc = 0;
  tDY->SetBranchAddress("diElMass_DileptonCalc", &diElMass_DileptonCalc);

  for(int ient = 0; ient < nEntries; ient++){
    tDY->GetEntry(ient);
    //std::cout<<"entry: "<<ient<<" and mass vector"<<diElMass_DileptonCalc->size()<<std::endl;

    if( diElMass_DileptonCalc->size()>0 && ( diElMass_DileptonCalc->at(0) <111) && ( diElMass_DileptonCalc->at(0) >76) ){ 
      //std::cout<<"diElMass is: "<<diElMass_DileptonCalc->at(0)<<std::endl;
      masshist->Fill(diElMass_DileptonCalc->at(0));
    }
  }
  masshist->Draw();
  //while (myReader.Next()){
  //std::cout<<"dEta is: "<<



}
