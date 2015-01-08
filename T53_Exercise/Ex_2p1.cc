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
  
  int nEntries = tDY->GetEntries();
  vector<double> diElMass_DileptonCalc;
  tDY->SetBranchAddress("diElMass_DileptonCalc", &diElMass_DileptonCalc);

  for(int ient = 0; ient < nEntries; ient++){
    tDY->GetEntry(ient);
    if( (diElMass_DileptonCalc.at(ient) <111) && diElMass_DileptonCalc.at(ient) >76){ 
      std::cout<<"diElMass is: "<<diElMass_DileptonCalc.at(ient)<<std::endl;
    }
  }

  //while (myReader.Next()){
  //std::cout<<"dEta is: "<<



}
