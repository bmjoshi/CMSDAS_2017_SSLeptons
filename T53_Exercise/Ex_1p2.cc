#include <iostream>
#include "TH1.h"
#include "TTree.h"


void Ex_1p2(){

  //load TFiles
  TFile* fDY = new TFile("ljmet_DY.root");
  TFile* fWZ = new TFile("ljmet_WZ.root");
  TFile* fWJets = new TFile("ljmet_WJets.root");
  TFile* fTT = new TFile("ljmet_TT.root");
  TFile* fTTZ = new TFile("ljmet_TTZ.root");

  //load the TTrees
  TTree* tDY = fDY->Get("ljmet_tree");
  TTree* tWZ = fWZ->Get("ljmet_tree");
  TTree* tWJets = fWJets->Get("ljmet_tree");
  TTree* tTT = fTT->Get("ljmet_tree");
  TTree* tTTZ = fTTZ->Get("ljmet_tree");

  //define our taget luminosity - we are picking 5 fb^{-1} to correspond to roughly half-way through the 2015 run
  float targetlumi = 5.0;

  //To simplify things a little I've placed the number of events ran over for each sample below:
  float nRunDY = pow( 10, 5);
  float nRunWZ = pow( 10, 5);
  float nRunWJets = pow( 10, 5);
  float nRunTT = pow( 10, 5);
  float nRunTTZ = pow( 10, 5);
  
  /* Now, to figure out how to scale things correctly we need to remember the following equation:

     N_events = lumi * cross section

     Hence, if we want to normalize to 5 fb^{-1} for each sample we want

     N_events / cross section = 5.0

     But, you will find that this is not always the case. NOTE: This equation is not true because we are forcing the events to be normalized to 5 inverse femtobarns.
     To make this equation true we will introduce a weight W_WZ for instance, which will make our equation true. What we really want then is to weight N_events as below:

     N_normalized = Weight_Sample * N_RunSample

     Where N_normalized is the normalized number of events we would see for that sample if we had collected 5^{-1] fb of data.

  */

  //To start, look up the cross sections for the various samples in MCFM and fill them out below:

  float xsecDY = ;
  float xsecWZ = ;
  float xsecWJets = ;
  float xsecTT = ;
  float xsecTTZ = ;

  //Here some math is done for you :)

  float weightDY = (targetlumi*xsecDY) / (nRunDY);
  float weightWZ = (targetlumi*xsecWZ) / (nRunWZ);
  float weightWJets = (targetlumi*xsecWJets) / (nRunWJets);
  float weightTT = (targetlumi*xsecTT) / (nRunTT);
  float weightTTZ = (targetlumi*xsecTTZ) / (nRunTTZ);

  //Now that we have the weights we can find out how many events passed out selection:

  int nSelDY = tDY->Draw("","76 < elDiMass < 106 ");
  int nSelWZ = tWZ->Draw("","76 < elDiMass < 106 ");
  int nSelWJets = tWJets->Draw("","76 < elDiMass < 106 ");
  int nSelTT = tTT->Draw("","76 < elDiMass < 106 ");
  int nSelTTZ = tTTZ->Draw("","76 < elDiMass < 106 ");

  //now weight them and print out the values

  float nNormDY = nSelDY * weightDY;
  float nNormWZ = nSelWZ * weightWZ;
  float nNormWJets = nSelWJets * weightWJets;
  float nNormTT = nSelTT * weightTT;
  float nNormTTZ = nSelTTZ * weightTTZ;
  

  std::cout<<"Number of events passing mass window cut from DY: "<<nNormDY<<std::endl;
  std::cout<<"Number of events passing mass window cut from WZ: "<<nNormWZ<<std::endl;
  std::cout<<"Number of events passing mass window cut from WJets: "<<nNormWJets<<std::endl;
  std::cout<<"Number of events passing mass window cut from TT: "<<nNormTT<<std::endl;
  std::cout<<"Number of events passing mass window cut from TTZ: "<<nNormTTZ<<std::endl;
}
