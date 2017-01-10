#include <iostream>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphAsymmErrors.h"
#include "ObjectID.C"
#include "TLorentzVector.h"
#include "TChain.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z

// Takes pointers to two int variables Ntot_<sample>, Nss_<sample> and a sample and increases them by the number of total events in mass window and number of events with same sign leptons in mass window respectively 
void runSample(int* Ntot_sample, int* Nss_sample, string sampleName = "ChargeMisID_Data_Run2016_Electrons_MVATightRC.root"){
  printf("Running over sample: %s\n", sampleName.c_str());
  TChain* t = new TChain("ChargeMisID");
  t->Add( ("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/"+sampleName ).c_str() );
  int nEntries = t->GetEntries();
  float elPts1, elEtas1, elPhis1, elEs1;
  float elPts2, elEtas2, elPhis2, elEs2;
  int elCharge1, elCharge2;
  int elTight1, elTight2;
  int lep1etabin, lep2etabin;

  //Set branch addresses
  t->SetBranchAddress("Lep1Pt", &elPts1);
  t->SetBranchAddress("Lep1Eta", &elEtas1);
  t->SetBranchAddress("Lep1Phi", &elPhis1);
  t->SetBranchAddress("Lep1E", &elEs1);
  t->SetBranchAddress("Lep1Charge", &elCharge1);
  t->SetBranchAddress("Lep1Tight", &elTight1);
  t->SetBranchAddress("Lep2Pt", &elPts2);
  t->SetBranchAddress("Lep2Eta", &elEtas2);
  t->SetBranchAddress("Lep2Phi", &elPhis2);
  t->SetBranchAddress("Lep2E", &elEs2);
  t->SetBranchAddress("Lep2Charge", &elCharge2);
  t->SetBranchAddress("Lep2Tight", &elTight2);

  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;

    // THE LINES BELOW DETERMINES THE QUALITY OF THE ELECTRONS
    //if (elTight1 != true) continue; 
    //if (elTight2 != true) continue;     

    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(elPts1, elEtas1, elPhis1, M_EL);
    v2.SetPtEtaPhiM(elPts2, elEtas2, elPhis2, M_EL);
    double mass = (v1+v2).M();

    if (mass > M_Z - dM && mass < M_Z + dM){
      (*Ntot_sample)++;
      //if (vEl.at(ui)->charge == vEl.at(uj)->charge){
      if ( elCharge1 == elCharge2){
	(*Nss_sample)++;
      }
    }	
  }
}

//The main function that will be called
void Ex_1p2(){
  /*
    THIS SCRIPT WILL LOOP OVER ALL EVENTS IN THE DIFFERENT SAMPLES AND PRINT OUT HOW MANY
    PASS THE MASS WINDOW CUT, AND ALSO HOW MANY PASS THE MASS WINDOW CUT AND ARE SAME SIGN
   */

  /*
  //load TFiles
  TFile* fDY = new TFile("/uscms_data/d3/clint/public/ljmet_tree_DY.root");
  TFile* fWZ = new TFile("/uscms_data/d3/clint/public/ljmet_tree_WZ.root");
  TFile* fWJets = new TFile("/uscms_data/d3/clint/public/ljmet_tree_WJets.root");
  TFile* fTT = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TT1.root");
  TFile* fTTZ = new TFile("/uscms_data/d3/clint/public/ljmet_tree_TTZ.root");

  //load the TTrees
  TTree* tDY = (TTree*)fDY->Get("ljmet");
  TTree* tWZ = (TTree*)fWZ->Get("ljmet");
  TTree* tWJets = (TTree*)fWJets->Get("ljmet");
  TTree* tTT = (TTree*)fTT->Get("ljmet");
  TTree* tTTZ = (TTree*)fTTZ->Get("ljmet");
*/

  //define our taget luminosity - we are picking 36.4 fb^{-1} to correspond to the total 2016 run
  float targetlumi = 36.4;

  //To simplify things a little I've placed the number of events ran over for each sample below:
  float nRunDY = 1366703;
  float nRunWZ = 237484;
  float nRunWJets = 3828404;
  float nRunTT = 2206600;
  float nRunTTZ = 249275;
  
  /* Now, to figure out how to scale things correctly we need to remember the following equation:

     N_events = lumi * cross section

     Hence, if we want to normalize to 5 fb^{-1} for each sample we want

     N_events / cross section = 5.0

     NOTE: This equation is not currently true because we are forcing the events to be normalized to 5 inverse femtobarns.
     To make this equation true we will introduce a weight W_WZ for instance, which will make our equation true. What we really want then is to weight N_events as below:

     N_normalized = Weight_Sample * N_RunSample

     Where N_normalized is the normalized number of events we would see for that sample if we had collected 5^{-1] fb of data.

  */

  //To start, look up the cross sections for the various samples in MCFM and fill them out below:

  float xsecDY = 1.0;
  float xsecWZ = 1.0;
  float xsecWJets = 1.0;
  float xsecTT = 1.0;
  float xsecTTZ = 1.0;

  //Here some math is done for you :)

  float weightDY = (targetlumi*xsecDY) / (nRunDY);
  float weightWZ = (targetlumi*xsecWZ) / (nRunWZ);
  float weightWJets = (targetlumi*xsecWJets) / (nRunWJets);
  float weightTT = (targetlumi*xsecTT) / (nRunTT);
  float weightTTZ = (targetlumi*xsecTTZ) / (nRunTTZ);

  //Now that we have the weights we can find out how many events passed out selection:

  int* Ntot_DY =new int(0);
  int* Ntot_WZ =new int(0);
  int* Ntot_WJets =new int(0);
  int* Ntot_TT =new int(0);
  int* Ntot_TTZ=new int(0);

  int* Nss_DY =new int(0);
  int* Nss_WZ =new int(0);
  int* Nss_WJets =new int(0);
  int* Nss_TT =new int(0);
  int* Nss_TTZ=new int(0);

  float nRunTest = 200000;
  float xsecTest = 1.0;
  float weightTest = (targetlumi*xsecTest) / (nRunTest);
  int* Ntot_Test=new int(0);
  int* Nss_Test=new int(0);


  //Now run over the samples and count up the number of (total or same signed) events and change the Ntot and Nss variable appropriately (see definition of runSample at top)
  runSample(Ntot_Test,Nss_Test);

  float nNormTest = *Ntot_Test * weightTest;
  float nSSNormTest = *Nss_Test * weightTest;
  
  std::cout<<"Number of events passing mass window cut from Test: "<<nNormTest<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from Test: "<<nSSNormTest<<std::endl;

  /*
  //now weight them and print out the values
  float nNormDY = Ntot_DY * weightDY;
  float nNormWZ = Ntot_WZ * weightWZ;
  float nNormWJets = Ntot_WJets * weightWJets;
  float nNormTT = Ntot_TT * weightTT;
  float nNormTTZ = Ntot_TTZ * weightTTZ;

  //now weight them and print out the values
  float nSSNormDY = Nss_DY * weightDY;
  float nSSNormWZ = Nss_WZ * weightWZ;
  float nSSNormWJets = Nss_WJets * weightWJets;
  float nSSNormTT = Nss_TT * weightTT;
  float nSSNormTTZ = Nss_TTZ * weightTTZ;
  
  std::cout<<"Number of events passing mass window cut from DY: "<<nNormDY<<std::endl;
  std::cout<<"Number of events passing mass window cut from WZ: "<<nNormWZ<<std::endl;
  std::cout<<"Number of events passing mass window cut from WJets: "<<nNormWJets<<std::endl;
  std::cout<<"Number of events passing mass window cut from TT: "<<nNormTT<<std::endl;
  std::cout<<"Number of events passing mass window cut from TTZ: "<<nNormTTZ<<std::endl;

  std::cout<<"Number of same-sign events passing mass window cut from DY: "<<nSSNormDY<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from WZ: "<<nSSNormWZ<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from WJets: "<<nSSNormWJets<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from TT: "<<nSSNormTT<<std::endl;
  std::cout<<"Number of same-sign events passing mass window cut from TTZ: "<<nSSNormTTZ<<std::endl;
*/

  /*ADD CODE TO PRINT OUT THE TOTAL PURITY YOU WILL HAVE TO DEFINE YOUR OWN VARIABLES AND PRINT THEM OUT. 
   YOU PROBABLY WANT ONE FOR TOTAL EVENTS AND ONE FOR NON DY EVENTS*/

}

