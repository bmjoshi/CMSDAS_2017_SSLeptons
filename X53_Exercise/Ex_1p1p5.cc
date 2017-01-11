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

int getEtaBin(float abseta){
  int bin=-1;
  if(abseta>2.0) bin=6;
  else if(abseta>1.556) bin=5;
  else if(abseta>1.442) bin=4;
  else if(abseta>0.8) bin=3;
  else if(abseta>0.4) bin=2;
  else if(abseta>=0.0) bin=1;

  return bin;
}

void Ex_1p1p5(){
  // This Macro is very similar to the previous one (Ex_1p1.cc) but finds the single lepton charge mis ID rate for high pt electrons (>100 GeV). 
  //Note that you could use Ex_1p1.cc to find the rate for high pt electrons by asking for both electrons to be high pt. You can even try this if you have time.
  // You will notice that the statistics for events that have both electrons at high pt is low, so instead, we will use what have learned about the charge misID rate of low pt electrons (from Ex_1p1.cc) and event that have one high and one low pt electron
  const Int_t NBINS = 6;
  double weightEtaLowPtFill[NBINS]={0.000249771, 0.000174497, 0.0011059, 0.0, 0.00788604, 0.00899985}; //weights array for you to replace by hand (after running and getting rates from low pt histograms in Ex_1p1)

  //load file and tree
  //TChain* t = new TChain("ljmet");
  TChain* t = new TChain("ChargeMisID");
  t->Add("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/ChargeMisID_Data_Run2016_Electrons_MVATightRC.root");
  //t->Add("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/ChargeMisID_Data_2016_Electrons_MVATightRC.root");
  //t->Add("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/DY.root");
  
  int nEntries = t->GetEntries();

  //kinematic variables
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
  //t->SetBranchAddress("Lep1Tight", &elTight1);
  t->SetBranchAddress("Lep2Pt", &elPts2);
  t->SetBranchAddress("Lep2Eta", &elEtas2);
  t->SetBranchAddress("Lep2Phi", &elPhis2);
  t->SetBranchAddress("Lep2E", &elEs2);
  t->SetBranchAddress("Lep2Charge", &elCharge2);
  //t->SetBranchAddress("Lep2Tight", &elTight2);

  // Histogram for charge misID rate measure per Eta and Pt binning for high pt (>100 GeV) electrons
  Double_t edges[NBINS+1] = {0.0, 0.4, 0.8, 1.442, 1.556, 2.0, 3.0};

  TH1F* ssEtaHighPt= new TH1F("ssEtaHighPt","#eta",NBINS, edges); //      (h3)
  TH1F* totEtaHighPt= new TH1F("totEtaHighPt","#eta",NBINS, edges); //     (h2)
  TH1F* weightEtaHighPt= new TH1F("weightEtaHighPt","#eta",NBINS, edges); //  (h1)

  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    if(ient % 1000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;

    // THE LINES BELOW DETERMINES THE QUALITY OF THE ELECTRONS // All electrons in sample are at least Loose electrons, the tight flag can be used to select only events with two (or at least one) tight electron
    //if ((elTight1 != true) || (elTight2 != true))) continue; //only consider events with two tight electrons

    // Fill histogram only if in same eta bin 
    // we do this so that we meassure the charge misID rate of a single lepton with respect to eta
    //if (getEtaBin( abs(elEtas1) )!= getEtaBin( abs(elEtas2)) ) continue;

    TLorentzVector v1, v2;
    v1.SetPtEtaPhiM(elPts1, elEtas1, elPhis1, M_EL);
    v2.SetPtEtaPhiM(elPts2, elEtas2, elPhis2, M_EL);
    double mass = (v1+v2).M();

    /*
    allmass->Fill(mass);

    totEtaHist->Fill(v1.Eta());
    totEtaHist->Fill(v2.Eta());
	  
    totPtHist->Fill(v1.Pt());
    totPtHist->Fill(v2.Pt());
    */

    float highPtEta;
    float lowPtEta;

    if (elPts1 >100 && elPts2 < 100) {// take events with one high and one low pt electron
      highPtEta=elEtas1;
      lowPtEta=elEtas2;
    }
    else if (elPts1 <100 && elPts2 > 100) {
      highPtEta=elEtas2;
      lowPtEta=elEtas1;
    }    
    else 
      continue;
    totEtaHighPt->Fill(highPtEta); // fill just once indicating one event with both leptons in same eta bin
    weightEtaHighPt->Fill(highPtEta, weightEtaLowPtFill[getEtaBin( abs(lowPtEta)) - 1] ); // fill a correction histogram weighted by charge misID rate of lower pt electron

    if (elCharge1 == elCharge2 ){
      ssEtaHighPt->Fill(highPtEta);
    }
  }

  double weightEtaHighPtFill[NBINS]={0};

  //single lepton charge mis ID rate for high pt (<100 GeV)
  for (int i=1; i<NBINS+1;i++ ){
    weightEtaHighPtFill[i-1] = (ssEtaHighPt->GetBinContent(i) - weightEtaHighPt->GetBinContent(i) )/totEtaHighPt->GetBinContent(i);
    printf("low pt single lepton charge Mis ID rate for eta bin %d is: %f\n", i, weightEtaHighPtFill[i-1]);
  }

  

  //Draw the histograms
  /*
  TCanvas c1;
  allmass->Draw();

  c1.Print("DiElectronInvariantMass_all.pdf");
  */
  /* ADD CODE TO PLOT MASS DISTRIBUTION FOR OPPOSITE SIGN AND SAME SIGN. 
     MOST OF IT IS THERE, YOU JUST NEED TO NORMALIZE OSMASS AND SSMASS BEFORE DRAWING THEM
  */ 
  /*
  TCanvas c2;
  osmass->Draw();

  c2.Print("DiElectronInvariantMass_os.pdf");

  TCanvas c3;
  ssmass->Draw();
  
  c3.Print("DiElectronInvariantMass_ss.pdf");  
  */

  //ADD CODE TO PRINT OUT CHARGE MISID RATE
  //Determine the charge misID rate with respect to abs(Eta) and Pt  
  // Bin abs(Eta) according to the bins given in getEtaBin (these are the bins used in the actual analysis)
  // and bin pt in a low pt (<100 GeV) bin and a high (>100GeV) bin
  // You will be asked for 12 rates 


  /* ADD CODE HERE TO MAKE PLOTS OF THE CHARGE MISID RATE VS PT, ETA. SOME OF IT IS THERE SO THAT WE ALL HAVE THE SAME OUTPUT FILES
     HINT1: MAKE USE OF TGRAPHASYMMERRORS, 
     HINT2: SET YOUR VS ETA BINS ACCORDING TO THE getEtaBin METHOD PROVIDED AT THE BOTTOM OF THIS FILE
   */
  //tcanvas for eta plot
  /*
  TCanvas c4;

  //add code!

  c4.Print("chargeMisID_vEta.pdf");

  //same as above but for pt
  TCanvas c5;

  c5.Print("chargeMisID_vPt.pdf");
  */
}
