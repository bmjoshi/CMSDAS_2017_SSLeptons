#include <iostream>
#include <vector>
#include "TTree.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TH1.h"
#include "THStack.h"
#include "TGraphAsymmErrors.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TColor.h"
#include "TLegend.h"
#include "TLatex.h"

const double M_EL = 0.000510998928; //Mass of electron in GeV
const double M_MU = 0.1056583715;   //Mass of muon in GeV
const double M_Z  = 91.1876;        //Mass of Z boson
const double dM   = 15;             //Size of window around Z

// Takes pointers to two int variables Ntot_<sample>, Nss_<sample> and a sample and increases them by the number of total events in mass window and number of events with same sign leptons in mass window respectively 
void runSample(TH1F* hist, string sampleName, string weightStr="noWeight", bool blinded=true){
  //Input samples should already be preselected for two tight same-sign dileptons
  printf("---------------------------------------------------------------\n");
  printf("Running over sample: %s with weight: %s\n", sampleName.c_str(), weightStr.c_str() );
  TChain* t = new TChain("tEvts_ssdl");
  t->Add( ("/eos/uscms/store/user/cmsdas/2017/long_exercises/Same-Sign-Dileptons/"+sampleName ).c_str() );
  int nEntries = t->GetEntries();
  //float lepPts1, lepEtas1, lepPhis1, lepEs1;
  //float lepPts2, lepEtas2, lepPhis2, lepEs2;
  float DilepMass, AssocMass, HT;
  int nJets;
  float NPWeight, ChargeMisIDWeight;
  float sampleWeight;

  //Set branch addresses
  /*
  t->SetBranchAddress("Lep1Pt", &lepPts1);
  t->SetBranchAddress("Lep1Eta", &lepEtas1);
  t->SetBranchAddress("Lep1Phi", &lepPhis1);
  t->SetBranchAddress("Lep1E", &lepEs1);
  t->SetBranchAddress("Lep2Pt", &lepPts2);
  t->SetBranchAddress("Lep2Eta", &lepEtas2);
  t->SetBranchAddress("Lep2Phi", &lepPhis2);
  t->SetBranchAddress("Lep2E", &lepEs2);
  */

  t->SetBranchAddress("DilepMass", &DilepMass);
  t->SetBranchAddress("AssocMass", &AssocMass);
  t->SetBranchAddress("cleanAK4HT", &HT);
  t->SetBranchAddress("nCleanAK4Jets", &nJets);
  t->SetBranchAddress("NPWeight", &NPWeight);
  t->SetBranchAddress("ChargeMisIDWeight", &ChargeMisIDWeight);

  for(int ient = 0; ient < nEntries; ient++){
    t->GetEntry(ient);
    if(ient % 100000 ==0) std::cout<<"Completed "<<ient<<" out of "<<nEntries<<" events"<<std::endl;
    
    if (blinded){
      if (HT>900) continue; //inverted HT cut to stay in control region, DO NOT REMOVE
    }
    if (AssocMass > M_Z - dM && AssocMass < M_Z + dM) continue; //veto associated Z-mass cut
    if (DilepMass > M_Z - dM && DilepMass < M_Z + dM) continue; //veto dilepton mass around Z-peak
    if (DilepMass < 20) continue; //require minimum Dilepton Mass of 20 GeV
    if (nJets < 2) continue; //require 2 or more jets

    //select appropriate weighting
    if (weightStr=="NPWeight"){
      sampleWeight=NPWeight;
    }
    else if (weightStr=="ChargeMisIDWeight"){
      sampleWeight=ChargeMisIDWeight;
    }
    else if (weightStr=="noWeight"){
      sampleWeight=1;
    }
    else{
      printf("WARNING: atypical weightStr, setting weight to 1\n");
      sampleWeight=1;
    }

    hist->Fill(HT, sampleWeight);

  }
}



//The main function that will be called
void Ex_3p1(){

  float targetlumi=36.4; //fb-1 - lumi for the full 2016 dataset
  float nRunX53L=293600; //number of left-handed X53 (mass=1000GeV) events that we ran over
  float nRunX53R=299000; //number of right-handed X53 (mass=1000GeV) events that we ran over
  float xsecX53=42.7; //xsec for a 1000 GeV X53 = 0.0427 pb-1
  float weightX53L = (targetlumi*xsecX53) / (nRunX53L);
  float weightX53R = (targetlumi*xsecX53) / (nRunX53R);

  float nRunTTW=252673.;
  float xsecTTW=204.3;
  float weightTTW=(targetlumi*xsecTTW) / (nRunTTW);

  float nRunTTZ=398600.;
  float xsecTTZ=252.9;
  float weightTTZ=(targetlumi*xsecTTZ) / (nRunTTZ);

  int nBins = 20;
  float xMin =0.;
  float xMax =2000.;

  //Create histograms
  //TH1F* TTbarHT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  THStack *hsHT = new THStack("hs", "HT (blinded)");
  TH1F* data_bkgnd_cm_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  TH1F* data_bkgnd_np_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  TH1F* data_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  TH1F* TTZ_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  TH1F* TTW_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  TH1F* mc_X53L_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  TH1F* mc_X53R_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);
  TH1F* h_err =  new TH1F("h_err","error",nBins,xMin,xMax);;

  TH1F* TTX_HT_h = new TH1F("firstStr","secondStr",nBins,xMin,xMax);

  //runSample(WZHT_h, "WZ_MuCBTightMiniIsoTight_ElMVATightRC_2016_B-H.root");
  // DO NOT CHANGE THIS SECTION OF THE CODE
  runSample(data_bkgnd_np_HT_h, "NonPromptData_MuCBTightMiniIsoTight_ElMVATightRC_2016B-H.root", "NPWeight");
  runSample(data_bkgnd_cm_HT_h, "ChargeMisID_MuCBTightMiniIsoTight_ElMVATightRC_2016B-H.root", "ChargeMisIDWeight");
  runSample(data_HT_h, "Data_MuCBTightMiniIsoTight_ElMVATightRC_2016_B-H.root", "noWeight");
  runSample(mc_X53L_HT_h, "X53X53m1000LH_Inc_MuCBTightMiniIsoTight_ElMVATightRC_2016_B-H.root", "noWeight", false);
  runSample(mc_X53R_HT_h, "X53X53m1000RH_Inc_MuCBTightMiniIsoTight_ElMVATightRC_2016_B-H.root", "noWeight", false);
  // IF YOU WOULD LIKE TO TRY TO ADD MORE SAMPLES, YOU CAN DO SO BELOW. BUT MAKE SURE NOT TO RUN THIS METHOD WITH THE FALSE FLAG, THIS WILL UNBLIND THE SAMPLE (technically, if it's a Monte Carlo sample, it's perfectly fine for you to unblind, and we have do so for the two monto carlo samples above, but when in doubt stay blinded)
  runSample(TTZ_HT_h, "TTZ_MuCBTightMiniIsoTight_ElMVATightRC_2016B-H.root", "noWeight");
  runSample(TTW_HT_h, "TTW_MuCBTightMiniIsoTight_ElMVATightRC_2016B-H.root", "noWeight");
  
  data_bkgnd_np_HT_h->SetFillColor(kRed);
  data_bkgnd_cm_HT_h->SetFillColor(kCyan);
  TTZ_HT_h->SetFillColor(kGreen+2);
  TTW_HT_h->SetFillColor(kGreen+2);
  TTX_HT_h->SetFillColor(kGreen+2);

  data_HT_h->SetMarkerStyle(21);
  mc_X53L_HT_h->SetLineColor(kBlue);
  mc_X53L_HT_h->SetLineStyle(2);
  mc_X53R_HT_h->SetLineColor(kBlue);
  //Scale MC by lumi
  // * SCALE THE SIGNAL MC HERE
  //mc_X53L_HT_h->Scale(weightX53L);
  //mc_X53R_HT_h->Scale(weightX53R);
  TTZ_HT_h->Scale(weightTTZ);
  TTW_HT_h->Scale(weightTTW);

  //Add the SM SSDL background together
  TTX_HT_h->Add(TTZ_HT_h);
  TTX_HT_h->Add(TTW_HT_h);

  //Add the SM SSDL to the stack hist
  hsHT->Add(TTX_HT_h); //Added the combine TT+X channel to stack
  hsHT->Add(data_bkgnd_cm_HT_h);
  hsHT->Add(data_bkgnd_np_HT_h);

  //Also add the same histos to h_err - so they have the same bin contents for plotting the error
  h_err->Add(TTX_HT_h); //Added the combine TT+X channel to stack
  h_err->Add(data_bkgnd_cm_HT_h);
  h_err->Add(data_bkgnd_np_HT_h);

  //Add errors for stack histo (add in quad)
  std::vector<float> errs;
  float sysErrorNP =0.0;  // * REPLACE WITH VALUE THAT REFLECT WHAT YOU THINK THE CHARGEMISID SYSTEMATIC ERROR SHOULD BE SET TO 
  float sysErrorCM =0.0; // * REPLACE WITH VALUE THAT REFLECT WHAT YOU THINK THE CHARGEMISID SYSTEMATIC ERROR SHOULD BE SET TO 
  for(unsigned int ibin=1; ibin<= data_bkgnd_np_HT_h->GetNbinsX(); ibin++){
    //nonprompt
    float etemp = pow(data_bkgnd_np_HT_h->GetBinError(ibin),2 ); //stat
    etemp = etemp + pow(sysErrorNP * data_bkgnd_np_HT_h->GetBinContent(ibin),2 ); //sys
    //chargeMisID
    etemp = etemp + pow(data_bkgnd_cm_HT_h->GetBinError(ibin),2 ); //stat
    etemp = etemp + pow(sysErrorCM * data_bkgnd_cm_HT_h->GetBinContent(ibin),2 ); //sys

    etemp = pow(etemp, 0.5);
    h_err->SetBinError(ibin,etemp);
    errs.push_back(etemp);
  }



  TCanvas c1;
  //c1.SetLogy();
  hsHT->Draw("hist");
  // * PLOT THE ERROR BELOW *
  h_err->SetFillStyle(3344);
  h_err->SetFillColor(1);
  h_err->Draw("SAME E2");
  data_HT_h->Draw("pe SAME");
  //UNCOMMENT THESE LINES TO PLOT THE SIGNAL MONTE CARLO
  //mc_X53L_HT_h->Draw("SAME");
  //mc_X53R_HT_h->Draw("SAME");

  //In general, when you present plots, you should report the luminosity and center of mass energy 
  TLatex* lumitex = new TLatex();
  lumitex->SetNDC();
  lumitex->SetTextSize(0.04);
  lumitex->DrawLatex(0.65,0.65,"36.4 fb^{-1} (13 TeV)");

  TLegend* leg = new TLegend(0.65,0.7,0.9,0.9);
  //leg->SetHeader("The Legend Title","C"); // option "C" allows to center the header
  leg->AddEntry(data_HT_h,"Data","pl");
  leg->AddEntry(data_bkgnd_np_HT_h,"Non-prompt","f");
  leg->AddEntry(data_bkgnd_cm_HT_h,"ChargeMisID","f");
  leg->AddEntry(TTX_HT_h,"TT+X","f");
  leg->AddEntry(h_err,"Bckgnd Error","f");
  //leg->AddEntry(TTZ_HT_h,"TT+Z","f");
  leg->AddEntry(mc_X53L_HT_h,"X53-L-M1000 MC","l");
  leg->AddEntry(mc_X53R_HT_h,"X53-R-M1000 MC","l");
  leg->Draw();

  c1.Print("blindedHT.pdf");

}

