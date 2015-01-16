#include <iostream>
#include "TFile.h"
#include "TH1.h"

void diff_plots(){

  TFile* fsig = new TFile("sig_HT.root");
  TFile* fbg  = new TFile("bg_nonPrompt.root");


  TH1F* sig_leadlepPt = fsig->Get("leadlep_pt");
  TH1F* bg_leadlepPt = fbg->Get("leadlep_pt");

  sig_leadlepPt->Scale(1/sig_leadlepPt->Integral());
  bg_leadlepPt->Scale(1/bg_leadlepPt->Integral());

  TH1F* sig_sublepPt = fsig->Get("sublep_pt");
  TH1F* bg_sublepPt = fbg->Get("sublep_pt");

  sig_sublepPt->Scale(1/sig_sublepPt->Integral());
  bg_sublepPt->Scale(1/bg_sublepPt->Integral());


  TCanvas c1;
  sig_leadlepPt->SetLineColor(kBlue);
  sig_leadlepPt->Draw();
  bg_leadlepPt->Draw("same");
  c1.Print("diff_leadlepPt.pdf");

  TCanvas c2;
  sig_sublepPt->SetLineColor(kBlue);
  sig_sublepPt->Draw();
  bg_sublepPt->Draw("same");
  c2.Print("diff_sublepPt.pdf");


  TH1F* sig_leadjetPt = fsig->Get("leadjet_pt");
  TH1F* bg_leadjetPt = fbg->Get("leadjet_pt");

  sig_leadjetPt->Scale(1/sig_leadjetPt->Integral());
  bg_leadjetPt->Scale(1/bg_leadjetPt->Integral());

  TH1F* sig_subjetPt = fsig->Get("subjet_pt");
  TH1F* bg_subjetPt = fbg->Get("subjet_pt");

  sig_subjetPt->Scale(1/sig_subjetPt->Integral());
  bg_subjetPt->Scale(1/bg_subjetPt->Integral());


  TCanvas c3;
  sig_leadjetPt->SetLineColor(kBlue);
  sig_leadjetPt->Draw();
  bg_leadjetPt->Draw("same");
  c3.Print("diff_leadjetPt.pdf");

  TCanvas c4;
  sig_subjetPt->SetLineColor(kBlue);
  sig_subjetPt->Draw();
  bg_subjetPt->Draw("same");
  c4.Print("diff_subjetPt.pdf");


  TH1F* sig_HT = fsig->Get("HT_h");
  TH1F* bg_HT = fbg->Get("HT_h");

  sig_HT->Scale(1/sig_HT->Integral());
  bg_HT->Scale(1/bg_HT->Integral());

  TH1F* sig_MET = fsig->Get("met_h");
  TH1F* bg_MET = fbg->Get("met_h");

  sig_MET->Scale(1/sig_MET->Integral());
  bg_MET->Scale(1/bg_MET->Integral());


  TCanvas c5;
  sig_HT->SetLineColor(kBlue);
  sig_HT->Draw();
  bg_HT->Draw("same");
  c5.Print("diff_HT.pdf");

  TCanvas c6;
  sig_MET->SetLineColor(kBlue);
  sig_MET->Draw();
  bg_MET->Draw("same");
  c6.Print("diff_MET.pdf");

  TH1F* sig_M = fsig->Get("dilep_mass");
  TH1F* bg_M =  fbg->Get("dilep_mass");

  sig_M->Scale(1 / sig_M->Integral());
  bg_M->Scale( 1 / bg_M->Integral());

  TCanvas c7;
  sig_M->SetLineColor(kBlue);
  sig_M->Draw();
  bg_M->Draw("same");
  c7.Print("diff_Mass.pdf");


  TH1F* sig_njet = fsig->Get("njets");
  TH1F* bg_njet =  fbg->Get("njets");

  sig_njet->Scale(1 / sig_njet->Integral());
  bg_njet->Scale( 1 / bg_njet->Integral());

  TCanvas c7;
  sig_njet->SetLineColor(kBlue);
  sig_njet->Draw();
  bg_njet->Draw("same");
  c7.Print("diff_Njets.pdf");

  TH1F* sig_TriMass = fsig->Get("TriMass");
  TH1F* bg_TriMass =  fbg->Get("TriMass");

  sig_TriMass->Scale(1 / sig_TriMass->Integral());
  bg_TriMass->Scale( 1 / bg_TriMass->Integral());

  TCanvas c8;
  sig_TriMass->SetLineColor(kBlue);
  sig_TriMass->Draw();
  bg_TriMass->Draw("same");
  c8.Print("diff_TriMass.pdf");



}
