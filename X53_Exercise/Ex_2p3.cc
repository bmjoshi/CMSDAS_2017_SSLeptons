#include <iostream>
#include <string>
#include "TFile.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "tdrstyle.C"

float weight(float fr){
  float weight = fr / (1 - fr);
  return weight;
};

void DrawAndSave(TH1F* h_pred, TH1F* h_obs, std::string pdfname, std::string flavor, std::string lepflavor){

  //rebin, currently bins are 5 GeV so rebin by 20 to make 100 GeV bins
  if(flavor.find("qcd")==std::string::npos){
    h_pred->Rebin(20);
    h_obs->Rebin(20);
    h_pred->Sumw2();
    h_obs->Sumw2();
  }
  else{};

  TCanvas* c = new TCanvas();
  c->SetLogy();
  std::string title = "Predicted vs. Observed H_{T} for "+flavor+" "+lepflavor+";H_{T};N_{Events}";
  h_pred->SetTitle(title.c_str());
  
  h_pred->SetLineColor(kRed+1);
  h_obs->SetMarkerStyle(20);
  h_obs->SetMarkerColor(kBlack);

  //make lumi weight
  float lumiweight = 36400 * (831.76 / 115091972);
  if(flavor.find("qcd")==std::string::npos){
    h_pred->Scale(lumiweight);
    h_obs->Scale(lumiweight);
  }
  
  //get fake rate weight
  float frweight=0;
  if(flavor=="light" && lepflavor=="Electrons") frweight = weight(0.218); 
  if(flavor=="charm"&& lepflavor=="Electrons") frweight = weight(0.161); 
  if(flavor=="bottom"&& lepflavor=="Electrons") frweight = weight(0.142); 
  if(flavor=="fake"&& lepflavor=="Electrons") frweight = weight(0.226); 
  if(flavor=="unmatched"&& lepflavor=="Electrons") frweight = weight(0.562); 
  if(flavor=="data"&& lepflavor=="Electrons") frweight = weight(0.206); 

  if(flavor=="light" && lepflavor=="Muons") frweight = weight(0.471); 
  if(flavor=="charm"&& lepflavor=="Muons") frweight = weight(0.420); 
  if(flavor=="bottom"&& lepflavor=="Muons") frweight = weight(0.353); 
  if(flavor=="fake"&& lepflavor=="Muons") frweight = weight(0.085); 
  if(flavor=="unmatched"&& lepflavor=="Muons") frweight = weight(0.014); 
  if(flavor=="data"&& lepflavor=="Muons") frweight = weight(0.427); 

  if(flavor.find("qcd")==std::string::npos) h_pred->Scale(frweight);

  //THESE LINES WILL MAKE A LATEX READY TABLE FOR YOU TO COMPARE THE NUMBER OF PREDICTED AND OBSERVED EVENTS BOTH FOR HT > 900 AND WITH NO HT CUT
  /*double obserr=0;
  float obs = h_obs->IntegralAndError(9,51,obserr);
  double prederr=0;
  float pred = h_pred->IntegralAndError(9,51,prederr);
  std::cout<<std::fixed<<setprecision(2)<<lepflavor<<"&"<<flavor<<"&$"<<obs<<"\\pm"<<obserr<<"$&$"<<pred<<"\\pm"<<prederr<<"$\\\\"<<std::endl;

  double totobserr=0;
  float totobs = h_obs->IntegralAndError(1,51,totobserr);
  double totprederr=0;
  float totpred = h_pred->IntegralAndError(1,51,totprederr);
  std::cout<<std::fixed<<setprecision(2)<<lepflavor<<"&"<<flavor<<"&$"<<totobs<<"\\pm"<<totobserr<<"$&$"<<totpred<<"\\pm"<<totprederr<<"$\\\\"<<std::endl;
  */

  h_pred->Draw("hist e");
  h_obs->Draw("pesame");

  TLegend* leg = new TLegend(0.4,0.7,0.9,0.9);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  std::string predleg = "Predicted Events for "+flavor+" "+lepflavor;
  std::string obsleg = "Observed Events for "+flavor+" "+lepflavor;
  leg->AddEntry(h_obs,obsleg.c_str(),"p");
  leg->AddEntry(h_pred,predleg.c_str(),"l");
  leg->Draw("same");
  c->Print(pdfname.c_str());
  delete c;

}

TH1F* getTotPredHist(TFile* f,std::string chan){


  if(chan=="el"){

    TH1F* h_light = (TH1F*) f->Get("h_el_pred_light");
    TH1F* h_charm = (TH1F*) f->Get("h_el_pred_charm");
    TH1F* h_bottom = (TH1F*) f->Get("h_el_pred_bottom");
    TH1F* h_fake = (TH1F*) f->Get("h_el_pred_fake");
    TH1F* h_unm = (TH1F*) f->Get("h_el_pred_unm");
    TH1F* h_tot = (TH1F*)h_light->Clone();
    h_tot->Add(h_charm);
    h_tot->Add(h_bottom);
    h_tot->Add(h_fake);
    h_tot->Add(h_unm);
    
    return h_tot;
  }

  else{
    TH1F* h_light = (TH1F*) f->Get("h_mu_pred_light");
    TH1F* h_charm = (TH1F*) f->Get("h_mu_pred_charm");
    TH1F* h_bottom = (TH1F*) f->Get("h_mu_pred_bottom");
    TH1F* h_fake = (TH1F*) f->Get("h_mu_pred_fake");
    TH1F* h_unm = (TH1F*) f->Get("h_mu_pred_unm");
    TH1F* h_tot = (TH1F*) h_light->Clone();
    h_tot->Add(h_charm);
    h_tot->Add(h_bottom);
    h_tot->Add(h_fake);
    h_tot->Add(h_unm);
    
    return h_tot;
  }


};


TH1F* getTotObsHist(TFile* f,std::string chan){


  if(chan=="el"){

    TH1F* h_light = (TH1F*) f->Get("h_el_obs_light");
    TH1F* h_charm = (TH1F*) f->Get("h_el_obs_charm");
    TH1F* h_bottom = (TH1F*) f->Get("h_el_obs_bottom");
    TH1F* h_fake = (TH1F*) f->Get("h_el_obs_fake");
    TH1F* h_unm = (TH1F*) f->Get("h_el_obs_unm");
    TH1F* h_tot = (TH1F*) h_light->Clone();
    h_tot->Add(h_charm);
    h_tot->Add(h_bottom);
    h_tot->Add(h_fake);
    h_tot->Add(h_unm);
    
    return h_tot;
  }

  else{
    TH1F* h_light = (TH1F*) f->Get("h_mu_obs_light");
    TH1F* h_charm = (TH1F*) f->Get("h_mu_obs_charm");
    TH1F* h_bottom = (TH1F*) f->Get("h_mu_obs_bottom");
    TH1F* h_fake = (TH1F*) f->Get("h_mu_obs_fake");
    TH1F* h_unm = (TH1F*) f->Get("h_mu_obs_unm");
    TH1F* h_tot = (TH1F*) h_light->Clone();
    h_tot->Add(h_charm);
    h_tot->Add(h_bottom);
    h_tot->Add(h_fake);
    h_tot->Add(h_unm);
    
    return h_tot;
  }


};


//THIS IS THE 'REAL' MACRO
void Ex_2p3(){

  gStyle->SetOptStat(kFALSE);

  setTDRStyle();
  
  TFile* f = new TFile("/uscms_data/d3/clint/public/SmartClosure_TTbar.root");
  
  DrawAndSave((TH1F*) f->Get("h_el_pred_light"), (TH1F*) f->Get("h_el_obs_light"), "Closure_Light_Electrons.pdf","light","Electrons");
  DrawAndSave((TH1F*) f->Get("h_el_pred_charm"), (TH1F*) f->Get("h_el_obs_charm"), "Closure_Charm_Electrons.pdf","charm","Electrons");
  DrawAndSave((TH1F*) f->Get("h_el_pred_bottom"), (TH1F*) f->Get("h_el_obs_bottom"), "Closure_Bottom_Electrons.pdf","bottom","Electrons");
  DrawAndSave((TH1F*) f->Get("h_el_pred_fake"), (TH1F*) f->Get("h_el_obs_fake"), "Closure_Fake_Electrons.pdf","fake","Electrons");
  DrawAndSave((TH1F*) f->Get("h_el_pred_unm"), (TH1F*) f->Get("h_el_obs_unm"), "Closure_UnMatched_Electrons.pdf","unmatched","Electrons");


  DrawAndSave((TH1F*) f->Get("h_mu_pred_light"), (TH1F*) f->Get("h_mu_obs_light"), "Closure_Light_Muons.pdf","light","Muons");
  DrawAndSave((TH1F*) f->Get("h_mu_pred_charm"), (TH1F*) f->Get("h_mu_obs_charm"), "Closure_Charm_Muons.pdf","charm","Muons");
  DrawAndSave((TH1F*) f->Get("h_mu_pred_bottom"), (TH1F*) f->Get("h_mu_obs_bottom"), "Closure_Bottom_Muons.pdf","bottom","Muons");
  DrawAndSave((TH1F*) f->Get("h_mu_pred_fake"), (TH1F*) f->Get("h_mu_obs_fake"), "Closure_Fake_Muons.pdf","fake","Muons");
  DrawAndSave((TH1F*) f->Get("h_mu_pred_unm"), (TH1F*) f->Get("h_mu_obs_unm"), "Closure_UnMatched_Muons.pdf","unmatched","Muons");

  f->Close();
  
  TFile* f1 = new TFile("/uscms_data/d3/clint/public/SmartClosure_TTbar.root");

  TH1F* h_pred_tot_el = getTotPredHist(f1,"el");
  TH1F* h_obs_tot_el = getTotObsHist(f1,"el");

  TH1F* h_pred_tot_mu = getTotPredHist(f1,"mu");
  TH1F* h_obs_tot_mu = getTotObsHist(f1,"mu");





}



