void HistPlotter () {
//load the files we are using
  TFile *fbkg1 = new TFile("bg_chargeMisID.root");
  TFile *fbkg2 = new TFile("bg_nonPrompt.root");
  TFile *fsig1 = new TFile("sig_HT.root");
//get the ttrees from the files
  //TTree *Tbkg1 = (TTree*)fbkg1->Get("T");
  //TTree *Tbkg2 = (TTree*)fbkg2->Get("T");
  //TTree *Tsig1 = (TTree*)fsig1->Get("T");
//initialize the histograms
  TH1F *histSSbkg1 = fbkg1->Get("ssHTHist");
  histSSbkg1->Scale(11.148);
  TH1F *histSSbkg2 = fbkg2->Get("HT_nonPrompt");
  histSSbkg2->Scale(11.148);
  TH1F *histSSsig1 = fsig1->Get("HT_h");
  histSSsig1->Scale( 0.0001666666667);
//assign the branches to the initialized histograms
  //Tbkg1->SetBranchAddress("ssHTHist",&histSSbkg1);
  //Tbkg2->SetBranchAddress("ssHTHist",&histSSbkg2);
  //Tsig1->SetBranchAddress("ssHTHist",&histSSsig1);
//get the entries
  // Tbkg1->GetEntry();
  //Tbkg2->GetEntry();
  //Tsig1->GetEntry();
//create the tgraphs for drawing errors (THStack can't do this)
  TGraphErrors *siggerr = new TGraphErrors(histSSsig1);
  TH1F *histSSbkg = (TH1F*)histSSbkg1->Clone("histSSbkg");
  histSSbkg->Add(histSSbkg2);
  TGraphErrors *bkggerr = new TGraphErrors(histSSbkg);
//stack the histograms for the background
  THStack *stackSSbkg = new THStack("stackSSbkg","Stacked background histograms");
  stackSSbkg->Add(histSSbkg1);
 stackSSbkg->Add(histSSbkg2);
//make things pretty
  histSSbkg1->SetLineColor(kBlue);
  histSSbkg1->SetFillColor(kBlue-2);
  histSSbkg1->SetLineWidth(2);
  histSSbkg2->SetLineColor(kRed);
  histSSbkg2->SetFillColor(kRed+2);
  histSSbkg2->SetLineWidth(2);
  histSSsig1->SetLineColor(kGreen+3);
  histSSsig1->SetFillStyle(0);
  histSSsig1->SetLineWidth(2);
  siggerr->SetFillStyle(3005);
  siggerr->SetFillColor(kGreen+1);
  bkggerr->SetFillStyle(3004);
  bkggerr->SetFillColor(kBlack);
//remove the stats box
  gStyle->SetOptStat(0);
  histSSsig1->SetStats(kFALSE);
//initialize the canvas
  TCanvas *c1 = new TCanvas("c1","SS Can",10,10,1200,1000);
  //draw everything
  stackSSbkg->Draw();
  stackSSbkg->GetYaxis()->SetRangeUser(0,500);
  stackSSbkg->Draw();
  histSSsig1->Draw("same");
  siggerr->Draw("sameP2");
  bkggerr->Draw("sameP2");
//setup and draw the legend
  leg = new TLegend(0.5,0.7,0.9,0.9);
  leg->SetHeader("Sample");
  leg->AddEntry(histSSsig1,"Signal","f");
  leg->AddEntry(histSSbkg1,"charge mis-ID","f");
  leg->AddEntry(histSSbkg2,"Non-Prompt 2","f");
  leg->Draw("same");
  c1.SetLogy();
  c1->Print("Final_HTPlot.pdf");
}
