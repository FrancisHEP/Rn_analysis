// Wenbo
// 2019-02-26
// The compile-able version of check_ba212.C
// 2019-02-12
// This is the better version of Plot_ba212.C
// But I just use it to check 24467&24468. 
// 2019-01-29
// update TBA[beta]:TBA[alpha] check.

#include "TLatex.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TF1.h"
#include "TH2.h"
#include "TH1.h"
#include "TChain.h"
#include "TString.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

int check_ba212_2446724468()
{
  gStyle->SetPalette(57); // kBird=57
  gStyle->SetNumberContours(100);
  TChain * chain = new TChain("rn220_BA");
  TString file24467 = "ba212/outTree_212_BiPo_24467.root";
  TString file24468 = "ba212/outTree_212_BiPo_24468.root";
  chain->Add(file24467);
  chain->Add(file24468);

  int Time = 6569642;
  int runNo, trigNo;
  double triggerTime;

  unsigned int nS1, nS2;
  int s1max, s2max;
  double qS1T[1000], qS1B[1000], xS1T[1000], yS1T[1000];
  double xS1B[1000], yS1B[1000];
  double qS1[1000], tS1[1000], pS1[1000], hS1[1000];
  double ttenS1[1000], wS1[1000], wtenS1[1000];
  double xS2PAF[1000], yS2PAF[1000];
  double qS1Veto[1000];
  int S1NPeaks[1000], nS1BSatur[1000], nS1TSatur[1000];

  double qS2T[1000], qS2B[1000], xS2T[1000], yS2T[1000];
  double xS2B[1000], yS2B[1000];
  double qS2[1000], tS2[1000], pS2[1000], hS2[1000];
  double ttenS2[1000], wS2[1000], wtenS2[1000];
  double xS2NN[1000], yS2NN[1000], chi2[1000];
  int nS2BSatur[1000], nS2TSatur[1000];
  int beta, alpha;

  chain->SetBranchAddress("trigNo", &trigNo);
  chain->SetBranchAddress("runNo", &runNo);
  chain->SetBranchAddress("triggerTime", &triggerTime);

  chain->SetBranchAddress("nS1", &nS1);
  chain->SetBranchAddress("s1max", &s1max);
  chain->SetBranchAddress("xS1B", xS1B);
  chain->SetBranchAddress("yS1B", yS1B);
  chain->SetBranchAddress("xS1T", xS1T);
  chain->SetBranchAddress("yS1T", yS1T);
  chain->SetBranchAddress("qS1T", qS1T);
  chain->SetBranchAddress("qS1B", qS1B);
  chain->SetBranchAddress("qS1", qS1);
  chain->SetBranchAddress("tS1", tS1);
  chain->SetBranchAddress("ttenS1", ttenS1);
  chain->SetBranchAddress("wS1", wS1);
  chain->SetBranchAddress("wtenS1", wtenS1);
  chain->SetBranchAddress("hS1", hS1);
  chain->SetBranchAddress("pS1", pS1);
  chain->SetBranchAddress("qS1Veto", qS1Veto);
  chain->SetBranchAddress("S1NPeaks", S1NPeaks);
  chain->SetBranchAddress("nS1BSatur", nS1BSatur);
  chain->SetBranchAddress("nS1TSatur", nS1TSatur);

  chain->SetBranchAddress("nS2", &nS2);
  chain->SetBranchAddress("s2max", &s2max);
  chain->SetBranchAddress("xS2T", xS2T);
  chain->SetBranchAddress("yS2T", yS2T);
  chain->SetBranchAddress("xS2B", xS2B);
  chain->SetBranchAddress("yS2B", yS2B);
  chain->SetBranchAddress("xS2NN", xS2NN);
  chain->SetBranchAddress("xS2PAF", xS2PAF);
  chain->SetBranchAddress("yS2NN", yS2NN);
  chain->SetBranchAddress("yS2PAF", yS2PAF);
  chain->SetBranchAddress("chi2", chi2);
  chain->SetBranchAddress("qS2T", qS2T);
  chain->SetBranchAddress("qS2B", qS2B);
  chain->SetBranchAddress("qS2", qS2);
  chain->SetBranchAddress("tS2", tS2);
  chain->SetBranchAddress("ttenS2", ttenS2);
  chain->SetBranchAddress("wS2", wS2);
  chain->SetBranchAddress("wtenS2", wtenS2);
  chain->SetBranchAddress("hS2", hS2);
  chain->SetBranchAddress("pS2", pS2);
  chain->SetBranchAddress("beta", &beta);
  chain->SetBranchAddress("alpha", &alpha);

  TCut FV = "tS2[s2max]-tS1[beta]>20000&&tS2[s2max]-tS1[beta]<350000";
  TCut R2 = "xS2PAF[s2max]*xS2PAF[s2max]+yS2PAF[s2max]*yS2PAF[s2max]<72000";

  // Depth:TBA cuts
  TCanvas *c1=new TCanvas("c1","c1",520,520);
  TH2F * h1 = new TH2F("h1","",100,-1,0,100,-400,0);
  chain->Draw("(tS1[beta]-tS2[s2max])/1000:(qS1T[beta]-qS1B[beta])/(qS1T[beta]+qS1B[beta])>>h1",FV,"colz");
  h1->SetTitle("Depth:TBA^{#beta}");
  h1->GetXaxis()->SetTitle("TBA^{#beta}");
  h1->GetYaxis()->SetTitle("Depth[us]");
  h1->SetMarkerColor(2);
  h1->Draw("colz");


  // fit and draw Depth:TBA
  TF1 *f1=new TF1("f1","[0]+[1]*x",-1,0);
  f1->SetParameters(20,400);
  h1->Fit("f1","R");

  double a = f1->GetParameter(1);
  double b = f1->GetParameter(0);

  TString Range1 = Form("%.2f*x+%.2f+50",a,b);
  TString Range2 = Form("%.2f*x+%.2f-50",a,b);
  
  TF1 *f4=new TF1("f4",Range1,-1,0);
  f4->SetLineStyle(2);
  f4->SetLineColor(4);
  f4->Draw("same");
  TF1 *f5=new TF1("f5",Range2,-1,0);
  f5->SetLineStyle(2);
  f5->SetLineColor(4);
  f5->Draw("same");

  TCut Cut3 = Form("(-tS2[s2max]+tS1[beta])/1000<%.2f*(qS1T[beta]-qS1B[beta])/(qS1T[beta]+qS2B[beta])+%.4f+50&&(-tS2[s2max]+tS1[beta])/1000>%.2f*(qS1T[beta]-qS1B[beta])/(qS1T[beta]+qS1B[beta])+%.2f-50",a,b,a,b);

  // S1&S2 with Z&R2 cuts
  TCanvas * c3=new TCanvas("c3","c3",640,640);
  c3->Divide(2,2);
  TH2F * h3 = new TH2F("h3","",100,0,350,100,0,100000);
  c3->cd(1);
  // qS1:dT
  chain->Draw("qS1[alpha]:(tS2[s2max]-tS1[beta])/1000>>h3",FV&&Cut3,"");
  h3->SetStats(0);
  h3->SetTitle("S1:Z");
  h3->SetMarkerColor(2);
  h3->SetMarkerStyle(22);
  h3->SetMarkerSize(1);
  h3->Draw("colz");
  TF1 *f3=new TF1("f3","[0]+[1]*x",0,350);
  f3->SetParameters(40000,50);
  f3->SetLineColor(4);
  f3->SetLineStyle(2);
  f3->Draw("same");

  TCut s1z = Form("qS1[alpha]>%.2f*(tS2[s2max]-tS1[beta])/1000+%.2f",f3->GetParameter(1),f3->GetParameter(0));
  TH2F * h4 = new TH2F("h4","",100,0,100000,100,0,100000);
  c3->cd(2);
  // qS1:r^2
  chain->Draw("qS1[alpha]:xS2PAF[s2max]^2+yS2PAF[s2max]^2>>h4",FV&&""&&Cut3&&""&&s1z,"");
  h4->SetStats(0);
  h4->SetTitle("S1:R^{2}");
  h4->SetMarkerColor(2);
  h4->SetMarkerStyle(22);
  h4->SetMarkerSize(1);
  h4->Draw("colz");
  TH2F * h5 = new TH2F("h5","",100,0,350,100,0,7);
  c3->cd(3);
  // log10qS2:dT
  chain->Draw("log10(qS2[s2max]):(tS2[s2max]-tS1[beta])/1000>>h5",FV&&Cut3&&s1z,"");
  h5->SetStats(0);
  h5->SetTitle("S2:Z");
  h5->SetMarkerColor(2);
  h5->SetMarkerStyle(22);
  h5->SetMarkerSize(1);
  h5->Draw("colz");
  TF1 *f6=new TF1("f6","[0]+[1]*x",0,350);
  f6->SetParameters(4,0);
  f6->SetLineColor(4);
  f6->SetLineStyle(2);
  f6->Draw("same");
  TCut s2z = Form("log10(qS2[s2max])>%.2f*(tS2[s2max]-tS1[beta])/1000+%.2f",f6->GetParameter(1),f6->GetParameter(0));
  TH2F * h6 = new TH2F("h6","",100,0,100000,100,0,7);
  c3->cd(4);
  // log10qS2:r^2
  chain->Draw("log10(qS2[s2max]):xS2PAF[s2max]^2+yS2PAF[s2max]^2>>h6",FV&&Cut3&&s1z&&s2z,"");
  h6->SetStats(0);
  h6->SetTitle("S2:R^{2}");
  h6->SetMarkerColor(2);
  h6->SetMarkerStyle(22);
  h6->SetMarkerSize(1);
  h6->Draw("colz");

  // updated TBA[alpha]:TBA[beta] 2019-01-29
  TCanvas * c_TBA = new TCanvas("c_TBA","c_TBA",520,520);
  TH2F * h_TBA = new TH2F("h_TBA","",100,-1,0.2,100,-1,0.2);
  chain->Draw("(qS1T[alpha]-qS1B[alpha])/(qS1T[alpha]+qS1B[alpha]):(qS1T[beta]-qS1B[beta])/(qS1T[beta]+qS1B[beta])>>h_TBA",FV,"colz");
  h_TBA->SetTitle("TBA^{#alpha}:TBA^{#beta}");
  h_TBA->GetXaxis()->SetTitle("TBA^{#beta}");
  h_TBA->GetYaxis()->SetTitle("TBA^{#alpha}");
  h_TBA->SetMarkerColor(2);
  chain->Draw("colz");

  // Decay time check the result
  TCanvas *c2=new TCanvas("c2","c2",520,520);
  TH1F * h2 = new TH1F("h2","",100,0,5);
  chain->Draw("(tS1[alpha]-tS1[beta])/1000>>h2",FV&&R2&&Cut3&&s1z&&s2z,"");//&&"(tS2[s2max]-tS1[beta])/1000<250","");
  h2->SetTitle("Decay time");
  h2->GetXaxis()->SetTitle("Dt[us]");
  h2->GetYaxis()->SetTitle("Count");
  int counts = h2->GetEntries();
  cout << "Counts=" << counts << endl;
  TF1 *f2=new TF1("f2","[0]*exp(-x/[1])",0.9,5); // changed from 0.3 to 0.6
  f2->SetParameters(counts/20,0.3*1.442695);
  h2->Fit("f2","R");
  TLatex *t1=new TLatex(1,5,Form("Dt=%.2f#pm%.2f#mus",f2->GetParameter(1),f2->GetParError(1)));
  t1->SetTextColor(2);
  t1->SetTextSize(0.04);
  t1->Draw("same");
  TLatex *t2=new TLatex(1,6,Form("%.4f #muBq",double(counts)/0.4885/Time*1.e6));
  t2->SetTextColor(2);
  t2->SetTextSize(0.04);
  t2->Draw("same");

  c1->SaveAs("Bi212c1.pdf");
  c1->SaveAs("Bi212c1.root");
  c2->SaveAs("Bi212c2.pdf");
  c2->SaveAs("Bi212c2.root");
  c3->SaveAs("Bi212c3.pdf");
  c3->SaveAs("Bi212c3.root");

  return counts;
}
