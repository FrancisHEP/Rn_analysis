/************************************************* 
Copyright:    PandaX-II Collaboration 
Author:       Andi Tan
Description:  This script aimed to draw Rn spectrum
              before and after Rn220 injection
Change Log:
              2017-11-21: Created src code.
**************************************************/

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include <TF1.h>
#include "TROOT.h"
#include "TBrowser.h"

void XY_Dist(int &RunNo)
{
  TChain *chain;
  chain = new TChain("signalTree");
  chain->Add(Form("ana2/run%d.root",RunNo));
  double eLife = 70;
  //qS1Top Correction on Z with respect to dt = 180us
  TF1 *f1 = new TF1("f1","pol4",0,360);
  double parT[5] = {19129.4,8.2409,-0.38895,0.00136134,-1.63133e-06};
  f1->SetParameters(parT);
  double mean_qS1T = f1->Eval(180);
  chain->SetAlias("qS1Tc",Form("qS1T[s1max]*%f/(%f+%f*(tS2[s2max]-tS1[s1max])/1000+%f*pow((tS2[s2max]-tS1[s1max])/1000,2)+%f*pow((tS2[s2max]-tS1[s1max])/10000,3)+%f*pow((tS2[s2max]-tS1[s1max])/10000,4))",mean_qS1T,parT[0],parT[1],parT[2],parT[3]*1e3,parT[4]*1e4));
  //qS1Bottom Correction on Z with respect to dt = 180us
  TF1 *f2 = new TF1("f2","pol5",0,360);
  double parB[6] = {26411,28.6145,0.336307,-0.00293665,1.53159e-05,-2.65163e-08};
  f2->SetParameters(parB);
  double mean_qS1B = f2->Eval(180);
  chain->SetAlias("qS1Bc",Form("qS1B[s1max]*%f/(%f+%f*(tS2[s2max]-tS1[s1max])/1000+%f*pow((tS2[s2max]-tS1[s1max])/1000,2)+%f*pow((tS2[s2max]-tS1[s1max])/10000,3)+%f*pow((tS2[s2max]-tS1[s1max])/10000,4)+%lf*pow((tS2[s2max]-tS1[s1max])/10000,5))",mean_qS1B,parB[0],parB[1],parB[2],parB[3]*1e3,parB[4]*1e4,parB[5]*1e5));
  //qS1Tc+qS1Bc correction on R2 with respect to R = 0
  TF1 *f3 = new TF1("f3","pol2",0,1200);
  double parR[3] = {53309,-3.88774,-0.000235078};
  f3->SetParameters(parR);
  chain->SetAlias("qS1c",Form("(qS1Tc+qS1Bc)*%f/(%f+%f*(yS2LRF[s2max]*yS2LRF[s2max]+xS2LRF[s2max]*xS2LRF[s2max])/100+%f*pow((yS2LRF[s2max]*yS2LRF[s2max]+xS2LRF[s2max]*xS2LRF[s2max])/100,2))",parR[0],parR[0],parR[1],parR[2]));
  chain->SetAlias("qS2c",Form("qS2[s2max]*exp((tS2[s2max]-tS1[s1max])/1000/%f)",eLife));
  //Define energy based on Rn220 alpha S1 at R = 0
  chain->SetAlias("E",Form("qS1c/%f*6.41*1.0056",parR[0]));
  
  /*
  TF1 *f1 = new TF1("f1","pol4",0,360);
  double parT[5] = {16819.9,15.2751,-0.446427,0.00166845,-2.15643e-06};
  f1->SetParameters(parT);
  double mean_qS1T = f1->Eval(180);
  chain->SetAlias("qS1Tc",Form("qS1T[s1max]*%f/(%f+%f*(tS2[s2max]-tS1[s1max])/1000+%f*pow((tS2[s2max]-tS1[s1max])/1000,2)+%f*pow((tS2[s2max]-tS1[s1max])/10000,3)+%f*pow((tS2[s2max]-tS1[s1max])/10000,4))",mean_qS1T,parT[0],parT[1],parT[2],parT[3]*1e3,parT[4]*1e4));
  double parB[5] = {23111.1,30.6084,0.0682929,0.000373337,-6.54284e-07};
  double mean_qS1T = f1->Eval(180);
  f1->SetParameters(parB);
  double mean_qS1B = f1->Eval(180);
  chain->SetAlias("qS1Bc",Form("qS1B[s1max]*%f/(%f+%f*(tS2[s2max]-tS1[s1max])/1000+%f*pow((tS2[s2max]-tS1[s1max])/1000,2)+%f*pow((tS2[s2max]-tS1[s1max])/10000,3)+%f*pow((tS2[s2max]-tS1[s1max])/10000,4))",mean_qS1B,parB[0],parB[1],parB[2],parB[3]*1e3,parB[4]*1e4));
  TF1 *f2 = new TF1("f2","pol2",0,120000);
  double parR[3] = {46653.3,-3.8,-9.7747e-04};
  f2->SetParameters(parR);
  chain->SetAlias("qS1c",Form("(qS1Tc+qS1Bc)*%f/(%f+%f*(yS2NN[s2max]*yS2NN[s2max]+xS2NN[s2max]*xS2NN[s2max])/100+%f*pow((yS2NN[s2max]*yS2NN[s2max]+xS2NN[s2max]*xS2NN[s2max])/100,2))",parR[0],parR[0],parR[1],parR[2]));
  chain->SetAlias("qS2c",Form("qS2[s2max]*exp((tS2[s2max]-tS1[s1max])/1000/%f)",eLife));
  chain->SetAlias("E",Form("qS1c/%f*5.59*1.012",parR[0]));
  
  */
  
  
  
  //TCanvas *c1 = new TCanvas("c1","c1",900,900);
  TCanvas *c1 = new TCanvas("c1","c1",1080,720);
  gStyle->SetOptStat(0);
  /*
  TCutG *cutg = new TCutG("mycut",12);
  cutg->SetVarX("y");
  cutg->SetVarY("x");
  float x[12] = {0};
  float y[12] = {0};
  for(int i = 0; i < 12; i++){
    x[i] = 335.43*TMath::Cos(TMath::Pi()/12+TMath::Pi()/6*i);
    y[i] = 335.43*TMath::Sin(TMath::Pi()/12+TMath::Pi()/6*i);
    cutg->SetPoint(i,x[i],y[i]);
  }
  cutg->SetPoint(12,x[0],y[0]);
  for(int i = 0; i < 12; i++){
    TH2F *hXY=new TH2F("hXY","",680,-340,340,680,-340,340);
    TCut FVCut=Form("(tS2[s2max]-tS1[s1max])/1000>%d&&(tS2[s2max]-tS1[s1max])/1000<=%d",i*30,(i+1)*30);
    chain->Draw("yS2LRF[s2max]:xS2LRF[s2max]>>hXY","s2max>=0&&s1max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<140&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=40&&hS1[s1max]/qS1[s1max]>0.16&&S1NPeaks[s1max]<=1&&qS2c>50&&log10(qS2c/qS1[s1max])>-0.4&&E>6.2&&E<6.4"&&FVCut,"");
    hXY->GetXaxis()->SetTitle("X [mm]");
    hXY->GetXaxis()->CenterTitle();
    hXY->GetYaxis()->SetTitle("Y [mm]");
    hXY->GetYaxis()->CenterTitle();
    hXY->SetMarkerStyle(20);
    hXY->SetMarkerColor(kBlue);
    cutg->Draw("same");
    c1->SaveAs(Form("plots/XY_dist_%d.pdf",i));
    delete hXY;
    }
  */
  TH2F *hVert=new TH2F("hVert","",300,4,10,300,-360,0);
  //chain->Draw("-(tS2[s2max]-tS1[s1max])/1000:E>>hVert","s2max>=0&&s1max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<140&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=40&&hS1[s1max]/qS1[s1max]>0.16&&S1NPeaks[s1max]<=1&&qS2c>0&&log10(qS2c/qS1[s1max])>-1","");
  chain->Draw("-(tS2[s2max]-tS1[s1max])/1000:E>>hVert","s2max>=0&&s1max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<140&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=40&&hS1[s1max]/qS1[s1max]>0.16&&S1NPeaks[s1max]<=1&&qS2c>50&&log10(qS2c/qS1[s1max])>-0.4","");
  hVert->GetXaxis()->SetTitle("E [MeV]");
  hVert->GetXaxis()->CenterTitle();
  hVert->GetYaxis()->SetTitle("-#Deltat [#mus]");
  hVert->GetYaxis()->CenterTitle();
  hVert->SetMarkerColor(kBlue);
  c1->SetGridx();
  /*
  c1->SaveAs("plots/Z_corr.pdf");
  c1->SetGridx(0);
  TH2F *hR=new TH2F("hR","",300,0,1200,300,4,10);
  chain->Draw("E:(yS2LRF[s2max]*yS2LRF[s2max]+xS2LRF[s2max]*xS2LRF[s2max])/100>>hR","s2max>=0&&s1max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<140&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=40&&hS1[s1max]/qS1[s1max]>0.16&&S1NPeaks[s1max]<=1&&qS2c>50&&log10(qS2c/qS1[s1max])>-0.4","");
  hR->GetXaxis()->SetTitle("R^{2} [cm^{2}]");
  hR->GetXaxis()->CenterTitle();
  hR->GetYaxis()->SetTitle("E [MeV]");
  hR->GetYaxis()->CenterTitle();
  hR->SetMarkerColor(kBlue);
  c1->SetGridy();
  c1->SaveAs("plots/R_corr.pdf");
  TH2F *hZX = new TH2F("hZX","",300,-340,340,300,-360,0);
  chain->Draw("-(tS2[s2max]-tS1[s1max])/1000:xS2LRF[s2max]>>hZX","s2max>=0&&s1max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<140&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=40&&hS1[s1max]/qS1[s1max]>0.16&&S1NPeaks[s1max]<=1&&qS2c>50&&log10(qS2c/qS1[s1max])>-0.4","");
  hZX->GetXaxis()->SetTitle("X [mm]");
  hZX->GetXaxis()->CenterTitle();
  hZX->GetYaxis()->SetTitle("-#Deltat [#mus]");
  hZX->GetYaxis()->CenterTitle();
  hZX->SetMarkerStyle(20);
  hZX->SetMarkerColor(kBlue);
  c1->SaveAs("plots/ZX_dist.pdf");
  TH2F *hZY = new TH2F("hZY","",300,-340,340,300,-360,0);
  chain->Draw("-(tS2[s2max]-tS1[s1max])/1000:yS2LRF[s2max]>>hZY","s2max>=0&&s1max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<140&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=40&&hS1[s1max]/qS1[s1max]>0.16&&S1NPeaks[s1max]<=1&&qS2c>50&&log10(qS2c/qS1[s1max])>-0.4","");
  hZY->GetXaxis()->SetTitle("Y [mm]");
  hZY->GetXaxis()->CenterTitle();
  hZY->GetYaxis()->SetTitle("-#Deltat [#mus]");
  hZY->GetYaxis()->CenterTitle();
  hZY->SetMarkerStyle(20);
  hZY->SetMarkerColor(kBlue);
  c1->SaveAs("plots/ZY_dist.pdf");
  TH2F *hAlpha = new TH2F("hAlpha","",400,1e4,9e4,400,-1,0.2);
  chain->Draw("(qS1T[s1max]-qS1B[s1max])/qS1[s1max]:qS1[s1max]>>hAlpha","s2max>=0&&s1max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<140&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=40&&hS1[s1max]/qS1[s1max]>0.16&&S1NPeaks[s1max]<=1&&qS2c>50&&log10(qS2c/qS1[s1max])>-0.4","");
  hAlpha->SetMarkerStyle(7);
  hAlpha->SetMarkerColor(kBlue);
  */

}
