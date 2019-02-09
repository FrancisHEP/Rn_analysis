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

void aa_Velocity(int& RunNo)
{
  TChain *ch;
  ch = new TChain("aa_tree");
  //ch->Add("tree_Th_aa_30SLPM.root");
  ch->Add(Form("tree_Th_aa_run%d.root",RunNo));
  int totalEvt;
  totalEvt = ch->GetEntries();
  cout << "totalEvents: " << totalEvt << endl;
  int runNo, fileNo, trigNo;
  double decayT;
  double xRn, yRn, xPo, yPo, dTRn, dTPo;
  double ERn, EPo;
  ch->SetBranchAddress("runNo",&runNo);
  ch->SetBranchAddress("fileNo",&fileNo);
  ch->SetBranchAddress("trigNo",&trigNo);
  ch->SetBranchAddress("decayT",&decayT);
  ch->SetBranchAddress("xRn",&xRn);
  ch->SetBranchAddress("yRn",&yRn);
  ch->SetBranchAddress("xPo",&xPo);
  ch->SetBranchAddress("yPo",&yPo);
  ch->SetBranchAddress("dTRn",&dTRn);
  ch->SetBranchAddress("dTPo",&dTPo);
  ch->SetBranchAddress("ERn",&ERn);
  ch->SetBranchAddress("EPo",&EPo);
  const int Nbin = 10;
  int nbin = Nbin;
  double Scale = 2.;//arraw_display_ratio
  TH2F *h_zy = new TH2F("h_zy","Event and Velosity Distribution in Z-Y plane",Nbin,-324,324,Nbin,0,600);
  h_zy->GetXaxis()->CenterTitle();
  h_zy->GetXaxis()->SetTitleOffset(1.2);
  h_zy->GetXaxis()->SetTitle("Y [mm] / #nu_{y} [mm/s]");
  h_zy->GetYaxis()->CenterTitle();
  h_zy->GetYaxis()->SetTitleOffset(1.4);
  h_zy->GetYaxis()->SetTitle("Z [mm] / #nu_{z} [mm/s]");
  int N_zy[Nbin][Nbin] = {0};
  double v_zy_y[Nbin][Nbin] = {0};
  double v_zy_z[Nbin][Nbin] = {0};
  TH2F *h_zx = new TH2F("h_zx","Event and Velosity Distribution in Z-X plane",Nbin,-324,324,Nbin,0,600);
  h_zx->GetXaxis()->CenterTitle();
  h_zx->GetXaxis()->SetTitleOffset(1.2);
  h_zx->GetXaxis()->SetTitle("X [mm] / #nu_{x} [mm/s]");
  h_zx->GetYaxis()->CenterTitle();
  h_zx->GetYaxis()->SetTitleOffset(1.4);
  h_zx->GetYaxis()->SetTitle("Z [mm] / #nu_{z} [mm/s]");
  int N_zx[Nbin][Nbin] = {0};
  double v_zx_x[Nbin][Nbin] = {0};
  double v_zx_z[Nbin][Nbin] = {0};
  TH2F *h_yx[12];
  for(int i = 0; i<12; i++){
    TString hpnt = Form("h_yx[%d]",i);
    h_yx[i] = new TH2F(hpnt,"",Nbin,-324,324,Nbin,-324,324);
    h_yx[i]->GetXaxis()->CenterTitle();
    h_yx[i]->GetXaxis()->SetTitleOffset(1.2);
    h_yx[i]->GetXaxis()->SetTitle("X [mm] / #nu_{x} [mm/s]");
    h_yx[i]->GetYaxis()->CenterTitle();
    h_yx[i]->GetYaxis()->SetTitleOffset(1.4);
    h_yx[i]->GetYaxis()->SetTitle("Y [mm] / #nu_{y} [mm/s]");
  }
  int N_yx[12][Nbin][Nbin] = {0};
  double v_yx_x[12][Nbin][Nbin] = {0};
  double v_yx_y[12][Nbin][Nbin] = {0};
  
  for(int j = 0 ; j < totalEvt; j++){
    ch->GetEntry(j);
    double v_x = (xPo - xRn)/decayT ;
    double v_y = (yPo - yRn)/decayT ;
    double v_z = (-dTPo+dTRn)/360.*600/decayT ;
    if(abs(v_x)>20||abs(v_y)>20||abs(v_z)>100) continue;
    v_x *= Scale;
    v_y *= Scale;
    v_z *= Scale;
    double zRn = (360-dTRn)/360.*600;
    int index_x = (xRn+324)/(648./Nbin);
    int index_y = (yRn+324)/(648./Nbin);
    int index_z = zRn/(600./Nbin);
    if(index_x==Nbin) index_x-=1;
    if(index_y==Nbin) index_y-=1;
    if(index_z==Nbin) index_z-=1;
    int layer_z = zRn/50;
    if(layer_z==12) layer_z-=1;
    N_zy[index_y][index_z] += 1;
    v_zy_y[index_y][index_z] += v_y;
    v_zy_z[index_y][index_z] += v_z;
    N_zx[index_x][index_z] += 1;
    v_zx_x[index_x][index_z] += v_x;
    v_zx_z[index_x][index_z] += v_z;
    N_yx[layer_z][index_x][index_y] += 1;
    v_yx_x[layer_z][index_x][index_y] += v_x;
    v_yx_y[layer_z][index_x][index_y] += v_y;
  }
  
  TCanvas *c1 = new TCanvas("c1","c1",2160,1080);
  gStyle->SetOptStat(0);
  gStyle->SetPalette(54);
  c1->Divide(2,1);
  c1->cd(1);
  ch->Draw("(360-dTRn)/360.*600:yRn>>h_zy","","colz");
  c1->cd(2);
  ch->Draw("(360-dTRn)/360.*600:xRn>>h_zx","","colz");
  TCanvas *c2 = new TCanvas("c2","c2",1440,1080);
  c2->Divide(4,3);
  for(int i = 0; i<12; i++){
    c2->cd(i+1);
    ch->Draw(Form("yRn:xRn>>h_yx[%d]",i),Form("dTRn>=%d*30&&dTRn<(%d+1)*30",11-i,11-i),"colz");
  }
  TArrow *arr_zy[Nbin][Nbin];
  TArrow *arr_zx[Nbin][Nbin];
  TArrow *arr_yx[12][Nbin][Nbin];
  for(int ii = 0; ii < nbin; ii++){
    for(int jj = 0; jj < nbin; jj++){
      if(N_zy[ii][jj]>16){
        c1->cd(1);
        double y_i = h_zy->GetXaxis()->GetBinCenter(ii+1);
        double z_i = h_zy->GetYaxis()->GetBinCenter(jj+1);
        arr_zy[ii][jj] = new TArrow(y_i,z_i,y_i+v_zy_y[ii][jj]/N_zy[ii][jj],z_i+v_zy_z[ii][jj]/N_zy[ii][jj],0.004,">");
        arr_zy[ii][jj]->SetAngle(40);
        arr_zy[ii][jj]->SetOption("|>");
        arr_zy[ii][jj]->SetLineColor(kRed);
        arr_zy[ii][jj]->SetLineWidth(1);
        arr_zy[ii][jj]->SetFillColor(kRed);
        arr_zy[ii][jj]->DrawArrow(y_i,z_i,y_i+v_zy_y[ii][jj]/N_zy[ii][jj],z_i+v_zy_z[ii][jj]/N_zy[ii][jj],0.004,"|>");
      }
      if(N_zx[ii][jj]>16){
        c1->cd(2);
        double x_i = h_zx->GetXaxis()->GetBinCenter(ii+1);
        double z_i = h_zx->GetYaxis()->GetBinCenter(jj+1);
        arr_zx[ii][jj] = new TArrow(x_i,z_i,x_i+v_zx_x[ii][jj]/N_zx[ii][jj],z_i+v_zx_z[ii][jj]/N_zx[ii][jj],0.001,">");
        arr_zx[ii][jj]->SetAngle(40);
        arr_zx[ii][jj]->SetOption("|>");
        arr_zx[ii][jj]->SetLineColor(kRed);
        arr_zx[ii][jj]->SetLineWidth(1);
        arr_zx[ii][jj]->SetFillColor(kRed);
        arr_zx[ii][jj]->DrawArrow(x_i,z_i,x_i+v_zx_x[ii][jj]/N_zx[ii][jj],z_i+v_zx_z[ii][jj]/N_zx[ii][jj],0.004,"|>");
      }
      for(int kk = 0; kk < 12; kk++)
      {
        if(N_yx[kk][ii][jj]>4){
          c2->cd(kk+1);
          double x_i = h_yx[kk]->GetXaxis()->GetBinCenter(ii+1);
          double y_i = h_yx[kk]->GetYaxis()->GetBinCenter(jj+1);
          arr_yx[kk][ii][jj] = new TArrow(x_i,y_i,x_i+v_yx_x[kk][ii][jj]/N_yx[kk][ii][jj],y_i+v_yx_y[kk][ii][jj]/N_yx[kk][ii][jj],0.001,">");
          arr_yx[kk][ii][jj]->SetAngle(40);
          arr_yx[kk][ii][jj]->SetOption("|>");
          arr_yx[kk][ii][jj]->SetLineColor(kRed);
          arr_yx[kk][ii][jj]->SetLineWidth(1);
          arr_yx[kk][ii][jj]->SetFillColor(kRed);
          arr_yx[kk][ii][jj]->DrawArrow(x_i,y_i,x_i+v_yx_x[kk][ii][jj]/N_yx[kk][ii][jj],y_i+v_yx_y[kk][ii][jj]/N_yx[kk][ii][jj],0.001,">");
        }
      }
    } 
  }
  //c1->Print(Form("/home/andy/Analysis/2017_Kr_Rn/plots/RnPo_Velocity_Z_%dSLPM.pdf",30));
  //  c2->Print(Form("/home/andy/Analysis/2017_Kr_Rn/plots/RnPo_Velocity_XY_%dSLPM.pdf",65));
}
