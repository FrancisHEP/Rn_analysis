#include <iostream>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#include "TROOT.h"
#include "TBrowser.h"
#include <fstream>
#include <vector>
#include <algorithm>


#include <TChain.h>
#include <TFileCollection.h>

// usage: 
// load ana2 file in directory ana2/
// ch->Add(Form("ana2/run%d.root",RunNo));
// and create aatree.root in directory aacoin_tree
// TFile *outFile = new TFile(Form("rootfile_aa_tree/tree_Th_aa_run%d.root",RunNo),"recreate");


void Rn220_Po216(int &RunNo)
{
  TChain *ch;
  ch = new TChain("signalTree");
  ch->Add(Form("ana2/run%d.root",RunNo));

  int totalEvt;
  int runNo, fileNo, trigNo;
  unsigned int time;
  double triggerTime;
  unsigned int nS1, nS2, SingleS2;
  int s1max, s2max;
  double tS1[1000], ttenS1[1000], wS1[1000], wtenS1[1000], fwhmS1[1000], hS1[1000],pS1[1000];
  double tS2[1000];
  double qS1[1000], qS1T[1000], qS1B[1000], qS2[1000];
  double xS2LRF[1000], yS2LRF[1000], LRFchi2[1000];

  ch->SetBranchAddress("runNo",&runNo);
  ch->SetBranchAddress("fileNo",&fileNo);
  ch->SetBranchAddress("trigNo",&trigNo);
  ch->SetBranchAddress("time",&time);
  ch->SetBranchAddress("triggerTime",&triggerTime);
  ch->SetBranchAddress("nS1", &nS1);
  ch->SetBranchAddress("s1max", &s1max);
  ch->SetBranchAddress("nS2", &nS2);
  ch->SetBranchAddress("s2max", &s2max);
  ch->SetBranchAddress("SingleS2", &SingleS2);
  ch->SetBranchAddress("qS1", qS1);
  ch->SetBranchAddress("qS1B", qS1B);
  ch->SetBranchAddress("qS1T", qS1T);
  ch->SetBranchAddress("tS1", tS1);
  ch->SetBranchAddress("ttenS1", ttenS1);
  ch->SetBranchAddress("wS1", wS1);
  ch->SetBranchAddress("wtenS1", wtenS1);
  ch->SetBranchAddress("fwhmS1", fwhmS1);
  ch->SetBranchAddress("hS1", hS1);
  ch->SetBranchAddress("pS1", pS1);

  ch->SetBranchAddress("xS2LRF",xS2LRF);
  ch->SetBranchAddress("yS2LRF",yS2LRF);
  ch->SetBranchAddress("LRFchi2",LRFchi2);
  ch->SetBranchAddress("tS2",tS2);
  ch->SetBranchAddress("qS2",qS2);

  totalEvt = ch->GetEntries();
  cout << "totalEvents: " << totalEvt << endl;
  
  int runno, fileno, evtno;
  double decayT;
  double xRn, yRn, xPo, yPo, dTRn, dTPo;
  double ERn, EPo;
  TFile *outFile = new TFile(Form("rootfile_aa_tree/tree_Th_aa_run%d.root",RunNo),"recreate");
  TTree *outTree = new TTree("aa_tree","radon220_AA_candidate_tree");
  outTree->Branch("runNo",&runno,"runNo/I");
  outTree->Branch("fileNo",&fileno,"fileNo/I");
  outTree->Branch("trigNo",&evtno,"trigNo/I");
  outTree->Branch("decayT",&decayT,"decayT/D");
  outTree->Branch("xRn",&xRn,"xRn/D");
  outTree->Branch("yRn",&yRn,"yRn/D");
  outTree->Branch("xPo",&xPo,"xPo/D");
  outTree->Branch("yPo",&yPo,"yPo/D");
  outTree->Branch("dTRn",&dTRn,"dTRn/D");
  outTree->Branch("dTPo",&dTPo,"dTPo/D");
  outTree->Branch("ERn",&ERn,"ERn/D");
  outTree->Branch("EPo",&EPo,"EPo/D");
  
  TF1 *f1 = new TF1("f1","pol4",0,360);
  double parT[5] = {19129.4,8.2409,-0.38895,0.00136134,-1.63133e-06};
  f1->SetParameters(parT);
  double mean_qS1T = f1->Eval(180);
  //qS1Bottom Correction on Z with respect to dt = 180us
  TF1 *f2 = new TF1("f2","pol5",0,360);
  double parB[6] = {26411,28.6145,0.336307,-0.00293665,1.53159e-05,-2.65163e-08};
  f2->SetParameters(parB);
  double mean_qS1B = f2->Eval(180);
  //qS1Tc+qS1Bc correction on R2 with respect to R = 0
  TF1 *f3 = new TF1("f3","pol2",0,1200);
  double parR[3] = {53309,-3.88774,-0.000235078};
  f3->SetParameters(parR);

  int N_T_half = 10;
  double dR_xy = 400;
  double qS1Tc, qS1Bc, qS1c, E;
  //totalEvt = 100;
  for (int ii = 0; ii < totalEvt; ii++){
    ch->GetEntry(ii);
    if((ii+1)%100000==0) cout << ii*1.00/totalEvt*100 <<"%"<< endl;
    if(nS1>1000||s1max<0||s2max<0) continue;
    if(hS1[s1max]<=2000||ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=140||ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<40||hS1[s1max]/qS1[s1max]<=0.16||qS2[s2max]<=50) continue;//pulse shape
    //if(tS2-tS1>270e3||tS2-tS1<50e3) continue;
    qS1Tc = qS1T[s1max]*mean_qS1T/(parT[0]+parT[1]*(tS2[s2max]-tS1[s1max])/1000+parT[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parT[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parT[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4));
    qS1Bc = qS1B[s1max]*mean_qS1B/(parB[0]+parB[1]*(tS2[s2max]-tS1[s1max])/1000+parB[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parB[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parB[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4)+parB[5]*1e5*pow((tS2[s2max]-tS1[s1max])/10000,5));
    qS1c = (qS1Tc+qS1Bc)*parR[0]/(parR[0]+parR[1]*(yS2LRF[s2max]*yS2LRF[s2max]+xS2LRF[s2max]*xS2LRF[s2max])/100+parR[2]*pow((yS2LRF[s2max]*yS2LRF[s2max]+xS2LRF[s2max]*xS2LRF[s2max])/100,2));
    E = qS1c/parR[0]*6.41*1.0056;
    //cout << trigNo << '\t' << qS1Tc <<  '\t' << qS1Bc <<  '\t' << E << endl;
    if(E<6.15||E>6.65) continue;
    double startT = time + triggerTime;
    runno = runNo;
    fileno = fileNo;
    evtno = trigNo;
    xRn = xS2LRF[s2max];
    yRn = yS2LRF[s2max];
    dTRn = (tS2[s2max]-tS1[s1max])/1000;
    ERn = E;
    for (int jj = ii+1; jj < totalEvt; jj++){
      ch->GetEntry(jj);
      if(nS1>1000||s1max<0||s2max<0) continue;
      if(time + triggerTime > startT + 0.15*N_T_half) break;
      //if(xTopLRF*xTopLRF+yTopLRF*yTopLRF>80000) continue;
      //if(tS2-tS1>270e3||tS2-tS1<50e3) continue;
      qS1Tc = qS1T[s1max]*mean_qS1T/(parT[0]+parT[1]*(tS2[s2max]-tS1[s1max])/1000+parT[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parT[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parT[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4));
      qS1Bc = qS1B[s1max]*mean_qS1B/(parB[0]+parB[1]*(tS2[s2max]-tS1[s1max])/1000+parB[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parB[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parB[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4)+parB[5]*1e5*pow((tS2[s2max]-tS1[s1max])/10000,5));
      qS1c = (qS1Tc+qS1Bc)*parR[0]/(parR[0]+parR[1]*(yS2LRF[s2max]*yS2LRF[s2max]+xS2LRF[s2max]*xS2LRF[s2max])/100+parR[2]*pow((yS2LRF[s2max]*yS2LRF[s2max]+xS2LRF[s2max]*xS2LRF[s2max])/100,2));
      E = qS1c/parR[0]*6.41*1.0056;
      if(E<6.65||E>7.5) continue;
      //cout << trigNo << '\t' << qS1Tc <<  '\t' << qS1Bc <<  '\t' << E << endl;
      if(pow(xS2LRF[s2max]-xRn,2)+pow(yS2LRF[s2max]-yRn,2) > 324*324./dR_xy)continue;
      decayT = time + triggerTime - startT;
      xPo = xS2LRF[s2max];
      yPo = yS2LRF[s2max];
      dTPo = (tS2[s2max]-tS1[s2max])/1000;
      if(pow(dTRn-dTPo,2) > 100)continue; //us^2
      EPo = E;
      outTree->Fill();
      break;
    }
  }
  outTree->Write();
  outFile->Close();
}
