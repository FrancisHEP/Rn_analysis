// Modified by Wenbo
// 2018-11-28
// select both RnPo evts and BiPo evts
// 2018-11-23
// use NN instead LRF for position reconstruction
// show progress

// usage: 
// load ana1 file
// ch->Add(Form("/store/px/data/udm/by-run/%d/new_correction/ana*",RunNo));
// and create aatree.root in directory aacoin_tree
// TFile *outFile = new TFile(Form("rootfile_aa_tree/tree_Th_aa_run%d.root",RunNo),"recreate");

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>
#include "TROOT.h"
#include "TBrowser.h"
#include <TChain.h>
#include <TFileCollection.h>

void coinSelector(int &RunNo)
{
  TChain *ch;
  ch = new TChain("signal_tree");
  ch->Add(Form("/store/px/data/udm/by-run/%d/new_correction/ana*",RunNo));

  int totalEvt;
  int runNo, fileNo, trigNo;
  unsigned int time;
  double triggerTime;
  unsigned int nS1, nS2, SingleS2;
  int s1max, s2max;
  double tS1[1000], ttenS1[1000], wS1[1000], wtenS1[1000], fwhmS1[1000], hS1[1000],pS1[1000];
  double tS2[1000];
  double qS1[1000], qS1T[1000], qS1B[1000], qS2[1000];
  double xS2PAF[1000], yS2PAF[1000], chi2_PAF[1000];
  double xS2NN[1000], yS2NN[1000];

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

  ch->SetBranchAddress("xS2PAF",xS2PAF);
  ch->SetBranchAddress("yS2PAF",yS2PAF);
  ch->SetBranchAddress("chi2_PAF",chi2_PAF);
  ch->SetBranchAddress("tS2",tS2);
  ch->SetBranchAddress("qS2",qS2);
  ch->SetBranchAddress("xS2NN",xS2NN);
  ch->SetBranchAddress("yS2NN",yS2NN);

  totalEvt = ch->GetEntries();
  cout << "totalEvents: " << totalEvt << endl;
  
  int runno, fileno, evtno;
  double decayT_AA;
  double xRn, yRn, xPo, yPo, dTRn, dTPo;
  double ERn, EPo;
  double decayT_BA;
  double xBi, yBi, xPo_BA, yPo_BA, dTBi, dTPo_BA;
  double EBi, EPo_BA;
  TFile *outFile = new TFile(Form("coin_tree/tree_Th_coin_run%d.root",RunNo),"recreate");
  TTree *outTree = new TTree("coin_tree","radon220_coinevts_candidate_tree");
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
  outTree->Branch("decayT_BA",&decayT_BA,"decayT_BA/D");
  outTree->Branch("xBi",&xBi,"xBi/D");
  outTree->Branch("yBi",&yBi,"yBi/D");
  outTree->Branch("xPo_BA",&xPo_BA,"xPo_BA/D");
  outTree->Branch("yPo_BA",&yPo_BA,"yPo_BA/D");
  outTree->Branch("dTBi",&dTBi,"dTBi/D");
  outTree->Branch("dTPo_BA",&dTPo_BA,"dTPo_BA/D");
  outTree->Branch("EBi",&EBi,"EBi/D");
  outTree->Branch("EPo_BA",&EPo_BA,"EPo_BA/D");

  
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


  // BiPo
  double ly(.0);
  double lowAlpha = 16000;
  ly = 43.0668*0.11; // LY in NEST sheet for 50keV at 315V/cm times PDE
  double lowBeta = 50*ly;
  ly = 28.3287*0.11; // LY in NEST sheet for 3.272MeV Q-value(Bi214) at 315V/cm times PDE
  double highBeta = 3272*ly;
  double t_low = 3000; //ns
  double t_high = 500000;  //ns
  ly = 43.0668*0.11; // LY in NEST sheet for 50keV at 315V/cm times PDE
  double lowBeta2 = 50*ly;
  ly = 28.5584*0.11; // LY in NEST sheet for 2.252MeV Q-value(Bi212) at 315V/cm times PDE
  double highBeta2 = 3272*ly;
  double t_low2 = 300; //ns
  double t_high2 = 3000;  //ns


  int N_T_half = 10;
  double dR_xy = 400;
  double qS1Tc, qS1Bc, qS1c, E;
  //totalEvt = 100;
  for (int ii = 0; ii < totalEvt; ii++){
    ch->GetEntry(ii);
    if((ii+1)%100000==0) cout << ii*1.00/totalEvt*100 <<"%"<< endl;
    if(nS1>1000||s1max<0||s2max<0) continue;
    if(hS1[s1max]<=2000||ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>=140||ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<40||hS1[s1max]/qS1[s1max]<=0.16||qS2[s2max]<=50) continue;//pulse shape
    /*
    // modified by Wenbo
    // 2018-11-23
    // nan values for xS2PAF 
    if(xS2PAF[s2max]!=xS2PAF[s2max]||yS2PAF[s2max]!=yS2PAF[s2max]) continue;
    */
    //if(tS2-tS1>270e3||tS2-tS1<50e3) continue;
    qS1Tc = qS1T[s1max]*mean_qS1T/(parT[0]+parT[1]*(tS2[s2max]-tS1[s1max])/1000+parT[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parT[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parT[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4));
    qS1Bc = qS1B[s1max]*mean_qS1B/(parB[0]+parB[1]*(tS2[s2max]-tS1[s1max])/1000+parB[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parB[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parB[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4)+parB[5]*1e5*pow((tS2[s2max]-tS1[s1max])/10000,5));
    qS1c = (qS1Tc+qS1Bc)*parR[0]/(parR[0]+parR[1]*(yS2PAF[s2max]*yS2PAF[s2max]+xS2PAF[s2max]*xS2PAF[s2max])/100+parR[2]*pow((yS2PAF[s2max]*yS2PAF[s2max]+xS2PAF[s2max]*xS2PAF[s2max])/100,2));
    E = qS1c/parR[0]*6.41*1.0056;
    //cout << trigNo << '\t' << qS1Tc <<  '\t' << qS1Bc <<  '\t' << E << endl;
    if(E<6.15||E>6.65) continue; // energy filter for Rn220
    double startT = time + triggerTime;
    runno = runNo;
    fileno = fileNo;
    evtno = trigNo;
    xRn = xS2NN[s2max];
    yRn = yS2NN[s2max];
    dTRn = (tS2[s2max]-tS1[s1max])/1000;
    ERn = E; // good 
    for (int jj = ii+1; jj < totalEvt; jj++){
      ch->GetEntry(jj);
      if(nS1>1000||s1max<0||s2max<0) continue;
      if(time + triggerTime > startT + 0.15*N_T_half) break;
      // // modified by Wenbo
      // // 2018-11-23
      // // nan values for xS2PAF 
      // if(xS2PAF[s2max]!=xS2PAF[s2max]||yS2PAF[s2max]!=yS2PAF[s2max]) continue;
      //if(xTopPAF*xTopPAF+yTopPAF*yTopPAF>80000) continue;
      //if(tS2-tS1>270e3||tS2-tS1<50e3) continue;
      qS1Tc = qS1T[s1max]*mean_qS1T/(parT[0]+parT[1]*(tS2[s2max]-tS1[s1max])/1000+parT[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parT[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parT[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4));
      qS1Bc = qS1B[s1max]*mean_qS1B/(parB[0]+parB[1]*(tS2[s2max]-tS1[s1max])/1000+parB[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parB[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parB[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4)+parB[5]*1e5*pow((tS2[s2max]-tS1[s1max])/10000,5));
      qS1c = (qS1Tc+qS1Bc)*parR[0]/(parR[0]+parR[1]*(yS2PAF[s2max]*yS2PAF[s2max]+xS2PAF[s2max]*xS2PAF[s2max])/100+parR[2]*pow((yS2PAF[s2max]*yS2PAF[s2max]+xS2PAF[s2max]*xS2PAF[s2max])/100,2));
      E = qS1c/parR[0]*6.41*1.0056;
      if(E<6.65||E>7.5) continue; // energy filter for Po216
      //cout << trigNo << '\t' << qS1Tc <<  '\t' << qS1Bc <<  '\t' << E << endl;
      if(pow(xS2PAF[s2max]-xRn,2)+pow(yS2PAF[s2max]-yRn,2) > 324*324./dR_xy)continue;
      decayT = time + triggerTime - startT;
      xPo = xS2NN[s2max]; if(xPo!=xPo) continue; // nan value filter Wenbo
      yPo = yS2NN[s2max]; if(yPo!=yPo) continue; // nan value filter Wenbo
      dTPo = (tS2[s2max]-tS1[s2max])/1000;
      if(pow(dTRn-dTPo,2) > 100)continue; //us^2
      EPo = E; if(EPo!=EPo) continue; // nan value
      outTree->Fill();
      break;
    }
  }
  for(int i=0; i<totalEvt; i++){
    ch->GetEntry(i);
    if((i+1)%100000==0) cout << i*1.00/totalEvt*100 <<"%"<< endl;
    beta = alpha = 0;
    if(nS2<1) continue;
    double qBeta(0);
    double qBeta2(0);
    for(int j = 0; j < nS1; j++){
      if(qS1[j]>lowAlpha && hS1[j]>2000 && ttenS1[j]+wtenS1[j]-pS1[j]<240 && ttenS1[j]+wtenS1[j]-pS1[j]>=90 && hS1[j]/qS1[j]>0.12){
        alpha = j;
        for(int k = 0; k < j; k++){
          // if(qS1[k]<highBeta && qS1[k]>lowBeta && tS1[alpha]-tS1[k]>t_low && tS1[alpha]-tS1[k]<t_high && tS1[k]<tS2[nS2-1] && qBeta<qS1[k]) {
          //   qBeta = qS1[k];
          //   beta = k;
          // }
          if(qS1[k]<highBeta2 && qS1[k]>lowBeta2 && tS1[alpha]-tS1[k]>t_low2 && tS1[alpha]-tS1[k]<t_high2 && tS1[k]<tS2[nS2-1] && qBeta2<qS1[k]) {
            qBeta2 = qS1[k];
            beta = k;
          }
        }
	//         if (qBeta != 0) outTree->Fill();
        if (qBeta2 != 0) outTree2->Fill();
      }
    }
  }
  outTree->Write();
  int outEvts = outTree->GetEntries();
  printf("outTree->GetEntries() = %d\n",outEvts);
  outFile->Close();
  printf("finished.\n");
}
