/************************************************* 
Copyright:PandaX-II Collaboration 
Author: Andi Tan
Date:2017-09-27
Description:  Looking up for the Bi-Po coincidence
              events in ana data in Run10. The aim
              of the study is to investigate the 
              alpha event on the electrodes.
**************************************************/

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

void BiPo_coin(int &RunNo)
{
  int totalEvt;

  TChain *chain;
  chain = new TChain("signal_tree");
  chain->Add(Form("/store/px/data/udm/by-run/%d/new_correction/ana*",RunNo));

  int runNo, trigNo;
  double triggerTime;

  unsigned int nS1, nS2;
  int s1max, s2max;
  double qS1T[1000], qS1B[1000], xS1T[1000], yS1T[1000];
  double xS1B[1000], yS1B[1000];
  double qS1[1000], tS1[1000], pS1[1000], hS1[1000], ;
  double ttenS1[1000], wS1[1000], wtenS1[1000];
  double xS2NN[1000], yS2NN[1000];
  double qS1Veto[1000];

  double qS2T[1000], qS2B[1000], xS2T[1000], yS2T[1000];
  double xS2B[1000], yS2B[1000];
  double qS2[1000], tS2[1000], pS2[1000], hS2[1000], ;
  double ttenS2[1000], wS2[1000], wtenS2[1000];
  double xS2NN[1000], yS2NN[1000], chi2[1000];

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

  chain->SetBranchAddress("nS2", &nS2);
  chain->SetBranchAddress("s2max", &s2max);
  chain->SetBranchAddress("xS2T", xS2T);
  chain->SetBranchAddress("yS2T", yS2T);
  chain->SetBranchAddress("xS2B", xS2B);
  chain->SetBranchAddress("yS2B", yS2B);
  chain->SetBranchAddress("xS2NN", xS2NN);
  chain->SetBranchAddress("yS2NN", yS2NN);
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

  totalEvt = chain->GetEntries();
  cout << "totalEvents: " << totalEvt<<endl;

  int beta, alpha;
  TFile *outFile = new TFile(Form("BA_tree/BACoin_run%d.root",RunNo),"recreate");
  TTree *outTree = new TTree("BiPo214","Bi-Po214_candidate");
  TTree *outTree2 = new TTree("BiPo212","Bi-Po212_candidate");
  outTree->Branch("trigNo", &trigNo, "trigNo/I");
  outTree->Branch("runNo", &runNo, "runNo/I");
  outTree->Branch("triggerTime", &triggerTime, "triggerTime/D");
  outTree->Branch("beta", &beta, "beta/I");
  outTree->Branch("alpha", &alpha, "alpha/I");
  outTree->Branch("nS1", &nS1, "nS1/I");
  outTree->Branch("s1max", &s1max, "s1max/I");
  outTree->Branch("xS1T", xS1T, "xS1T[nS1]/D");
  outTree->Branch("xS1B", xS1B, "xS1B[nS1]/D");
  outTree->Branch("yS1T", yS1T, "yS1T[nS1]/D");
  outTree->Branch("yS1B", yS1B, "yS1B[nS1]/D");
  outTree->Branch("qS1T", qS1T, "qS1T[nS1]/D");
  outTree->Branch("qS1B", qS1B, "qS1B[nS1]/D");
  outTree->Branch("qS1", qS1, "qS1[nS1]/D");
  outTree->Branch("tS1", tS1, "tS1[nS1]/D");
  outTree->Branch("ttenS1", ttenS1, "ttenS1[nS1]/D");
  outTree->Branch("wS1", wS1, "wS1[nS1]/D");
  outTree->Branch("wtenS1", wtenS1, "wtenS1[nS1]/D");
  outTree->Branch("hS1", hS1, "hS1[nS1]/D");
  outTree->Branch("pS1", pS1, "pS1[nS1]/D");
  outTree->Branch("qS1Veto",qS1Veto,"qS1Veto[nS1]/D");
  outTree->Branch("nS2", &nS2, "nS2/I");
  outTree->Branch("s2max", &s2max, "s2max/I");
  outTree->Branch("xS2T", xS2T, "xS2T[nS2]/D");
  outTree->Branch("yS2T", yS2T, "yS2T[nS2]/D");
  outTree->Branch("xS2B", xS2B, "xS2B[nS2]/D");
  outTree->Branch("yS2B", yS2B, "yS2B[nS2]/D");
  outTree->Branch("xS2NN", xS2NN, "xS2NN[nS2]/D");
  outTree->Branch("yS2NN", yS2NN, "yS2NN[nS2]/D");
  outTree->Branch("chi2", chi2, "chi2[nS2]/D");
  outTree->Branch("tS2", tS2, "tS2[nS2]/D");
  outTree->Branch("ttenS2", ttenS2, "ttenS2[nS2]/D");
  outTree->Branch("qS2T", qS2T, "qS2T[nS2]/D");
  outTree->Branch("qS2B", qS2B, "qS2B[nS2]/D");
  outTree->Branch("qS2", qS2, "qS2[nS2]/D");
  outTree->Branch("wS2", wS2, "wS2[nS2]/D");
  outTree->Branch("wtenS2", wtenS2, "wtenS2[nS2]/D");
  outTree->Branch("hS2", hS2, "hS2[nS2]/D");
  outTree->Branch("pS2", pS2, "pS2[nS2]/D");

  outTree2->Branch("trigNo", &trigNo, "trigNo/I");
  outTree2->Branch("runNo", &runNo, "runNo/I");
  outTree2->Branch("triggerTime", &triggerTime, "triggerTime/D");
  outTree2->Branch("beta", &beta, "beta/I");
  outTree2->Branch("alpha", &alpha, "alpha/I");
  outTree2->Branch("nS1", &nS1, "nS1/I");
  outTree2->Branch("s1max", &s1max, "s1max/I");
  outTree2->Branch("xS1T", xS1T, "xS1T[nS1]/D");
  outTree2->Branch("xS1B", xS1B, "xS1B[nS1]/D");
  outTree2->Branch("yS1T", yS1T, "yS1T[nS1]/D");
  outTree2->Branch("yS1B", yS1B, "yS1B[nS1]/D");
  outTree2->Branch("qS1T", qS1T, "qS1T[nS1]/D");
  outTree2->Branch("qS1B", qS1B, "qS1B[nS1]/D");
  outTree2->Branch("qS1", qS1, "qS1[nS1]/D");
  outTree2->Branch("tS1", tS1, "tS1[nS1]/D");
  outTree2->Branch("ttenS1", ttenS1, "ttenS1[nS1]/D");
  outTree2->Branch("wS1", wS1, "wS1[nS1]/D");
  outTree2->Branch("wtenS1", wtenS1, "wtenS1[nS1]/D");
  outTree2->Branch("hS1", hS1, "hS1[nS1]/D");
  outTree2->Branch("pS1", pS1, "pS1[nS1]/D");
  outTree2->Branch("qS1Veto",qS1Veto,"qS1Veto[nS1]/D");
  outTree2->Branch("nS2", &nS2, "nS2/I");
  outTree2->Branch("s2max", &s2max, "s2max/I");
  outTree2->Branch("xS2T", xS2T, "xS2T[nS2]/D");
  outTree2->Branch("yS2T", yS2T, "yS2T[nS2]/D");
  outTree2->Branch("xS2B", xS2B, "xS2B[nS2]/D");
  outTree2->Branch("yS2B", yS2B, "yS2B[nS2]/D");
  outTree2->Branch("xS2NN", xS2NN, "xS2NN[nS2]/D");
  outTree2->Branch("yS2NN", yS2NN, "yS2NN[nS2]/D");
  outTree2->Branch("chi2", chi2, "chi2[nS2]/D");
  outTree2->Branch("tS2", tS2, "tS2[nS2]/D");
  outTree2->Branch("ttenS2", ttenS2, "ttenS2[nS2]/D");
  outTree2->Branch("qS2T", qS2T, "qS2T[nS2]/D");
  outTree2->Branch("qS2B", qS2B, "qS2B[nS2]/D");
  outTree2->Branch("qS2", qS2, "qS2[nS2]/D");
  outTree2->Branch("wS2", wS2, "wS2[nS2]/D");
  outTree2->Branch("wtenS2", wtenS2, "wtenS2[nS2]/D");
  outTree2->Branch("hS2", hS2, "hS2[nS2]/D");
  outTree2->Branch("pS2", pS2, "pS2[nS2]/D");
  
  double qS1c, Ec;
  outTree->Branch("qS1c", &qS1c, "qS1c/D");
  outTree->Branch("Ec", &Ec, "Ec/D");
  outTree2->Branch("qS1c", &qS1c, "qS1c/D");
  outTree2->Branch("Ec", &Ec, "Ec/D");
  
  
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
  
  double ly(.0);
  double lowAlpha = 16000;
  ly = 43.0668*0.11; // LY in NEST sheet for 50keV at 315V/cm times PDE
  double lowBeta = 50*ly;
  ly = 28.3287*0.11; // LY in NEST sheet for 3.272MeV Q-value(Bi214) at 315V/cm times PDE
  double highBeta = 3272*ly;
  double t_low = 3000; //ns
  double t_high = 500000; //ns
  ly = 43.0668*0.11; // LY in NEST sheet for 50keV at 315V/cm times PDE
  double lowBeta2 = 50*ly;
  ly = 28.5584*0.11; // LY in NEST sheet for 2.252MeV Q-value(Bi212) at 315V/cm times PDE
  double highBeta2 = 3272*ly;
  double t_low2 = 300; //ns
  double t_high2 = 3000; //ns

  // totalEvt = 50000;
  for(int i=0; i<totalEvt; i++){
    chain->GetEntry(i);
    if(i%10000==0) cout << i*1.00/totalEvt*100 <<"%"<< endl;
    beta = alpha = 0;
    if(nS2<1) continue;
    double qBeta(0);
    double qBeta2(0);
    double qS1Tc, qS1Bc; // energy cut
    for(int j = 0; j < nS1; j++){
      // addition filter for alpha; E
      if (s1max<0) continue;
      qS1Tc = qS1T[s1max]*mean_qS1T/(parT[0]+parT[1]*(tS2[s2max]-tS1[s1max])/1000+parT[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parT[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parT[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4));
      qS1Bc = qS1B[s1max]*mean_qS1B/(parB[0]+parB[1]*(tS2[s2max]-tS1[s1max])/1000+parB[2]*pow((tS2[s2max]-tS1[s1max])/1000,2)+parB[3]*1e3*pow((tS2[s2max]-tS1[s1max])/10000,3)+parB[4]*1e4*pow((tS2[s2max]-tS1[s1max])/10000,4)+parB[5]*1e5*pow((tS2[s2max]-tS1[s1max])/10000,5));
      qS1c = (qS1Tc+qS1Bc)*parR[0]/(parR[0]+parR[1]*(yS2NN[s2max]*yS2NN[s2max]+xS2NN[s2max]*xS2NN[s2max])/100+parR[2]*pow((yS2NN[s2max]*yS2NN[s2max]+xS2NN[s2max]*xS2NN[s2max])/100,2));
      Ec = qS1c/parR[0]*6.41*1.0056;
      if (qS1c<0) continue;
      if (Ec<0) continue;
      // Energy region is [8.6, 9.6]MeV, according to Ec spectrum. 
      bool Egood =Ec>8.6&&Ec<9.6; 
      // Egood = 1;
      if(Egood && qS1[j]>lowAlpha && hS1[j]>2000 && ttenS1[j]+wtenS1[j]-pS1[j]<240 && ttenS1[j]+wtenS1[j]-pS1[j]>=90 && hS1[j]/qS1[j]>0.12){
      alpha = j;
      for(int k = 0; k < j; k++){
	// BiPo214 is commented
	// if(qS1[k]<highBeta && qS1[k]>lowBeta && tS1[alpha]-tS1[k]>t_low && tS1[alpha]-tS1[k]<t_high && tS1[k]<tS2[nS2-1] && qBeta<qS1[k]) {
	//   qBeta = qS1[k];
	//   beta = k;
	// }
	if(qS1[k]<highBeta2 && qS1[k]>lowBeta2 && tS1[alpha]-tS1[k]>t_low2 && tS1[alpha]-tS1[k]<t_high2 && tS1[k]<tS2[nS2-1] && qBeta2<qS1[k]) {
	  qBeta2 = qS1[k];
	  beta = k;
	}
      }
      // if (qBeta != 0) outTree->Fill();
      if (qBeta2 != 0) outTree2->Fill();
      }
    }
  }  
  // outTree->Write();
  outTree2->Write();
  outFile->Close();
  printf("100\%\n");
}
