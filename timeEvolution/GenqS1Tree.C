/************************************************* 
Copyright: PandaX-II Collaboration 
Author: Wenbo Ma
Date: 2018-11-30
Description: 
Generate qS1Tree, which only contain qS1 and qS1c. 
qS1c is the qS1 correction for high-energy alpha events, such as Rn220, Po216, Bi212, Po212, Tl208.
This tree is for automatic qS1 spectrum normalization and alpha rate calculation in the next program, CalAPeaksRate.C
Usage:
Read ana1 data on bl;
Write qS1Tree in directory qS1_tree/
**************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>

void GenqS1Tree(int & runNo)
{
  TChain * chain;
  chain = new TChain("signal_tree");
  chain->Add(Form("/store/px/data/udm/by-run/%d/new_correction/ana*",runNo));

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

  int totalEvt = chain->GetEntries();
  cout << "totalEvents: " << totalEvt<<endl;

  double o_tS1,o_tS2,o_qS1,o_qS1c,o_xS2,o_yS2;
  TTree * otree = new TTree("qS1Tree","qS1Tree for high energy alpha events");
  otree->Branch("tS1", &o_tS1, "tS1/D");
  otree->Branch("tS2", &o_tS2, "tS2/D");
  otree->Branch("qS1", &o_qS1, "qS1/D");
  otree->Branch("qS1c", &o_qS1c, "qS1c/D");
  otree->Branch("xS2", &o_xS2, "xS2/D");
  otree->Branch("yS2", &o_yS2, "yS2/D");
  
  for (int i=0;i<totalEvt;i++){
    chain->GetEntry(i);
    if(i%(totalEvt/999)==0) cout << i*1.00/totalEvt*100 << "%" << "\r" ;
    fflush(stdout);
    bool flag;
    flag = s1max>=0&&s2max>=0; if(!flag)continue;
    flag = qS1[s1max]>=3&&qS2B[s2max]>30; if(!flag)continue;
    flag = tS2[s2max]-tS1[s1max]<250e3&&tS2[s2max]-tS1[s1max]>50e3&&xS2NN[s2max]*yS2NN[s2max]+yS2NN[s2max]*yS2NN[s2max]<6e4; if(!flag)continue;
    o_tS1 = tS1[s1max];
    o_tS2 = tS2[s2max];
    o_qS1 = qS1[s1max];
    o_qS1c = qS1[s1max]*exp((tS1[s1max]-tS2[s2max]+150000)/1200000)/(1-1E-6*(xS2NN[s2max]*xS2NN[s2max]+yS2NN[s2max]*yS2NN[s2max]));
    o_xS2 = xS2NN[s2max];
    o_yS2 = yS2NN[s2max];
    otree->Fill();
  }

  TFile * ofile = new TFile(Form("qS1_tree/qS1_tree_run%d.root",runNo),"recreate");
  otree->Write();
  ofile->Close();

  printf("finished\n");

}
