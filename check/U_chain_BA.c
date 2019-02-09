//////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                              //
//      Andy Tan Revised on Jan. 31st 2016                                                      //
//      For the calculation of Rn222 Level in Run8;                                             //
//      Finding Bi-214 and Po-214 conincidence                                                  //
//                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////



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


void U_Chain_BA()
{
    
  int totalEvt;
    

    
    TChain *chain;
    TFileCollection *fc;
    
    chain = new TChain("signal_tree");
	ifstream inf("./good.lst");
        for(std::string line; std::getline(inf,line);) chain->Add(line.c_str());

    int runNo, trigNo;
    double triggerTime;
    
    unsigned int nS1, nS2;
    int s1max, s2max;
    double qS1T[10000], qS1B[10000], xS1T[10000], yS1T[10000];
    double xS1B[10000], yS1B[10000];
    double qS1[10000], tS1[10000], pS1[10000], hS1[10000];
    double ttenS1[10000], wS1[10000], wtenS1[10000];
    double xS2NN[10000], yS2NN[10000];
    double qS1Veto[10000];
    int S1NPeaks[10000], nS1BSatur[10000], nS1TSatur[10000];
    
    double qS2T[10000], qS2B[10000], xS2T[10000], yS2T[10000];
    double xS2B[10000], yS2B[10000];
    double qS2[10000], tS2[10000], pS2[10000], hS2[10000], ;
    double ttenS2[10000], wS2[10000], wtenS2[10000];
    double xS2NN[10000], yS2NN[10000], chi2[10000];
    double xS2PAF[10000], yS2PAF[10000];
    int nS2BSatur[10000], nS2TSatur[10000];


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
    chain->SetBranchAddress("nS2BSatur", nS2BSatur);
    chain->SetBranchAddress("nS2TSatur", nS2TSatur);


    
  totalEvt = chain->GetEntries();
  cout << "totalEvents: " << totalEvt<<endl;

  int beta, alpha;
  TFile *outFile = new TFile("outTree_222_BA.root","recreate");
  TTree *outTree = new TTree("rn222_BA","radon222_candidate");
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
    outTree->Branch("S1NPeaks",S1NPeaks,"S1NPeaks[nS1]/I");
    outTree->Branch("nS1TSatur",nS1TSatur,"nS1TSatur[nS1]/I");
    outTree->Branch("nS1BSatur",nS1BSatur,"nS1BSatur[nS1]/I");
    
    
    outTree->Branch("nS2", &nS2, "nS2/I");
    outTree->Branch("s2max", &s2max, "s2max/I");
    outTree->Branch("xS2T", xS2T, "xS2T[nS2]/D");
    outTree->Branch("yS2T", yS2T, "yS2T[nS2]/D");
    outTree->Branch("xS2B", xS2B, "xS2B[nS2]/D");
    outTree->Branch("yS2B", yS2B, "yS2B[nS2]/D");
    outTree->Branch("xS2NN", xS2NN, "xS2NN[nS2]/D");
    outTree->Branch("xS2PAF", xS2PAF, "xS2PAF[nS2]/D");
    outTree->Branch("yS2NN", yS2NN, "yS2NN[nS2]/D");
    outTree->Branch("yS2PAF", yS2PAF, "yS2PAF[nS2]/D");
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
    outTree->Branch("nS2TSatur",nS2TSatur,"nS1TSatur[nS2]/I");
    outTree->Branch("nS2BSatur",nS2BSatur,"nS1BSatur[nS2]/I");



  for(int i=0; i<totalEvt; i++){
    chain->GetEntry(i);
	double ly_Alpha;
	ly_Alpha = 45424.2/5580;
	double ly_Beta_50 = 4.315*1.29/1.199;		// 164kev, the ly=4.315 then times the ratio of Energy and Beta/Gamma!
	double ly_Beta_4000 = 4.315*0.745/1.052;	// just the same, using 164kev


    double lowAlpha = 4000*ly_Alpha;
    double lowBeta = 50*ly_Beta_50;
    double highBeta = 4000*ly_Beta_4000;
    double t_low = 3000; //ns
    double t_high = 500000;  //ns
    if(i%1000000==0) cout << i*1.00/totalEvt*100 <<"%"<< endl;
    beta = alpha = 0;
    if(nS2<1) continue;
      int Beta[100] = {0}, Alpha[100] = {0};
      int nBeta(0), nAlpha(0);
      for(int j=0; j<nS1; j++){
          if(qS1[j]<highBeta&&qS1[j]>lowBeta&&tS1[j]>490000&&tS1[j]<510000){
            Beta[nBeta] = j;
            nBeta++;
          }
          if(qS1[j]>lowAlpha&&hS1[j]>2000&&ttenS1[j]+wtenS1[j]-pS1[j]<140&&ttenS1[j]+wtenS1[j]-pS1[j]>=90&&hS1[j]/qS1[j]>0.16){
              Alpha[nAlpha] = j;
              nAlpha++;
          }
      }
    if(nBeta*nAlpha==0) continue;
    if(tS1[Alpha[nAlpha-1]]<tS1[Beta[0]]) continue;
    for(int j=0; j<nAlpha; j++){
      for(int k=0; k<nBeta; k++){
        if(tS1[Alpha[j]]-tS1[Beta[k]]>t_low&&tS1[Alpha[j]]-tS1[Beta[k]]<t_high&&tS1[Beta[k]]<tS2[nS2-1]){
          beta = Beta[k];
          alpha = Alpha[j];
          //cout << eventNumber << "\t" << qBeta << "\t" << tBeta << "\t" << qAlpha << "\t" << tAlpha << "\t" << tAlpha-tBeta << endl;
          outTree->Fill();
        }
      }
    }
  }
  outTree->Write();
  outFile->Close();
}
