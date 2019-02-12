//////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                              //
//      Andy Tan Revised on Jan. 30th 2016                                                      //
//      For the calculation of Rn-220 Level in Run8                                             //
//      By Rn220-Po216 alpha-alpha coincidence                                                  //
//                                                                                              //
//////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

void Th_Chain_AA(int & runNo)
{
  int totalEvt;
  TChain * chain;
  chain = new TChain("signal_tree");
  chain->Add(Form("/store/px/data/udm/by-run/%d/new_correction/ana*",runNo));
  // int totalEvt;
  // TChain *chain;
  // TFileCollection *fc;
  // chain = new TChain("signal_tree");
  // ifstream inf("./good.lst");
  // for(std::string line; std::getline(inf,line);) chain->Add(line.c_str());
  int runNo, trigNo;
  double triggerTime;
    
  unsigned int nS1, nS2;
  int s1max, s2max;
  double qS1T[10000], qS1B[10000], xS1T[10000], yS1T[10000];
  double xS1B[10000], yS1B[10000];
  double qS1[10000], tS1[10000], pS1[10000], hS1[10000], ;
  double ttenS1[10000], wS1[10000], wtenS1[10000];
  double xS2NN[10000], yS2NN[10000];
    
  double qS2T[10000], qS2B[10000], xS2T[10000], yS2T[10000];
  double xS2B[10000], yS2B[10000];
  double qS2[10000], tS2[10000], pS2[10000], hS2[10000], ;
  double ttenS2[10000], wS2[10000], wtenS2[10000];
  double xS2NN[10000], yS2NN[10000], chi2[10000];
  double xS2PAF[10000], yS2PAF[10000];
    
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
    
  totalEvt = chain->GetEntries();
  cout << "totalEvents: " << totalEvt<<endl;

  int alpha1, alpha2;
  double deltaT;
  int nS1_pre, nS2_pre, s1max_pre, s2max_pre;
  double xS2T_pre[1000], yS2T_pre[1000], xS2B_pre[1000], yS2B_pre[1000], xS2NN_pre[1000], yS2NN_pre[1000],xS2PAF_pre[1000],yS2PAF_pre[1000];
  double chi2_pre[1000];
  double tS1_pre[1000], ttenS1_pre[1000], qS1T_pre[1000], qS1B_pre[1000], qS1_pre[1000], wS1_pre[1000], wtenS1_pre[1000], hS1_pre[1000], pS1_pre[1000];
  double tS2_pre[1000], ttenS2_pre[1000], qS2T_pre[1000], qS2B_pre[1000], qS2_pre[1000], wS2_pre[1000], wtenS2_pre[1000], hS2_pre[1000], pS2_pre[1000];

  TTree *outTree = new TTree("rn220_AA","radon220_AA_candidate_tree");
  outTree->Branch("trigNo", &trigNo, "trigNo/I");
  outTree->Branch("runNo", &runNo, "runNo/I");
  outTree->Branch("triggerTime", &triggerTime, "triggerTime/D");
  outTree->Branch("deltaT", &deltaT, "deltaT/D");
  outTree->Branch("alpha1", &alpha1, "alpha1/I");
  outTree->Branch("alpha2", &alpha2, "alpha2/I");

  outTree->Branch("nS1_pre", &nS1_pre, "nS1_pre/I");
  outTree->Branch("s1max_pre", &s1max_pre, "s1max_pre/I");
  outTree->Branch("tS1_pre", tS1_pre, "tS1_pre[nS1_pre]/D");
  outTree->Branch("ttenS1_pre", ttenS1_pre, "ttenS1_pre[nS1_pre]/D");
  outTree->Branch("qS1T_pre", qS1T_pre, "qS1T_pre[nS1_pre]/D");
  outTree->Branch("qS1B_pre", qS1B_pre, "qS1B_pre[nS1_pre]/D");
  outTree->Branch("qS1_pre", qS1_pre, "qS1_pre[nS1_pre]/D");
  outTree->Branch("wS1_pre", wS1_pre, "wS1_pre[nS1_pre]/D");
  outTree->Branch("wtenS1_pre", wtenS1_pre, "wtenS1_pre[nS1_pre]/D");
  outTree->Branch("hS1_pre", hS1_pre, "hS1_pre[nS1_pre]/D");
  outTree->Branch("pS1_pre", pS1_pre, "pS1_pre[nS1_pre]/D");

  outTree->Branch("nS2_pre", &nS2_pre, "nS2_pre/I");
  outTree->Branch("s2max_pre", &s2max_pre, "s2max_pre/I");
  outTree->Branch("xS2T_pre", xS2T_pre, "xS2T_pre[nS2_pre]/D");
  outTree->Branch("yS2T_pre", yS2T_pre, "yS2T_pre[nS2_pre]/D");
  outTree->Branch("xS2B_pre", xS2B_pre, "xS2B_pre[nS2_pre]/D");
  outTree->Branch("yS2B_pre", yS2B_pre, "yS2B_pre[nS2_pre]/D");
  outTree->Branch("xS2NN_pre", xS2NN_pre, "xS2NN_pre[nS2_pre]/D");
  outTree->Branch("xS2PAF_pre", xS2PAF_pre, "xS2PAF_pre[nS2_pre]/D");
  outTree->Branch("yS2NN_pre", yS2NN_pre, "yS2NN_pre[nS2_pre]/D");
  outTree->Branch("yS2PAF_pre", yS2PAF_pre, "yS2PAF_pre[nS2_pre]/D");
  outTree->Branch("chi2_pre", chi2_pre, "chi2_pre[nS2_pre]/D");
  outTree->Branch("tS2_pre", tS2_pre, "tS2_pre[nS2_pre]/D");
  outTree->Branch("ttenS2_pre", ttenS2_pre, "ttenS2_pre[nS2_pre]/D");
  outTree->Branch("qS2T_pre", qS2T_pre, "qS2T_pre[nS2_pre]/D");
  outTree->Branch("qS2B_pre", qS2B_pre, "qS2B_pre[nS2_pre]/D");
  outTree->Branch("qS2_pre", qS2_pre, "qS2_pre[nS2_pre]/D");
  outTree->Branch("wS2_pre", wS2_pre, "wS2_pre[nS2_pre]/D");
  outTree->Branch("wtenS2_pre", wtenS2_pre, "wtenS2_pre[nS2_pre]/D");
  outTree->Branch("hS2_pre", hS2_pre, "hS2_pre[nS2_pre]/D");
  outTree->Branch("pS2_pre", pS2_pre, "pS2_pre[nS2_pre]/D");

  outTree->Branch("nS1", &nS1, "nS1/I");
  outTree->Branch("s1max", &s1max, "s1max/I");
  outTree->Branch("tS1", tS1, "tS1[nS1]/D");
  outTree->Branch("ttenS1", ttenS1, "ttenS1[nS1]/D");
  outTree->Branch("qS1T", qS1T, "qS1T[nS1]/D");
  outTree->Branch("qS1B", qS1B, "qS1B[nS1]/D");
  outTree->Branch("qS1", qS1, "qS1[nS1]/D");
  outTree->Branch("wS1", wS1, "wS1[nS1]/D");
  outTree->Branch("wtenS1", wtenS1, "wtenS1[nS1]/D");
  outTree->Branch("hS1", hS1, "hS1[nS1]/D");
  outTree->Branch("pS1", pS1, "pS1[nS1]/D");

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

  int Alpha1[1000]={0}, Alpha2[1000]={0}, Alpha1_pre[1000]={0};
  int nAlpha1(0), nAlpha2(0), nAlpha1_pre(0);
  double triggerTime_pre(.0);
  for(int i=0; i<totalEvt; i++){
    //    for(int i=0; i<3000; i++){
    chain->GetEntry(i);
    //        double ly(.0);
    //        ly = 4.4/1.1732;
    double ly_Alpha;
    ly_Alpha = 45424.2/5580;
    double lowAlpha1 = 3200*ly_Alpha;
    double lowAlpha2 = 3600*ly_Alpha;
    double highAlpha1 = 10000*ly_Alpha;
    double highAlpha2 = 12000*ly_Alpha;
        
    if((i+1)%10000==0) cout << i*1.00/totalEvt*100 <<"%"<< endl;
    if(nS2<1) continue;
    if(nS1<1) continue;
    nAlpha1 = nAlpha2 = 0;
        
    for(int j=0; j<nS1; j++){
      if(qS1[j]>lowAlpha1&&tS1[j]<tS2[s2max]&&hS1[j]>2000&&ttenS1[j]+wtenS1[j]-pS1[j]<140&&ttenS1[j]+wtenS1[j]-pS1[j]>=90&&hS1[j]/qS1[j]>0.16){
	Alpha1[nAlpha1] = j;
	nAlpha1++;
      }
      if(qS1[j]>lowAlpha2&&tS1[j]<tS2[s2max]&&hS1[j]>2000&&ttenS1[j]+wtenS1[j]-pS1[j]<140&&ttenS1[j]+wtenS1[j]-pS1[j]>=90&&hS1[j]/qS1[j]>0.16){
	Alpha2[nAlpha2] = j;
	nAlpha2++;
      }
    }
        
    if(nAlpha1_pre*nAlpha2!=0){
      for(int j=0; j<nAlpha2; j++){
	for(int k=0; k<nAlpha1_pre; k++){
	  //                    if(qS1_pre[k]>=qS1[j]) continue;
	  alpha1 = Alpha1_pre[k];
	  alpha2 = Alpha2[j];
	  deltaT = 1.e9*(triggerTime - triggerTime_pre)+tS1[alpha2]-tS1_pre[alpha1];
	  outTree->Fill();
	}
      }
    }
    if(nAlpha1>0){
      nAlpha1_pre = nAlpha1;
      triggerTime_pre = triggerTime;
      for(int z=0; z<nAlpha1; z++){
	Alpha1_pre[z] = Alpha1[z];
      }
      nS1_pre = nS1;
      s1max_pre = s1max;
      for(int z=0; z<nS1; z++){
	tS1_pre[z] = tS1[z];
	ttenS1_pre[z] = ttenS1[z];
	qS1T_pre[z] = qS1T[z];
	qS1B_pre[z] = qS1B[z];
	qS1_pre[z] = qS1[z];
	wS1_pre[z] = wS1[z];
	wtenS1_pre[z] = wtenS1[z];
	hS1_pre[z] = hS1[z];
	pS1_pre[z] = pS1[z];

      }
      nS2_pre = nS2;
      s2max_pre = s2max;
      for(int z=0; z<nS2; z++){
	tS2_pre[z] = tS2[z];
	ttenS2_pre[z] = ttenS2[z];
	qS2T_pre[z] = qS2T[z];
	qS2B_pre[z] = qS2B[z];
	qS2_pre[z] = qS2[z];
	wS2_pre[z] = wS2[z];
	wtenS2_pre[z] = wtenS2[z];
	hS2_pre[z] = hS2[z];
	pS2_pre[z] = pS2[z];
	xS2T_pre[z] = xS2T[z];
	yS2T_pre[z] = yS2T[z];
	xS2B_pre[z] = xS2B[z];
	yS2B_pre[z] = yS2B[z];
	xS2NN_pre[z] = xS2NN[z];
	xS2PAF_pre[z] = xS2PAF[z];
	yS2NN_pre[z] = yS2NN[z];
	yS2PAF_pre[z] = yS2PAF[z];
      }

    }
    nAlpha1_pre = nAlpha1;
  }
  TFile * outFile = new TFile(Form("aa220/outTree_220_216_%d.root",runNo),"recreate");
  outTree->Write();
  outFile->Close();
}
