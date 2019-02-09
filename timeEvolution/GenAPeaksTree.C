/************************************************* 
Copyright: PandaX-II Collaboration 
Author: Wenbo Ma
Date: 2018-12-1
Description: 
Calculate rate of alpha peaks from qS1Tree.
**************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "SetRunList.h"
using namespace std;

double * CalAPeaksCount(int & runNo);

void GenAPeaksTree(TString runlist)
{

  // vector<double> v;
  // CalAPeaksCount(runNo,v);
  // for (int i=0;i<v.size();++i) printf("%f ",v[i]);
  // printf("\n");

  int runNo;
  double Rn222, Po218, Rn220, Po216, Po212;
  TTree * tree;
  tree = new TTree("apeaks_tree","alpha peaks tree");
  tree->Branch("runNo", &runNo, "runNo/I");
  tree->Branch("Rn222", &Rn222, "Rn222/D");
  tree->Branch("Po218", &Po218, "Po218/D");
  tree->Branch("Rn220", &Rn220, "Rn220/D");
  tree->Branch("Po216", &Po216, "Po216/D");
  tree->Branch("Po212", &Po212, "Po212/D");

  vector<int> runvec; SetRunList(runlist, runvec);
  for (int i=0;i<runvec.size();++i) {
    double * daughters;
    daughters = CalAPeaksCount(runvec[i]);
    runNo = runvec[i];
    Rn222 = daughters[0];
    Po218 = daughters[1];
    Rn220 = daughters[3];
    Po216 = daughters[4];
    Po212 = daughters[5];
    tree->Fill();
  }
  
  TFile * ofile = new TFile("APeaksTreeTest.root","recreate");
  tree->Write();
  ofile->Close();
  
  printf("GenAPeaksTree finished. \n");

}

double * CalAPeaksCount(int & runNo)
{
  TChain * chain;
  chain = new TChain("qS1Tree");
  chain->Add(Form("qS1_tree/qS1_tree_run%d.root",runNo));
  // printf("run%d ", runNo); fflush(stdout);

  // draw normalized qS1 spectrum and fit
  // Po210 5.31 Rn222 5.59 Po218 6.00 [MeV]
  // Bi212 6.09 Rn220 6.29 Po216 6.78 Po212 8.78 [MeV]
  double peakvalue = 52850;
  double Rn222, Po218; // 
  double Rn220, Po216, Po212;
  qS1Tree->SetAlias("EqS1","qS1c*6.288/52850");
  // 2018
  // Rn222 = qS1Tree->GetEntries("EqS1>5.47-0.10*2&&EqS1<5.47+0.10*2");
  // Po218 = qS1Tree->GetEntries("EqS1>6.02-0.10*2&&EqS1<6.02+0.10*2"); // 6.22
  // Rn220 = qS1Tree->GetEntries("EqS1>6.39-0.09*2&&EqS1<6.39+0.09*2"); // 6.21
  // Po216 = qS1Tree->GetEntries("EqS1>6.92-0.09*2&&EqS1<6.92+0.09*2");
  // Po212 = qS1Tree->GetEntries("EqS1>9.23-0.28*2&&EqS1<9.23+0.28*2");
  // 2017
  Rn222 = qS1Tree->GetEntries("EqS1>5.35-0.118 *2&&EqS1<5.35+0.118 *2");
  Po218 = qS1Tree->GetEntries("EqS1>5.93-0.0935*2&&EqS1<5.93+0.0935*2"); // 6.22
  Rn220 = qS1Tree->GetEntries("EqS1>6.29-0.0912*2&&EqS1<6.29+0.0912*2"); // 6.21
  Po216 = qS1Tree->GetEntries("EqS1>6.81-0.0951*2&&EqS1<6.81+0.0951*2");
  Po212 = qS1Tree->GetEntries("EqS1>9.18-0.296 *2&&EqS1<9.18+0.296 *2");

  // printf("Rn222  Po218  Rn220  Po216  Po212\n");
  // %-12f: at least 12 wide && left-justified
  printf("%6d %-12f %-12f %-12f %-12f %-12f \n", runNo, Rn222, Po218, Rn220, Po216, Po212);
  
  double arr[5];
  arr[0] = Rn222 ;
  arr[1] = Po218 ;
  arr[2] = Rn220 ;
  arr[3] = Po216 ;
  arr[4] = Po212 ;
  
  // printf("CalAPeaksCount for run%d finished.\n", runNo);
  //    */  
  delete chain;
  return arr;
}
