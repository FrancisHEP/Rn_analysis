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

const double ksqrt2pi = sqrt(2*3.14159265359);

vector<double> CalAPeaksCount(int & runNo) 
{
  TChain * chain;
  chain = new TChain("qS1Tree");
  chain->Add(Form("qS1_tree/qS1_tree_run%d.root",runNo));
  printf("CalAPeaksCount for run%d started...\n", runNo);

  // draw normalized qS1 spectrum and fit
  // Po210 5.31 Rn222 5.59 Po218 6.00 [MeV]
  // Bi212 6.09 Rn220 6.29 Po216 6.78 Po212 8.78 [MeV]
  double peakvalue = 52850;
  double Rn222, Po218; // 
  double Rn220, Po216, Po212;
  qS1Tree->SetAlias("EqS1","qS1c*6.288/52850");
  Rn222 = qS1Tree->GetEntries("EqS1>5.47-0.10*2&&EqS1<5.47+0.10*2");
  Po218 = qS1Tree->GetEntries("EqS1>6.02-0.10*2&&EqS1<6.02+0.10*2");
  Rn220 = qS1Tree->GetEntries("EqS1>6.39-0.09*2&&EqS1<6.39+0.09*2");
  Po216 = qS1Tree->GetEntries("EqS1>6.92-0.09*2&&EqS1<6.92+0.09*2");
  Po212 = qS1Tree->GetEntries("EqS1>9.23-0.28*2&&EqS1<9.23+0.28*2");

  // printf("Rn222  Po218  Rn220  Po216  Po212\n");
  // printf("%f %f %f %f %f \n", Rn222, Po218, Rn220, Po216, Po212);
  
  vector<double> daughter_vec;
  daughter_vec.clear();
  daughter_vec.push_back(Rn222);
  daughter_vec.push_back(Po218);
  daughter_vec.push_back(Rn220);
  daughter_vec.push_back(Po216);
  daughter_vec.push_back(Po212);
  for (int i=0;i<daughter_vec.size();++i) printf("%f ", daughter_vec[i]);
  
  printf("CalAPeaksCount for run%d finished.\n", runNo);
  //    */  
  return daughter_vec;
}
