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

void CalAPeaksRate(int & runNo, vector<double> & daughter_vec) 
{
  TChain * chain;
  chain = new TChain("qS1Tree");
  chain->Add(Form("qS1_tree/qS1_tree_run%d.root",runNo));
  printf("CalAPeaksRate for run%d started...\n", runNo);


  // calculate norm factor by iterations
  // h1, h2 and h3 gives three peak values for Rn220 6288keV
  /*
  TH1F * h1 = new TH1F("h1","",200,30e3,100e3);
  qS1Tree->Draw("qS1c>>h1");
  int peakbin;
  double peakvalue, peaklow, peakup;
  peakbin = h1->GetMaximumBin();
  peakvalue = h1->GetBinCenter(peakbin);
  peaklow = peakvalue - 5e3;
  peakup  = peakvalue + 5e3;
  TH1F * h2 = new TH1F("h2","",100,peaklow,peakup);
  qS1Tree->Draw("qS1c>>h2");
  peakbin = h2->GetMaximumBin();
  peakvalue = h2->GetBinCenter(peakbin);
  // the gaus fit below is useless
  // peaklow = peakvalue - 2e3; // printf("peaklow = %f\n", peaklow);
  // peakup  = peakvalue + 2e3; // printf("peakup  = %f\n", peakup );
  // TF1 * f1 = new TF1("f1","gaus",peaklow,peakup);
  // h2->Fit("f1","R");
  // peakvalue = f1->GetParameter(1);
  printf("peakvalue = %f\n", peakvalue);
  delete h1;
  delete h2;
  // delete f1;
  */
  
  // draw normalized qS1 spectrum and fit
  // Rn222 5.29 Po210 5.40 Po218 6.00 [MeV]
  // Bi212 6.09 Rn220 6.29 Po216 6.78 Po212 8.78 [MeV]
  double peakvalue = 52850;
  double Rn222, Po210; // Po218 is ignored
  double Bi212, Rn220, Po216, Po212;
  TH1F * h = new TH1F("h","",300,4.000,10.000);
  qS1Tree->Draw(Form("qS1c*6.288/%f>>h",peakvalue));
  const double kbinwidth = (10.0-5)/500;
  // 5 MeV region
  TF1 * g1 = new TF1("g1","gaus",5.2,5.6); // Rn222 5.59
  TF1 * g2 = new TF1("g2","gaus",5.8,6.2); // Po218 6.00
  TF1 * G1 = new TF1("G1","g1+g2",5.0,6.2); // Rn222 + Po210
  double par_G1[6]; 
  h->Fit(g1,"RQ+"); // Q: Quiet mode (minimum printing) 
  h->Fit(g2,"RQ+"); // R: Use the range specified in the function range
  g1->GetParameters(&par_G1[0]); // address of array begin
  g2->GetParameters(&par_G1[3]); // address of array[3] begin
  G1->SetParameters(par_G1); // address of array begin
  h->Fit(G1,"R+");
  G1->GetParameters(&par_G1[0]); // address of array begin
  Rn222 = par_G1[0]*par_G1[2]*ksqrt2pi/kbinwidth; // printf("Rn222 = %f\n", Rn222);
  Po210 = par_G1[3]*par_G1[5]*ksqrt2pi/kbinwidth;
  // 6 MeV region 
  TF1 * g3 = new TF1("g3", "gaus", 5.8, 6.09); // Bi212 6.09
  TF1 * g4 = new TF1("g4", "gaus", 6.1, 6.6); // Rn220 6.29
  TF1 * g5 = new TF1("g5", "gaus", 6.6, 7.2); // Po216 6.78 (shifted)
  TF1 * G2 = new TF1("G2", "g3+g4+g5", 5.8, 7.2); // all
  double par_G2[9];
  h->Fit(g3,"RQ+");
  h->Fit(g4,"RQ+");
  h->Fit(g5,"RQ+");
  g3->GetParameters(&par_G2[0]); 
  g4->GetParameters(&par_G2[3]); 
  g5->GetParameters(&par_G2[6]); 
  //  par_G2[0] = par_G2[3]/10; par_G2[1] = 6.14; par_G2[2] = par_G2[5]/10;
  G2->SetParameters(par_G2); 
  h->Fit(G2,"R+");
  G2->GetParameters(&par_G2[0]);
  Bi212 = par_G2[0]*par_G2[2]*ksqrt2pi/kbinwidth;
  Rn220 = par_G2[3]*par_G2[5]*ksqrt2pi/kbinwidth;
  Po216 = par_G2[6]*par_G2[8]*ksqrt2pi/kbinwidth;
  // 8 MeV region
  TF1 * g6 = new TF1("g6", "gaus", 8.6, 9.8); // Po212 8784 (shifted)
  double par_g6[3];
  h->Fit(g6,"R+");
  g6->GetParameters(&par_g6[0]);
  Po212 = par_g6[0]*par_g6[2]*ksqrt2pi/kbinwidth;

  printf("Rn222  Po210  Bi212  Rn220  Po216  Po212\n");
  printf("%f %f %f %f %f %f \n", Rn222, Po210, Bi212, Rn220, Po216, Po212);
  
  // printf("%f \n", Rn222);
  vector<double> vec;
  vec.push_back(Rn222);
  vec.push_back(Po210);
  vec.push_back(Bi212);
  vec.push_back(Rn220);
  vec.push_back(Po216);
  vec.push_back(Po212);
  // printf("%f \n", vec[0]);
  // daughter_arr[0] = Rn222;
  daughter_vec = vec;

  printf("CalAPeaksRate for run%d finished.\n", runNo);
  //    */  
}
