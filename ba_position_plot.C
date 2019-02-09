 /************************************************* 
Copyright:    PandaX-II Collaboration 
Author:       Wenbo Ma
Description:  This script aimed to draw position
              dependence of BiPo coincidence events.
Change Log:
              2018-11-29 created
**************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
using namespace std;
void SetRunList(TString path_file);

void ba_position_plot(TString runlistfile) // (int& RunNo)
{
  // initializing...
  TChain *ch;
  ch = new TChain("BiPo212"); 
  // ch->Add(Form("tree_Th_aa_run%d.root",RunNo));
  vector<int> runvec; SetRunList(runlistfile, runvec);
  for (int i=0;i<runvec.size();++i) {
    ch->Add(Form("/home/mawenbo/radonRun/Rn_analysis/BA_tree/BACoin_run%d.root",
		 runvec[i]));
    printf("loading...%d\n",runvec[i]);
  }

  // modify bad PMTs. 2018-12-15
  ch->SetAlias("PMTcut1","xS2NN[alpha]<-160||xS2NN[alpha]>-70||yS2NN[alpha]<-320||yS2NN[alpha]>-260||tS2[alpha]-tS1[alpha]>20e3");
  ch->SetAlias("PMTcut2","xS2NN[alpha]<-140||xS2NN[alpha]>-90||yS2NN[alpha]<20||yS2NN[alpha]>80||tS2[alpha]-tS1[alpha]>20e3");

  // printf bulk && cathode evts
  cout << "all BiPo212 cnt: "<<BiPo212->GetEntries("(PMTcut1)&&(PMTcut2)") << endl;
  cout << "bulk cnt: "<<BiPo212->GetEntries("(((qS1T[alpha]-qS1B[alpha])/qS1[alpha]>=-0.5&&qS1[alpha]>55e3)||qS1[alpha]>62e3)&&(PMTcut1)&&(PMTcut2)") << endl;
  cout << "cathode cnt: "<<BiPo212->GetEntries("((qS1T[alpha]-qS1B[alpha])/qS1[alpha]<-0.5&&qS1[alpha]<62e3)&&(PMTcut1)&&(PMTcut2)") << endl;

  // check qS1 & energy spectrum
  TCanvas * c1 = new TCanvas("c1","E spectrum check",1200,800);
  BiPo212->SetAlias("driftTimeCut","tS2[alpha]-tS1[alpha]>0 && tS2[alpha]-tS1[alpha]<400e3");
  c1->Divide(1,3);
  c1->cd(1);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.05);
  //  gPad->SetLogy();
  BiPo212->Draw("qS1[s1max]>>h1(100,0,120000)");
  h1->GetXaxis()->SetTitle("qS1[s1max]");
  c1->cd(2);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.05);
  //  gPad->SetLogy();
  BiPo212->Draw("qS1c>>h2(100,0,120000)","driftTimeCut","");
  h2->GetXaxis()->SetTitle("qS1c");
  c1->cd(3);
  gPad->SetLeftMargin(0.08);
  gPad->SetRightMargin(0.05);
  //  gPad->SetLogy();
  BiPo212->Draw("Ec>>h3(200,0,12)","driftTimeCut",""); // MeV
  h3->GetXaxis()->SetTitle("Ec");

  // check TB ratio
  TCanvas * c2 = new TCanvas("c2","TB ratio check",1400,600);
  BiPo212->SetAlias("TBratio","(qS1T[alpha]-qS1B[alpha])/qS1[alpha]");
  c2->Divide(2,1);
  c2->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  BiPo212->Draw("TBratio:qS1[alpha]>>h21","driftTimeCut","colz");
  h21->GetXaxis()->SetTitle("qS1[alpha]"); 
  h21->GetYaxis()->SetTitle("TB ratio");
  c2->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  BiPo212->Draw("tS2[alpha]-tS1[alpha]:TBratio>>h22","driftTimeCut","colz");
  h22->GetXaxis()->SetTitle("TBratio");
  h22->GetYaxis()->SetTitle("dt"); 

  // draw z-r^2
  TCanvas *c3 = new TCanvas("c3","z v.s. r-square",1400,600);
  BiPo212->SetAlias("x","xS2NN[alpha]");
  BiPo212->SetAlias("y","yS2NN[alpha]");
  BiPo212->SetAlias("dt","tS1[alpha]-tS2[alpha]");
  c3->Divide(2,1);
  c3->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  BiPo212->Draw("dt:(x^2+y^2)>>h31","driftTimeCut&&(PMTcut1)&&(PMTcut2)","colz");
  h31->GetXaxis()->SetRangeUser(-10e3,120e3);
  h31->GetXaxis()->SetTitle("x^2+y^2");
  h31->GetYaxis()->SetTitle("tS2-tS1"); 
  c3->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetRightMargin(0.15);
  BiPo212->Draw("y:x>>h32","driftTimeCut&&(PMTcut1)&&(PMTcut2)","colz");
  h32->GetXaxis()->SetRangeUser(-340,340);
  h32->GetYaxis()->SetRangeUser(-340,340);
  h32->GetXaxis()->SetTitle("x"); 
  h32->GetYaxis()->SetTitle("y");
  // h32->SetMaximum(340);
  // h32->SetMinimum(-340);

  
}

void SetRunList(TString str_path_file, vector<int> & runlistvec)
{
  // runlistvec.clear();
  ifstream path_file(str_path_file);
  if (path_file.good()) {
    int current_number = 0;
    while (path_file >> current_number){
      runlistvec.push_back(current_number);
    }
    path_file.close();
  } else {
    cout << "Error! No run list found";
    _exit(0);
  }
}
