/*
2019-02-27 Wenbo
It could plot the evolution of ba212 at a given runlist, 
and evaluate the radioactive level during that period.

You input a special run list of interest that includes [ runNo ], later you read 
both run-query big table that includes [ runNo | startTime | duration ],
and access the directory of counts data, e.g. coin_analysis/ba212_counts/

inputs:
runlist of interest: coin_analysis/lst_files/2018_physical_lists/2018_Kr_before.txt
run-query big table: coin_analysis/lst_files/2018_Kr_Rn_mid.txt
counts data directory: coin_analysis/ba212_counts/

As for the batch script, see coin_analysis/plotevolution_ba212_batch.sh
 */

#include "TStyle.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TDatime.h"
#include "TGraph.h"
#include "TString.h"
#include <iostream>
#include <fstream>
using namespace std;

void plotevolution_ba212(TString myrunlist)
{
  // read the big table about run-query webpage
  TString bigtable = "lst_files/2018_Kr_Rn_mid.txt"; // input run-query big table
  ifstream fin1(bigtable);
  int fin1line = 0;
  int runNo[1000];
  double startTime[1000];
  double duration[1000];
  int year[1000], month[1000], day[1000], hour[1000], min[1000], sec[1000];
  TDatime date[1000];
  if (fin1.good()) {
    int i = 0;
    double t1,t2,t3,t4,t5,t6,t7;
    // the format of big table is:
    // runNo | startTime............... | duration
    // runNo | month | day | hour | min | hour | min 
    while(fin1 >>
	  t1 >> t2 >>  t3 >> t4 >>  t5 >> t6 >>  t7)
      {
	runNo[i]=t1; 
	year[i]=2018;
	month[i]=t2;
	day[i]=t3;
	hour[i]=t4;
	min[i]=t6;
	sec[i]=0;
	date[i]=TDatime(year[i],month[i],day[i],hour[i],min[i],sec[i]);
	startTime[i]=date[i].Convert(); // [sec] very long very large number! started from 1970
	duration[i] = t6/24. + t7/60./24.; // [day]
	duration[i] *= 86400; // [sec]
	if (duration[i]==0) printf("runNo%d duration %f\n",runNo[i],duration[i]);
	// std::cout << runNo[i] << year[i] << day[i] <<endl; // check value assigning
	++i;
      }
    fin1.close();
    fin1line = i;
  } else {
    printf("Error! fin1 not found\n");
    _exit(0);
  }

  // match input runlist of interest with the big table.
  // TString myrunlist = "lst_files/2018_physical_lsts/one_file.txt";
  ifstream fin2(myrunlist);
  int fin2line;
  int myrunno[1000];
  double mytime[1000];
  double myduration[1000];
  double mycounts[1000]; // used in the mext part
  if (fin2.good()) {
    int runno_t;
    int i = 0;
    while (fin2 >> runno_t){
      myrunno[i] = runno_t;
      for (int j=0;j<fin1line;++j) {
	if (runno_t==runNo[j]) {
	  // std::cout << "L76 successfully == \n"; // check matching
	  mytime[i] = startTime[j];
	  myduration[i] = duration[j];
	}
      }
      ++i;
    }
    fin2.close();
    fin2line = i;
  } else {
    std::cout << "Error! No runlist was found";
    _exit(0);
  }

  // look into counts data and search for counts of a given run
  // counts data are in coin_analysis/ba212_counts/
  for (int i=0;i<fin2line;++i) {
    ifstream fin3(Form("ba212_counts/ba212_%d_counts.txt",myrunno[i]));
    int count_t;
    if (!fin3.good()) {printf("fin3 no file found\n"); _exit(0);}
    while (fin3 >> count_t) {mycounts[i] = count_t; fin3.close();}
  }
  
  // Finally, myrunno[], mytime[], myduration[], and mycounts[] 
  // are successfully prepared. 
  // And myrate[], myerror[] can be calculated.
  // At last we plot it into a TGraph, and get a summary about how much 
  // radioisotopes inside it.
  // 2019-02-27 Wenbo
  int totrunno = fin2line;
  int totcount = 0;
  double totduration = 0;
  for (int i=0;i<totrunno;++i) {
    totcount += mycounts[i];
    totduration += myduration[i];
  }
  double BqNum = double(totcount)/totduration; // get the physical result about radioactivity.
  printf("total count: %d \n", totcount);
  printf("total duration: %f sec \n", totduration);
  printf("radioactivity: %f Bq \n", BqNum);

  double nullarr[1000];
  double myrate[1000];
  double myerror[1000];
  for (int i=0;i<totrunno;++i) { 
    mytime[i] += myduration[i]/2; // move time to central of duration
    myrate[i] = mycounts[i]/myduration[i];
    myerror[i] = sqrt(mycounts[i])/myduration[i];
  }

  // we check the result before plot it.
  // std::cout << totrunno << endl;
  // for (int i=0;i<totrunno;++i) {
  //   std::cout<<myrunno[i]<<" "<<mytime[i]<<" "<<mycounts[i]<<" "<<myrate[i]<<endl;
  // }

  // Plot and save it.
  TCanvas * c1 = new TCanvas("c1","c1",1200,600);
  c1->cd();
  // c1->SetLogy();
  c1->SetTickx();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.09);
  gStyle->SetOptFit(1011);
  TGraph * g1;
  for (int i=0;i<totrunno;++i) myduration[i] /= 2.; // use duration/2 as xerr
  { // check inf and nan values
    // please keep this part remain.
    // 02-28 wenbo
    for (int i=0;i<totrunno;++i) {
      if (mycounts[i]!=mycounts[i]) mycounts[i]=0;
      if (myduration[i]!=myduration[i]) myduration[i]=0;
      if (myerror[i]!=myerror[i]) myerror[i]=0;
    }
    for (int i=0;i<totrunno;++i) {
      printf("%d ",i+1);
      printf("%f ",mytime[i]);
      printf("%d ",mycounts[i]);
      printf("%f ",myduration[i]);
      printf("%f\n",myerror[i]);
    }
  }
  g1 = new TGraphErrors(totrunno,mytime,myrate,myduration,myerror);  // x y xerr yerr
  g1->SetTitle("");
  g1->SetMarkerStyle(20);
  //  g1->SetMarkerSize(0.5);
  //  g1->SetMarkerColor(kRed+10);
  //  g1->SetLineWidth(1.5);
  //  g1->SetLineColor(kRed);
  // int xmin = TDatime(2018,09,10,00,00,00).Convert();
  // int xmax = TDatime(2018,10,31,00,00,00).Convert();
  // g1->GetXaxis()->SetLimits(xmin,xmax);
  g1->GetXaxis()->SetLabelOffset(0.05);
  g1->GetXaxis()->SetLabelFont(132);
  g1->GetXaxis()->SetTimeDisplay(1);
  g1->GetXaxis()->SetTimeFormat("#splitline{%b.%d}{%Y}");
  g1->GetXaxis()->SetTimeOffset(0,"cst");
  g1->GetXaxis()->SetTickLength(0.02);
  // double ymin = 0;
  // double ymax = 155; 
  // g1->GetYaxis()->SetRangeUser(ymin,ymax);
  g1->GetYaxis()->SetTitleFont(132);
  g1->GetYaxis()->SetTitleSize(0.06);
  g1->GetYaxis()->SetTitleOffset(1);
  g1->GetYaxis()->CenterTitle();
  g1->GetYaxis()->SetTitle("Activity [Bq]");
  g1->GetYaxis()->SetTickLength(0.01);
  //  g1->SetFillColor(38);
  g1->Draw("AP");
  c1->SaveAs("plotevolution_result.C");
  // c1->SaveAs("evo_ba212_2017.root");

  //


  // printf("sum: %d[evts]\n",sum);
  // printf("duration_sum: %d[h]\n",duration_sum/3600);
}
