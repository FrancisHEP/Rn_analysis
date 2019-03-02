/*
  03-02 getCounts
  03-01 TDatimeAutofit

  2019-02-27 
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

  Wenbo
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

// I wrote it to fit evolution precisely. 03-01
void TDatimeAutofit(TCanvas * c1, TGraphErrors * g0, TGraph * gra, double fitmin, double fitmax);
// get counts in a give region from TGraph
double * getCounts(TGraphErrors * g0, TDatime tmin, TDatime tmax);
double * getCounts(TGraphErrors * g0, double xmin, double xmax);

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
    exit(0);
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
    exit(0);
  }

  // look into counts data and search for counts of a given run
  // counts data are in coin_analysis/ba212_counts/
  for (int i=0;i<fin2line;++i) {
    ifstream fin3(Form("ba212_counts/ba212_%d_counts.txt",myrunno[i]));
    int count_t;
    if (!fin3.good()) {printf("fin3 no file found\n"); exit(0);}
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
  TCanvas * c1 = new TCanvas("c1","c1",1000,400);
  c1->cd();
  // c1->SetLogy();
  c1->SetTickx();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.15);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.09);
  gStyle->SetOptFit(1011);
  TGraphErrors * g1;
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
      printf("%f ",mycounts[i]);
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
  double xmin, xmax;
  //  xmin = TDatime(2018,09,10,00,00,00).Convert();
  //  xmax = TDatime(2018,10,31,00,00,00).Convert();
  // g1->GetXaxis()->SetLimits(xmin,xmax);
  g1->GetXaxis()->SetLabelOffset(0.04); // [0, 0.1]
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
  g1->GetYaxis()->SetTitleOffset(0.7); // [0.5, 1]
  g1->GetYaxis()->CenterTitle(); 
  g1->GetYaxis()->SetTitle("Activity [Bq]");
  g1->GetYaxis()->SetTickLength(0.01);
  //  g1->SetFillColor(38);
  g1->Draw("AP");
  c1->SaveAs("plotevolution_result.C");
  // c1->SaveAs("evo_ba212_2017.root");

  TString analyze_cmd = "studyRnOct";

  if (analyze_cmd=="studyRnOct") {
    // lines below are of special fits and physical analysis.
    // analyze Rn decline 1 week stopping injection; 1 month stopping injection
    // analyze Rn bkg before Rn injection and lone time after Rn injection
    // 2019-03-01 Wenbo
    xmin = TDatime(2018,10,8,22,00,00).Convert();
    xmax = TDatime(2018,10,19,18,00,00).Convert();
    TGraphErrors * g_region1 = (TGraphErrors *)g1->Clone(); 
    g_region1->GetXaxis()->SetRangeUser(xmin,xmax);
    /* g_region1->Draw("ap"); */
    /* TF1 * f_region1 = new TF1("f_region1","[0]+[1]*exp(-(x-[2])/[3])",xmin,xmax); */
    /* f_region1->SetParameters(1e-4,1e-4,xmin,86400); */
    /* g_region1->Fit("f_region1","R"); // R: fix in given region */
    TCanvas * c_region1_autofit = new TCanvas("c_region1_autofit","c_region1_autofit",1200,400);
    TGraph * g_region1_autofit;
    TDatimeAutofit(c_region1_autofit,g_region1,g_region1_autofit,xmin,xmax);

    xmin = xmin; //TDatime(2018,10,11,22,00,00).Convert();
    xmax = TDatime(2018,11,6,0,00,00).Convert();
    TGraphErrors * g_region2 = (TGraphErrors *)g1->Clone(); 
    g_region2->GetXaxis()->SetRangeUser(xmin,xmax); 
    TCanvas * c_region2_autofit = new TCanvas("c_region2_autofit","c_region2_autofit",1200,400);
    TGraph * g_region2_autofit;
    TDatimeAutofit(c_region2_autofit,g_region2,g_region2_autofit,xmin,xmax);

    double * counts_arr = getCounts(g1,TDatime(2018,8,9,0,0,0).Convert(),TDatime(2018,9,18,0,0,0).Convert());
    printf("%20f %20f %f \n", counts_arr[0], counts_arr[1], counts_arr[2]);
    counts_arr = getCounts(g1,TDatime(2018,10,18,0,0,0).Convert(),xmax);
    printf("%20f %20f %f \n", counts_arr[0], counts_arr[1], counts_arr[2]);
    counts_arr = getCounts(g1,TDatime(2018,1,1,0,0,0).Convert(),TDatime(2018,3,1,0,0,0).Convert());
    printf("%20f %20f %f \n", counts_arr[0], counts_arr[1], counts_arr[2]);

  }


  // printf("sum: %d[evts]\n",sum);
  // printf("duration_sum: %d[h]\n",duration_sum/3600);
}


// use TDatime::Convert() to get fitmin, fitmax value. 
// please prepare for these parameters beforehand:
// pointer to 1200*400 canvas;
// TGraphErrors with values;
// null pointer to TGraph;
// proper fit region [fitmin, fitmax];
void TDatimeAutofit(TCanvas * c1, TGraphErrors * g0, TGraph * gra, double fitmin, double fitmax)
{
  int n = g0->GetN();
  double * x = g0->GetX(); double * y = g0->GetY();
  double * ex = g0->GetEX(); double * ey = g0->GetEY();
  gra = new TGraph(n,x,y);
  gra->SetMarkerStyle(20);
  double x0 = x[n-1]; 
  for (int i=0;i<n;++i) {
    gra->SetPoint(i,(x[i]-x0)/86400,y[i]); // caution: x[i] is from future to past!!!
  }
  fitmin-=x0; fitmax-=x0; fitmin/=86400; fitmax/=86400;
  gra->GetXaxis()->SetRangeUser(fitmin,fitmax);
  TF1 * f1 = new TF1("f1","[0]+[1]*exp(-(x-[2])/[3])",fitmin,fitmax);
  f1->SetParameters(1e-4,1e-4,fitmin,1);
  gra->Fit("f1","R");

  c1->Divide(2,1);
  c1->cd(1);
  gPad->SetLogy();
  g0->Draw("ap");
  c1->cd(2);
  gPad->SetLogy();
  gra->GetXaxis()->SetTitle("time[day]");
  gra->Draw("ap");
}

double * getCounts(TGraphErrors * g0, TDatime tmin, TDatime tmax) 
{
  /* prepare these parameters beforehand;
     g0;
     TDatime tmin = TDatime(2018,10,18,0,00,00);
     TDatime tmax = TDatime(2018,11,6,0,00,00);
     return array[duration, counts, rate]
  */
  double xmin = tmin.Convert();
  double xmax = tmax.Convert();
  int n = g0->GetN();
  // x y ex ey
  // time rate duration/2 err_rate
  double * x = g0->GetX(); double * y = g0->GetY();
  double * ex = g0->GetEX(); double * ey = g0->GetEY();
  double totduration = 0;
  double totcount = 0;
  for (int i=0;i<n;++i) {
    if (x[i]<xmax&&x[i]>xmin) {
      totduration += ex[i]*2;
      totcount += y[i]*ex[i]*2;
    }
  }
  double arr[3];
  arr[0] = totduration;
  arr[1] = totcount;
  arr[2] = totcount/totduration;
  return arr;
}
double * getCounts(TGraphErrors * g0, double xmin, double xmax) 
{
  /* prepare these parameters beforehand;
     g0;
     xmin;
     xmax;
     return array[duration, counts, rate]
  */
  int n = g0->GetN();
  // x y ex ey
  // time rate duration/2 err_rate
  double * x = g0->GetX(); double * y = g0->GetY();
  double * ex = g0->GetEX(); double * ey = g0->GetEY();
  double totduration = 0;
  double totcount = 0;
  for (int i=0;i<n;++i) {
    if (x[i]<xmax&&x[i]>xmin) {
      totduration += ex[i]*2;
      totcount += y[i]*ex[i]*2;
      printf("%f %f \n",totduration,totcount);
    }
  }
  double arr[3];
  arr[0] = totduration;
  arr[1] = totcount;
  arr[2] = totcount/totduration;
  return arr;
}
