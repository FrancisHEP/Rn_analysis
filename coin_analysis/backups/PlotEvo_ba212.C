#include <fstream>
#include <iostream>
#include <vector>
#include <map>
using namespace std;

void PlotEvo_ba212(const char * inputfilename1, const char * inputfilename2)
{
  // const char * inputfilename1 = "2017_daqinfo.txt"
  // const char * inputfilename2 = "2017_ba212_rate.txt"
  ifstream fin1(inputfilename1);
  int fin1line = 0;
  int runno[1000];
  double time[1000];
  double duration[1000];
  int year[1000], month[1000], day[1000], hour[1000], min[1000], sec[1000];
  TDatime date[1000];
  if (fin1.good()) {
    int i = 0;
    double t1,t2,t3,t4,t5,t6,t7,t8,t9,t10;
    while(fin1>>t1>>t2>>t3>>t4>>t5>>t6>>t7>>t8>>t9>>t10)
      {
	runno[i]=t1;
	year[i]=t2;
	month[i]=t3;
	day[i]=t4;
	hour[i]=t5;
	min[i]=t6;
	sec[i]=t7;
	date[i]=TDatime(year[i],month[i],day[i],hour[i],min[i],sec[i]);
	time[i]=date[i].Convert();
	duration[i]=t8/24.+t9/60./24.+t10/60./60./24.; // [day]
	duration[i] *= 86400; // [sec]
	duration[i] /= 2; // use for error plot
	time[i] += duration[i]; if (duration[i]==0) printf("%d\n",duration[i]);
	++i;
      }
    fin1.close();
    fin1line = i;
  } else {
    printf("Error! fin1 not found\n");
    _exit(0);
  }
  ifstream fin2(inputfilename2);
  int fin2line = 0;
  int runno2[1000];
  double ba212[1000];
  double ba212err[1000];
  if (fin2.good()) {
    int i = 0;
    double t1,t2;
    while(fin2>>t1>>t2)
      {
	runno2[i]=t1;
	ba212err[i]=sqrt(t2)/duration[i];
	ba212[i]=t2/duration[i];
	i++;
      }
    fin2.close();
    fin2line = i;
  } else {
    printf("Error! fin2 not found\n");
    _exit(0);
  }
  if (fin1line!=fin2line) {
    cout << "line num doesn't match!" << endl;
    _exit(0);
  }

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
  double nullarr[1000];
  g1 = new TGraphErrors(i,time,ba212,duration,nullarr);  // x y xerr yerr
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
  // c1->SaveAs("evo_ba212_2017.pdf");
  // c1->SaveAs("evo_ba212_2017.root");

}
