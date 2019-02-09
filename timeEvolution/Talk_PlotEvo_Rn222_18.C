#include <fstream>
#include <iostream>
#include <vector>
#include <map>
using namespace std;
// daqfile, rnEvofile
// void PlotRnEvo2018()
void Talk_PlotEvo_Rn222(const char * inputfilename1, const char * inputfilename2)
{
  // const char * inputfilename1 = "evo_lists/2018_inj_1st.daqinfo"; // daq info
  // const char * inputfilename2 = "APeaks_2018_inj_1st.lst"; // a peaks info 
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
  double Rn222[1000], Po218[1000], Rn220[1000], Po216[1000], Po212[1000];
  double Rn222err[1000], Po218err[1000], Rn220err[1000], Po216err[1000], Po212err[1000];
  if (fin2.good()) {
    int i = 0;
    double t1,t2,t3,t4,t5,t6;
    while(fin2>>t1>>t2>>t3>>t4>>t5>>t6)
      {
	runno2[i]=t1;
	Rn222err[i]=sqrt(t2)/duration[i];
	Po218err[i]=sqrt(t3)/duration[i];
	Rn220err[i]=sqrt(t4)/duration[i];
	Po216err[i]=sqrt(t5)/duration[i];
	Po212err[i]=sqrt(t6)/duration[i];
	Rn222[i]=t2/duration[i];
	Po218[i]=t3/duration[i];
	Rn220[i]=t4/duration[i];
	Po216[i]=t5/duration[i];
	Po212[i]=t6/duration[i];
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

  double nullarr[1000];
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  c1->cd();
  c1->SetTickx();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.13);
  c1->SetLeftMargin(0.12);
  c1->SetRightMargin(0.09);
  gStyle->SetOptFit(1011);
  TGraph * g1;
  g1 = new TGraphErrors(i,time,Rn222,duration,Rn222err);
  g1->SetTitle("");
  g1->SetMarkerStyle(24);
  g1->GetXaxis()->SetLabelOffset(0.05);
  g1->GetXaxis()->SetLabelFont(132);
  g1->GetXaxis()->SetLabelSize(0.04);
  g1->GetXaxis()->SetTimeDisplay(1);
  g1->GetXaxis()->SetTimeFormat("#splitline{%b.%d}{%Y}");
  g1->GetXaxis()->SetTimeOffset(0,"cst");
  g1->GetXaxis()->SetTickLength(0.02);
  g1->GetYaxis()->SetTitleFont(132);
  g1->GetYaxis()->SetTitleSize(0.06);
  g1->GetYaxis()->SetLabelSize(0.04);
  g1->GetYaxis()->SetTitleOffset(0.8);
  g1->GetYaxis()->CenterTitle();
  g1->GetYaxis()->SetTitle("Activity [Bq]");
  g1->GetYaxis()->SetTickLength(0.01);
  g1->Draw("AP");
  //  g1->GetYaxis()->SetRangeUser(1e-2,3e-1);
  //  c1->SetLogy();

  int t0 = time[0]; printf("t0 %d\n",t0);
  for (int i=0;i<45;++i) {time[i]-=t0; time[i]/=86400; duration[i]/=86400;}
  TCanvas * c2 = new TCanvas("c2","c2",800,600);  
  c2->cd();
  c2->SetTickx();
  c2->SetTopMargin(0.05);
  c2->SetBottomMargin(0.13);
  c2->SetLeftMargin(0.12);
  c2->SetRightMargin(0.09);
  gStyle->SetOptFit(1011);
  TGraph * g2 = new TGraphErrors(i,time,Rn222,duration,Rn222err);
  g2->SetTitle("");
  g2->SetMarkerStyle(24);
  g2->GetXaxis()->SetTickLength(0.02);
  g2->GetYaxis()->SetTickLength(0.01);
  g2->Draw("AP");
  TF1 * f1 = new TF1("f1","[0]*exp(-(x-0)/[1])+[2]",0,20);
  f1->SetParameters(0.22,4,0.04);
  g2->Draw("AP");
  g2->Fit("f1","RE");

  // c1->cd();
  // TF1 * f2 = new TF1("f2","[0]*exp(-(x-1539220339)/[1])+[2]",0,30);
  // double * pars = f1->GetParameters();
  // f2->SetParameter(0,f1->GetParameter(0));
  // f2->SetParameter(1,f1->GetParameter(1)*86400);
  // f2->SetParameter(2,f1->GetParameter(2));
  // // TCanvas * c3 = new TCanvas("c3","c3",800,600);  
  // // c3->cd();
  // f2->Draw("");
  
  
  /*
  // TGraph * g2 = new TGraph(i,time,Po212);
  // g2->Draw("");
  // TF1 * f1 = new TF1("name","[0]+[1]*exp(-x*[2])",time[0],time[50]);
  // f1->SetParameters(0.05,2.2e3,1.4e-6);
  // g2->Fit(f1,"RE");
*/

}

