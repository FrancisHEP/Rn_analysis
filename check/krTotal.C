// this script is to trace a signal
#include "/home/zhangdan/rootheader.h"
#include <cstdio>
#include <map>
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <fstream>

//#include <bambooshoot2/BambooShoot.h>
//#include <bambooshoot2/Reader.h>
//#include <bambooshoot2/Event.h>
//#include <bambooshoot2/ClusterData.h>
//#include <bambooshoot2/RawData.h>
//#include <bambooshoot2/HitData.h>
//#include <bambooshoot2/SignalData.h>

#include <TSQLServer.h>
#include <TSQLRow.h>
#include <TSQLResult.h>
#include <TCanvas.h>
#include <TH2Poly.h>
#include <TEllipse.h>
#include <TLine.h>
#include <TLatex.h>
#include <TMath.h>
#include <TBranch.h>
#include <TCut.h>
#include <TDatime.h>

using namespace std;
//using namespace BambooShoot;
int timeBKG = 78.117;
vector<int> fileTime;
std::map<int, string> fileAna3;

void LoadFileTime(int runNumber){
	
//	vector<string> anafiles;
	stringstream ss;
	ss<<runNumber;
	string line;
	string filelist = "/home/zhangdan/ana3/"+ss.str()+"/ana3_"+ss.str()+".lst";
	ifstream fin(filelist.c_str());
	while (getline(fin,line)){
	//	anafiles.push_back(line);
	//	cout<<line.c_str()<<endl;	
		istringstream sin(line);
		vector<string> fields;
		string field;
		while (getline(sin, field,'_')){
		
			fields.push_back(field);
		
		}
		int fileNo = atoi(TString(fields[3]));
		fileAna3[fileNo] = line; 
		int date = atoi(TString(fields[4]));
//		int year = date/10000;
//		int month = int(1.*(date-year*10000)/100);
//		int day = date%100;

		int timeDay = atoi(TString(fields[5]));
//		int hour = int(1.*timeDay/10000);
//		int min = int(1.*(timeDay-hour*1e4)/100);
//		int sec = timeDay-hour*10000-min*100;
		TDatime timeFile;
		timeFile.Set(date,timeDay);
		int timeF = timeFile.Convert();
		//cout<<fileNo<<endl;
		fileTime.push_back(timeF);
		//cout<<year<<month<<day;
		std::sort (fileTime.begin(),fileTime.end());
	}
	


}

vector<string> loadAna3(int runNumber){
	vector<string> anafiles;
	stringstream ss;
	ss<<runNumber;
	string line;
	string filelist = "/home/zhangdan/ana3/"+ss.str()+"/ana3_"+ss.str()+".lst";
	ifstream fin(filelist.c_str());
	while (getline(fin,line)){
		anafiles.push_back(line);
		//cout<<line.c_str()<<endl;	
	
	}
	return anafiles;

}



vector<string> loadAnaTotal(int runNumber){
	vector<string> anafiles;
	stringstream ss;
	ss<<runNumber;
	string line;
//	string filelist = "/nas/a/zhangdan/ana3/"+ss.str()+"/ana3_run"+ss.str()+".root";
	string filelist = "/nas/a/zhangdan/ana3/"+ss.str()+"/ana3_total_run"+ss.str()+"_run10paf.root";
	//string filelist = "/home/zhangdan/ana3/"+ss.str()+"/ana3_run"+ss.str()+".root";
	anafiles.push_back(filelist);
//	anafiles.push_back("/home/zhangdan/ana3/ana3_total.root");
//	anafiles.push_back("/home/zhangdan/ana3/23997/ana3_total_run23997_run10paf.root");

	
	
	return anafiles;

}


int main(int argc, char* argv[])
{
//argv[1] runNumber;
//argv[2] fileNumber;
  //if(argc!=4) {
  //  print_usage(argv[0]);
  //  return 1;
  //}
	int runNumber = atoi(TString(argv[1]));
	double duration = atof(TString(argv[2]));
	LoadFileTime(runNumber);
	TCut total = "tS2-tS1>20e3&&tS2-tS1<350e3&&xTopNN*xTopNN+yTopNN*yTopNN<6e4";
	gStyle->SetOptStat(0);
	TFile *bkgFile = new TFile("/home/zhangdan/kr83mCalibration/bkg.root","READ");
	TTree *bkgTree = (TTree*)bkgFile->Get("out_tree");
	//TH1F *hEBKGhigh = new TH1F("hEBKGhigh","hEBKGhigh",2000,1,3000);

//  int eventNumber;
//  double tS1,tS2;
//  double qS1,qS2;
//  double xx, yy;
//  UInt_t time1,time2;
//	vector<double> Time;
//	vector<double> Kr83;
	vector<string> anafiles = loadAnaTotal(runNumber);
	int Nana = anafiles.size();
	int time0=0;
	string ana1 = anafiles[0];
	int bkgEvents = bkgTree->Draw("13.6e-3*(qS1/0.108+qS2/22.3/0.617)>>hEBKGhigh0(200,0,3000)",total&&"13.6e-3*(qS1/0.108+qS2/22.3/0.617)>0&&13.6e-3*(qS1/0.108+qS2/22.3/0.617)<3000","colz");
	TH1F *hEBKGhigh =(TH1F*)gDirectory->Get("hEBKGhigh0");
	hEBKGhigh->Sumw2();

	TFile *anaFile1 = new TFile(ana1.c_str(),"READ");
  TTree *anaTree1 = (TTree*)anaFile1->Get("out_tree");
	int kr83 = anaTree1->Draw("13.6e-3*(qS1/0.108+qS2/22.3/0.617)>>hEhigh(200,0,3000)",total&&"13.6e-3*(qS1/0.108+qS2/22.3/0.617)>0"&&"13.6e-3*(qS1/0.108+qS2/22.3/0.617)<3000");	

	//TH1F *hEhigh=new TH1F("hEhigh","hEhigh",2000,1,3000);
	
	double scale = 1.*duration/timeBKG;

	cout<<"scale"<<scale<<endl;
	TH1F* hEhigh=(TH1F*)gDirectory->Get("hEhigh");
	hEhigh->Sumw2();

	TCanvas *cs1s2 = new TCanvas("cs1s2","cs1s2",800,600);


	gPad->SetLeftMargin(0.1);
	
	double bkgRb = 1.*scale*hEBKGhigh->Integral();
	hEBKGhigh->Scale(scale);
	double sigRb = hEhigh->Integral();
	double errbkgRb = TMath::Sqrt(hEBKGhigh->Integral())*scale;
	double errsigRb = TMath::Sqrt(sigRb);
	cout<<"sig <10keV\t"<<sigRb<<endl;
	cout<<"total bkg<10keV\t"<<bkgRb<<endl;
	hEhigh->Draw("ehist");
	hEBKGhigh->Draw("esame");
	hEhigh->SetTitle(Form("run%d",runNumber));
	hEhigh->GetXaxis()->SetTitle("E/keV");
	hEhigh->GetYaxis()->SetTitle("Events");
	hEhigh->GetXaxis()->SetTitleSize(0.05);
	hEhigh->GetYaxis()->SetTitleSize(0.05);
	hEhigh->GetXaxis()->SetTitleOffset(0.95);
	hEhigh->GetXaxis()->SetLabelSize(0.05);
	hEhigh->GetYaxis()->SetLabelSize(0.05);
	hEhigh->GetXaxis()->SetNdivisions(6,5,0,kTRUE);
	hEhigh->GetYaxis()->SetNdivisions(9,5,0,kTRUE);
	cs1s2->Update();
	
	double y_max = gPad->GetUymax()*1.2;
	double y_min = 0;
	hEhigh->GetYaxis()->SetRangeUser(y_min,y_max);
	hEhigh->GetXaxis()->CenterTitle();
	hEhigh->GetYaxis()->CenterTitle();
	hEBKGhigh->SetLineWidth(3);
	hEBKGhigh->SetLineColor(2);	
	hEhigh->SetLineWidth(3);
	TH1F *hDiff=(TH1F*) hEhigh->Clone("hDiff");
	hDiff->Add(hEBKGhigh,-1);
	hDiff->SetLineColor(5);
	hDiff->Draw("same");


	TLegend *lg3=new TLegend(0.6,0.7,0.9,0.9);
	lg3->AddEntry(hEBKGhigh,"scaled bkg","l");
	lg3->AddEntry(hEhigh,"signal","l");
	lg3->AddEntry(hDiff,"subtraction","l");
	lg3->Draw();

	cs1s2->Print(Form("/nas/a/zhangdan/check/energyhigh_%d.png",runNumber,runNumber));
	return 0;	
}
