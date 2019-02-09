#include <iostream>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#include "TROOT.h"
#include "TBrowser.h"
#include <fstream>
#include <vector>
#include <algorithm>


void Single_Alpha(){

	int runNo = 0;
     unsigned int nS1, nS2;
     double qS1[1000], tS1[1000];
     double qS2[1000], tS2[1000];
     double qS1T[1000], qS1B[1000];
     double qS2T[1000], qS2B[1000];
     double xS2NN[1000], yS2NN[1000];
     double xS2PAF[1000], yS2PAF[1000];
     double xS2T[1000], yS2T[1000];
     double qS1Veto[1000];
     double qS1Max;
     double tS1Max;
     int s1max, s2max;

	double pS1[1000], hS1[1000];
	double pS2[1000], hS2[1000];
	double ttenS1[1000], wS1[1000], wtenS1[1000];
	double ttenS2[1000], wS2[1000], wtenS2[1000];
	
	TChain ch("signal_tree");
	ch.Add("rootfile/ana_run*.root");
//	ifstream inf("./good.lst");
//        for(std::string line; std::getline(inf,line);) ch.Add(line.c_str());

	int time= 4.06260000000000000e+06;

     ch.SetBranchAddress("runNo", &runNo);
     ch.SetBranchAddress("nS1", &nS1);
     ch.SetBranchAddress("qS1",qS1);
     ch.SetBranchAddress("tS1",tS1);
     ch.SetBranchAddress("qS1T",qS1T);
     ch.SetBranchAddress("qS1B",qS1B);
     ch.SetBranchAddress("qS2T",qS2T);
     ch.SetBranchAddress("qS2B",qS2B);
     ch.SetBranchAddress("nS2", &nS2);
     ch.SetBranchAddress("tS2", tS2);
     ch.SetBranchAddress("qS2", qS2);
     ch.SetBranchAddress("xS2NN",xS2NN);
     ch.SetBranchAddress("xS2PAF",xS2PAF);
     ch.SetBranchAddress("yS2NN",yS2NN);
     ch.SetBranchAddress("yS2PAF",yS2PAF);
     ch.SetBranchAddress("s1max",&s1max);
     ch.SetBranchAddress("s2max",&s2max);
     ch.SetBranchAddress("tS1Max",&tS1Max);
     ch.SetBranchAddress("qS1Max",&qS1Max);
     ch.SetBranchAddress("xS2T",xS2T);
     ch.SetBranchAddress("yS2T",yS2T);
     ch.SetBranchAddress("qS1Veto",qS1Veto);

    ch.SetBranchAddress("ttenS1", ttenS1);
    ch.SetBranchAddress("wS1", wS1);
    ch.SetBranchAddress("wtenS1", wtenS1);
    ch.SetBranchAddress("hS1", hS1);
    ch.SetBranchAddress("pS1", pS1);
    ch.SetBranchAddress("ttenS2", ttenS2);
    ch.SetBranchAddress("wS2", wS2);
    ch.SetBranchAddress("wtenS2", wtenS2);
    ch.SetBranchAddress("hS2", hS2);
    ch.SetBranchAddress("pS2", pS2);

	double totalEvt = ch.GetEntries();
	cout << "totalEvt = " << totalEvt << endl;

	signal_tree->SetAlias("qS2cr",Form("(qS2[s2max]*exp((tS2[s2max]-tS1Max)/400000)/(1+3.7E-8*(xS2NN[s2max]^2+yS2NN[s2max]^2)-2.39E-11*(xS2NN[s2max]^2+yS2NN[s2max]^2)^2))"));
	signal_tree->SetAlias("qS1Maxcr",Form("qS1Max*exp((tS1Max-tS2[s2max]+150000)/1200000)/(1-1E-6*(xS2NN[s2max]^2+yS2NN[s2max]^2))"));

	TCut FV = "tS2[s2max]-tS1Max>25000&&tS2[s2max]-tS1Max<340000&&xS2PAF[s2max]^2+yS2PAF[s2max]^2<72000";
	TCut Pulse_Shape = "s1max>=0&&s2max>=0&&hS1[s1max]>2000&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]>90&&ttenS1[s1max]+wtenS1[s1max]-pS1[s1max]<250&&hS1[s1max]/qS1Max>0.1&&nS2<15&&S1NPeaks[s1max]<=1";

	TCut Ra_Po_cut = "log10(qS2cr/qS1Maxcr)>0&&log10(qS2cr/qS1Maxcr)<1&&log10(qS1Maxcr)>4.6&&log10(qS1Maxcr)<4.74"; 



	// Draw Log(S2/S1):Log(S1) Pic
	TCanvas *c1=new TCanvas("c1","c1",520,520);
    TH2F *h1 = new TH2F("h1","",200,4,5.5,200,-4,4);
	signal_tree->Draw("log10(qS2cr/qS1Maxcr):log10(qS1Maxcr)>>h1",Pulse_Shape&&""&&FV,"colz");
	h1->SetTitle("Alpha Events Distribution");
	h1->GetXaxis()->SetTitle("log10(qS1Maxcr)");
	h1->GetYaxis()->SetTitle("log10(qS2cr/qS1Maxcr)");
	
	TLine *l1=new TLine(4.6,0,4.6,1);
	l1->SetLineColor(2);
	l1->SetLineWidth(2);
	l1->Draw("same");
        TLine *l2=new TLine(4.74,0,4.74,1);
        l2->SetLineColor(2);
        l2->SetLineWidth(2);
        l2->Draw("same");
        TLine *l3=new TLine(4.6,0,4.74,0);
        l3->SetLineColor(2);
        l3->SetLineWidth(2);
        l3->Draw("same");
        TLine *l4=new TLine(4.6,1,4.74,1);
        l4->SetLineColor(2);
        l4->SetLineWidth(2);
        l4->Draw("same");

	// Draw the S1 spectrum of Rn222 and Po218 and Bi212;
	TCanvas *c2=new TCanvas("c2","c2",520,520);
	gStyle->SetOptFit(1011);	
	TH1F *h2 = new TH1F("h2","",100,10000,100000);
	signal_tree->Draw("qS1Maxcr>>h2",Pulse_Shape&&""&&FV,"");
	h2->SetTitle("qS1 Spectrum");
	h2->GetXaxis()->SetTitle("qS1Maxcr(pe)");
	h2->GetYaxis()->SetTitle("count");

	double Counts = h2->GetEntries();
	
	TF1 *f1=new TF1("f1","gaus(0)",43000,48000);
	f1->SetParameters(Counts/11,46000,1000);
	TF1 *f2=new TF1("f2","gaus(0)",48000,54000);
	f2->SetParameters(Counts/15,51000,1000);
	TF1 *f3=new TF1("f3","f1+f2",43000,54000);
	h2->Fit("f3","R");

	double p0_1 = f3->GetParameter(0);
	double p2_1 = f3->GetParameter(2);
	double p0_2 = f3->GetParameter(3);
	double p2_2 = f3->GetParameter(5);

	cout << "Mean 222Rn= " << f3->GetParameter(2) << endl;
        cout << "Mean 218Po= " << f3->GetParameter(4) << endl;

	double Count_1 = p0_1*p2_1*pow(2*3.1415926,0.5)/300;
	double Count_2 = p0_2*p2_2*pow(2*3.1415926,0.5)/300;
	cout << "Count_222Rn= " << Count_1 << endl;
	cout << "Count_218Po= " << Count_2 << endl;

	// Draw Energy spectrum based on 222Rn is 5.59MeV
	TCanvas *c3=new TCanvas("c3","c3",520,520);

	double Mean_1 = f3->GetParameter(1);
	TH1F *h3 = new TH1F("h3","",100,4,10);
	TString energy = Form("qS1Maxcr/%.3f*5.59>>h3",Mean_1);
	
	signal_tree->Draw(energy,Pulse_Shape&&""&&FV,"");
	h3->SetStats(0);
	double Count = h3->GetEntries();
	h3->SetTitle("Alpha Spectrum");
	h3->GetXaxis()->SetTitle("Energy(MeV)");
	h3->GetYaxis()->SetTitle("Counts");
	h3->GetYaxis()->SetRangeUser(0,10000);

	//////////////// 222Rn Chain
	// 222Rn
	TF1 *ff1=new TF1("ff1","gaus(0)",5.2,5.8);
	ff1->SetParameters(Count/12,5.59,0.15);
	ff1->SetLineColor(2);
	
	TLine *ll1=new TLine(5.59,0,5.59,10000);
	ll1->SetLineColor(2);
	ll1->SetLineStyle(9);
	ll1->SetLineWidth(2);
	ll1->Draw("same");

	TLatex *tt1=new TLatex(5.56,9000,"^{222}Rn");
	tt1->SetTextColor(2);
	tt1->SetTextAngle(90);
	tt1->SetTextSize(0.03);
	tt1->Draw("same");

	// 218Po + 212Bi
	TF1 *ff2=new TF1("ff2","gaus(0)",5.8,6.4);
	ff2->SetParameters(Count/15,6.16,0.15);
	ff2->SetLineColor(3);

        TLine *ll2=new TLine(6.16,0,6.16,10000);
        ll2->SetLineColor(3);
        ll2->SetLineStyle(9);
        ll2->SetLineWidth(2);
        ll2->Draw("same");

        TLatex *tt2=new TLatex(6.13,8000,"^{218}Po+^{212}Bi");
        tt2->SetTextColor(3);
        tt2->SetTextAngle(90);
        tt2->SetTextSize(0.03);
        tt2->Draw("same");

	// 214Po
	TF1 *ff3=new TF1("ff3","gaus(0)",7.4,8.6);
	ff3->SetParameters(Count/150,7.84,0.4);
	ff3->SetLineColor(4);

        TLine *ll3=new TLine(7.84,0,7.84,10000);
        ll3->SetLineColor(4);
        ll3->SetLineStyle(9);
        ll3->SetLineWidth(2);
        ll3->Draw("same");

        TLatex *tt3=new TLatex(7.81,9000,"^{214}Po");
        tt3->SetTextColor(4);
        tt3->SetTextAngle(90);
        tt3->SetTextSize(0.03);
        tt3->Draw("same");

	// 210Po
	TF1 *ff4=new TF1("ff4","gaus(0)",4.8,5.6);
	ff4->SetParameters(Count/80,5.3,0.3);
	ff4->SetLineColor(5);

        TLine *ll4=new TLine(5.3,0,5.3,10000);
        ll4->SetLineColor(5);
        ll4->SetLineStyle(9);
        ll4->SetLineWidth(2);
        ll4->Draw("same");

        TLatex *tt4=new TLatex(5.27,9000,"^{210}Po");
        tt4->SetTextColor(5);
        tt4->SetTextAngle(90);
        tt4->SetTextSize(0.03);
        tt4->Draw("same");

	//////////////// 224Ra(220Rn) Chain
	//220Rn
	TF1 *ff5=new TF1("ff5","gaus(0)",6.0,6.8);
	ff5->SetParameters(Count/1200,6.41,0.3);
	ff5->SetLineColor(6);

        TLine *ll5=new TLine(6.41,0,6.41,10000);
        ll5->SetLineColor(6);
        ll5->SetLineStyle(9);
        ll5->SetLineWidth(2);
        ll5->Draw("same");

        TLatex *tt5=new TLatex(6.38,9000,"^{220}Rn");
        tt5->SetTextColor(6);
        tt5->SetTextAngle(90);
        tt5->SetTextSize(0.03);
        tt5->Draw("same");

	//216Po
	TF1 *ff6=new TF1("ff6","gaus(0)",6.4,7.4);
	ff6->SetParameters(Count/1200,6.91,0.1);
	ff6->SetLineColor(7);

        TLine *ll6=new TLine(6.91,0,6.91,10000);
        ll6->SetLineColor(7);
        ll6->SetLineStyle(9);
        ll6->SetLineWidth(2);
        ll6->Draw("same");

        TLatex *tt6=new TLatex(6.88,9000,"^{216}Po");
        tt6->SetTextColor(7);
        tt6->SetTextAngle(90);
        tt6->SetTextSize(0.03);
        tt6->Draw("same");

	// 212Po
	TF1 *ff7=new TF1("ff7","gaus(0)",8.4,9.2);
//	ff7->SetParLimits(0,0,50);
//	ff7->SetParLimits(1,8.5,9.3);
//	ff7->SetParLimits(2,0.1,0.5);
	ff7->SetParameters(Count/3000,8.82,0.1);
	ff7->SetLineColor(8);

        TLine *ll7=new TLine(8.82,0,8.82,10000);
        ll7->SetLineColor(8);
        ll7->SetLineStyle(9);
        ll7->SetLineWidth(2);
        ll7->Draw("same");

        TLatex *tt7=new TLatex(8.79,9000,"^{212}Po");
        tt7->SetTextColor(8);
        tt7->SetTextAngle(90);
        tt7->SetTextSize(0.03);
        tt7->Draw("same");
	
	TF1 *ffE=new TF1("ffE","ff1+ff2+ff3+ff5+ff6+ff7",4.8,9.2);
	h3->Fit("ffE","R");


	// Fit 216Po by energy cut alone 
	TCanvas *c4=new TCanvas("c4","c4",520,520);
	gStyle->SetOptFit(1011);
     TH1F *h4 = new TH1F("h4","",100,6.7,7.3);
        TString energy_1 = Form("qS1Maxcr/%.3f*5.59>>h4",Mean_1);

        signal_tree->Draw(energy_1,Pulse_Shape&&""&&FV,"");
	h4->SetStats(0);
	h4->SetTitle("^{216}Po Fit by Energy Selection");
	h4->GetXaxis()->SetTitle("Energy(MeV)");
	h4->GetYaxis()->SetTitle("Counts");

	TF1 *fff1 =new TF1("fff1","gaus(0)",6.7,7.3);
	fff1->SetParameters(Count/300,7.0,0.11);
	h4->Fit("fff1","R");
	
	double Count_216Po = (fff1->GetParameter(0))*(fff1->GetParameter(2))*pow(2*3.1415926,0.5)/0.06;
        cout << "Count_216Po= " << Count_216Po << endl;
        cout << "Rate_216Po= " << Count_216Po/time*1000 << endl;
	cout << "Mean 216Po= " << fff1->GetParameter(1) << endl;

	// Fit 212Po with 214Po by Energy Cut alone
	TCanvas *c5=new TCanvas("c5","c5",520,520);
	gStyle->SetOptFit(1011);
     TH1F *h5 = new TH1F("h5","",50,7,10);
	TString energy_2 = Form("qS1Maxcr/%.3f*5.59>>h5",Mean_1);

        signal_tree->Draw(energy_2,Pulse_Shape&&""&&FV,"");
        h5->SetStats(0);
        h5->SetTitle("^{214}Po with ^{212}Po Fit by Energy Selection");
        h5->GetXaxis()->SetTitle("Energy(MeV)");
        h5->GetYaxis()->SetTitle("Counts");
	
	TF1 *fff2=new TF1("fff2","gaus(0)",7.3,8.7);
	TF1 *fff3=new TF1("fff3","gaus(0)",8.7,9.7);
 
	TF1 *fff4=new TF1("fff4","fff2+fff3",7.3,9.6);
	fff4->SetParameters(Count/110,8.1,0.3,Count/1500,9.2,0.3);
	fff4->SetParLimits(5,0.1,0.4);
//	fff4->SetParLimits(4,9.3,9.5);

	h5->Fit("fff4","R");

        double Count_214Po = (fff4->GetParameter(0))*(fff4->GetParameter(2))*pow(2*3.1415926,0.5)/0.06;
        cout << "Count_214Po= " << Count_214Po << endl;
        cout << "Rate_214Po= " << Count_214Po/time*1000 << endl;
	cout << "Mean 214Po= " << fff4->GetParameter(1) << endl;

	double Count_212Po = (fff4->GetParameter(3))*(fff4->GetParameter(5))*pow(2*3.1415926,0.5)/0.06;
        cout << "Count_212Po= " << Count_212Po << endl;
        cout << "Rate_212Po= " << Count_212Po/time*1000 << endl;
	cout << "Mean 212Po= " << fff4->GetParameter(4) << endl;

	// Fit 220Rn with the Gaussian of 218Po of Energy Cuts (also with 212Bi)
	TCanvas *c6=new TCanvas("c6","c6",520,520);
        gStyle->SetOptFit(1011);
     TH1F *h6 = new TH1F("h6","",20,5.8,7);
        TString energy_3 = Form("qS1Maxcr/%.3f*5.59>>h6",Mean_1);

        signal_tree->Draw(energy_3,Pulse_Shape&&""&&FV,"");
        h6->SetStats(0);
        h6->SetTitle("^{218}Po with ^{220}Rn with ^{212}Bi Fit by Energy Selection");
        h6->GetXaxis()->SetTitle("Energy(MeV)");
        h6->GetYaxis()->SetTitle("Counts");

        TF1 *fff5=new TF1("fff5","gaus(0)",6,6.4);
        TF1 *fff6=new TF1("fff6","gaus(0)",6.3,6.6);
	TF1 *ffff1=new TF1("ffff1","gaus(0)",6.6,7.3);

	TF1 *fff7=new TF1("fff7","fff5+fff6+ffff1",6,7.3);
//	fff7->SetParameters(Count/12,6.16,0.1,Count/1100,6.4,0.1,Count/1300,6.21,0.1);
	fff7->SetParameters(Count/12,6.16,0.1,Count/1000,6.4,0.1,Count/1200,7,0.2);
	fff7->SetParLimits(3,0,400);
	fff7->SetParLimits(4,6.3,6.5);
	fff7->SetParLimits(5,0,0.4);
//	fff7->SetParLimits(6,0,10);
//	fff7->SetParLimits(7,6.18,6.3);
//	fff7->SetParLimits(8,0.1,0.4);
	h6->Fit("fff7","R");

        double Count_218Po = (fff7->GetParameter(0))*(ffE->GetParameter(2))*pow(2*3.1415926,0.5)/0.06;
        cout << "Count_218Po= " << Count_218Po << endl;
        cout << "Rate_218Po= " << Count_218Po/time*1000 << endl;
	cout << "Mean 218Po= " << fff7->GetParameter(1) <<endl;

        double Count_220Rn = (fff7->GetParameter(3))*(ffE->GetParameter(5))*pow(2*3.1415926,0.5)/0.06;
        cout << "Count_220Rn= " << Count_220Rn << endl;
        cout << "Rate_220Rn= " << Count_220Rn/time*1000 << endl;
	cout << "Mean 220Rn= " << fff7->GetParameter(4) << endl;

        double Count_212Bi = (fff7->GetParameter(6))*(ffE->GetParameter(8))*pow(2*3.1415926,0.5)/0.06;
        cout << "Count_212Bi= " << Count_212Bi << endl;
        cout << "Rate_212Bi= " << Count_212Bi/time*1000 << endl;
        cout << "Mean 212Bi= " << fff7->GetParameter(7) << endl;

	// Fit All Peaks by using 216Po and 212Po and 218Po and 220Rn Gaussian functions
	TCanvas *c7=new TCanvas("c7","c7",520,520);
	     TH1F *h7 = new TH1F("h7","",100,4,10);
        TString energy_all = Form("qS1Maxcr/%.3f*5.59>>h7",Mean_1);

        signal_tree->Draw(energy_all,Pulse_Shape&&""&&FV,"");
        h7->SetStats(0);
        h7->SetTitle("Alpha Spectrum");
        h7->GetXaxis()->SetTitle("Energy(MeV)");
        h7->GetYaxis()->SetTitle("Counts");
	h7->GetYaxis()->SetRangeUser(0,10000);

        TF1 *ff1=new TF1("ff1","gaus(0)",5.2,5.9);
        ff1->SetLineColor(2);

	ff1->SetParameters(Count/35,5.59,0.15);
	h7->Fit("ff1","R");

        double Count_222Rn = (ff1->GetParameter(0))*(ff1->GetParameter(2))*pow(2*3.1415926,0.5)/0.06;
        cout << "Count_222Rn= " << Count_222Rn << endl;
        cout << "Rate_222Rn= " << Count_222Rn/time*1000 << endl;
	cout << "Mean 222Rn= " << ff1->GetParameter(1) << endl;

	// Draw the Finnal Picture of each Gaussian 
	
	fff1->Draw("same");
	fff1->SetLineColor(7);

        TF1 *fff2=new TF1("fff2","gaus(0)",7.3,8.8);
        fff2->SetParameters(fff4->GetParameter(0),fff4->GetParameter(1),fff4->GetParameter(2));
	fff2->Draw("same");
	fff2->SetLineColor(4);
        TF1 *fff3=new TF1("fff3","gaus(0)",8.4,9.6);
        fff3->SetParameters(fff4->GetParameter(3),fff4->GetParameter(4),fff4->GetParameter(5));
	fff3->Draw("same");
	fff3->SetLineColor(8);

	TF1 *fff5=new TF1("fff5","gaus(0)",5.9,6.4);
	fff5->SetParameters(fff7->GetParameter(0),fff7->GetParameter(1),fff7->GetParameter(2));
        fff5->Draw("same");
        fff5->SetLineColor(3);
        TF1 *fff6=new TF1("fff6","gaus(0)",6.1,6.7);
        fff6->SetParameters(fff7->GetParameter(3),fff7->GetParameter(4),fff7->GetParameter(5));
        fff6->Draw("same");
        fff6->SetLineColor(6);

        ll1->Draw("same");
        tt1->Draw("same");

        ll2->Draw("same");
        tt2->Draw("same");

        ll3->Draw("same");
        tt3->Draw("same");

	ll4->Draw("same");
        tt4->Draw("same");

        ll5->Draw("same");
        tt5->Draw("same");

        ll6->Draw("same");
        tt6->Draw("same");

        ll7->Draw("same");
        tt7->Draw("same");
TFile * out = new TFile("SinAlp_out.root","recreate");
h1->Write();
h2->Write();
h3->Write();
h4->Write();
h5->Write();
h6->Write();
h7->Write();



}
