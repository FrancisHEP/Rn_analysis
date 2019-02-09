#include <iostream>
#include <TString.h>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#include <TChain.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCut.h>
#include "TROOT.h"

#include "TROOT.h"
#include "TBrowser.h"
#include <fstream>
#include <vector>
#include <algorithm>

#include <TChain.h>
#include <TFileCollection.h>


void Pic_Rn220(int runnumber){
	gStyle->SetOptFit(1111);
	    TChain *chain;
    chain = new TChain("rn220_AA");
	chain->Add(Form("outTree_220_216_%d.root",runnumber));

	int Time = 6569642;
    int runNo, trigNo;
    double triggerTime;

    unsigned int nS1, nS2;
    int s1max, s2max;
    double qS1T[1000], qS1B[1000], xS1T[1000], yS1T[1000];
    double xS1B[1000], yS1B[1000];
    double qS1[1000], tS1[1000], pS1[1000], hS1[1000] ;
    double ttenS1[1000], wS1[1000], wtenS1[1000];
    double xS2PAF[1000], yS2PAF[1000];
    double qS1Veto[1000];
    int S1NPeaks[1000], nS1BSatur[1000], nS1TSatur[1000];

    double qS2T[1000], qS2B[1000], xS2T[1000], yS2T[1000];
    double xS2B[1000], yS2B[1000];
    double qS2[1000], tS2[1000], pS2[1000], hS2[1000] ;
    double ttenS2[1000], wS2[1000], wtenS2[1000];
    double xS2NN[1000], yS2NN[1000], chi2[1000];
    int nS2BSatur[1000], nS2TSatur[1000];
	int alpha1, alpha2;
	double deltaT;

    int nS1_pre, nS2_pre, s1max_pre, s2max_pre;
    double xS2T_pre[1000], yS2T_pre[1000], xS2B_pre[1000], yS2B_pre[1000], xS2NN_pre[1000], yS2NN_pre[1000],xS2PAF_pre[1000], yS2PAF_pre[1000];
    double chi2_pre[1000];
    double tS1_pre[1000], ttenS1_pre[1000], qS1T_pre[1000], qS1B_pre[1000], qS1_pre[1000], wS1_pre[1000], wtenS1_pre[1000], hS1_pre[1000], pS1_pre[1000];
    double tS2_pre[1000], ttenS2_pre[1000], qS2T_pre[1000], qS2B_pre[1000], qS2_pre[1000], wS2_pre[1000], wtenS2_pre[1000], hS2_pre[1000], pS2_pre[1000];
    chain->SetBranchAddress("trigNo", &trigNo);
    chain->SetBranchAddress("runNo", &runNo);
    chain->SetBranchAddress("triggerTime", &triggerTime);
    chain->SetBranchAddress("nS1", &nS1);
    chain->SetBranchAddress("s1max", &s1max);
    chain->SetBranchAddress("qS1T", qS1T);
    chain->SetBranchAddress("qS1B", qS1B);
    chain->SetBranchAddress("qS1", qS1);
    chain->SetBranchAddress("tS1", tS1);
    chain->SetBranchAddress("ttenS1", ttenS1);
    chain->SetBranchAddress("wS1", wS1);
    chain->SetBranchAddress("wtenS1", wtenS1);
    chain->SetBranchAddress("hS1", hS1);
    chain->SetBranchAddress("pS1", pS1);
    chain->SetBranchAddress("nS2", &nS2);
    chain->SetBranchAddress("s2max", &s2max);
    chain->SetBranchAddress("xS2NN", xS2NN);
    chain->SetBranchAddress("xS2PAF", xS2PAF);
    chain->SetBranchAddress("yS2NN", yS2NN);
    chain->SetBranchAddress("yS2PAF", yS2PAF);
    chain->SetBranchAddress("chi2", chi2);
    chain->SetBranchAddress("qS2T", qS2T);
    chain->SetBranchAddress("qS2B", qS2B);
    chain->SetBranchAddress("qS2", qS2);
    chain->SetBranchAddress("tS2", tS2);
    chain->SetBranchAddress("ttenS2", ttenS2);
    chain->SetBranchAddress("wS2", wS2);
    chain->SetBranchAddress("wtenS2", wtenS2);
    chain->SetBranchAddress("hS2", hS2);
    chain->SetBranchAddress("pS2", pS2);
    chain->SetBranchAddress("alpha1", &alpha1);
    chain->SetBranchAddress("alpha2", &alpha2);
    chain->SetBranchAddress("deltaT", &deltaT);

    chain->SetBranchAddress("nS1_pre", &nS1_pre);
    chain->SetBranchAddress("s1max_pre", &s1max_pre);
    chain->SetBranchAddress("qS1T_pre", qS1T_pre);
    chain->SetBranchAddress("qS1B_pre", qS1B_pre);
    chain->SetBranchAddress("qS1_pre", qS1_pre);
    chain->SetBranchAddress("tS1_pre", tS1_pre);
    chain->SetBranchAddress("ttenS1_pre", ttenS1_pre);
    chain->SetBranchAddress("wS1_pre", wS1_pre);
    chain->SetBranchAddress("wtenS1_pre", wtenS1_pre);
    chain->SetBranchAddress("hS1_pre", hS1_pre);
    chain->SetBranchAddress("pS1_pre", pS1_pre);
    chain->SetBranchAddress("nS2_pre", &nS2_pre);
    chain->SetBranchAddress("s2max_pre", &s2max_pre);
    chain->SetBranchAddress("xS2NN_pre", xS2NN_pre);
    chain->SetBranchAddress("xS2PAF_pre", xS2PAF_pre);
    chain->SetBranchAddress("yS2NN_pre", yS2NN_pre);
    chain->SetBranchAddress("yS2PAF_pre", yS2PAF_pre);
    chain->SetBranchAddress("chi2_pre", chi2_pre);
    chain->SetBranchAddress("qS2T_pre", qS2T_pre);
    chain->SetBranchAddress("qS2B_pre", qS2B_pre);
    chain->SetBranchAddress("qS2_pre", qS2_pre);
    chain->SetBranchAddress("tS2_pre", tS2_pre);
    chain->SetBranchAddress("ttenS2_pre", ttenS2_pre);
    chain->SetBranchAddress("wS2_pre", wS2_pre);
    chain->SetBranchAddress("wtenS2_pre", wtenS2_pre);
    chain->SetBranchAddress("hS2_pre", hS2_pre);
    chain->SetBranchAddress("pS2_pre", pS2_pre);


	TCut FV = "tS2_pre[s2max_pre]-tS1_pre[alpha1]>20000&&tS2_pre[s2max_pre]-tS1_pre[alpha1]<350000&&tS2[s2max]-tS1[alpha2]>20000&&tS2[s2max]-tS1[alpha2]<350000";
	TCut R2 = "xS2PAF[s2max_pre]*xS2PAF[s2max_pre]+yS2PAF[s2max_pre]*yS2PAF[s2max_pre]<72000&&xS2PAF[s2max]*xS2PAF[s2max]+yS2PAF[s2max]*yS2PAF[s2max]<72000";
	// Depth:TBA cuts

	TCanvas *c1=new TCanvas("c1","c1",520,520);
    TH2F * h1 = new TH2F("h1","",100,-1,0,100,-400,0);
	chain->Draw("(tS1_pre[alpha1]-tS2_pre[s2max_pre])/1000:(qS1T_pre[alpha1]-qS1B_pre[alpha1])/(qS1T_pre[alpha1]+qS1B_pre[alpha1])>>h1",FV,"");
	h1->SetTitle("Depth:TBA^{#alpha}");
	h1->GetXaxis()->SetTitle("TBA^{#alpha}");
	h1->GetYaxis()->SetTitle("Depth[us]");
	h1->SetMarkerColor(2);
	h1->Draw();
    TH2F * h11 = new TH2F("h11","",100,-1,0,100,-400,0);
	chain->Draw("(tS1[alpha2]-tS2[s2max])/1000:(qS1T[alpha2]-qS1B[alpha2])/(qS1T[alpha2]+qS1B[alpha2])>>h11",FV,"same");
	h11->SetMarkerColor(4);
	h1->Draw();
	h11->Draw("same");


	TF1 *f1=new TF1("f1","[0]+[1]*x",-0.66,-0.2);
	f1->SetParameters(50,400);
	h1->Fit("f1","R");

        double a = f1->GetParameter(1);
        double b = f1->GetParameter(0);

        TString Range1 = Form("%.2f*x+%.2f+50",a,b);
        TString Range2 = Form("%.2f*x+%.2f-50",a,b);

        TF1 *f4=new TF1("f4",Range1,-1,0);
        f4->SetLineStyle(2);
        f4->SetLineColor(4);
        f4->Draw("same");
        TF1 *f5=new TF1("f5",Range2,-1,0);
        f5->SetLineStyle(2);
        f5->SetLineColor(4);
        f5->Draw("same");

	delete f1;

	// delata Z cuts
	TCanvas *c3=new TCanvas("c3","c3",520,520);
    TH1F * h3 = new TH1F("h3","",100,-100,100);
	TCut Cut3 = Form("(-tS2_pre[s2max_pre]+tS1_pre[alpha1])/1000<%.2f*(qS1T_pre[alpha1]-qS1B_pre[alpha1])/(qS1T_pre[alpha1]+qS1B_pre[alpha1])+%.2f+50&&(-tS2_pre[s2max_pre]+tS1_pre[alpha1])/1000>%.2f*(qS1T_pre[alpha1]-qS1B_pre[alpha1])/(qS1T_pre[alpha1]+qS1B_pre[alpha1])+%.2f-50&&(-tS2[s2max]+tS1[alpha2])/1000<%.2f*(qS1T[alpha2]-qS1B[alpha2])/(qS1T[alpha2]+qS1B[alpha2])+%.2f+50&&(-tS2[s2max]+tS1[alpha2])/1000>%.2f*(qS1T[alpha2]-qS1B[alpha2])/(qS1T[alpha2]+qS1B[alpha2])+%.2f-50",a,b,a,b,a,b,a,b);

	chain->Draw("((-tS2_pre[s2max_pre]+tS1_pre[alpha1])/1000)-((-tS2[s2max]+tS1[alpha2])/1000)>>h3",FV&&Cut3,"");
	h3->SetTitle("Delta Z");
        h3->GetXaxis()->SetTitle("DZ[us]");
        h3->GetYaxis()->SetTitle("Count");
	double ymax = 500;

	double min = -10;
	double max = 10;
	TLine *l1=new TLine(min,0,min,ymax);
	l1->SetLineStyle(2);
	l1->SetLineColor(2);
	l1->Draw("same");
        TLine *l2=new TLine(max,0,max,ymax);
        l2->SetLineStyle(2);
        l2->SetLineColor(2);
        l2->Draw("same");


	// Position distribution of alpha1 and alpah2
	TCanvas *c4=new TCanvas("c4","c4",520,520);
    TH2F * h4 = new TH2F("h4","",1000,0,100000,100,-360,0);
	TCut Zcut = Form("((-tS2_pre[s2max_pre]+tS1_pre[alpha1])/1000)-((-tS2[s2max]+tS1[alpha2])/1000)>%.2f&&((-tS2_pre[s2max_pre]+tS1_pre[alpha1])/1000)-((-tS2[s2max]+tS1[alpha2])/1000)<%.2f",min,max);
	chain->Draw("(-tS2_pre[s2max_pre]+tS1_pre[alpha1])/1000:xS2PAF_pre[s2max_pre]^2+yS2PAF_pre[s2max_pre]^2>>h4",FV&&Cut3&&Zcut,"");
	//h4->SetStats(0);
	h4->SetTitle("Position Distribution");
        h4->GetXaxis()->SetTitle("R^{2}");
        h4->GetYaxis()->SetTitle("Dt[us]");
	h4->SetMarkerColor(2);
	h4->SetMarkerStyle(4);
	h4->Draw();
    TH2F * h42 = new TH2F("h42","",1000,0,100000,100,-360,0);    
    h42->SetMarkerStyle(6);
	h42->SetMarkerColor(4);
	chain->Draw("(-tS2[s2max]+tS1[alpha2])/1000:xS2PAF[s2max]^2+yS2PAF[s2max]^2>>h42",FV&&Cut3&&Zcut,"same");
	h4->Draw();
	h42->Draw("same");
	//TLegend *lgh4=new TLegend()
	c4->Print(Form("zR2_run%d",runnumber));


	//xy distribution
	TCanvas *cDis=new TCanvas("cDis","cDis",700,700);
	TH2F *hxy=new TH2F("hxy","hxy",100,-350,350,100,-350,350);
	chain->Draw("yS2PAF[s2max]:xS2PAF[s2max]>>hxy",FV&&Cut3&&Zcut&&R2,"colz");
	hxy->SetTitle(Form("x-y distribution run%d",runnumber));
	hxy->GetXaxis()->SetTitle("x/mm");
	hxy->GetYaxis()->SetTitle("y/mm");
	hxy->GetYaxis()->SetTitleOffset(1.3);
	hxy->Draw("colz");
	cDis->Print(Form("xy_run%d",runnumber));
	

	
	// Decay time check the result
	TCanvas *c2=new TCanvas("c2","c2",520,520);
    TH1F * h2 = new TH1F("h2","",100,0,1200);
	chain->Draw("deltaT/1000000>>h2",FV&&Cut3&&Zcut&&R2,"");
	h2->SetTitle("Decay time");
	h2->GetXaxis()->SetTitle("Dt[ms]");
	h2->GetYaxis()->SetTitle("Count");
	int counts = h2->GetEntries();
	cout << "Count= " << counts << endl;
	

	TF1 *f2=new TF1("f2","[0]*exp(-x/[1])",40,1200);
	f2->SetParameters(counts/20,202);
	h2->Fit("f2","R");

	TLatex *t1=new TLatex(600,20,Form("Dt=%.2f#pm%.2fms",f2->GetParameter(1),f2->GetParError(1)));
	t1->SetTextColor(2);
	t1->SetTextSize(0.04);
	t1->Draw("same");
	TLatex *t2=new TLatex(600,30,Form("%.2fmBq in FV",counts/0.9942/Time*1000));
        t2->SetTextColor(2);
        t2->SetTextSize(0.04);
        t2->Draw("same");
c1->SaveAs("Rn220c1.pdf");
c1->SaveAs("Rn220c1.root");
c2->SaveAs("Rn220c2.pdf");
c2->SaveAs("Rn220c2.root");
c3->SaveAs("Rn220c3.pdf");
c3->SaveAs("Rn220c3.root");
c4->SaveAs("Rn220c4.pdf");
c4->SaveAs("Rn220c4.root");


}
