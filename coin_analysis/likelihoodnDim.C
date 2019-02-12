#include<algorithm>
#include <math.h>
// #include "Math/Minimizer.h"
// #include "Math/Factory.h"
// #include "Math/Functor.h"
// #include "Math/PdfFuncMathCore.h"
//#include <Math/ProbFuncMathCore.h>
//#include "TMath.h"
//#include "TFile.h"
//#include "rootheader.h"

double IngetralHist(TH1F* h1, double x_Min,double x_Max)
{


  TAxis *xaxis = h1->GetXaxis();
  double xmax = xaxis->GetXmax();
  double xmin = xaxis->GetXmin();
  if (x_Min > x_Max){return 0;}

  double x_low = 0;
  double x_up = 0;
  if (x_Min < xmin){x_low = xmin;}
  else if (x_Max > xmax){x_up = xmax;}
  else {
    x_low = x_Min;
    x_up = x_Max;
	
  }
  int n_low = xaxis->FindBin(x_low)-1;
  if (n_low<0) n_low=0;
  int n_up = xaxis->FindBin(x_up)-1;
  double sum=h1->Integral(n_low,n_up);
  return sum;

}

vector<double> Mu_limit(TGraph *g1,double alpha_c)
{
	
  vector<double> x_limit;

  int N = g1->GetN();
  if(N<2){
    cout<<"N<2! "<<endl;
    exit(1);
  }
  double x1;
  double x2;
  double y1 = y1 = g1->GetY()[0];
  double y2;
  double gotone = 0;//when gotone = 0, it's finding lower bound; when gotone = 1, it's finding upper bound; when gotone =2, x_limit is found already
  double x_interval = 0.001;
  double x_low=g1->GetX()[0];
  double x_up =g1->GetX()[N-1];
  if (y1 > alpha_c){gotone = 1;}//check whether the lowest accessible point is within the confidence interval or not.
  int j=0;
  double y_low,y_up;
  while (!(gotone==2)&(j<N-1)){
    x1 = g1->GetX()[j];
    x2 = g1->GetX()[j+1];
    y1 = g1->GetY()[j];
    y2 = g1->GetY()[j+1];
    if(( (y1-alpha_c)*(y2-alpha_c)<0) & (gotone==0)){
      cout<<"lower bound found"<<endl;
      y_low = y1;
      for (double xm =x1;xm<x2;xm+=x_interval){
	//y_low = y1;
	y_up = g1->Eval(xm+x_interval);
	if ((y_low-alpha_c)*(y_up-alpha_c)< 0 ){
	  x_low = xm;
	  gotone = 1;
	  xm = x2+0.01;
	}
	y_low = y_up;
      }
			
    }
		
    if(((y1-alpha_c)*(y2-alpha_c)<0) & (gotone==1)){
      cout<<"upper bound found"<<endl;	
      y_low = y1;
				
      for (double xm1 =x1;xm1<x2;xm1+=x_interval){
	//y_low = y1;
	y_up = g1->Eval(xm1+x_interval);
	if ((y_low-alpha_c)*(y_up-alpha_c)< 0 ){
	  gotone = 2;
	  x_up = xm1;
	  xm1 = x2+0.01;
	}
	y_low = y_up;
      }
    }
    j++;
  }
  x_limit.push_back(x_low);
  x_limit.push_back(x_up);
  return x_limit;
}

double likelihoods(const double *xx1){
  const double mu = xx1[0];
  const double delta =xx1[1];
  const double data = int(xx1[2]);
  const double sig_exp = xx1[3];
  const double bkg_exp = xx1[4];
  const double bkg_exp_error = xx1[5];
  double l =-2*TMath::Log(ROOT::Math::poisson_pdf(data,(mu*sig_exp+bkg_exp*(1+delta)))*ROOT::Math::gaussian_pdf(delta,bkg_exp_error,0)); 
  //double l =1./(ROOT::Math::poisson_pdf(data,(mu*sig_exp+bkg_exp))); 
  return l;
}

double NumericalMinimization(vector<double> xx, vector<double> delta_domain, int fixmu)
{
  const char * minName = "Minuit";
  const char *algoName = "Migrad"  ;
  // create minimizer giving a name and a name (optionally) for the specific
  // algorithm
  // possible choices are:
  //     minName                  algoName
  // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
  //  Minuit2                     Fumili2
  //  Fumili
  //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
  //                              BFGS2, SteepestDescent
  //  GSLMultiFit
  //   GSLSimAn
  //   Genetic
  ROOT::Math::Minimizer* min = 
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  const int npar = 6; // ???
  ROOT::Math::Functor f(&likelihoods,npar); 
  min->SetFunction(f);
  min->SetErrorDef(0.5);
  // set tolerance , etc...
  min->SetMaxFunctionCalls(1000000); // for Minuit/Minuit2
  min->SetMaxIterations(10000);  // for GSL
  min->SetTolerance(0.001);
  //min->SetPrintLevel(2);
  min->SetPrintLevel(0);
  // create funciton wrapper for minmizer
  // a IMultiGenFunction type
  //ROOT::Math::Functor f(&likelihood_binned,npar);
  //ROOT::Math::Functor f(&likelihood_1D,npar); //for both unbinned and binned TH1 
  double step[npar];
  for(int j=0; j< npar; j++){
    step[j]= 0.001;
  }
  // starting point
  double variable[npar];
  //variable[0] = mu;
  for(int j=0; j<npar; j++){
    variable[j] = xx[j]; 
  }
  double mu = xx[0];
  string names[npar];
  names[0]= "mu";
  names[1]="delta";
  names[2] = "data";
  names[3]="sig_exp";
  names[4]="bkg_exp";
  names[5]="bkg_exp_error";
  for(int j=0; j<npar; j++){
    min->SetVariable(j,names[j],variable[j], step[j]);
  }
  for (int j=2;j<npar;j++){
    min->FixVariable(j);
  }
  double delta_min = delta_domain[0];
  double delta_max = delta_domain[1];
  min->SetVariableLimits(1,delta_min,delta_max);
  //min->SetVariable(1,names[1],0, step[1]);
  min->SetVariableValue(1,0);
  //min->FixVariable(1);
  if(fixmu) {
    min->FixVariable(0);  //for numeritaor of  L(d|mu,bhat_mu)/L(d|muhat,bhat)
    //if(mu==0){////no need to fit for eps_s
    //	min->FixVariable(npar-1);
    //}
  }
  else {
    min->SetVariableLimits(0,0,mu); //for denominator
    min->SetVariableValue(0,mu/2);
    if(mu==0){
      min->FixVariable(0);
    }
  }
  //for(int j=1; j<npar-1; j++){
  //	min->SetVariableLimits(j,max(0.,1-5*relBErr[j-1]),1+5*relBErr[j-1]);
  //}
  //min->SetVariableLimits(npar-1,max(0.,1-5*relSErr),1+5*relSErr);
  
  //cout<<"soldier m0"<<endl;
  // do the minimization
  min->Minimize();
  //cout<<"soldier  m1"<<endl;
  if(min->Status() !=0){
    //cout<<"min status "<< min->Status() <<endl;
    //min->SetPrintLevel(2);
    if(min->X()[0]<0.01){
      min->SetVariableLimits(0,0,0.01); //for denominator
    }
    min->SetStrategy(2);
    min->Minimize();
    //     if(min->Status() !=0){
    //       //cout<<"min failed to converge with stra=2 status="<<  min->Status() <<endl;
    //       status = min->Status();
    //     }
  }
  //	int status = min->Status();
  const double *xs = min->X();
  double l = likelihoods(xs); 
  return l;
}

void likelihoodnDim(){
  //	unsigned int Ndata = 15367;
  unsigned int Ndata = 1;

  double N_bkg_exp = 0.54; //  bkg level before injection
  // 7000h
  // same cut signal bkg
  //double N_bkg_exp_error = TMath::Sqrt(N_bkg_exp);
  double N_bkg_exp_error = 0.12;
  //double N_bkg_exp_error = 200;
  double N_bkg_exp_error_percent = N_bkg_exp_error/N_bkg_exp;
  double N_sig_exp = 1;
  double alpha = 0.1;
  vector<double> mu_test;
  vector<double> Cls_Ndata;
  vector<double> q_mu_Ndata;
	
  TRandom3 rnd(0);
  gRandom->SetSeed(0);

  double mu_temp, bkg_temp, sig_temp, data_temp, q_mu_data,cls_Ndata;
  vector<double> confidenceI;

  TFile file1("fq_mu.root","recreate");
  double q_mu,q_mu_bkg,l_mu,l_mu_bkg,l_mu_c,mu_c;
  double q_mu_max = 0;

  int N = 20;//number of mu calculated

  double mu_min = 1;  // 0 [count]
  double mu_max = 20;  // 3 sigma [count]
  double interval_mu = (mu_max-mu_min)/(N-1);
	
  int NMC = 20000;//number of simulated experiments for each mu
  double qu_temp[500]={0};
  double norm = NMC;
  double bin_max=2*mu_max;
  double bin_min=0;
  int n_bin=500;
  double bin_w = (bin_max-bin_min)/n_bin;
	
  int nPar = 6;
  vector<double> xxL;
  xxL.push_back(0);
  xxL.push_back(0);
  xxL.push_back(Ndata);
  xxL.push_back(N_sig_exp);
  xxL.push_back(N_bkg_exp);
  xxL.push_back(N_bkg_exp_error_percent);

  double delta_min = -1,delta_max=3;
  vector<double> delta_domain;
  delta_domain.push_back(delta_min);
  delta_domain.push_back(delta_max);


  double lde,lno;

  //	cout<<"soldier   0"<<endl;
  //	ROOT::Math::Minimizer *min = NumericalMinimization(xxL,delta_domain,true);
  double q1;
  TString pdffilename = "PDFcomparison.pdf";
  TCanvas *c0 = new TCanvas("c0","c0",800,600);
  c0->Print(pdffilename+"[");
  for(int k=0;k<N;k++){
    mu_temp = 1.*k*interval_mu +mu_min;
    //double mu_temp_up = mu_temp+0.001;
    mu_test.push_back(mu_temp);
		
    xxL[0] = mu_temp;
    xxL[2] = Ndata;
    //double soldier = xxL[0];
    //cout<<soldier<<endl;
    lde =  NumericalMinimization(xxL,delta_domain,true);
    lno =  NumericalMinimization(xxL,delta_domain,false);
    q_mu_data = lde - lno;

		
    TH1F* mu_f_q_mu=new TH1F("mu_tenp","h1",n_bin,bin_min,bin_max);
    TH1F* mu_f_q_mu_bkg=new TH1F("bkg_temp","h1",n_bin,bin_min,bin_max);
    for (int l = 0;l<NMC;l++){
			
      bkg_temp =rnd.Gaus(N_bkg_exp,N_bkg_exp_error);
      sig_temp = 1;
      if (bkg_temp<0){bkg_temp=0;}
      //bkg_temp = 1;	
      //pdf for signal hypothesis
      data_temp = rnd.Poisson(mu_temp*sig_temp+bkg_temp);
      xxL[2] = data_temp;
      xxL[1] = (bkg_temp-N_bkg_exp)/N_bkg_exp;
      lde = NumericalMinimization(xxL,delta_domain,true);
      lno =  NumericalMinimization(xxL,delta_domain,false);
      q_mu = lde-lno;
      mu_f_q_mu->Fill(q_mu);

			
      bkg_temp =rnd.Gaus(N_bkg_exp,N_bkg_exp_error);
      sig_temp = 1;
      if (bkg_temp<0){bkg_temp=0;}
      //bkg_temp = 1;	
      //pdf for bkg hypothesis
      data_temp = rnd.Poisson(bkg_temp);
      xxL[2] = data_temp;
      xxL[1] = (bkg_temp-N_bkg_exp)/N_bkg_exp;
      lde = NumericalMinimization(xxL,delta_domain,true);
      lno =  NumericalMinimization(xxL,delta_domain,false);
      q_mu = lde-lno;
      mu_f_q_mu_bkg->Fill(q_mu);

			
    }
				
    mu_f_q_mu->Scale(1./norm);
    mu_f_q_mu_bkg->Scale(1./norm);
    double norm_test = mu_f_q_mu->Integral();
    double norm_test_bkg = mu_f_q_mu_bkg->Integral();
    cout<<"norm_test:\t"<<norm_test<<"; \tnorm_test_bkg:\t"<<norm_test_bkg<<endl;
    mu_f_q_mu->Draw("hist");
    mu_f_q_mu_bkg->Draw("same l");
    //mu_f_q_mu_bkg->Draw("hist");
    //mu_f_q_mu->Draw("same");
    mu_f_q_mu_bkg->SetLineColor(2);
    TLegend *lg0 = new TLegend(0.7,0.9,0.9,0.98);
    lg0->AddEntry(mu_f_q_mu,"f(q_{#mu}|#mu)","pl");
    lg0->AddEntry(mu_f_q_mu_bkg,"f(q_{#mu}|0)","pl");
    lg0->Draw();
    c0->Print(pdffilename);

    double clsPlusb = IngetralHist(mu_f_q_mu,q_mu_data,bin_max);
    double clb = IngetralHist(mu_f_q_mu_bkg,q_mu_data,bin_max);
    cls_Ndata = clsPlusb/clb;
		
		
		
		
    if (clb==0)cls_Ndata=-1;
    Cls_Ndata.push_back(cls_Ndata);
    q_mu_Ndata.push_back(q_mu_data);
    lg0->Delete();	
    mu_f_q_mu->Delete();
    mu_f_q_mu_bkg->Delete();
    cout<<"  mu:"<<mu_temp<<"q_mu_obs:"<<q_mu_data<<"  CLs+b:"<<clsPlusb<<"  CLb:"<<clb<<"  CLs:"<<cls_Ndata<<endl;		
	
  }
  c0->Print(pdffilename+"]");
  TCanvas *c1 = new TCanvas("c1","multipads",900,700);
  //c1->SetLogy();	
  TGraph *g2 = new TGraph(N,&mu_test[0],&Cls_Ndata[0]);
  confidenceI = Mu_limit(g2,alpha);
  cout<<confidenceI[1]<<endl;
  g2->SetLineColor(2);
  g2->GetXaxis()->SetTitle("#mu");
  g2->GetYaxis()->SetTitle("CL_{s}");
  g2->Draw("apl");
  c1->Print("test.png");



	
}
