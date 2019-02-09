// int kqS1Low = 40e3;
// int kqS1High = 100e3;

void alpha_peaks_plot(int & RunNo)
{
  TChain * ch = new TChain("signal_tree");
  ch->Add(Form("/store/px/data/udm/by-run/%d/new_correction/ana*",RunNo));
  // qS1 vs z 
  TCanvas * c1 = new TCanvas("alpha_band_after_tuning","",1);
  ch->Draw("tS1[s1max]-tS2[s2max]:qS1[s1max]*exp((tS1[s1max]-tS2[s2max]+150000)/1200000)/(1-1E-6*(xS2NN[s2max]^2+yS2NN[s2max]^2))>>h0(300,0,120e3,300,-270e3,-30e3)","s1max>=0&&s2max>=0&&qS1[s1max]>=3&&qS2B[s2max]>30&&tS2[s2max]-tS1[s1max]<250e3&&tS2[s2max]-tS1[s1max]>50e3&&xS2NN[s2max]^2+yS2NN[s2max]^2<6e4","colz");
  // qS1 histogram
  TCanvas * c2 = new TCanvas("alpha_spectrum","",1);
  ch->Draw("qS1[s1max]*exp((tS1[s1max]-tS2[s2max]+150000)/1200000)/(1-1E-6*(xS2NN[s2max]^2+yS2NN[s2max]^2))>>h1(100, 40e3, 60e3)","s1max>=0&&s2max>=0&&qS1[s1max]>=3&&qS2B[s2max]>30&&tS2[s2max]-tS1[s1max]<250e3&&tS2[s2max]-tS1[s1max]>50e3&&xS2NN[s2max]^2+yS2NN[s2max]^2<6e4","");
  printf("hint:\n");
  printf("fastGausFit(h1,44e3,48e3,48e3,54e3);\n");
  printf("gausCalculator(constant,sigma,time[h]);\n");
}

void fastGausFit(TH1F * h, int a, int b)
{
  h->Fit("gaus","+","",a,b);
}

void fastGausFit(TH1F * h, int a, int b, int c, int d)
{
  h->Fit("gaus","+","",a,b);
  h->Fit("gaus","+","",c,d);
}

double gausCalculator(double constant, double sigma)
{
  return constant*sigma*sqrt(TMath::Pi()*2);
}

double gausCalculator(double constant, double sigma, double time)
{
  // scale is the bin width 
  // e.g. 100./(60e3-40e3)
  double scale = 100./(60e3-40e3);
  time = time*3600;
  return constant*sigma*sqrt(TMath::Pi()*2)*scale/time;
}
