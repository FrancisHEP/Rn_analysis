void plot_evolution_exposure_elifetime_Run10Run11(const char * inputfilename1, const char * inputfilename2)
{
  ifstream fin1(inputfilename1);
  int runno[1000];
  int year[1000], month[1000], day[1000], hour[1000], min[1000], sec[1000];
  double livetime[1000];
  double exposure[1000], totalexposure[1000];
  double totalexposure_error[1000] = {0};
  double time_error[1000] = {0};
  TDatime date[1000];
  double time[1000];
  double t1,t2,t3,t4,t5,t6,t7;
  double t8,t9,t10;
  int i(0);
  while(fin1>>t1>>t2>>t3>>t4>>t5>>t6>>t7>>t8>>t9>>t10)
  {
    runno[i]=t1;
    year[i]=t2;
    month[i]=t3;
    day[i]=t4;
    hour[i]=t5;
    min[i]=t6;
    sec[i]=t7;
    livetime[i]=t8/24.+t9/60./24.+t10/60./60./24.; //unit: day
    exposure[i]=361.5/1000.*livetime[i]; //unit: ton*day
    if(i==0) totalexposure[i] = 328.9/1000.*79.6+exposure[i];
    else totalexposure[i]=totalexposure[i-1]+exposure[i];
    date[i]=TDatime(year[i],month[i],day[i],hour[i],min[i],sec[i]);
    time[i]=date[i].Convert();
    //cout<<runno[i]<<"  "<<time[i]<<"  "<<totalexposure[i]<<endl;
    i++;
  }
  cout<<runno[i-1]<<"  "<<time[i-1]<<"  "<<totalexposure[i-1]<<endl;
  
  double leftmax = 155.;
  double rightmax = 1550.;
  double scale = leftmax/rightmax; 

  ifstream fin2(inputfilename2);
  double elifetime[1000];
  double elifetime_error[1000] = {0};
  double time2[1000];
  double time2_error[1000] = {0};
  int j(0);
  while(fin2>>t1>>t2>>t3>>t4>>t5>>t6>>t7>>t8>>t9)
  {
    runno[j]=t1;
    year[j]=t2;
    month[j]=t3;
    day[j]=t4;
    hour[j]=t5;
    min[j]=t6;
    sec[j]=t7;
    elifetime[j]=t8; //unit: us or ns
    elifetime_error[j]=t9; //unit: us or ns
    if (elifetime[j]>10000)
    {
      elifetime[j]/=1000.; //unit: us
      elifetime_error[j]/=1000.; //unit: us
    }
    elifetime[j]*=scale; //scale e-lifetime(right) to the left coordinates
    elifetime_error[j]*=scale; //scale e-lifetime(right) to the left coordinates

    time2[j]=TDatime(year[j],month[j],day[j],hour[j],min[j],sec[j]).Convert();
    //cout<<runno[j]<<"  "<<time2[j]<<"  "<<elifetime[j]<<endl;
    j++;
  }

  TCanvas * c1 = new TCanvas("c1","c1",1200,600);
  c1->cd();
  c1->SetTickx();
  //c1->SetGridx();
  //c1->SetGridy();
  c1->SetTopMargin(0.05);
  c1->SetBottomMargin(0.13);
  c1->SetLeftMargin(0.09);
  c1->SetRightMargin(0.09);
  gStyle->SetOptFit(1011);
  TGraph * g1 = new TGraphErrors(i,time,totalexposure,time_error,totalexposure_error);
  g1->SetTitle("");
  g1->SetMarkerStyle(25);
  g1->SetMarkerSize(0.5);
  g1->SetMarkerColor(kAzure+10);
  g1->SetLineWidth(1.5);
  g1->SetLineColor(kAzure+10);

  int xmin = TDatime(2017,03,01,00,00,00).Convert();
  int xmax = TDatime(2018,09,10,00,00,00).Convert();
  g1->GetXaxis()->SetLimits(xmin,xmax);
  g1->GetXaxis()->SetLabelOffset(0.04);
  g1->GetXaxis()->SetLabelFont(132);
  g1->GetXaxis()->SetTimeDisplay(1);
  g1->GetXaxis()->SetTimeFormat("#splitline{%b.%d}{%Y}");
  //g1->GetXaxis()->SetTimeFormat("%b.%d");
  g1->GetXaxis()->SetTimeOffset(0,"cst");
  g1->GetXaxis()->SetTickLength(0.02);
  
  double ymin = 0;
  double ymax = 155; 
  g1->GetYaxis()->SetRangeUser(ymin,ymax);
  g1->GetYaxis()->SetTitleFont(132);
  g1->GetYaxis()->SetTitleSize(0.06);
  g1->GetYaxis()->SetTitleOffset(0.65);
  g1->GetYaxis()->CenterTitle();
  g1->GetYaxis()->SetTitle("Exposure [ton-day]");
  g1->GetYaxis()->SetTickLength(0.01);
  g1->Draw("APL");

  TGraph * g2 = new TGraphErrors(j,time2,elifetime,time2_error,elifetime_error);
  g2->SetMarkerStyle(20);
  g2->SetMarkerSize(0.7);
  g2->SetMarkerColor(kRed);
  g2->SetLineStyle(1);
  g2->SetLineWidth(2);
  g2->SetLineColor(kRed);
  g2->Draw("P");

  //draw an axis on the right side
  //TGaxis *axis = new TGaxis(c1->GetUxmax(),c1->GetUymin(),
                            //c1->GetUxmax(),c1->GetUymax(),0,rightmax,510,"+L");
  TGaxis *axis = new TGaxis(xmax,ymin,xmax,ymax,0,rightmax,510,"+L");
  axis->SetLineColor(kRed);
  axis->SetLabelColor(kRed);
  //axis->SetTickSize(0.01);
  axis->SetTitleFont(132);
  axis->SetTitleSize(0.06);
  axis->SetTitleColor(kRed);
  axis->SetTitleOffset(0.7);
  axis->CenterTitle();
  axis->SetTitle("e^{-} Lifetime [#mus]");
  axis->Draw();


  TLine *ly_Run9 = new TLine(TDatime(2016,09,22,21,02,23).Convert(),ymin,TDatime(2016,09,22,21,02,23).Convert(),ymax);
  TLine *ly_Run10 = new TLine(TDatime(2017,07,16,10,42,13).Convert(),ymin,TDatime(2017,07,16,10,42,13).Convert(),ymax);
  ly_Run9->SetLineStyle(1);
  ly_Run9->SetLineWidth(1);
  ly_Run9->SetLineColor(1);
  ly_Run9->Draw("same");
  ly_Run10->SetLineStyle(1);
  ly_Run10->SetLineWidth(1);
  ly_Run10->SetLineColor(1);
  ly_Run10->Draw("same");

  TLine *lx_Run10 = new TLine(xmin,54,TDatime(2017,07,16,10,42,13).Convert(),54);
  lx_Run10->SetLineStyle(2);
  lx_Run10->SetLineWidth(1);
  //lx_Run10->SetLineColor(kAzure+10);
  lx_Run10->SetLineColor(kBlack);
  lx_Run10->Draw("same");

  //TLine *lx_Run11 = new TLine(xmin,139.165,TDatime(2018,07,25,08,09,01).Convert(),139.165);
  TLine *lx_Run11 = new TLine(xmin,146.288,TDatime(2018,09,04,08,30,11).Convert(),146.288); //final exposure line
  lx_Run11->SetLineStyle(2);
  lx_Run11->SetLineWidth(1);
  lx_Run11->SetLineColor(kBlack);
  lx_Run11->Draw("same");


  //run condition 
  //Run10+Run11 periods: Run10DM(19304,20152),Run11DM(20200,20401),Run11AmBeCo(20402,20904),Run11DM(20905,21199),AmBe/ZLE(21200,21302),Run11DM(21303,21778),Rn220(21779,22102),Run11DM(22103,22279),AmBe/ZLE(22280,22473),Run11DM(22474,23546)
  int n = 2;
  double yMax[2] = {ymax,ymax};
  double yMin[2] = {ymin,ymin};
  int x[2] = {0};
  //int const nprd = 14;
  //int T_start[nprd] = {TDatime(2017,04,22,10,15,25).Convert(),TDatime(2017,07,17,17,58,06).Convert(),TDatime(2017,08,06,15,40,03).Convert(),TDatime(2017,08,27,11,42,27).Convert(),TDatime(2017,09,16,08,36,15).Convert(),TDatime(2017,09,24,11,47,42).Convert(),TDatime(2017,11,21,11,00,54).Convert(),TDatime(2017,12,28,12,51,09).Convert(),TDatime(2018,01,19,12,32,21).Convert(),TDatime(2018,02,13,10,00,31).Convert()};
  //int T_end[nprd]  =  {TDatime(2017,07,16,10,42,13).Convert(),TDatime(2017,08,06,15,37,51).Convert(),TDatime(2017,08,27,11,30,14).Convert(),TDatime(2017,09,16,08,35,41).Convert(),TDatime(2017,09,24,11,46,38).Convert(),TDatime(2017,11,21,10,58,23).Convert(),TDatime(2017,12,28,12,50,39).Convert(),TDatime(2018,01,19,12,31,52).Convert(),TDatime(2018,02,13,09,59,12).Convert(),TDatime(2018,07,25,08,09,01).Convert()};
  int const nprd = 4;
  int T_start[nprd] = {TDatime(2017,08,06,15,40,03).Convert(), TDatime(2017,09,16,08,36,15).Convert(), TDatime(2017,11,21,11,00,54).Convert(), TDatime(2018,01,19,12,32,21).Convert()};
  int T_end[nprd]  =  {TDatime(2017,08,27,11,30,14).Convert(), TDatime(2017,09,24,11,46,38).Convert(), TDatime(2017,12,28,12,50,39).Convert(), TDatime(2018,02,13,09,59,12).Convert()};
  TGraph *grshade[nprd];
  for (int jj = 0; jj<nprd;jj++){
    x[0] = T_start[jj];
    x[1] = T_end[jj];
    grshade[jj] = new TGraph(2*n);
    for (int k=0; k<n ;k++) {
      grshade[jj]->SetPoint(k,x[k],yMax[k]);
      grshade[jj]->SetPoint(n+k,x[n-k-1],yMin[n-k-1]);
    }
    //if (jj%2 ==0) grshade[jj]->SetFillStyle(3004);
    //else grshade[jj]->SetFillStyle(3022);
    grshade[jj]->SetFillStyle(3004);
    grshade[jj]->SetFillColor(8);
    grshade[jj]->Draw("fsame");
  }


 
/*
  TString sa,sb,sc;
  sa.Form("#tau1 = %.2f day",-1.0/para[1]/86400.);
  sb.Form("#tau2 = %.2f day",-1.0/parb[1]/86400.);
  sc.Form("#tau3 = %.2f day",-1.0/parc[1]/86400.);
  TLegend *leg1=new TLegend(0.4,0.5,0.6,0.7);
  leg1->SetTextFont(42);
  leg1->AddEntry(fa,sa,"l");
  leg1->AddEntry(fb,sb,"l");
  leg1->AddEntry(fc,sc,"l");
  leg1->Draw("same");
*/

  c1->Print("./figures/evolution_exposure_elifetime_Run10Run11_final.png");

}
