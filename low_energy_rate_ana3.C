void low_energy_rate_ana3(int runNo)
{
  TChain ch("out_tree");
  ch.Add(Form("ana3/ana3_run%d.root",runNo));
  // ch.SetAlias("QualityCut","1==1");
  ch.SetAlias("QualityCut","S1NoiseFilter&&S1Coincidence&&S2Asy&&S1Pattern&&TimeDiff>0.01&&qS1V==0&&nPostCutS1<=2&&log10(qS2/qS1)>1.074-0.595*exp(-qS1/6.38)");

  int low_energy_entries = ch.GetEntries("xTopNN^2+yTopNN^2<72000&&tS2-tS1>20e3&&tS2-tS1<350e3&&qS1>3&&qS1<100&&qS2Top+qS2Bottom>100&&qS2<10000&&QualityCut");
  printf("qS1[3,100] entries is %d\n",low_energy_entries);
  ch.Draw("log10(qS2/qS1):qS1>>(100,0,100,80,0,4)","xTopNN^2+yTopNN^2<72000&&tS2-tS1>20e3&&tS2-tS1<350e3&&qS1>3&&qS1<100&&qS2Top+qS2Bottom>100&&qS2<10000&&QualityCut","colz");

}
