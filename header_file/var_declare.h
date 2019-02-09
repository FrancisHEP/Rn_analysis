//Variable Declaration
  //For hitData
  int runNumber, eventNumber;
  int channelId;
  int nhits;
  const int maxhit = 9999;
  double prePedstal[maxhit];
  double postPedstal[maxhit];
  double rmsPreBaseline[maxhit];
  double rmsPostBaseline[maxhit];
  int startTime[maxhit];
  int peakTime[maxhit];
  double height[maxhit];
  int width[maxhit];
  int hitType[maxhit];
  int pmtId,pmtNo[20000];

  //For input file, anaData
  int runNo, trigNo, fileNo;
  double triggerTime;
  unsigned int nS1, nS2;
  int s1max, s2max;
  double qS1T[1000], qS1B[1000], xS1T[1000], yS1T[1000];
  double xS1B[1000], yS1B[1000];
  double qS1[1000], tS1[1000], pS1[1000], hS1[1000];
  double ttenS1[1000], wS1[1000], wtenS1[1000];
  double qS1Veto[1000];
  unsigned int S1NPeaks[1000], nS1BSatur[1000], nS1TSatur[1000];
  double qS2T[1000], qS2B[1000], xS2T[1000], yS2T[1000];
  double xS2B[1000], yS2B[1000];
  double qS2[1000], tS2[1000], pS2[1000], hS2[1000];
  double ttenS2[1000], wS2[1000], wtenS2[1000];
  double xS2NN[1000], yS2NN[1000], chi2[1000];
  unsigned int nS2BSatur[1000], nS2TSatur[1000];

  //For additional ones in output file
  int beta, alpha;
  double xS2LRF[1000], yS2LRF[1000], LRFchi2[1000];
