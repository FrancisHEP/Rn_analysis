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
    //General
    int runNo, trigNo, fileNo;
    double triggerTime, Time2PreviousTrig, TotQPreviousTrig, TotQThisTrig;
    unsigned int time, nS1, nS2, SingleS2;
    int s1max, s2max;
    //for S1s
    double xS1B[9999], yS1B[9999], xS1T[9999], yS1T[9999];
    double qS1[9999], qS1T[9999], qS1B[9999];
    double tS1[9999], ttenS1[9999], wS1[9999], wtenS1[9999], fwhmS1[9999], hS1[9999], pS1[9999];
    double qS1Veto[9999];
    unsigned int S1Density[9999], S1NPeaks[9999], nS1BSatur[9999], nS1TSatur[9999];
    double S1PS40[9999], S1PS60[9999], S1PS100[9999], S1PSS40[9999], S1PSS60[9999], S1PSS100[9999], S1Decay10[9999], S1Decay30[9999], S1Decay60[9999];
    //for S2s
    double xS2B[9999], yS2B[9999], xS2T[9999], yS2T[9999], xS2NN[9999], yS2NN[9999], chi2[9999];
    double qS2[9999], qS2T[9999], qS2B[9999];
    double tS2[9999], ttenS2[9999], wS2[9999], wtenS2[9999], fwhmS2[9999], hS2[9999], pS2[9999];
    unsigned int nS2BSatur[9999], nS2TSatur[9999];
  
    //For additional ones in output file ana2
    double xS2LRF[9999], yS2LRF[9999], LRFchi2[9999];
