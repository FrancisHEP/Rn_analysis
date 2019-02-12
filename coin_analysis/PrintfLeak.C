// printf count and duration info
// at a given runNo interval
// const char * daqfile = "2017_daqinfo.txt"
// const char * inputfilename2 = "2017_lowEnergy_rate.txt"
void PrintfLeak(TString daqfile, TString ratefile, int runNo_low, int runNo_up)
{

  int sum = 0;
  double duration_sum = 0;
  // int runNo_low = 0; 
  // int runNo_up  = 21778;

  ifstream fin1(daqfile);
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
	time[i] += duration[i]; if (duration[i]==0) printf("%d\n",duration[i]);
	
	if (runno[i]>=runNo_low&&runno[i]<=runNo_up) duration_sum += duration[i];
	
	++i;
      }
    fin1.close();
    fin1line = i;
  } else {
    printf("Error! fin1 not found\n");
    _exit(0);
  }

  ifstream path_file(ratefile);
  if (path_file.good()) {
    int runNo;
    int count;
    while (path_file >> runNo >> count){
      if (runNo>=runNo_low&&runNo<=runNo_up) sum+= count;
    }
    path_file.close();
  } else {
    cout << "Error! No run list found";
    _exit(0);
  }
  
  printf("sum: %d[evts]\n",sum);
  printf("duration_sum: %d[h]\n",duration_sum/3600);
}
