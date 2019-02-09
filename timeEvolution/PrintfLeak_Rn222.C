#include <fstream>
#include <iostream>
#include <vector>
#include <map>
using namespace std;
// const char * inputfilename1 = "evo_lists/2018_inj_1st.daqinfo"; // daq info
// const char * inputfilename2 = "APeaks_2018_inj_1st.lst"; // a peaks info 
void PrintfLeak_Rn222(const char * inputfilename1, const char * inputfilename2, int runNo_low, int runNo_up)
{
  int sum = 0;
  double duration_sum = 0;

  ifstream fin1(inputfilename1);
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
  ifstream fin2(inputfilename2);
  int fin2line = 0;
  int runno2[1000];
  double Rn222[1000], Po218[1000], Rn220[1000], Po216[1000], Po212[1000];
  double Rn222err[1000], Po218err[1000], Rn220err[1000], Po216err[1000], Po212err[1000];
  if (fin2.good()) {
    int i = 0;
    double t1,t2,t3,t4,t5,t6;
    while(fin2>>t1>>t2>>t3>>t4>>t5>>t6)
      {
	runno2[i]=t1;
	Rn222err[i]=sqrt(t2)/duration[i];
	Po218err[i]=sqrt(t3)/duration[i];
	Rn220err[i]=sqrt(t4)/duration[i];
	Po216err[i]=sqrt(t5)/duration[i];
	Po212err[i]=sqrt(t6)/duration[i];
	Rn222[i]=t2;
	Po218[i]=t3;
	Rn220[i]=t4;
	Po216[i]=t5;
	Po212[i]=t6;
	if (runno2[i]>=runNo_low&&runno2[i]<=runNo_up) sum+=Rn222[i];
	i++;
      }
    fin2.close();
    fin2line = i;
  } else {
    printf("Error! fin2 not found\n");
    _exit(0);
  }
  if (fin1line!=fin2line) {
    cout << "line num doesn't match!" << endl;
    _exit(0);
  }
  printf("sum: %d[evts]\n",sum);
  printf("duration_sum: %d[h]\n",duration_sum/3600);

  
}
