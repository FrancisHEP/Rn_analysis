/************************************************* 
Copyright: PandaX-II Collaboration 
Author: Wenbo Ma
Date: 2018-12-12
Description: 
Print a low-energy rate list to text file.
**************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include "SetRunList.h"
#include "CalRate_ba212.C"
using namespace std;

//
void PrintList_ba212(TString runlist)
{
  vector<int> runvec; SetRunList(runlist, runvec);
  for (int i=0;i<runvec.size();++i) {
    int rate = CalRate_ba212(runvec[i]);
  }
  printf("PrintList_ba212 finished. \n");
}

