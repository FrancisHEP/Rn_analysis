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
#include "CalRate_lowEnergy.C"
using namespace std;

//
void PrintList_lowEnergy(TString runlist)
{
  vector<int> runvec; SetRunList(runlist, runvec);
  for (int i=0;i<runvec.size();++i) {
    int rate = CalRate_lowEnergy(runvec[i]);
  }
  printf("PrintList_lowEnergy finished. \n");
}

