// Wenbo 
// A formal one that read a txt file.
#ifndef SETRUNLIST_H
#define SETRUNLIST_H

// SetRunList("/store/bl2/pandax/zhouxp/PandaX-II/leakBDT/Run11Rn_run.lst", runlistvec);
void SetRunList(TString str_path_file, vector<int> & runlistvec)
{
  // runlistvec.clear();
  ifstream path_file(str_path_file);
  if (path_file.good()) {
    int current_number = 0;
    while (path_file >> current_number){
      runlistvec.push_back(current_number);
    }
    path_file.close();
  } else {
    cout << "Error! No run list found";
    _exit(0);
  }
}

#endif 
