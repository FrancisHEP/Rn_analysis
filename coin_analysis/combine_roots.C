// combine BiPo212 root files 24467 and 24468
// Wenbo Ma
void combine_roots()
{
  vector<int> runNo = {24467, 24468};
  TChain * chain = new TChain("rn220_BA");
  for (int i=0; i<runNo.size(); ++i) {
    TString pathfile = Form("ba212/outTree_212_BiPo_%d.root",runNo[i]);
    chain->Add(pathfile);
  }
  TFile * file = new TFile("ba212/BiPo212_24467-24468.root","recreate");
  chain->Write();
  file->Close();
  
}
