{
  int RunNo = 24327;
  TChain *chain;
  chain = new TChain("signal_tree");
  chain->Add(Form("/store/px/data/udm/by-run/%d/new_correction/ana*",RunNo));  
  chain->Draw("qS1");
}


