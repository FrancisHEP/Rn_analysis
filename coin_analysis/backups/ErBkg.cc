#include <TTree.h>

#include "stats-common.hh"

ValueMap ErBkg(require_ana3_t, TTree& out_tree)
{
  double counts = out_tree.GetEntries("S1Coincidence&&S1NoiseFilter&&S1Pattern&&S2Asy&&TimeDiff>0.01&&qS1>3&&qS1<45&&qS2Top+qS2Bottom>100&&qS2<10000&&qS1V==0&&nPostCutS1<=2&&nS2==1&&log10(qS2/qS1)>1.074-0.595*exp(-qS1/6.38)&&tS2-tS1>20e3&&tS2-tS1<340e3&&xTopNN*xTopNN+yTopNN*yTopNN<60000");

  return {
    { "low_energy_er_counts", counts }
  };
}
