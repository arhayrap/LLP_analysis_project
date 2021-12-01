#ifndef DEF_WplusHToSS_SToTauTau_ms15_pl1000
#define DEF_WplusHToSS_SToTauTau_ms15_pl1000

#include "RazorAnalyzer.h"

class WplusHToSS_SToTauTau_ms15_pl1000: public RazorAnalyzer {
    public: 
        WplusHToSS_SToTauTau_ms15_pl1000(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
