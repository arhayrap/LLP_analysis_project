#ifndef DEF_WplusHToSS_STodd_ms15_pl1000
#define DEF_WplusHToSS_STodd_ms15_pl1000

#include "RazorAnalyzer.h"

class WplusHToSS_STodd_ms15_pl1000: public RazorAnalyzer {
    public: 
        WplusHToSS_STodd_ms15_pl1000(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
