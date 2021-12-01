#ifndef DEF_GMSB_STau_mass247_ctau3000
#define DEF_GMSB_STau_mass247_ctau3000

#include "RazorAnalyzer.h"

class GMSB_STau_mass247_ctau3000: public RazorAnalyzer {
    public: 
        GMSB_STau_mass247_ctau3000(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
