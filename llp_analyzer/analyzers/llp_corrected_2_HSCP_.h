#ifndef DEF_llp_corrected_2_HSCP_
#define DEF_llp_corrected_2_HSCP_

#include "RazorAnalyzer.h"

class llp_corrected_2_HSCP_: public RazorAnalyzer {
    public: 
        llp_corrected_2_HSCP_(TTree *tree=0): RazorAnalyzer(tree) { }
        void Analyze(bool isData, int option, string outputFileName, string label);
};

#endif
