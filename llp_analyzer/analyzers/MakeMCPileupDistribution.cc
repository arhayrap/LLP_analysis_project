#include "MakeMCPileupDistribution.h"
#include "JetCorrectorParameters.h"

//C++ includes
#include <sys/stat.h>

//ROOT includes
#include "TH1F.h"

using namespace std;


void NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b)<hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }
}


void MakeMCPileupDistribution::Analyze(bool isData, int option, string outputFileName, string label)
{
  //initialization: create one TTree for each analysis box
  std::cout << "Initializing..." << std::endl;
  if (outputFileName.empty()){
    cout << "MakeMCPileupDistribution: Output filename not specified!" << endl << "Using default output name MCPileupDistribution.root" << endl;
    outputFileName = "MCPileupDistribution.root";
  }
  //---------------------------

  TFile outFile(outputFileName.c_str(), "RECREATE");

  string Label = label;
  if (label != "") Label = "_"+label;
  TH1F* histPUMean =  new TH1F( ("PUMean"+Label).c_str(),";nPUMean;Number of Events", 200, -0.5, 199.5);

  // TH1F* histPUMean =  new TH1F( ("PUMean"+Label).c_str(),";nPUMean;Number of Events", 200, 0, 200);

  //begin loop
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;


  //cout << "Total Events: " << fChain->GetEntries() << "\n";
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //begin event
    if(jentry % 1000000 == 0) cout << "Processing entry " << jentry << endl;
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    //find in-time pileup
    float intime_PUmean = 0;
    for (int i=0; i<nBunchXing; ++i) {
      if (BunchXing[i] == 0) {
	       intime_PUmean = nPUmean[i];
      }
    }

    //fill normalization histogram
    histPUMean->Fill(intime_PUmean);

  }//end of event loop

  //Normalize the histograms
  NormalizeHist(histPUMean);

  cout << "Writing output ..." << endl;
  outFile.WriteTObject(histPUMean, ("PUMean"+Label).c_str(), "WriteDelete");
  outFile.Close();
}
