#ifndef slantBin_hh
#define slantBin_hh
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <list>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TMath.h>
#include "TBenchmark.h"
#include <TColor.h>
#include <TStyle.h>

#include "TProfile.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoManager.h"
#include "TPaveText.h"
#include "TDatabasePDG.h"
#include "TVirtualGeoTrack.h"
#include "TVector3.h"
#include "TApplication.h"
#include "TGeoCompositeShape.h"
#include "TGraphErrors.h"
#include "TRandom.h"

class slantBin
{
	public:
		slantBin();
		int main();
}

#endif
;
