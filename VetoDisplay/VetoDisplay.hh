#ifndef VetoDisplay_hh
#define VetoDisplay_hh
#include "VetoDisplay.hh"
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

class VetoDisplay
{
	public:
		VetoDisplay();
		int coloring(int qdc,int mode);
		int isNextTo(int panel1, int panel2);
		int isLayerHit(int panel1, int panel2);
		std::vector<Double_t> hitLocation(int panel1, int panel2);
		void hitFinder(char* innum1, char* innum2);
		void DrawEvent(Int_t qdcVals[], Int_t numberOfPanelsHit, Int_t totalQDC, Int_t runNumber, Int_t eventCount);
		void fillPlots(Int_t qdcvals[], Int_t totalQDC, Int_t numberOfPanelsHit, Int_t ievent);
		void drawPlots();
		void printPlots();
}

#endif
;
