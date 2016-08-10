#define SlantDepth_hh
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

#include "TStyle.h"
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

Double_t shooter(Double_t r);
void SlantDepth();
void SlantDepth(char* innum1, char* innum2);
Double_t SlantDepth(Double_t phiin, Double_t thetain);
