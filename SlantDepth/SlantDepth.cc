// SlantDepth.cc
// topographical data
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "SlantDepth.hh"
#include "../VetoDisplay/VetoDisplay.hh"

using namespace std;

Double_t slantPhi = 0.0;
Double_t slantTheta = 0.0;
static Double_t PI = 3.14159265;
static Double_t EPSILON = 1.0;
static Double_t MJDDEPTH = 1478.0; //in meters
Double_t r = 1478.0; //initial guess
Double_t deltar = 10; //how much r is incremented by
Double_t diffOld = 0.0;
Double_t diffNew = 0.0;
int counter = 0;
int counter2 = 0;
//int i = 0;
Double_t thetaStep = 0;
//Double_t deltatheta = log(91-thetaStep)/4;
Double_t deltatheta = 1;

TFile *f = new TFile("../Data/Surf4850_2.85_out.root");
TH2F *h = (TH2F*)f->Get("landscape2");
	
Int_t binnum = h->FindBin(0,0,0); //location of experiment
Double_t heightAt0 = h->GetBinContent(binnum); 
Int_t nbinx = h->GetNbinsX();
Int_t nbiny = h->GetNbinsY();

// main function -------------------------------------------------------
void SlantDepth() {
	TList* hList = new TList();
	TH2F *slantHist = new TH2F("slantHist","Slant Depth",360,0,360, 90,0,90);
	TH2F *slantLowTheta = new TH2F("slantLowTheta","Slant Depth low Theta",720,0,360, 120,0,60);
	TH2F *slantHighTheta = new TH2F("slantHighTheta","Slant Depth High Theta",720,0,360, 60,60,90);
	TH2F *lookupHist = new TH2F("lookupHist","lookupHist",nbinx,0,360, nbiny,0,90);
	
	//filling 0-90 theta histo -----------------------------------------
	for(int i = 0; i < 360; i++) {
		while(thetaStep < 90) {
			slantPhi = (i) * (PI/180);
			slantTheta = thetaStep * (PI/180);
			r = 1478.0;
			deltar = 10;
			diffNew = 0.0;
			diffOld = 0.0;
			counter = 0;
			slantHist->Fill(i,thetaStep,shooter(r));
			//deltatheta = log(91-thetaStep)/4;
			deltatheta = 1;
			thetaStep += deltatheta;
			//cout << thetaStep << " ";
			counter2++;
			if(counter2 == 140) {
				break;
			}
			//thetaStep++;
		}
		cout << i << endl;
		thetaStep = 0;
		counter2 = 0;
	}

	//filling 0-60 theta histo -----------------------------------------
	for(double counteri1 = 0; counteri1 < 360; counteri1+=0.5) {
		for(double counterj1 = 0; counterj1 < 60; counterj1+=0.5) {
			slantPhi = (counteri1) * (PI/180);
			slantTheta = counterj1 * (PI/180);
			r = 1478.0;
			deltar = 10;
			diffNew = 0.0;
			diffOld = 0.0;
			counter = 0;
			slantLowTheta->Fill(counteri1,counterj1,shooter(r));
		}
		cout << counteri1 << endl;
	}
	
	
	
	//filling 60-90 theta histo ----------------------------------------
	for(double counteri2 = 0; counteri2 < 360; counteri2+=0.5) {
		for(double counterj2 = 60; counterj2 < 90; counterj2+=0.5) {
			slantPhi = (counteri2) * (PI/180);
			slantTheta = counterj2 * (PI/180);
			r = 1478.0;
			deltar = 10;
			diffNew = 0.0;
			diffOld = 0.0;
			counter = 0;
			slantHighTheta->Fill(counteri2,counterj2,shooter(r));
		}
		cout << counteri2 << endl;
	}
	//editing,writing, and drawing histos ------------------------------
	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	
	c1->cd();
		slantHist->SetStats(0);
		slantHist->SetXTitle("Phi");
		slantHist->SetYTitle("Theta");
		c1->SetLogz();
		slantHist->Draw("colz");
	
	c2->cd();
		slantLowTheta->SetStats(0);
		slantLowTheta->SetXTitle("Phi");
		slantLowTheta->SetYTitle("Theta");
		//c2->SetLogz();
		slantLowTheta->Draw("colz");
		
	c3->cd();
		slantHighTheta->SetStats(0);
		slantHighTheta->SetXTitle("Phi");
		slantHighTheta->SetYTitle("Theta");
		c3->SetLogz();
		slantHighTheta->Draw("colz");
		
	
	TFile *f2 = new TFile("output/slantHist.root","NEW");
	slantHist->Write();
	slantLowTheta->Write();
	slantHighTheta->Write();
	
	char c1print[150];
    sprintf(c1print,"output/slantHist.png");
    c1->Print(c1print,"png");
    
    char c2print[150];
    sprintf(c2print,"output/slantLowTheta.png");
    c2->Print(c2print,"png");
    
    char c3print[150];
    sprintf(c3print,"output/slantHighTheta.png");
    c3->Print(c3print,"png");
    
    f2->Close();
	f->Close();
}

//command line input ---------------------------------------------------
void SlantDepth(char* innum1, char* innum2) {
	bool validNum1 = false;
	bool validNum2 = false;
	
	stringstream myStream(innum1);
	if (myStream >> slantPhi) {
		slantPhi = atoi(innum1);
		validNum1 = true;
	}
	else {
		cout << innum1 << " is not a valid number. Please try again. " << endl;
	}
	
	stringstream myStream2(innum2);
	if (myStream2 >> slantTheta) {
		slantTheta = atoi(innum2);
		validNum2 = true;
	}
	else {
		cout << innum2 << " is not a valid number. Please try again. " << endl;
	}
	
	
	if(validNum1 == true && validNum2 == true) {
		slantPhi = (slantPhi) * (PI/180); //cos and sin need to be input in rads
		slantTheta = slantTheta * (PI/180);
		cout << "Slant Depth: " << shooter(r) << endl;
	}
}

//implementation for VetoDisplay
Double_t SlantDepth(Double_t phiin, Double_t thetain) {
	slantPhi = (phiin) * (PI/180);
	slantTheta = thetain * (PI/180);
	r = 1478.0;
	deltar = 10;
	diffNew = 0.0;
	diffOld = 0.0;
	counter = 0;
	return shooter(r);
}

//recursive method -----------------------------------------------------
Double_t shooter(Double_t r) {
	Double_t xAtrTip = cos(slantPhi) * sin(slantTheta) * r;
	Double_t yAtrTip = sin(slantPhi) * sin(slantTheta) * r;
	Double_t zAtrTip = (heightAt0 - MJDDEPTH) + r*cos(slantTheta);
	
	Int_t binnum = h->FindBin(xAtrTip,yAtrTip,0);
	Double_t zAtSurf = h->GetBinContent(binnum);
	diffNew = zAtrTip - zAtSurf;
	
	int binx,biny,binz;
	h->GetBinXYZ(binnum,binx,biny,binz);
	
	Double_t binWidthx = h->GetXaxis()->GetBinWidth(0);
	Double_t binWidthy = h->GetYaxis()->GetBinWidth(0);
	Double_t binCenterx = h->GetXaxis()->GetBinCenter(binx);
	Double_t binCentery = h->GetYaxis()->GetBinCenter(biny);
	Double_t distToEdgex = binWidthx/2 - abs(binCenterx - xAtrTip);
	Double_t distToEdgey = binWidthy/2 - abs(binCentery - yAtrTip);
	
	/*
	cout << "xAtrTip: " << xAtrTip << endl;
	cout << "yAtrTip: " << yAtrTip << endl;
	cout << "zAtrTip: " << zAtrTip << endl;
	cout << "zAtSurf: " << zAtSurf << endl;
	cout << "DiffNew: " << diffNew << endl;
	cout << "DiffOld: " << diffOld << endl;
	cout << "DeltaR: " << deltar << endl;
	cout << "Dist to edge x: " << distToEdgex << " Of bin: " << binx << endl;
	cout << "Dist to edge y: " << distToEdgey << " Of bin: " << biny << endl;
	cout << "Bin width x: " << binWidthx << endl;
	cout << "Bin width y: " << binWidthy << endl;
	cout << "Bin center x: " << binCenterx << endl;
	cout << "Bin center y: " << binCentery << endl;
	cout << "r: " << r << endl;
	*/
	
	//fixing errors near bin edges -------------------------------------
	if(abs(diffNew) < 10) { //when close to surface
		Double_t slope = 0.0;
		Int_t zNextBin = 0;
		Double_t line = 0;
		
		if( distToEdgex < (binWidthx /6) ) { //and close to bin edge in x
			if( binCenterx > xAtrTip ) { //if left edge
				zNextBin = h->GetBinContent(binx-1);
				slope = (zAtSurf - zNextBin) / (binCenterx - h->GetXaxis()->GetBinCenter(binx-1)); //y2-y1 / x2-x1
				//cout << "left edge" << endl;
				line = (slope / binWidthx) * (binCenterx - xAtrTip);
				diffNew = zAtrTip - (zAtSurf-line);
				//cout << "line: " << line << endl;
				//cout << "new diffNew: " << diffNew << endl;
			}
			if( binCenterx < xAtrTip ) { //if right edge
				zNextBin = h->GetBinContent(binx+1);
				slope = (zNextBin - zAtSurf) / (h->GetXaxis()->GetBinCenter(binx+1) - binCenterx); //y2-y1 / x2-x1
				//cout << "right edge" << endl;
				line = (slope / binWidthx) * (xAtrTip - binCenterx);
				diffNew = zAtrTip - (zAtSurf+line);
				//cout << "line: " << line << endl;
				//cout << "new diffNew: " << diffNew << endl;
			}
		}
		if( distToEdgey < (binWidthy /6) ) { //and close to bin edge in y
			if( binCentery > yAtrTip ) { //if lower edge
				zNextBin = h->GetBinContent(biny-1);
				slope = (zAtSurf - zNextBin) / (binCentery - h->GetYaxis()->GetBinCenter(biny-1)); //y2-y1 / x2-x1
				//cout << "lower edge" << endl;
				line = (slope / binWidthy) * (binCentery - yAtrTip);
				diffNew = zAtrTip - (zAtSurf-line);
				//cout << "line: " << line << endl;
				//cout << "new diffNew: " << diffNew << endl;
			}
			if( binCentery < yAtrTip ) { //if upper edge
				zNextBin = h->GetBinContent(biny+1);
				slope = (zNextBin - zAtSurf) / (h->GetYaxis()->GetBinCenter(biny+1) - binCentery); //y2-y1 / x2-x1
				//cout << "upper edge" << endl;
				line = (slope / binWidthy) * (yAtrTip - binCentery);
				diffNew = zAtrTip - (zAtSurf+line);
				//cout << "line: " << line << endl;
				//cout << "new diffNew: " << diffNew << endl;
			}
		}
		//cout << "slope: " << slope << endl;
	}
	
	//exit condition ---------------------------------------------------
	if(counter == 10) {
		//return r;
	}
	
	if(deltar < 1) { 
		//cout << "something may have went wrong.." << endl;
		//cout << "counter: " << counter << endl;
		//cout << "i: " << i << endl;
		//cout << "j: " << thetaStep << endl << endl;
		return r;
	}
	
	if(abs(diffNew) < EPSILON) { //exit condition
	//	cout << "times it took: " << counter << endl;
		return r;
	}
	// -----------------------------------------------------------------
	
	if (diffNew * diffOld < 0) { //if sign flipped on new difference
		deltar = deltar / 2;
	}
	if(diffNew > 0) {
		r = r - deltar;
	}
	if(diffNew < 0) {
		r = r + deltar;
	}
	//cout << endl;
	diffOld = diffNew;
	counter++;
	shooter(r);
}
// ---------------------------------------------------------------------
 //can run stand alone, but must comment out call at the end of VetoDisplay
int main(int argc, char* argv[]) {
	if(argc == 3) {
		SlantDepth(argv[1],argv[2]);
	}
	else {
		SlantDepth();
	}
	return 0;
}
