// SlantDepth.cc
// topographical data
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "SlantDepth.hh"

using namespace std;

Double_t phi = 0.0;
Double_t theta = 0.0;
static Double_t PI = 3.14159265;
static Double_t EPSILON = 1.0;
static Double_t MJDDEPTH = 1478; //in meters
Double_t r = 1478; //initial guess
Double_t deltar = 739; //how much r is incremented by
Double_t diffOld = 0.0;
Double_t diffNew = 0.0;
int counter = 0;

TFile *f = new TFile("../Data/Surf4850_2.85_out.root");
TH2F *h = (TH2F*)f->Get("landscape2");
	
Int_t binnum = h->FindBin(0,0,0); //location of experiment
Double_t heightAt0 = h->GetBinContent(binnum); 

//user input -----------------------------------------------------------
void SlantDepth() {
	cout << "Enter your phi: ";
	cin >> phi;
	cout << "Enter your theta: ";
	cin >> theta;
	
	phi = phi * (PI/180); //cos and sin need to be input in rads
	theta = theta * (PI/180);
	cout << "Slant Depth: " << shooter(r) << endl;
}

//command line input ---------------------------------------------------
void SlantDepth(char* innum1, char* innum2) {
	bool validNum1 = false;
	bool validNum2 = false;
	
	stringstream myStream(innum1);
	if (myStream >> phi) {
		phi = atoi(innum1);
		validNum1 = true;
	}
	else {
		cout << innum1 << " is not a valid number. Please try again. " << endl;
	}
	
	stringstream myStream2(innum2);
	if (myStream2 >> theta) {
		theta = atoi(innum2);
		validNum2 = true;
	}
	else {
		cout << innum2 << " is not a valid number. Please try again. " << endl;
	}
	
	
	if(validNum1 == true && validNum2 == true) {
		phi = phi * (PI/180); //cos and sin need to be input in rads
		theta = theta * (PI/180);
		cout << "Slant Depth: " << shooter(r) << endl;
	}
}

//recursive method -----------------------------------------------------
int shooter(int r) {
	Double_t xAtrTip = sin(phi) * sin(theta) * r;
	Double_t yAtrTip = cos(phi) * sin(theta) * r;
	Double_t zAtrTip = (heightAt0 - MJDDEPTH) + sqrt( pow(r,2) - pow(r*sin(theta),2) );
	
	Int_t binnum = h->FindBin(xAtrTip,yAtrTip,0);
	Double_t zAtSurf = h->GetBinContent(binnum);
	diffNew = zAtrTip - zAtSurf;
	
	cout << "xAtrTip: " << xAtrTip << endl;
	cout << "yAtrTip: " << yAtrTip << endl;
	cout << "zAtrTip: " << zAtrTip << endl;
	cout << "zAtSurf: " << zAtSurf << endl;
	cout << "DiffNew: " << diffNew << endl;
	cout << "DiffOld: " << diffOld << endl;
	cout << "DeltaR: " << deltar << endl << endl;
	cout << "r: " << r << endl;
	
	if(deltar < 1) { //for high thetas that go off the chart, look into better way
		cout << "something may have went wrong.." << endl;
		cout << "counter: " << counter << endl;
		return r;
	}
	
	if(abs(diffNew) < EPSILON) { //exit condition
		cout << "times it took: " << counter << endl;
		return r;
	}
	
	if (diffNew * diffOld < 0) { //if sign flipped on new difference
		deltar = deltar / 2;
	}
	if(diffNew > 0) {
		r = r - deltar;
	}
	if(diffNew < 0) {
		r = r + deltar;
	}
	
	diffOld = diffNew;
	counter++;
	shooter(r);
}
// ---------------------------------------------------------------------

int main(int argc, char* argv[]) {
	if(argc == 3) {
		SlantDepth(argv[1],argv[2]);
	}
	else {
		SlantDepth();
	}
	return 0;
}
