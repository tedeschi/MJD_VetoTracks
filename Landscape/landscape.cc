// landscape.cc
// topographical data
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "landscape.hh"
#include "TApplication.h"
#include <math.h> 
#include "TStyle.h"

using namespace std;

void landscape() {
	TFile *f = new TFile("../Data/Surf4850_2.85_out.root");
	TH2F *h = (TH2F*)f->Get("landscape2");
	h->SetDirectory(0);
	
	Int_t binx = h->GetNbinsX();
	Int_t biny = h->GetNbinsY();
	Double_t min = h->GetMinimum();
	Double_t max = h->GetMaximum();
	Double_t xminbin = h->FindFirstBinAbove(0.0 , 1); //1=x-axis, 2=y-axis
	Double_t xmaxbin = h->FindLastBinAbove(0.0 , 1);
	Double_t yminbin = h->FindFirstBinAbove(0.0, 2);
	Double_t ymaxbin = h->FindLastBinAbove(0.0 , 2);
	Int_t binnum = h->FindBin(0,0,0); //location of experiment
	Double_t depthAt0 = h->GetBinContent(binnum); 
	Double_t theta = 0.0;
	Double_t phi = 0.0;
	Double_t r = 0.0;
	Double_t surfDist = 0.0;
	Double_t PI = 3.14159265;
	Double_t depthAtxy = 0;
	Double_t slantDepth = 0;
	Double_t binWidthx = 0;
	Double_t binWidthy = 0;
	Int_t MJDbinx,MJDbiny,MJDbinz; //location of experiment bin numbers
	Double_t xminRange = h->GetXaxis()->GetXmin();
	Double_t xmaxRange = h->GetXaxis()->GetXmax();
	Double_t yminRange = h->GetYaxis()->GetXmin();
	Double_t ymaxRange = h->GetYaxis()->GetXmax();
	TH2F *thetaPhiHist = new TH2F("thetaPhiHist","Theta Phi",binx,0,95,biny,0,365);
	TH2F *slantHist = new TH2F("slantHist","Slant depth",binx,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),biny,h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
	TH2F *cosThetaPhiHist;
	
	cout << "number of bins x: " << binx << endl;
	cout << "number of bins y: " << biny << endl;
	cout << "min heigth: " << min << endl;
	cout << "max heigth: " << max << endl;
	cout << "depth at (0,0): " << depthAt0 << endl;
	cout << "Bin number at (0,0): " << binnum << endl;
	cout << "X min bin: " << xminbin << endl;
	cout << "X max bin: " << xmaxbin << endl;
	cout << "Y min bin: " << yminbin << endl;
	cout << "Y max bin: " << ymaxbin << endl;
	cout << "X min range: " << xminRange << endl;
	cout << "X max range: " << xmaxRange << endl;
	cout << "Y min range: " << yminRange << endl;
	cout << "Y max range: " << ymaxRange << endl;
	
	h->GetBinXYZ(binnum,MJDbinx,MJDbiny,MJDbinz);
	cout << "MJDbinx: " << MJDbinx << " MJDbiny: " << MJDbiny << " MJDbinz: " << MJDbinz << endl;
	binWidthx = h->GetXaxis()->GetBinWidth(0);
	binWidthy = h->GetYaxis()->GetBinWidth(0);
	cout << "bin widthX: " << binWidthx << endl;
	cout << "bin widthY: " << binWidthy << endl;
	
	/*
	//some manual testing
	double xcoord = xmaxbin; //2702
	double ycoord = ymaxbin; //1624
	surfDist = sqrt((((xcoord-MJDbinx)*binWidthx)*((xcoord-MJDbinx)*binWidthx))+(((ycoord-MJDbiny)*binWidthy)*((ycoord-MJDbiny)*binWidthy)));
	cout << "surfdist: " << surfDist << endl;
	double k = h->GetBinContent(xcoord,ycoord);
	cout << "depth at coords: " << k << endl;
	r = sqrt((((xcoord-MJDbinx)*binWidthx)*((xcoord-MJDbinx)*binWidthx)) + (((ycoord-MJDbiny)*binWidthy)*((ycoord-MJDbiny)*binWidthy)) + ((k - (depthAt0 - 1478))*(k - (depthAt0 - 1478))));
	cout << "r: " << r << endl;
	cout << "theta: " << asin(surfDist/r) * (180/PI) << endl;
	phi = atan((ycoord-1624)/(xcoord-2702)) * (180/PI);
	if(xcoord >= MJDbinx && ycoord >= MJDbiny) { //quadrant 1
		phi = 90-phi;
	}
	if(xcoord >= MJDbinx && ycoord < MJDbiny) { //quadrant 2
		phi = abs(phi)+90;
	}
	if(xcoord < MJDbinx && ycoord < MJDbiny) { //quadrant 3
		phi = (90-phi)+180;
	}
	if(xcoord < MJDbinx && ycoord >= MJDbiny) { //quadrant 4
		phi = abs(phi)+270;
	}
	cout << "phi: " << phi << endl;
	
	//cout << h->GetBinContent(h->FindBin(h->GetXaxis()->GetXmax()-67,h->GetYaxis()->GetXmax()-140));
	
	int i,o,p;
	h->GetBinXYZ(h->FindBin(-10000,-50000),i,o,p);
	cout << "i: " << i << " o: " << o << endl;
	*/
	
	getchar();
	for(double i = 0; i < binx; i++) {
		for(double j = 0; j < biny; j++) { //(x,y,z) = (i,j,k)
			double k = h->GetBinContent(i,j);
			binnum = h->GetBin(i,j); //global bin number
			
			//we subtract from MJDbin to get coordinates relative to MJD, and multiply by binWidth to get units in meters rather than bins
			phi = atan((j-MJDbiny)/(i-MJDbinx)) * (180/PI);
			surfDist = sqrt((((i-MJDbinx)*binWidthx)*((i-MJDbinx)*binWidthx))+(((j-MJDbiny)*binWidthy)*((j-MJDbiny)*binWidthy)));
			r = sqrt((((i-MJDbinx)*binWidthx)*((i-MJDbinx)*binWidthx)) + (((j-MJDbiny)*binWidthy)*((j-MJDbiny)*binWidthy)) + ((k - (depthAt0 - 1478))*(k - (depthAt0 - 1478))));
			theta = asin(surfDist/r) * (180/PI);
			
			if(i >= MJDbinx && j >= MJDbiny) { //quadrant 1
				phi = 90-phi;
			}
			if(i >= MJDbinx && j < MJDbiny) { //quadrant 2
				phi = abs(phi)+90;
			}
			if(i < MJDbinx && j < MJDbiny) { //quadrant 3
				phi = (90-phi)+180;
			}
			if(i < MJDbinx && j >= MJDbiny) { //quadrant 4
				phi = abs(phi)+270;
			}
			thetaPhiHist->Fill(theta,phi,r); 
			slantHist->SetBinContent(binnum,r);
		}
		cout << "i: " << i << endl;
	}
	
	TApplication *myapp=new TApplication("myapp",0,0);
	TCanvas *c1 = new TCanvas("c1","c1",1200,1200);
	TCanvas *c2 = new TCanvas("c2","c2",1200,1200);
	TCanvas *axisProj = new TCanvas("axisProj","X and Y axis projections",1200,1200);
	
	//some coloring
	const Int_t rgbn = 5;
    const Int_t contn = 999;
    Double_t stops[rgbn] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[rgbn]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[rgbn] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[rgbn]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(rgbn, stops, red, green, blue, contn);
    gStyle->SetNumberContours(contn);
	
		c1->cd();
	c1->SetLogz();
	thetaPhiHist->SetXTitle("Theta");
	thetaPhiHist->SetYTitle("Phi");
	thetaPhiHist->SetStats(0);
	thetaPhiHist->GetYaxis()->SetTitleOffset(1.3);
	thetaPhiHist->Draw("colz");
	
		c2->cd();
	//c2->SetLogz();
	slantHist->SetStats(0);
	slantHist->SetXTitle("Meters");
	slantHist->SetYTitle("Meters");
	slantHist->Draw("colz");
	
	axisProj->Divide(2,2);
		axisProj->cd(1);
	TH1D *proj1 = thetaPhiHist->ProjectionX();
	proj1->SetTitle("Theta Phi X Projection");
	proj1->Draw();
		axisProj->cd(2);
	TH1D *proj2 = thetaPhiHist->ProjectionY();
	proj2->SetTitle("Theta Phi Y Projection");
	proj2->Draw();
		axisProj->cd(3);
	TH1D *proj3 = slantHist->ProjectionX();
	proj3->SetTitle("Slant Depth X Projection");
	proj3->Draw();
		axisProj->cd(4);
	TH1D *proj4 = slantHist->ProjectionY();
	proj4->SetTitle("Slant Depth Y Projection");
	proj4->Draw();
	
	//print out the canvases
	cout << "Making pictures..." << endl;
	char thetaPhiHistcan[150];
    sprintf(thetaPhiHistcan,"thetaPhiTop.png");
    c1->Print(thetaPhiHistcan,"png");
    char slantcan[150];
    sprintf(slantcan,"slantDepthTop.png");
    c2->Print(slantcan,"png");
    cout << "Done! Use ctrl+c to exit." << endl;
	
	myapp->Run();

	f->Close();
}

int main(int argc, char* argv[]) {
	landscape();
	return 0;
}
