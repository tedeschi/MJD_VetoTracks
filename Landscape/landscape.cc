// landscape.cc
// topographical data
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "landscape.hh"
#include "TApplication.h"
#include <math.h> 

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
	Int_t MJDbinx,MJDbiny,MJDbinz; //location of experiment bin numbers
	Double_t xminRange = h->GetXaxis()->GetXmin();
	Double_t xmaxRange = h->GetXaxis()->GetXmax();
	Double_t yminRange = h->GetYaxis()->GetXmin();
	Double_t ymaxRange = h->GetYaxis()->GetXmax();
	TH2F *thetaPhiHist = new TH2F("thetaPhiHist","Theta Phi",binx,0,90,biny,0,365);
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
	cout << "bin width: " << h->GetXaxis()->GetBinWidth(0) << endl;
	
	
	//some manual testing
	double xcoord = MJDbinx-5; //2702
	double ycoord = MJDbiny+5; //1624
	surfDist = sqrt(((xcoord-MJDbinx)*(xcoord-MJDbinx))+((ycoord-MJDbiny)*(ycoord-MJDbiny)));
	cout << "surfdist: " << surfDist << endl;
	double k = h->GetBinContent(xcoord,ycoord);
	cout << "depth at coords: " << k << endl;
	r = sqrt(((xcoord-MJDbinx)*(xcoord-MJDbinx)) + ((ycoord-MJDbiny)*(ycoord-MJDbiny)) + ((k - (depthAt0 - 1478))*(k - (depthAt0 - 1478))));
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
		phi = phi+180;
	}
	if(xcoord < MJDbinx && ycoord >= MJDbiny) { //quadrant 4
		phi = abs(phi)+270;
	}
	cout << "phi: " << phi << endl;
	
	
	getchar();
	for(double i = xminbin; i < xmaxbin; i++) {
		for(double j = yminbin; j < ymaxbin; j++) { //(x,y,z) = (i,j,k)
			double k = h->GetBinContent(i,j);
			binnum = h->GetBin(i,j,k); //global bin number
			
			phi = atan((j-MJDbiny)/(i-MJDbinx)) * (180/PI);
			surfDist = sqrt(((i-MJDbinx)*(i-MJDbinx))+((j-MJDbiny)*(j-MJDbiny)));
			r = sqrt(((i-MJDbinx)*(i-MJDbinx)) + ((j-MJDbiny)*(j-MJDbiny)) + ((k - (depthAt0 - 1478))*(k - (depthAt0 - 1478))));
			theta = asin(surfDist/r) * (180/PI);
			
			if(i >= MJDbinx && j >= MJDbiny) { //quadrant 1
				phi = 90-phi;
			}
			if(i >= MJDbinx && j < MJDbiny) { //quadrant 2
				phi = abs(phi)+90;
			}
			if(i < MJDbinx && j < MJDbiny) { //quadrant 3
				phi = phi+180;
			}
			if(i < MJDbinx && j >= MJDbiny) { //quadrant 4
				phi = abs(phi)+270;
			}
			
			thetaPhiHist->Fill(theta,phi,r); 
			slantHist->SetBinContent(binnum,r);
		}
		cout << "i: " << i << endl;
	}
	
	binnum = thetaPhiHist->FindBin(0.282895,315);
	cout << thetaPhiHist->GetBinContent(binnum) << endl;
	
	
	TApplication *myapp=new TApplication("myapp",0,0);
	TCanvas *c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(2,1);
	c1->cd(1);
	thetaPhiHist->SetXTitle("Theta");
	thetaPhiHist->SetYTitle("Phi");
	thetaPhiHist->Draw("colz");
	c1->cd(2);
	slantHist->Draw("colz");
	myapp->Run();

	f->Close();
}

int main(int argc, char* argv[]) {
	landscape();
	return 0;
}
