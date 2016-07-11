// landscape.cc
// topographical data
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "landscape.hh"

using namespace std;

void landscape() {
	TFile *f = new TFile("../Data/Surf4850_2.85_out.root");
	TH2F *h = (TH2F*)f->Get("landscape2");
	TList* hList = new TList();
	
	static Double_t PI = 3.14159265;
	static Double_t MJDDEPTH = 1478; //in meters
	Int_t nbinx = h->GetNbinsX();
	Int_t nbiny = h->GetNbinsY();
	Double_t minHeight = h->GetMinimum();
	Double_t maxHeight = h->GetMaximum();
	Int_t binnum = h->FindBin(0,0,0); //location of experiment
	Double_t heightAt0 = h->GetBinContent(binnum); 
	Double_t theta = 0.0;
	Double_t phi = 0.0;
	Double_t r = 0.0;
	Double_t surfDist = 0.0;
	Double_t slantDepth = 0;
	Double_t binWidthx = h->GetXaxis()->GetBinWidth(0);
	Double_t binWidthy = h->GetYaxis()->GetBinWidth(0);
	Int_t MJDbinx,MJDbiny,MJDbinz; //location of experiment bin numbers
	Double_t xminRange = h->GetXaxis()->GetXmin();
	Double_t xmaxRange = h->GetXaxis()->GetXmax();
	Double_t yminRange = h->GetYaxis()->GetXmin();
	Double_t ymaxRange = h->GetYaxis()->GetXmax();
	//TH2F *thetaPhiHist = new TH2F("thetaPhiHist","Calculated polar angles for topology",180,0,360,nbinx,0,93);
	TH2F *thetaPhiHist = new TH2F("thetaPhiHist","Calculated polar angles for topology",nbinx,0,93,nbiny,0,360);
	TH2F *cosThetaPhiHist = new TH2F("cosThetaPhiHist","Calculated polar angles for topology",nbinx,0,1,nbiny,0,365);
	TH2F *slantHist = new TH2F("slantHist","Calculated Slant depth for topology",nbinx,h->GetXaxis()->GetXmin(),h->GetXaxis()->GetXmax(),nbiny,h->GetYaxis()->GetXmin(),h->GetYaxis()->GetXmax());
	
	cout << "number of bins x: " << nbinx << endl;
	cout << "number of bins y: " << nbiny << endl;
	cout << "min heigth: " << minHeight << endl;
	cout << "max heigth: " << maxHeight << endl;
	cout << "depth at (0,0): " << heightAt0 << endl;
	cout << "Bin number at (0,0): " << binnum << endl;
	cout << "X min range: " << xminRange << endl;
	cout << "X max range: " << xmaxRange << endl;
	cout << "Y min range: " << yminRange << endl;
	cout << "Y max range: " << ymaxRange << endl;
	cout << "bin widthX: " << binWidthx << endl;
	cout << "bin widthY: " << binWidthy << endl;
	
	h->GetBinXYZ(binnum,MJDbinx,MJDbiny,MJDbinz);
	cout << "MJDbinx: " << MJDbinx << " MJDbiny: " << MJDbiny << " MJDbinz: " << MJDbinz << endl;
	
	cout << "MJD bin center x: " << h->GetXaxis()->GetBinCenter(MJDbinx) << endl;
	cout << "MJD bin center y: " << h->GetYaxis()->GetBinCenter(MJDbiny) << endl;
	
	//main loops -------------------------------------------------------
	cout << "Press enter to draw graphs" << endl;
	getchar();
	
	for(Int_t i = 0; i < nbinx; i++) {
		for(Int_t j = 0; j < nbiny; j++) {
			Int_t xBinAtI,yBinAtI,zBinAtI;
			
			binnum = h->GetBin(i,j);
			h->GetBinXYZ(binnum,xBinAtI,yBinAtI,zBinAtI); //global bin number
			
			Double_t xAtI = h->GetXaxis()->GetBinCenter(xBinAtI);
			Double_t yAtI = h->GetYaxis()->GetBinCenter(yBinAtI);
			Double_t heightAtI = h->GetBinContent(i,j);
			
			phi = atan(yAtI/xAtI) * (180/PI);
			surfDist = sqrt( pow(xAtI,2) + pow(yAtI,2) );
			r = sqrt( pow(xAtI,2) + pow(yAtI,2) + pow(heightAtI-(heightAt0-MJDDEPTH),2) );
			theta = asin(surfDist/r);
			
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
			
			thetaPhiHist->Fill(theta*(180/PI),phi,r);
			//thetaPhiHist->Fill(phi,theta*(180/PI),r);
			theta = (cos(theta));
			cosThetaPhiHist->Fill(theta,phi,r);
			slantHist->SetBinContent(binnum,r);
		}
		cout << "i: " << i << " / " << nbinx-1 << endl;
	}
	
	//canvas declarations ----------------------------------------------
	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	TCanvas *axisProj = new TCanvas("axisProj","X and Y axis projections",1200,800);
	TCanvas *top = new TCanvas("top","topology",1200,800);
	
	//some coloring ----------------------------------------------------
	const Int_t rgbn = 5;
    const Int_t contn = 999;
    Double_t stops[rgbn] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[rgbn]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[rgbn] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[rgbn]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(rgbn, stops, red, green, blue, contn);
    gStyle->SetNumberContours(contn);
	
	//edit and draw graphs ---------------------------------------------
		c1->cd();
	c1->SetLogz();
	thetaPhiHist->SetStats(0);
	thetaPhiHist->SetXTitle("Theta");
	thetaPhiHist->SetYTitle("Phi");
	thetaPhiHist->GetYaxis()->SetTitleOffset(1.3);
	thetaPhiHist->Draw("colz");
	//thetaPhiHist->Draw("polsurf");
	
		c2->cd();
	//c2->SetLogz();
	slantHist->SetStats(0);
	slantHist->SetXTitle("Meters");
	slantHist->SetYTitle("Meters");
	slantHist->GetYaxis()->SetTitleOffset(1.6);
	slantHist->GetXaxis()->SetNoExponent(true);
	slantHist->GetYaxis()->SetNoExponent(true);
	slantHist->Draw("colz");
	
		c3->cd();
	c3->SetLogz();
	cosThetaPhiHist->SetStats(0);
	cosThetaPhiHist->SetXTitle("Cos(Theta)");
	cosThetaPhiHist->SetYTitle("Phi");
	cosThetaPhiHist->GetYaxis()->SetTitleOffset(1.3);
	cosThetaPhiHist->Draw("colz");
	
	axisProj->Divide(2,2);
		axisProj->cd(1);
	TH1D *proj1 = thetaPhiHist->ProjectionX();
	proj1->SetTitle("Polar angles X Projection");
	proj1->SetYTitle("Total Slant Depth");
	proj1->GetYaxis()->SetTitleOffset(1.3);
	proj1->Draw();
		axisProj->cd(2);
	TH1D *proj2 = thetaPhiHist->ProjectionY();
	proj2->SetTitle("Polar angles Y Projection");
	proj2->SetYTitle("Total Slant Depth");
	proj2->GetYaxis()->SetTitleOffset(1.3);
	proj2->Draw();
		axisProj->cd(3);
	TH1D *proj3 = slantHist->ProjectionX();
	proj3->SetTitle("Slant Depth X Projection");
	proj3->SetYTitle("Total Slant Depth");
	proj3->GetYaxis()->SetTitleOffset(1.3);
	proj3->Draw();
		axisProj->cd(4);
	TH1D *proj4 = slantHist->ProjectionY();
	proj4->SetTitle("Slant Depth Y Projection");
	proj4->SetYTitle("Total Slant Depth");
	proj4->GetYaxis()->SetTitleOffset(1.3);
	proj4->Draw();
	
		top->cd();
	h->SetStats(0);
	h->SetXTitle("Meters");
	h->SetYTitle("Meters");
	h->GetYaxis()->SetTitleOffset(1.3);
	h->Draw("colz");

	//print out the canvases -------------------------------------------
	cout << "Making pictures..." << endl;
	char thetaPhiHistcan[150];
    sprintf(thetaPhiHistcan,"output/polarAnglesTop.png");
    c1->Print(thetaPhiHistcan,"png");
    
    char slantcan[150];
    sprintf(slantcan,"output/slantDepthTop.png");
    c2->Print(slantcan,"png");
    
    char cosThetaPhiHistcan[150];
    sprintf(cosThetaPhiHistcan,"output/polarAnglesTop2.png");
    c3->Print(cosThetaPhiHistcan,"png");
    
    char axisProjcan[150];
    sprintf(axisProjcan,"output/axisProjections.png");
    axisProj->Print(axisProjcan,"png");
    
    char topcan[150];
    sprintf(topcan,"output/Topology.png");
    top->Print(topcan,"png");
    
	// write objects to file -------------------------------------------
	cout << "Writing objects to root file..." << endl;
	TFile *f2 = new TFile("output/TopologyCalc.root","NEW");
	thetaPhiHist->Write();
	cosThetaPhiHist->Write();
	slantHist->Write();
	proj1->Write();
	proj2->Write();
	proj3->Write();
	proj4->Write();
	h->Write();
	
	f2->Close();
	f->Close();
	cout << "Done!" << endl;
}

int main(int argc, char* argv[]) {
	landscape();
	return 0;
}
