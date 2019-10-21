#include "SlantDepthGrapher.hh"
using namespace std;

int main() {
	
	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	TCanvas *c2 = new TCanvas("c2","c2",1200,800);
	TCanvas *c3 = new TCanvas("c3","c3",1200,800);
	
	TFile *f = new TFile("../SlantDepth/output/slantHist.root");
	TH2F *slantHist = (TH2F*)f->Get("slantHist");
	TH2F *slantLowTheta = (TH2F*)f->Get("slantLowTheta");
	TH2F *slantHighTheta = (TH2F*)f->Get("slantHighTheta");
	
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
		
	char c1print[150];
    sprintf(c1print,"../SlantDepth/output/slantHist.png");
    c1->Print(c1print,"png");
    
    char c2print[150];
    sprintf(c2print,"../SlantDepth/output/slantLowTheta.png");
    c2->Print(c2print,"png");
    
    char c3print[150];
    sprintf(c3print,"../SlantDepth/output/slantHighTheta.png");
    c3->Print(c3print,"png");
    
    f->Close();
    return 0;
}
