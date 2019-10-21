// grapher.cc
// QDC and angle error checking
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "grapher.hh"

void grapher() {
	Double_t means[6];
	Double_t devs[6];
	
	TFile *f = new TFile("../VetoDisplay/output/changeMyName.root");
	
	TH1D *h2 = (TH1D*)f->Get("totalqdcangle2");
	TH1D *h3 = (TH1D*)f->Get("totalqdcangle3");
	TH1D *h4 = (TH1D*)f->Get("totalqdcangle4");
	TH1D *h5 = (TH1D*)f->Get("totalqdcangle5");
	TH1D *h6 = (TH1D*)f->Get("totalqdcangle6");
	TH1D *h7 = (TH1D*)f->Get("totalqdcangle7");
	
	//y axis -----------------------------------------------------------
	means[0] = h2->GetMean(1); 
	means[1] = h3->GetMean(1);
	means[2] = h4->GetMean(1);
	means[3] = h5->GetMean(1);
	means[4] = h6->GetMean(1);
	means[5] = h7->GetMean(1);
	
	/*
	Int_t binmax = h2->GetMaximumBin();
	Double_t x = h2->GetXaxis()->GetBinCenter(binmax);
	means[0] = x;
	
	binmax = h3->GetMaximumBin();
	x = h3->GetXaxis()->GetBinCenter(binmax);
	means[0] = x;
	
	binmax = h4->GetMaximumBin();
	x = h4->GetXaxis()->GetBinCenter(binmax);
	means[0] = x;

	binmax = h5->GetMaximumBin();
	x = h5->GetXaxis()->GetBinCenter(binmax);
	means[0] = x;
	
	binmax = h6->GetMaximumBin();
	x = h6->GetXaxis()->GetBinCenter(binmax);
	means[0] = x;
	
	binmax = h7->GetMaximumBin();
	x = h7->GetXaxis()->GetBinCenter(binmax);
	means[0] = x;
	*/
	
	//Ey ---------------------------------------------------------------
	devs[0] = h2->GetStdDev(1); 
	devs[1] = h3->GetStdDev(1);
	devs[2] = h4->GetStdDev(1);
	devs[3] = h5->GetStdDev(1);
	devs[4] = h6->GetStdDev(1);
	devs[5] = h7->GetStdDev(1);
	
	//manual x axis and Ex ---------------------------------------------
	Double_t xVals[6] = {4.5, 11.0, 15.0, 20.0, 27.5, 38.5}; //center values for each graph
	
	Double_t exVals[6] = {16.0, 16.4, 17.1, 16.0, 16.0, 16.0};

	TCanvas *c1 = new TCanvas("c1","A Simple Graph with error bars",200,10,700,500);
	TCanvas *c2 = new TCanvas("c2","TotalQDCinThetaRange",1200, 800);
	//c1->SetGrid();
	   
	TGraphErrors *gr = new TGraphErrors(6,xVals,means,exVals,devs);
	TF1 *myfit2 = new TF1("myfit2","[0]*2/cos(x*(pi/180))",0,45);
	
	c1->cd();
	myfit2->SetParameter(0,0);
	gr->GetXaxis()->SetRangeUser(0,60);
	gr->Fit("myfit2");
	gr->SetTitle("QDC and angle variablility;Degrees theta;Total QDC top panels");
	gr->GetYaxis()->SetTitleOffset(1.4);
	gr->Draw("AP*");
	
	c2->Divide(3,2);
	c2->cd(1);
	h2->Draw("bar");
	c2->cd(2);
	h3->Draw("bar");
	c2->cd(3);
	h4->Draw("bar");
	c2->cd(4);
	h5->Draw("bar");
	c2->cd(5);
	h6->Draw("bar");
	c2->cd(6);
	h7->Draw("bar");
	
	char errCan[150];
		sprintf(errCan,"./errors.pdf");
		c1->Print(errCan,"pdf");
		
	char angleCan[150];
		sprintf(angleCan,"./angles.pdf");
		c2->Print(angleCan,"pdf");
		
	f->Close();
}

int main() {
	grapher();
	return 0;
}
