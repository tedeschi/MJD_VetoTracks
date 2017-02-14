//
// macro to plot hists from file
// 
// David J Tedeschi, University of South Carolina
//
// don't close file at end - else hists disappear (root directory thing)
//--------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>


void hplot(){

	// open the file
	TFile *f = new TFile("AccCalc.root");

 	//get hists
	TH2D *hTopXY[144];
	TH2D *hBotXY[144];
	char *histname = new char[20];
	for (int i=0; i<144; i++){
		sprintf(histname, "hTopXY_%d",i);
		hTopXY[i] = (TH2D*)f->Get(histname);
		sprintf(histname, "hBotXY_%d",i);
		hBotXY[i] = (TH2D*)f->Get(histname);
	}


	gStyle->SetOptStat(0);
	
	//-----------------------------------------------------------------------------------------
  	//create slant depth TCanvasses
	int aCounter = 33;
	/*
  	TCanvas *can1 = new TCanvas("can1","slant depth 1",0,0,800,800);
	can1->Divide(6,6,0,0);
	for (Int_t i=0;i<36;i++){
		can1->cd(1+i);
		hSlantDepth[i]->Draw();
	}
	
	for (Int_t i=1;i<=36;i++){
		can1->cd(i);
		hSlantDepth[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}

  	TCanvas *can2 = new TCanvas("can2","slant depth 2 ",0,50,800,800);
	can2->Divide(6,6,0,0);
	aCounter = 69;
	for (Int_t i=36;i<72;i++){
		can2->cd(1+i-36);
		hSlantDepth[i]->Draw();
	}
	
	for (Int_t i=37;i<=72;i++){
		can2->cd(i-36);
		hSlantDepth[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}

  	TCanvas *can3 = new TCanvas("can3","slant depth 3",0,100,800,800);
	can3->Divide(6,6,0,0);
	aCounter = 105;
	for (Int_t i=72;i<108;i++){
		can3->cd(1+i-72);
		hSlantDepth[i]->Draw();
	}
	
	for (Int_t i=73;i<=108;i++){
		can3->cd(i-72);
		hSlantDepth[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}

  	TCanvas *can4 = new TCanvas("can4","slant depth 4",0,150,800,800);
	can4->Divide(6,6,0,0);
	aCounter = 141;
	for (Int_t i=108;i<144;i++){
		can4->cd(1+i-108);
		hSlantDepth[i]->Draw();
	}
	
	for (Int_t i=109;i<=144;i++){
		can4->cd(i-108);
		hSlantDepth[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	*/
	//-----------------------------------------------------------------------------------------
 	//create XY coverage TCanvasses
  	TCanvas *Tcan1 = new TCanvas("Tcan1","XY 1",100,0,800,800);
	Tcan1->Divide(6,6,0,0);
	/*for (Int_t i=0;i<36;i++){
		Tcan1->cd(1+i);
		hTopXY[i]->SetFillColor(2);
		hTopXY[i]->Draw("box");
		hBotXY[i]->SetFillColor(3);
		hBotXY[i]->Draw("same");

	}
	*/
	aCounter = 33;
	for(Int_t i = 1; i <= 36; i++) {
		Tcan1->cd(i);
		hTopXY[aCounter]->SetFillColor(2);
		hTopXY[aCounter]->Draw("box");
		hBotXY[aCounter]->SetFillColor(3);
		hBotXY[aCounter]->Draw("same");
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	
  	TCanvas *Tcan2 = new TCanvas("Tcan2","XY 2 ",100,50,800,800);
	Tcan2->Divide(6,6,0,0);
	/*for (Int_t i=36;i<72;i++){
		Tcan2->cd(1+i-36);
		hTopXY[i]->SetFillColor(2);
		hTopXY[i]->Draw("box");
		hBotXY[i]->SetFillColor(3);
		hBotXY[i]->Draw("same");
	}
	*/
	aCounter = 69;
	for(Int_t i = 37; i <= 72; i++) {
		Tcan2->cd(i-36);
		hTopXY[aCounter]->SetFillColor(2);
		hTopXY[aCounter]->Draw("box");
		hBotXY[aCounter]->SetFillColor(3);
		hBotXY[aCounter]->Draw("same");
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	
  	TCanvas *Tcan3 = new TCanvas("Tcan3","XY 3",100,100,800,800);
	Tcan3->Divide(6,6,0,0);
	
	/*for (Int_t i=72;i<108;i++){
		Tcan3->cd(1+i-72);
		hTopXY[i]->SetFillColor(2);
		hTopXY[i]->Draw("box");
		hBotXY[i]->SetFillColor(3);
		hBotXY[i]->Draw("same");

	}
	*/
	aCounter = 105;
	for(Int_t i = 73; i <= 108; i++) {
		Tcan3->cd(i-72);
		hTopXY[aCounter]->SetFillColor(2);
		hTopXY[aCounter]->Draw("box");
		hBotXY[aCounter]->SetFillColor(3);
		hBotXY[aCounter]->Draw("same");
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	
  	TCanvas *Tcan4 = new TCanvas("Tcan4","XY 4",100,150,800,800);
	Tcan4->Divide(6,6,0,0);
	/*for (Int_t i=108;i<144;i++){
		Tcan4->cd(1+i-108);
		hTopXY[i]->SetFillColor(2);
		hTopXY[i]->Draw("box");
		hBotXY[i]->SetFillColor(3);
		hBotXY[i]->Draw("same");
	}
	*/
	aCounter = 141;
	for(Int_t i = 109; i <= 144; i++) {
		Tcan4->cd(i-108);
		hTopXY[aCounter]->SetFillColor(2);
		hTopXY[aCounter]->Draw("box");
		hBotXY[aCounter]->SetFillColor(3);
		hBotXY[aCounter]->Draw("same");
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	//-----------------------------------------------------------------------------------------
  	//create slant depth TCanvasses
  	/*
  	TCanvas *Acan1 = new TCanvas("Acan1","ThetaPhi 1",200,0,800,800);
	Acan1->Divide(6,6,0,0);
	aCounter = 33;
	for (Int_t i=0;i<36;i++){
		Acan1->cd(1+i);
		hThetaPhi[i]->Draw();
	}
	
	for(Int_t i = 1; i <= 36; i++) {
		Acan1->cd(i);
		hThetaPhi[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	
  	TCanvas *Acan2 = new TCanvas("Acan2","ThetaPhi 2 ",200,50,800,800);
	Acan2->Divide(6,6,0,0);
	aCounter = 69;
	for (Int_t i=36;i<72;i++){
		Acan2->cd(1+i-36);
		hThetaPhi[i]->Draw();
	}
	
	for(Int_t i = 37; i <= 72; i++) {
		Acan2->cd(i-36);
		hThetaPhi[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	
  	TCanvas *Acan3 = new TCanvas("Acan3","ThetaPhi 3",200,100,800,800);
	Acan3->Divide(6,6,0,0);
	aCounter = 105;
	for (Int_t i=72;i<108;i++){
		Acan3->cd(1+i-72);
		hThetaPhi[i]->Draw();
	}
	
	for(Int_t i = 73; i <= 108; i++) {
		Acan3->cd(i-72);
		hThetaPhi[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	
  	TCanvas *Acan4 = new TCanvas("Acan4","ThetaPhi 4",200,150,800,800);
	Acan4->Divide(6,6,0,0);
	aCounter = 141;
	for (Int_t i=108;i<144;i++){
		Acan4->cd(1+i-108);
		hThetaPhi[i]->Draw();
	}
	
	for(Int_t i = 109; i <= 144; i++) {
		Acan4->cd(i-108);
		hThetaPhi[aCounter]->Draw();
		if(i%3==0) {
			aCounter = aCounter-6;
		}
		aCounter++;
	}
	*/
}
