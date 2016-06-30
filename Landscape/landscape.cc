// landscape.cc
// topographical data
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "landscape.hh"

using namespace std;

void landscape() {
	TFile *f = new TFile("../Data/Surf4850_2.85_out.root");
	TH2F *h = (TH2F*)f->Get("landscape2");
	h->SetDirectory(0);
	Int_t binx = h->GetNbinsX();
	Int_t biny = h->GetNbinsY();
	Double_t min = h->GetMinimum();
	Double_t max = h->GetMaximum();
	Int_t binnum = h->FindBin(0,0);
	Double_t heigthAt0 = h->GetBinContent(binnum);
	
	cout << "number of bins x: " << binx << endl;
	cout << "number of bins y: " << biny << endl;
	cout << "min heigth: " << min << endl;
	cout << "max heigth: " << max << endl;
	cout << "heigth at (0,0): " << heigthAt0 << endl;
	
	//TCanvas *can = new TCanvas("can","landscape",0,0,800,800);
//	can->cd();
//	h->Draw();
//	can->Update();
	
//	char landcan[150];
//		sprintf(landcan,"landcan.pdf");
//		can->Print(landcan,"pdf");
//	
	f->Close();
	delete f;
}

int main(int argc, char* argv[]) {
	landscape();
	return 0;
}
