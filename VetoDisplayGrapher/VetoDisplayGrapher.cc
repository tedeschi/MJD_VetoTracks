#include "VetoDisplayGrapher.hh"
using namespace std;

int main() {

	TFile *f = new TFile("../VetoDisplay/vetoPlots.root");

	TCanvas *graphs = new TCanvas("graphs","Analysis graphs", 900, 0, 900, 600);
	TCanvas *QDCcanvas = new TCanvas("QDC","QDC", 900, 0, 900, 600);
	TCanvas *projcan = new TCanvas("Projections","Panels hit vs. Multiplicity", 900, 0, 700, 700);
	TCanvas *qdctotalcan = new TCanvas("qdctotalcan", "QDC Totals", 900, 0, 1401, 400);
	TCanvas *qdctotaln = new TCanvas("qdctotaln", "QDC totals for each n", 1000, 500);
	TCanvas *thetaPhican = new TCanvas("thetaPhican","Theta and Phi of particle's tracks", 900, 0, 900, 600);
	TCanvas *costhetaPhican = new TCanvas("costhetaPhican","Cos(Theta) and Phi of particle's tracks", 900, 0, 900, 600);
	TCanvas *QDCanglecan = new TCanvas("QDCanglecan"," ", 1200, 800);
	TCanvas *QDCanglecan2 = new TCanvas("QDCanglecan2"," ", 1200, 800);
	TCanvas *QDCanglecan3 = new TCanvas("QDCanglecan3"," ", 1200, 800);
	TCanvas *QDCslantcan = new TCanvas("QDCslantcan"," ", 1200, 800);
	TCanvas *inThrucan = new TCanvas("inThrucan"," ", 1200, 800);
	TCanvas *totalQDCanglecan = new TCanvas("totalQDCanglecan"," ", 1200, 800);
	TCanvas *totalQDCanglecan2 = new TCanvas("totalQDCanglecan2"," ", 1200, 800);
	TCanvas *totalQDCanglecan3 = new TCanvas("totalQDCanglecan3"," ", 1200, 800);
	TCanvas *totalQDCslantcan = new TCanvas("totalQDCslantcan"," ", 1600, 700);
	TCanvas *angleOverview = new TCanvas("angleOverview", " ", 1200, 800);

	TH1F *graph1 = (TH1F*)f->Get("graph1");
	TH1F *graph2 = (TH1F*)f->Get("graph2");
	TH2F *graph3 = (TH2F*)f->Get("graph3");
	TH2F *graph4 = (TH2F*)f->Get("graph4");
	TH2F *graph5 = (TH2F*)f->Get("graph5");
	TH2F *graph6 = (TH2F*)f->Get("graph6");
	TH2F *graph7 = (TH2F*)f->Get("graph7");
	TH1F *graph8 = (TH1F*)f->Get("graph8");
	TH1F *graphn2 = (TH1F*)f->Get("graphn2");
	TH1F *graphn3 = (TH1F*)f->Get("graphn3");
	TH1F *graphn4 = (TH1F*)f->Get("graphn4");
	TH1F *graphn5 = (TH1F*)f->Get("graphn5");
	TH1F *graphn6 = (TH1F*)f->Get("graphn6");
	TH1F *graphn7 = (TH1F*)f->Get("graphn7");
	TH1F *graphn8 = (TH1F*)f->Get("graphn8");
	TH1F *graphn9 = (TH1F*)f->Get("graphn9");
	TH2F *thetaPhi = (TH2F*)f->Get("thetaPhi");
	TH2F *costhetaPhi = (TH2F*)f->Get("costhetaPhi");
	TH2D *QDCangle = (TH2D*)f->Get("QDCangle");
	TH2D *QDCangleTwo = (TH2D*)f->Get("QDCangleTwo");
	TH2D *QDCangleThree = (TH2D*)f->Get("QDCangleThree");
	TH2D *QDCangleFour = (TH2D*)f->Get("QDCangleFour");
	TH2D *QDCangleFive = (TH2D*)f->Get("QDCangleFive");
	TH2D *QDCangleSix = (TH2D*)f->Get("QDCangleSix");
	TH2F *QDCslant = (TH2F*)f->Get("QDCslant");
	TH2F *thetaSlant = (TH2F*)f->Get("thetaSlant");
	TH2F *inThruHist = (TH2F*)f->Get("inThruHist");
	TH2D *totalQDCangle = (TH2D*)f->Get("totalQDCangle");
	TH2D *totalQDCangleTwo = (TH2D*)f->Get("totalQDCangleTwo");
	TH2D *totalQDCangleThree = (TH2D*)f->Get("totalQDCangleThree");
	TH2D *totalQDCangleFour = (TH2D*)f->Get("totalQDCangleFour");
	TH2D *totalQDCangleFive = (TH2D*)f->Get("totalQDCangleFive");
	TH2D *totalQDCangleSix = (TH2D*)f->Get("totalQDCangleSix");
	TH2D *totalQDCangleSeven = (TH2D*)f->Get("totalQDCangleSeven");
	//TH1D *totalQDCangleprojy2 = totalQDCangleTwo->ProjectionY("totalqdcangle2",1,180);
	//TH1D *totalQDCangleprojy3 = totalQDCangleThree->ProjectionY("totalqdcangle3",1,180);
	//TH1D *totalQDCangleprojy4 = totalQDCangleFour->ProjectionY("totalqdcangle4",1,180);
	//TH1D *totalQDCangleprojy5 = totalQDCangleFive->ProjectionY("totalqdcangle5",1,180);
	//TH1D *totalQDCangleprojy6 = totalQDCangleSix->ProjectionY("totalqdcangle6",1,180);
	//TH1D *totalQDCangleprojy7 = totalQDCangleSeven->ProjectionY("totalqdcangle7",1,180);
	TH2F *totalQDCslant = (TH2F*)f->Get("totalQDCslant");
	TH1F *allSlantsHist = (TH1F*)f->Get("allSlantsHist");
	TH2F *northWestTop = (TH2F*)f->Get("northWestTop");
	TH2F *northEastTop = (TH2F*)f->Get("northEastTop");
	TH2F *southWestTop = (TH2F*)f->Get("southWestTop");
	TH2F *southEastTop = (TH2F*)f->Get("southEastTop");
	
	//some coloring for "colz" option ----------------------------------
	const Int_t rgbn = 5;
    const Int_t contn = 999;
    Double_t stops[rgbn] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[rgbn]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[rgbn] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[rgbn]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(rgbn, stops, red, green, blue, contn);
    gStyle->SetNumberContours(contn);
    
	//Drawing graphs canvas
	graphs->Divide(3,2);
	
	graphs->cd(1);
		graph1->SetXTitle("Multiplicity");
		graph1->SetYTitle("Count");
		graph1->GetYaxis()->SetTitleOffset(1.45);
		graph1->SetFillColor(1);
		graph1->Draw("bar");
	
	graphs->cd(2);
		graph2->SetXTitle("Panel number");
		graph2->SetYTitle("Count");
		graph2->GetYaxis()->SetTitleOffset(1.45);
		graph2->SetFillColor(1);
		graph2->Draw("bar");
	
	graphs->cd(3);
		graph4->SetXTitle("Multiplicity");
		graph4->SetYTitle("Total QDC");
		graph4->SetMarkerStyle(kFullDotSmall);
		graph4->GetYaxis()->SetTitleOffset(1.6);
		graph4->SetStats(0);
		graph4->Draw("colz");
	
	graphs->cd(4);
		graph3->SetYTitle("Panel number");
		graph3->SetXTitle("QDC");
		graph3->SetMarkerStyle(kFullDotSmall);
		graph3->SetStats(0);
		graph3->Draw("colz");
    
    graphs->cd(5);
		graph5->SetYTitle("Panel number");
		graph5->SetXTitle("QDC");
		graph5->SetMarkerStyle(kFullDotSmall);
		graph5->SetStats(0);
		graph5->Draw("colz");
    
    graphs->cd(6);
		graph6->SetYTitle("Panel number");
		graph6->SetXTitle("QDC");
		graph6->SetMarkerStyle(kFullDotSmall);
		graph6->SetStats(0);
		graph6->Draw("colz");
	
	//Drawing qdccanvas canvas
	QDCcanvas->Divide(3,2);
	
	QDCcanvas->cd(1);
		graph3->Draw("colz");
	
	TH1D *qdcpanelproj = graph3->ProjectionX();
	QDCcanvas->cd(4);
		qdcpanelproj->SetTitle("QDC distribution for all hits");
		qdcpanelproj->SetYTitle("Count");
		qdcpanelproj->GetYaxis()->SetTitleOffset(1.45);
		qdcpanelproj->Draw("bar");	
	
	QDCcanvas->cd(2);
		graph5->Draw("colz");
	
	TH1D *qdcpanelprojn2 = graph5->ProjectionX();
	QDCcanvas->cd(5);
		qdcpanelprojn2->SetTitle("QDC distribution for 2 panels hit");
		qdcpanelprojn2->SetYTitle("Count");
		qdcpanelprojn2->GetYaxis()->SetTitleOffset(1.45);
		qdcpanelprojn2->Draw("bar");
	
	QDCcanvas->cd(3);
		graph6->Draw("colz");
	
	TH1D *qdcpanelprojn4 = graph6->ProjectionX();
	QDCcanvas->cd(6);
		qdcpanelprojn4->SetTitle("QDC distribution for 4 panels hit");
		qdcpanelprojn4->SetYTitle("Count");
		qdcpanelprojn4->GetYaxis()->SetTitleOffset(1.45);
		qdcpanelprojn4->Draw("bar");
	
	//Drawing projcan canvas
	projcan->Divide(2,2);
	
	projcan->cd(1);
		graph7->SetXTitle("Multiplicity");
		graph7->SetYTitle("Panel number");
		graph7->SetMarkerStyle(kFullDotSmall);
		graph7->SetStats(0);
		graph7->Draw("colz");
	
	TH1D *countingPanelProj = graph7->ProjectionY();
	projcan->cd(2);
		countingPanelProj->SetTitle("Times each panel hit");
		countingPanelProj->SetYTitle("Count");
		countingPanelProj->Draw("hbar");
	
	TH1D *countingMultProj = graph7->ProjectionX();
	projcan->cd(3);
		countingMultProj->SetTitle("Multiplicity");
		countingMultProj->SetYTitle("Count");
		countingMultProj->GetYaxis()->SetTitleOffset(1.6);
		countingMultProj->Draw("bar");
	
	projcan->cd(4);
		graph8->GetYaxis()->SetTitleOffset(1.6);
		graph8->Draw("bar");
	
	//Drawing qdctotalcan canvas
	qdctotalcan->Divide(3,1);
	
	qdctotalcan->cd(1);
		graph8->SetXTitle("QDC total");
		graph8->SetYTitle("Count");
		graph8->Draw("bar");
	
	qdctotalcan->cd(2);
		graphn2->Draw("bar");
	
	qdctotalcan->cd(3);
		graphn4->Draw("bar");
	
	//Drawing qdctotaln canvas
	qdctotaln->Divide(4,2);
	int n2max = graphn2->GetMaximum(); //set each max to the first graph's max
	
	qdctotaln->cd(1);
		graphn2->SetXTitle("Total QDC");
		graphn2->SetYTitle("Count");
		graphn2->GetYaxis()->SetTitleOffset(1.45);
		graphn2->Draw("bar");
    
    qdctotaln->cd(2);
		graphn3->SetMaximum(n2max);
		graphn3->SetXTitle("Total QDC");
		graphn3->SetYTitle("Count");
		graphn3->GetYaxis()->SetTitleOffset(1.45);
		graphn3->Draw("bar");
    
    qdctotaln->cd(3);
		graphn4->SetMaximum(n2max);
		graphn4->SetXTitle("Total QDC");
		graphn4->SetYTitle("Count");
		graphn4->GetYaxis()->SetTitleOffset(1.45);
		graphn4->Draw("bar");
    
    qdctotaln->cd(4);
		graphn5->SetMaximum(n2max);
		graphn5->SetXTitle("Total QDC");
		graphn5->SetYTitle("Count");
		graphn5->GetYaxis()->SetTitleOffset(1.45);
		graphn5->Draw("bar");
    
    qdctotaln->cd(5);
		graphn6->SetMaximum(n2max);
		graphn6->SetXTitle("Total QDC");
		graphn6->SetYTitle("Count");
		graphn6->GetYaxis()->SetTitleOffset(1.45);
		graphn6->Draw("bar");
    
    qdctotaln->cd(6);
		graphn7->SetMaximum(n2max);
		graphn7->SetXTitle("Total QDC");
		graphn7->SetYTitle("Count");
		graphn7->GetYaxis()->SetTitleOffset(1.45);
		graphn7->Draw("bar");
    
    qdctotaln->cd(7);
		graphn8->SetMaximum(n2max);
		graphn8->SetXTitle("Total QDC");
		graphn8->SetYTitle("Count");
		graphn8->GetYaxis()->SetTitleOffset(1.45);
		graphn8->Draw("bar");
    
    qdctotaln->cd(8);
		graphn9->SetMaximum(n2max);
		graphn9->SetXTitle("Total QDC");
		graphn9->SetYTitle("Count");
		graphn9->GetYaxis()->SetTitleOffset(1.45);
		graphn9->Draw("bar");

	//drawing thetaPhican
	thetaPhican->Divide(2,2);
	
	thetaPhican->cd(1);
		thetaPhi->SetStats(0);
		thetaPhi->SetXTitle("Degrees(Phi)");
		thetaPhi->SetYTitle("Degrees(Theta)");
		thetaPhi->Draw("colz");
	
	TH1D *thetaPhiprojy = thetaPhi->ProjectionY();
	thetaPhican->cd(2);
		thetaPhiprojy->SetTitle("Theta");
		thetaPhiprojy->SetXTitle("Degrees");
		thetaPhiprojy->SetYTitle("Count");
		thetaPhiprojy->Draw("hbar");
	
	TH1D *thetaPhiprojx = thetaPhi->ProjectionX();
	thetaPhican->cd(3);
		thetaPhiprojx->SetTitle("Phi");
		thetaPhiprojx->SetXTitle("Degrees");
		thetaPhiprojx->SetYTitle("Count");
		thetaPhiprojx->Draw("bar");
	
	//drawing costhetaPhican
	costhetaPhican->Divide(2,2);
	
	costhetaPhican->cd(1);
		costhetaPhi->SetStats(0);
		costhetaPhi->SetXTitle("Degrees(Phi)");
		costhetaPhi->SetYTitle("Cos(Theta)");
		costhetaPhi->Draw("colz");
	
	TH1D *costhetaPhiprojy = costhetaPhi->ProjectionY();
	costhetaPhican->cd(2);
		costhetaPhiprojy->SetTitle("Cos(Theta)");
		costhetaPhiprojy->SetXTitle("cos(theta)");
		costhetaPhiprojy->SetYTitle("Count");
		costhetaPhiprojy->Draw("hbar");
	
	TH1D *costhetaPhiprojx = costhetaPhi->ProjectionX();
	costhetaPhican->cd(3);
		costhetaPhiprojx->SetTitle("Phi");
		costhetaPhiprojx->SetXTitle("Degrees");
		costhetaPhiprojx->SetYTitle("Count");
		costhetaPhiprojx->Draw("bar");
	
	//QDC angle canvas
	QDCanglecan->Divide(2,2);
	
	TF1 *myfit = new TF1("myfit","[0]*(1/cos(x*(pi/180)))",0,45);
	QDCanglecan->cd(1);
		QDCangle->SetXTitle("Degrees Theta");
		QDCangle->SetYTitle("QDC per panel");
		QDCangle->GetYaxis()->SetTitleOffset(1.3);
		myfit->SetParameter(0,1000);
		QDCangle->Fit("myfit");
		QDCangle->Draw("colz");
	
	QDCanglecan->cd(2);
		TH1D *QDCangleprojy = QDCangle->ProjectionY();
		QDCangleprojy->SetTitle("Projection Y");
		QDCangleprojy->SetYTitle("Count");
		QDCangleprojy->Fit("landau");
		QDCangleprojy->Draw("bar");
	
	QDCanglecan->cd(3);
		TH1D *QDCangleprojx = QDCangle->ProjectionX();
		QDCangleprojx->SetTitle("Projection X");
		QDCangleprojx->SetYTitle("Count");
		QDCangleprojx->Draw("bar");
	
	TProfile *QDCangleprofilex = QDCangle->ProfileX();
	QDCanglecan->cd(4);
		QDCangleprofilex->SetTitle("Profile X");
		QDCangleprofilex->SetYTitle("QDC per panel");
		QDCangleprofilex->GetYaxis()->SetTitleOffset(1.3);
		QDCangleprofilex->GetYaxis()->SetRangeUser(0.,5000.);
		QDCangleprofilex->Fit("myfit");

	//qdcanglecan2
	QDCanglecan2->Divide(3,2);
	
	QDCanglecan2->cd(1);
		QDCangle->Draw("colz");
	
	QDCanglecan2->cd(4);
		QDCangleprojy->Draw("bar");
	
	QDCanglecan2->cd(2);
		QDCangleTwo->Draw("colz");
	
	TH1D *QDCangleprojyTwo = QDCangleTwo->ProjectionY("qdcangle2",1,180);
	QDCanglecan2->cd(5);
		QDCangleprojyTwo->SetTitle("Projection Y");
		QDCangleprojyTwo->SetYTitle("Count");
		QDCangleprojyTwo->Draw("bar");
	
	QDCanglecan2->cd(3);
		QDCangleThree->Draw("colz");

	TH1D *QDCangleprojyThree = QDCangleThree->ProjectionY("qdcangle3",1,180);
	QDCanglecan2->cd(6);
		QDCangleprojyThree->SetTitle("Projection Y");
		QDCangleprojyThree->SetYTitle("Count");
		QDCangleprojyThree->Draw("bar");
	
	//qdcanglecan3
	QDCanglecan3->Divide(3,2);
	
	QDCanglecan3->cd(1);
		QDCangleFour->Draw("colz");
	
	TH1D *QDCangleprojyFour = QDCangleFour->ProjectionY("qdcangle4",1,180);
	QDCanglecan3->cd(4);
		QDCangleprojyFour->SetTitle("Projection Y");
		QDCangleprojyFour->SetYTitle("Count");
		QDCangleprojyFour->Draw("bar");
	
	QDCanglecan3->cd(2);
		QDCangleFive->Draw("colz");
	
	TH1D *QDCangleprojyFive = QDCangleFive->ProjectionY("qdcangle5",1,180);
	QDCanglecan3->cd(5);
		QDCangleprojyFive->SetTitle("Projection Y");
		QDCangleprojyFive->SetYTitle("Count");
		QDCangleprojyFive->Draw("bar");
	
	QDCanglecan3->cd(3);
		QDCangleSix->Draw("colz");
	
	TH1D *QDCangleprojySix = QDCangleSix->ProjectionY("qdcangle6",1,180);
	QDCanglecan3->cd(6);
		QDCangleprojySix->SetTitle("Projection Y");
		QDCangleprojySix->SetYTitle("Count");
		QDCangleprojySix->Draw("bar");
	
	//drawing QDCslantcan
	QDCslantcan->Divide(1,1);
	
	QDCslantcan->cd(1);
		allSlantsHist->SetXTitle("Slant Depth(meters)");
		allSlantsHist->SetYTitle("Count");
		allSlantsHist->Draw("bar");
	/*
	QDCslantcan->cd(1);
		QDCslant->SetYTitle("QDC per panel");
		QDCslant->GetYaxis()->SetTitleOffset(1.45);
		QDCslant->SetXTitle("Slant Depth(meters)");
		QDCslant->Draw("colz");
	
	QDCslantcan->cd(2);
		thetaSlant->SetYTitle("Degrees theta");
		thetaSlant->SetXTitle("Slant Depth(meters)");
		thetaSlant->Draw("colz");
	
	TProfile *QDCslantcanprofilex = QDCslant->ProfileX();
	QDCslantcan->cd(3);
		QDCslantcanprofilex->SetTitle("Profile X");
		QDCslantcanprofilex->SetYTitle("QDC per panel");
		QDCslantcanprofilex->GetYaxis()->SetTitleOffset(1.3);
		QDCslantcanprofilex->GetYaxis()->SetRangeUser(0.,5000.);
		QDCslantcanprofilex->Fit("pol1");
	*/
	//drawing inThrucan
	inThrucan->Divide(2,2);
	
	inThrucan->cd(1);
		inThruHist->SetXTitle("Inches");
		inThruHist->SetYTitle("QDC per panel");
		inThruHist->Draw("colz");
	
	TProfile *inThruHistprofilex = inThruHist->ProfileX();
	inThrucan->cd(2);
		//gStyle->SetOptStat(0);
		//gStyle->SetOptFit(1);
		inThruHistprofilex->SetTitle("Profile X");
		inThruHistprofilex->SetYTitle("QDC per panel");
		inThruHistprofilex->GetYaxis()->SetTitleOffset(1.3);
		inThruHistprofilex->GetYaxis()->SetRangeUser(0.,5000.);
		inThruHistprofilex->Fit("pol1");
	
	TH1D *inThruprojy = inThruHist->ProjectionY("inThruHistprojy",1,180);
	inThrucan->cd(3);
		inThruprojy->SetTitle("Projection Y");
		inThruprojy->SetYTitle("Count");
		inThruprojy->Fit("landau");
		inThruprojy->Draw("bar");
	
	//totalQDCanglecan
	totalQDCanglecan->Divide(2,2);
	/*
	TF1 *myfit2 = new TF1("myfit2","[0]*2/cos(x*(pi/180))",0,90);
	totalQDCanglecan->cd(1);
		totalQDCangle->SetXTitle("Degrees theta");
		totalQDCangle->SetYTitle("Total QDC");
		totalQDCangle->SetStats(0);
	*/
	TF1 *myfit2 = new TF1("myfit2","[0]*2/cos(x*(pi/180))",0,45);
	totalQDCanglecan->cd(1);
		totalQDCangle->SetXTitle("Degrees theta");
		totalQDCangle->SetYTitle("Total QDC");
		totalQDCangle->GetYaxis()->SetTitleOffset(1.3);
		myfit2->SetParameter(0,2000);
		totalQDCangle->Fit("myfit2");
		totalQDCangle->Draw("colz");
	/*
	TProfile *totalQDCangleprofilex = totalQDCangle->ProfileX();
	totalQDCanglecan->cd(2);
		//gStyle->SetOptStat(0);
		//gStyle->SetOptFit(1);
		totalQDCangleprofilex->SetTitle("Profile X");
		totalQDCangleprofilex->SetYTitle("total QDC");
		totalQDCangleprofilex->GetYaxis()->SetTitleOffset(1.3);
		totalQDCangleprofilex->GetYaxis()->SetRangeUser(0.,12000.);
		totalQDCangleprofilex->Fit("pol1");
	*/
	TH1D *totalQDCangleprojy = totalQDCangle->ProjectionY("totalqdcangle",1,180);
	totalQDCanglecan->cd(3);
		totalQDCangleprojy->SetTitle("Projection Y");
		totalQDCangleprojy->SetYTitle("Count");
		totalQDCangleprojy->Fit("landau");
		totalQDCangleprojy->Draw("bar");
	
	//totalQDCanglecan2
	TFile *f2 = new TFile("../VetoDisplay/output/totalQDCangles.root","RECREATE");
	
	totalQDCanglecan2->Divide(3,2);
	
	totalQDCanglecan2->cd(1);
		totalQDCangleTwo->Draw("colz");
	
	TH1D *totalQDCangleprojy2 = totalQDCangleTwo->ProjectionY("totalqdcangle2",1,180);
	totalQDCanglecan2->cd(4);
		totalQDCangleprojy2->SetTitle("Projection Y");
		totalQDCangleprojy2->SetYTitle("Count");
		totalQDCangleprojy2->Draw("bar");
		totalQDCangleprojy2->Write();
		
	
	totalQDCanglecan2->cd(2);
		totalQDCangleThree->Draw("colz");
	
	TH1D *totalQDCangleprojy3 = totalQDCangleThree->ProjectionY("totalqdcangle3",1,180);
	totalQDCanglecan2->cd(5);
		totalQDCangleprojy3->SetTitle("Projection Y");
		totalQDCangleprojy3->SetYTitle("Count");
		totalQDCangleprojy3->Draw("bar");
		totalQDCangleprojy3->Write();

	totalQDCanglecan2->cd(3);
		totalQDCangleFour->Draw("colz");
	
	TH1D *totalQDCangleprojy4 = totalQDCangleFour->ProjectionY("totalqdcangle4",1,180);
	totalQDCanglecan2->cd(6);
		totalQDCangleprojy4->SetTitle("Projection Y");
		totalQDCangleprojy4->SetYTitle("Count");
		totalQDCangleprojy4->Draw("bar");
		totalQDCangleprojy4->Write();
	
	//totalQDCanglecan3
	totalQDCanglecan3->Divide(3,2);
	
	totalQDCanglecan3->cd(1);
		totalQDCangleFive->Draw("colz");
	
	TH1D *totalQDCangleprojy5 = totalQDCangleFive->ProjectionY("totalqdcangle5",1,180);
	totalQDCanglecan3->cd(4);
		totalQDCangleprojy5->SetTitle("Projection Y");
		totalQDCangleprojy5->SetYTitle("Count");
		totalQDCangleprojy5->Draw("bar");
		totalQDCangleprojy5->Write();
	
	totalQDCanglecan3->cd(2);
		totalQDCangleSix->Draw("colz");
	
	TH1D *totalQDCangleprojy6 = totalQDCangleSix->ProjectionY("totalqdcangle6",1,180);
	totalQDCanglecan3->cd(5);
		totalQDCangleprojy6->SetTitle("Projection Y");
		totalQDCangleprojy6->SetYTitle("Count");
		totalQDCangleprojy6->Draw("bar");
		totalQDCangleprojy6->Write();
	
	totalQDCanglecan3->cd(3);
		totalQDCangleSeven->Draw("colz");
	
	TH1D *totalQDCangleprojy7 = totalQDCangleSeven->ProjectionY("totalqdcangle7",1,180);
	totalQDCanglecan3->cd(6);
		totalQDCangleprojy7->SetTitle("Projection Y");
		totalQDCangleprojy7->SetYTitle("Count");
		totalQDCangleprojy7->Draw("bar");
		totalQDCangleprojy7->Write();
		
	f2->Close();
	
	//totalQDCslantcan
	totalQDCslantcan->Divide(2,1);
	
	totalQDCslantcan->cd(1);
		totalQDCslant->SetXTitle("Slant Depth");
		totalQDCslant->SetYTitle("total QDC");
		totalQDCslant->GetYaxis()->SetTitleOffset(1.3);
		totalQDCslant->Draw("colz");
	
	totalQDCslantcan->cd(2);
		totalQDCslant->FitSlicesY();
		TH1D *totalQDCslant_0 = (TH1D*)gDirectory->Get("totalQDCslant_0");
		totalQDCslant_0->Draw();
		
	//angleOverview canvas
	angleOverview->Divide(2,2);
	
	angleOverview->cd(1);
	//northEastTop->GetYaxis()->SetLimits(7,13);
	northEastTop->SetStats(0);
	northEastTop->Draw("box");
	angleOverview->DrawFrame(-7, -12.5, -0.5, -6.5,"North East Top");
	northEastTop->Draw("samebox");
	//ReverseXAxis(northEastTop);
	//ReverseYAxis(northEastTop);
	
	angleOverview->cd(2);
	//northWestTop->GetYaxis()->SetLimits(7,13);
	northWestTop->SetStats(0);
	northWestTop->Draw("box");
	angleOverview->DrawFrame(-6.5, -12.5, -0.5, -6.5,"North West Top");
	northWestTop->Draw("samebox");
	//ReverseXAxis(northWestTop);
	//ReverseYAxis(northWestTop);
	
	angleOverview->cd(3);
	//southEastTop->GetYaxis()->SetLimits(7,13);
	//southEastTop->GetYaxis()->SetRange(6.5,12.5);
	southEastTop->SetStats(0);
	southEastTop->Draw("box");
	angleOverview->DrawFrame(-6.5, -12.5, -0.5, -6.5,"South East Top");
	southEastTop->Draw("samebox");
	//ReverseXAxis(southEastTop);
	//ReverseYAxis(southEastTop);
	
	angleOverview->cd(4);
	//southWestTop->GetYaxis()->SetLimits(7,13);
	southWestTop->SetStats(0);
	southWestTop->Draw();
	angleOverview->DrawFrame(-6.5, -12.5, -0.5, -6.5,"South West Top");
	southWestTop->Draw("samebox");
	//ReverseXAxis(southWestTop);
	//ReverseYAxis(southWestTop);

	char graphscan[150];
		sprintf(graphscan,"../VetoDisplay/output/plots/graphs.pdf");
		graphs->Print(graphscan,"pdf");
		
	char qdccan[150];
		sprintf(qdccan,"../VetoDisplay/output/plots/qdcDisplay.pdf");
		QDCcanvas->Print(qdccan,"pdf");
	
	char projcanvas[150];
		sprintf(projcanvas,"../VetoDisplay/output/plots/Projections.pdf");
		projcan->Print(projcanvas,"pdf");
	
	char qdctotcan[150];
		sprintf(qdctotcan,"../VetoDisplay/output/plots/QDCtotal.pdf");
		qdctotalcan->Print(qdctotcan,"pdf");
		
	char qdctotncan[150];
        sprintf(qdctotncan,"../VetoDisplay/output/plots/QDCtotaln.pdf");
        qdctotaln->Print(qdctotncan,"pdf");
        
    char thetaPhicanprint[150];
        sprintf(thetaPhicanprint,"../VetoDisplay/output/plots/thetaPhi.pdf");
        thetaPhican->Print(thetaPhicanprint,"pdf");
        
    char costhetaPhicanprint[150];
        sprintf(costhetaPhicanprint,"../VetoDisplay/output/plots/cosThetaPhi.pdf");
        costhetaPhican->Print(costhetaPhicanprint,"pdf");

	char QDCanglecanprint[150];
		sprintf(QDCanglecanprint,"../VetoDisplay/output/plots/QDCangle.pdf");
		QDCanglecan->Print(QDCanglecanprint,"pdf");
		
	char QDCslantcanprint[150];
		sprintf(QDCslantcanprint,"../VetoDisplay/output/plots/QDCslant.pdf");
		QDCslantcan->Print(QDCslantcanprint,"pdf");
		
	char QDCanglecan2print[150];
		sprintf(QDCanglecan2print,"../VetoDisplay/output/plots/QDCangle2.pdf");
		QDCanglecan2->Print(QDCanglecan2print,"pdf");
		
	char QDCanglecan3print[150];
		sprintf(QDCanglecan3print,"../VetoDisplay/output/plots/QDCangle3.pdf");
		QDCanglecan3->Print(QDCanglecan3print,"pdf");
		
	char inThrucanprint[150];
		sprintf(inThrucanprint,"../VetoDisplay/output/plots/inThru.pdf");
		inThrucan->Print(inThrucanprint,"pdf");
		
	char totalQDCanglecanprint[150];
		sprintf(totalQDCanglecanprint,"../VetoDisplay/output/plots/totalQDCangle.pdf");
		totalQDCanglecan->Print(totalQDCanglecanprint,"pdf");
		
	char totalQDCanglecanprint2[150];
		sprintf(totalQDCanglecanprint2,"../VetoDisplay/output/plots/totalQDCangle2.pdf");
		totalQDCanglecan2->Print(totalQDCanglecanprint2,"pdf");
		
	char totalQDCanglecanprint3[150];
		sprintf(totalQDCanglecanprint3,"../VetoDisplay/output/plots/totalQDCangle3.pdf");
		totalQDCanglecan3->Print(totalQDCanglecanprint3,"pdf");
		
	char totalQDCslantcanprint[150];
		sprintf(totalQDCslantcanprint,"../VetoDisplay/output/plots/totalQDCslant.pdf");
		totalQDCslantcan->Print(totalQDCslantcanprint,"pdf");
		
	char angleOverviewprint[150];
		sprintf(angleOverviewprint,"../VetoDisplay/output/plots/angleOverview.pdf");
		angleOverview->Print(angleOverviewprint,"pdf");

	f->Close();
	return 0;
}
