// accCalc.cc
// Sample tracks through "detector squares"
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "accCalc.hh"
using namespace std;
static int howManyIters = 1000000;

map<string,vector<vector<double>>> doTheParse(vector<vector<string>> theParse);
bool passThroughSquare(double p1x, double p1y, double p1z, double p2x, 
		double p2y, double p2z, double sx1, double sx2, double sy1, double sy2, double sz);

void accCalc() {
	vector<vector<string>> bottom; //the non-parsed vectors read in
	vector<vector<string>> top;
	map<string,vector<vector<double>>> bottomCoords; //mapped detector numbers with their parsed coordinates
	map<string,vector<vector<double>>> topCoords;
	
	ifstream infile;
	infile.open("DetCoordsUpdated2.txt");
	string line,cell;
	int counter = 0;
	
	//for simplicities sake, we read everything into vectors without parsing
	while(getline(infile,line)) { 
		stringstream lineStream(line);
		
		vector<string> aCell;
		while(getline(lineStream,cell,'\t')) {
			if(cell[0] == '"') { //get rid of quotation marks to make parsing easier
				cell = cell.substr(1,cell.size()-2);
			}
			if(cell[cell.length()-1] == '"') {
				cell = cell.substr(0,cell.size()-1);
			}
			aCell.push_back(cell);
			
		}
		counter++;
		
		if(counter <= 36) {
			bottom.push_back(aCell);
		}
		if(counter > 37 && counter <= 41) {
			top.push_back(aCell);
		}
		//push back (to the) future values here
		
		// ------------------------------
	}
	infile.close();
	
	//now to parse
	bottomCoords = doTheParse(bottom);
	topCoords = doTheParse(top);
	
	//an example on how to iterate through the map
	
	/*map<string, vector<vector<double>>>::iterator it;
	for (it = bottomCoords.begin(); it != bottomCoords.end(); ++it) {
		vector<vector<double>> tmpvcr = (it->second);
		cout << it->first << endl;
		for(int i = 0; i < tmpvcr.size(); i++) {
			for(int j = 0; j < tmpvcr[i].size(); j++) {
				cout << tmpvcr[i][j] << " ";	
			}
		cout << endl;
		}
	}*/
	

//
// histogram definitions
//  
	TFile *f = new TFile("../SlantDepth/output/slantHist.root");
	TFile *f2 = new TFile("AccCalc.root","RECREATE");
	//TH2D *h = (TH2D*)_file0->Get("");
	
	TH2D *hTopXY[144]; 
	TH2D *hBotXY[144];
	char *histname = new char[20];

	//TH3D *upperPart = new TH3D("upperPart","upperPart", 1000, -550, 550, 1000, -550, 550, 1000, -200, 300);
	TH1D *cosThetaPlot = new TH1D("cosThetaPlot","cosThetaUniformDist", 1000, -1, 1);
	TH2D *xyPlot = new TH2D("xyPlot","xyUniformDist", 1000, -550, 550, 1000, -550, 550);
	TH2D *xyUpperPlot = new TH2D("xyUpperPlot","xyUpperPlot", 1000, -550, 550, 1000, -550, 550);
	TH2D *thetaPhiUpperPlot = new TH2D("thetaPhiUpperPlot","thetaPhiUpperPlot", 1000, 0, 1, 1000, 0, 360);
	
	for (int i=0; i<144; i++) {
		sprintf(histname, "hTopXY_%d",i);
		hTopXY[i] = new TH2D(histname,histname, 1000, -550, 550, 1000, -550, 550);
		sprintf(histname, "hBotXY_%d",i);
		hBotXY[i] = new TH2D(histname,histname, 1000, -550, 550, 1000, -550, 550);
	}

	double howManyHits[145] = {0};
	double howManyTries[145] = {0};
	double detAcc[145] = {0};
	string theDetNums[145];
	int theCounter = 0;
	//now comes the permuation dance
	// outer loop for top x,y,z
	map<string, vector<vector<double>>>::iterator it;
	for (it = topCoords.begin(); it != topCoords.end(); ++it) {
		vector<vector<double>> tmpvcr = (it->second); //a vector (each of the corners) of vectors (the x,y,z coordinates) for top squares
		
		double xtlo,ytlo;
		double xthi,ythi;	
		double zTop = tmpvcr[0][2];	
		if(tmpvcr[0][0] < tmpvcr[2][0]) {
			xtlo = tmpvcr[0][0];
			xthi = tmpvcr[2][0];
		}
		else {
			xtlo = tmpvcr[2][0];
			xthi = tmpvcr[0][0];
		}

		if(tmpvcr[0][1] < tmpvcr[2][1]) {
			ytlo = tmpvcr[0][1];
			ythi = tmpvcr[2][1];
		}
		else {
			ytlo = tmpvcr[2][1];
			ythi = tmpvcr[0][1];
		}
		//next loop for botton x,y,z
		vector<vector<double>> tmpvcr2;
		map<string, vector<vector<double>>>::iterator it2;
		for (it2 = bottomCoords.begin(); it2 != bottomCoords.end(); ++it2) {
			tmpvcr2 = (it2->second);

			double xblo,yblo;
			double xbhi,ybhi;
			double zBot= tmpvcr2[0][2];		
			if(tmpvcr2[0][0] < tmpvcr2[2][0]) {
				xblo = tmpvcr2[0][0];
				xbhi = tmpvcr2[2][0];
			}
			else {
				xblo = tmpvcr2[2][0];
				xbhi = tmpvcr2[0][0];
			}
			
			if(tmpvcr2[0][1] < tmpvcr2[2][1]) {
				yblo = tmpvcr2[0][1];
				ybhi = tmpvcr2[2][1];
			}
			else {
				yblo = tmpvcr2[2][1];
				ybhi = tmpvcr2[0][1];
			}
			
			theCounter++;
			theDetNums[theCounter] = (it->first);
			theDetNums[theCounter] += " and ";
			theDetNums[theCounter] += (it2->first);
			// innermost loop iterates within the top/bot combo
			for (int i=0;i<howManyIters;i++){
				//acceptance -------------------------------------------
				//bottom circle
				double theLowerX = gRandom->Uniform(-500.0, 500.0); //make a box of radius 5 meters
				double theLowerY = gRandom->Uniform(-500.0, 500.0);
				double theLowerZ = -200.0;
				double outsideBound = pow(theLowerX,2) + pow(theLowerY,2); //eq for the circle inside box
				while(outsideBound >= pow(500,2)) { //if the point falls outside the circle, select a new one
					theLowerX = gRandom->Uniform(-500.0, 500.0);  
					theLowerY = gRandom->Uniform(-500.0, 500.0);
					outsideBound = pow(theLowerX,2) + pow(theLowerY,2);
				}

				//upper half sphere
				double theUpperPhi = gRandom->Uniform(0, 360);
				double theUpperTheta = acos(gRandom->Uniform(0, 1));
				double theUpperX = 500 * cos(theUpperPhi) * sin(theUpperTheta);
				double theUpperY = 500 * sin(theUpperPhi) * sin(theUpperTheta);
				double theUpperZ = 500 * cos(theUpperTheta) + theLowerZ; //polar coordinates from origin
			//cout << theUpperTheta << endl;
				
				cosThetaPlot->Fill(theUpperTheta);
				xyPlot->Fill(theLowerX, theLowerY);
				thetaPhiUpperPlot->Fill(cos(theUpperTheta), theUpperPhi);
				xyUpperPlot->Fill(theUpperX, theUpperY);
				//upperPart->Fill(theUpperX, theUpperY, theUpperZ);
				
				//check if the line went through the dets.
				if(passThroughSquare(theUpperX, theUpperY, theUpperZ, theLowerX, theLowerY, theLowerZ, xtlo, xthi, ytlo, ythi, zTop) == true 
					&& passThroughSquare(theUpperX, theUpperY, theUpperZ, theLowerX, theLowerY, theLowerZ, xblo, xbhi, yblo, ybhi, zBot) == true) {
					howManyHits[theCounter]++;
					hTopXY[theCounter-1]->Fill(theUpperX, theUpperY);
					hBotXY[theCounter-1]->Fill(theLowerX, theLowerY);
					/*
					cout << "its a hit!" << endl;
					cout << "TheUpperX: " << theUpperX << " TheUpperY: " << theUpperY << " TheUpperZ: " << theUpperZ << endl;
					cout << "TheLowerX: " << theLowerX << " TheLowerY: " << theLowerY << " TheLowerZ: " << theLowerZ << endl;
					cout << xtlo <<" "<< xthi<<" " << ytlo<<" " << ythi<<" " << zTop << endl;
					cout << xblo<<" " << xbhi<<" " << yblo<<" " << ybhi<<" " << zBot << endl << endl;
					*/
				}
				howManyTries[theCounter]++;
				//cout << passThroughSquare(2,2,2,-2,-2,-2,-1,1,-1,1,0);
				// -----------------------------------------------------
			}
		}
	}

	f2->Write();		
	f2->Close();
	
	double totalAcc = 0;
	for(int i = 1; i <= 144; i++) {
		detAcc[i] = (howManyHits[i]/howManyTries[i])*(2*M_PI*pow(500.0,2)*M_PI); //(2*M_PI*pow(500.0,2)*M_PI) ~= 5,000,000
		
		cout << "Det nums: " << theDetNums[i] << " How many Hits: " << howManyHits[i] << " How many tries: " << howManyTries[i] << " Acc: " << detAcc[i] << endl;
		totalAcc += detAcc[i];
	}
	cout << "Total Acc: " << totalAcc << endl;
}

//takes in two 3d points and x,y,z coords of square
bool passThroughSquare(double p1x, double p1y, double p1z, double p2x, 
		double p2y, double p2z, double sx1, double sx2, double sy1, double sy2, double sz) {
	
	//get parametric eq for line
	double parx = p1x - p2x;
	double pary = p1y - p2y;
	double parz = p1z - p2z;
	//(p1x + parx*t, p1y + pary*t, p1z + parz*t)
	
	double t = (sz - p1z)/parz; //find where the line is when z=zsquare
	double parxAtsz = p1x + (parx*t);
	double paryAtsz = p1y + (pary*t);
	//cout << "p1x:" << p1x << cout << " p2x:" << p2x << " sx1:" << sx1 << " sx2:" << sx2 << endl;
	
	if(parxAtsz >= sx1 && parxAtsz <= sx2) { //if the line falls in the square
		if(paryAtsz >= sy1 && paryAtsz <= sy2) {
			return true;
		}
	}
	
	return false;
}

map<string,vector<vector<double>>> doTheParse(vector<vector<string>> theParse) {
	map<string,vector<vector<double>>> theReturnMap;
	
	for(int i = 0; i < theParse.size(); i++) {
		vector<vector<double>> theCoords;
		
		for(int j = 1; j < theParse[i].size(); j++) {
			string split;
			size_t pos = 0;
			vector<double> theCoords2;
			
			while((pos = theParse[i][j].find(',')) != string::npos) {
				split = theParse[i][j].substr(0,pos);
				theParse[i][j].erase(0, pos + 1);
				//cout << split << endl;
				
				theCoords2.push_back(stod(split));
			}
			theCoords2.push_back(stod(theParse[i][j]));
			theCoords.push_back(theCoords2);
		}
		theReturnMap[theParse[i][0]] = theCoords;
		//cout << endl;
	}
	
	return theReturnMap;
}

int main() {
	accCalc();
	return 0;
}
