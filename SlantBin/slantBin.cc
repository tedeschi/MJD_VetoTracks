// slantBin.cc
// Sample tracks through "detector squares"
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "slantBin.hh"
using namespace std;

map<string,vector<vector<double>>> doTheParse(vector<vector<string>> theParse);

void slantBin() {
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
	/*
	map<string, vector<vector<double>>>::iterator it;
	for (it = bottomCoords.begin(); it != bottomCoords.end(); ++it) {
		vector<vector<double>> tmpvcr = (it->second);
		cout << it->first << endl;
		for(int i = 0; i < tmpvcr.size(); i++) {
			for(int j = 0; j < tmpvcr[i].size(); j++) {
				cout << tmpvcr[i][j] << " ";	
			}
		cout << endl;
		}
	}
	*/
  
	TFile *f = new TFile("../SlantDepth/output/slantHist.root");
	TH2F *h = (TH2F*)f->Get("slantLowTheta");
	
	//now comes the permuation dance
	map<string, vector<vector<double>>>::iterator it;
	for (it = topCoords.begin(); it != topCoords.end(); ++it) {
		vector<vector<double>> tmpvcr = (it->second); //a vector (each of the corners) of vectors (the x,y,z coordinates) for top squares
		
		vector<vector<double>> tmpvcr2;
		map<string, vector<vector<double>>>::iterator it2;
		for (it2 = bottomCoords.begin(); it2 != bottomCoords.end(); ++it2) {
			tmpvcr2 = (it2->second);
			
			
			//needs to be done 1000 times
			double xTopRand; //goal is to generate random numbers between tmpvcr[0][0] and tmpvcr[2][0]
			double yTopRand; //goal is to generate random numbers between tmpvcr[0][1] and tmpvcr[2][1]
			double zTop; //just takes tmpvcr[0][2]
			double xBotRand; //goal is to generate random numbers between tmpvcr2[0][0] and tmpvcr2[2][0]
			double yBotRand; //goal is to generate random numbers between tmpvcr2[0][1] and tmpvcr2[2][1]
			double zBot; //just takes tmpvcr[0][2]
			double lo;
			double hi;
			
			if(tmpvcr[0][0] > tmpvcr[2][0]) {
				lo = tmpvcr[0][0];
				hi = tmpvcr[2][0];
			}
			else {
				lo = tmpvcr[2][0];
				hi = tmpvcr[0][0];
			}
			
			//xTopRand = rand()%(hi-lo+1)+lo;
			
			//once we have random numbers we can calculate phi and theta, and interpolate the slantDepth
			TVector3 r1;
			TVector3 r2;
			TVector3 t1;
			double phi;
			double theta;
			r1.SetXYZ(xTopRand,yTopRand,zTop);
			r2.SetXYZ(xBotRand,yBotRand,zBot);
			t1 = r2-r1;
	
			phi = t1.Phi()*(180/M_PI); //180/pi to get deg rather than rad
			theta = t1.Theta()*(180/M_PI);
			
			double slantDepth = h->Interpolate(phi, theta);
			
			//cout << it->first << endl;
			//cout << tmpvcr[0][0] << " " << tmpvcr[2][0] << endl;
			//cout << xTopRand << endl;
		}
		
		
	}
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
	slantBin();
	return 0;
}
