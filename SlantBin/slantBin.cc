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
		//push back future values here
		
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
	
	TH2F *h2 = (TH2F*)f->Get("slantLowTheta");
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
