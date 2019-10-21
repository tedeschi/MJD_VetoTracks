// slantBin.cc
// Sample tracks through "detector squares"
// Bradley McClain, University of South Carolina

//-----------------------------------------------------------------------
#include "slantBin.hh"
using namespace std;

map<string,vector<vector<double>>> doTheParse(vector<vector<string>> theParse);
static list<string> insertionOrder; //this is needed because maps DO NOT keep insertion order

void slantBin() {
	vector<vector<string>> bottom; //the non-parsed vectors read in
	vector<vector<string>> top;
	map<string,vector<vector<double>>> bottomCoords; //mapped detector numbers with their parsed coordinates
	map<string,vector<vector<double>>> topCoords;
	
	ifstream infile;
	infile.open("DetCoordsUpdated3.txt");
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
	insertionOrder.pop_back(); //I know, it's ugly. But it works.
	insertionOrder.pop_back();
	insertionOrder.pop_back();
	insertionOrder.pop_back();
	
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
	TH2F *h = (TH2F*)f->Get("slantLowTheta");

	TFile *f2 = new TFile("slantBin.root","RECREATE");
	ofstream outputfile;
	outputfile.open("slantBin.txt");

	TH1D *hSlantDepth[145];
	TH2D *hTopXY[145]; 
	TH2D *hBotXY[145]; 
	TH2D *hThetaPhi[145];
	char *histname = new char[20];

	for (int i=1; i<=144; i++){
		sprintf(histname, "hSlantDepth_%d",i);
//		hSlantDepth[i] = new TH1D(histname,histname, 1000, 1400., 2800.);
		hSlantDepth[i] = new TH1D(histname,histname, 1000, 4., 8.);
		sprintf(histname, "hTopXY_%d",i);
		hTopXY[i] = new TH2D(histname,histname, 240, -120., 120., 240, -120. , 120.);
		sprintf(histname, "hBotXY_%d",i);
		hBotXY[i] = new TH2D(histname,histname, 240, -120., 120., 240, -120. , 120.);
		sprintf(histname, "hThetaPhi_%d",i);
		hThetaPhi[i] = new TH2D(histname,histname, 100, 0, 360., 100, 0. , 90.);
	}
	gRandom->SetSeed(0); // uses time for seed
	int idet=1;
	
	//now comes the permuation dance
	// outer loop for top x,y,z
	map<string, vector<vector<double>>>::iterator it;
	for (it = topCoords.begin(); it != topCoords.end(); ++it) {
		vector<vector<double>> tmpvcr = (it->second); //a vector (each of the corners) of vectors (the x,y,z coordinates) for top squares

		// check top coords
		//cout << tmpvcr[0][0]<< " " << tmpvcr[2][0] << " " << tmpvcr[0][2] << endl;
		//cout << tmpvcr[0][1]<< " " << tmpvcr[2][1] << " " << tmpvcr[0][2] << endl;
		double xtlo,ytlo;
		double xthi,ythi;
		double xtmid,ytmid;			
		if(tmpvcr[0][0] < tmpvcr[2][0]) {
			xtlo = tmpvcr[0][0];
			xthi = tmpvcr[2][0];
		}
		else {
			xtlo = tmpvcr[2][0];
			xthi = tmpvcr[0][0];
		}
		xtmid = xtlo + 0.5*(xthi-xtlo);

		if(tmpvcr[0][1] < tmpvcr[2][1]) {
			ytlo = tmpvcr[0][1];
			ythi = tmpvcr[2][1];
		}
		else {
			ytlo = tmpvcr[2][1];
			ythi = tmpvcr[0][1];
		}
		ytmid = ytlo + 0.5*(ythi-ytlo);

		
		//next loop for botton x,y,z
		//map<string, vector<vector<double>>>::iterator it2;
		list<string>::iterator it2;
		for (it2 = insertionOrder.begin(); it2 != insertionOrder.end(); ++it2) {
			string theIteration = (string)*it2;
			vector<vector<double>> tmpvcr2 = bottomCoords.find(theIteration)->second;
			
			//check bottom coords
			//cout << tmpvcr2[0][0] << " " << tmpvcr2[2][0] << " " << tmpvcr2[0][2] << endl;
			//cout << tmpvcr2[0][1] << " " << tmpvcr2[2][1] << " " << tmpvcr2[0][2] << endl;
			double xblo,yblo;
			double xbhi,ybhi;			
			double xbmid,ybmid;			
			if(tmpvcr2[0][0] < tmpvcr2[2][0]) {
				xblo = tmpvcr2[0][0];
				xbhi = tmpvcr2[2][0];
			}
			else {
				xblo = tmpvcr2[2][0];
				xbhi = tmpvcr2[0][0];
			}
			xbmid = xblo + 0.5*(xbhi-xblo);
			
			if(tmpvcr2[0][1] < tmpvcr2[2][1]) {
				yblo = tmpvcr2[0][1];
				ybhi = tmpvcr2[2][1];
			}
			else {
				yblo = tmpvcr2[2][1];
				ybhi = tmpvcr2[0][1];
			}
			ybmid = yblo + 0.5*(ybhi-yblo);

			// innermost loop iterates within the top/bot combo
			for (int i=0;i<10000;i++){
				double xTopRand=gRandom->Uniform(xtlo,xthi);        //goal is to generate random numbers between tmpvcr[0][0] and tmpvcr[2][0]
				double yTopRand=gRandom->Uniform(ytlo,ythi);        //goal is to generate random numbers between tmpvcr[0][1] and tmpvcr[2][1]
				double zTop = tmpvcr[0][2];

				double xBotRand=gRandom->Uniform(xblo,xbhi);      //goal is to generate random numbers between tmpvcr2[0][0] and tmpvcr2[2][0]
				double yBotRand=gRandom->Uniform(yblo,ybhi);      //goal is to generate random numbers between tmpvcr2[0][1] and tmpvcr2[2][1]
				double zBot= tmpvcr2[0][2];                                         //just takes tmpvcr2[0][2]

				hTopXY[idet]->Fill(xTopRand,yTopRand);
				hBotXY[idet]->Fill(xBotRand,yBotRand);
		
				//calculate phi and theta, and interpolate the slantDepth
				TVector3 r1;
				TVector3 r2;
				TVector3 t1;
				double phi;
				double theta;
				r2.SetXYZ(xTopRand,yTopRand,zTop);         
				r1.SetXYZ(xBotRand,yBotRand,zBot);
				t1 = r2-r1;
				phi = t1.Phi()*(180/M_PI); //180/pi to get deg rather than rad
				if(phi < 0) {
					phi += 360; //makes negative phi's the correct positive
				}
				theta = t1.Theta()*(180/M_PI);
				hThetaPhi[idet]->Fill(phi,theta);

				double slantDepth = h->Interpolate(phi, theta);
				double slantDepth_kmwe = slantDepth*2.86/1000.;
				//weight contribution by cross section (use Meei-Hime fit)
				double weight=8.6E-6*exp(-slantDepth_kmwe/0.45)+0.44E-6*exp(-slantDepth_kmwe/0.87);
				hSlantDepth[idet]->Fill(slantDepth_kmwe,weight);
			}

			TVector3 r1mid;
			TVector3 r2mid;
			TVector3 t1mid;
			r2mid.SetXYZ(xtmid,ytmid,tmpvcr[0][2]);         
			r1mid.SetXYZ(xbmid,ybmid,tmpvcr2[0][2]);
			t1mid = r2mid-r1mid;
			double phimid;
			double thetamid;

			phimid = t1mid.Phi()*(180/M_PI); //180/pi to get deg rather than rad
			if(phimid < 0) {
				phimid += 360; //makes negative phi's the correct positive
			}
			thetamid = t1mid.Theta()*(180/M_PI);
			
			double slantDepthmid = h->Interpolate(phimid, thetamid);

			//cout << "det=" << idet << " " << "theta=" << thetamid << " " <<"phi=" << phimid << " " <<"depth=" << slantDepthmid << endl;
			cout  << hSlantDepth[idet]->GetMean()<< " "  << hSlantDepth[idet]->GetStdDev() << endl;
			outputfile << idet << " " << hSlantDepth[idet]->GetMean()<< " "  << hSlantDepth[idet]->GetStdDev() << endl;
			idet++;

		}
	}

	f2->Write();		
	f2->Close();
	outputfile.close();
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
		insertionOrder.push_back(theParse[i][0]);
		//cout << endl;
	}
	
	return theReturnMap;
}

int main() {
	slantBin();
	return 0;
}
