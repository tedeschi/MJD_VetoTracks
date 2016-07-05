// VetoAnalysis.cc
// macro to analyze Majorana Data
// David J Tedeschi, Clint Wiseman and Bradley McClain, University of South Carolina
//
// .x VetoDisplay.C++
// to clear output all output subdirectories:
// cd to the main output directory "assembly_output" and use the command:
// find . ! -type d -exec rm '{}' \;
// (careful, dangourous command if used from any other directory)
//-----------------------------------------------------------------------
#include "VetoDisplay.hh"

using namespace std;

string inputFile = "../Data/vetoskims/vList_DS1_1-6.txt"; //EDIT THIS FOR INPUT FILE

int isNextTo(int panel1, int panel2);
int isLayerHit(int panel1, int panel2);
std::vector<Double_t> hitLocation(int panel1, int panel2);

int coloring(int qdc,int mode){

	int v[15]={0};

	// x^2 coloring (separates low-qdc events better) 
	// 15.556 x^2 + 500, x=1,15
	if (mode == 1)
		v[0]=500,  v[1]=562,  v[2]=640,   v[3]=749,   v[4]=889,   v[5]=1060,  v[6]=1262,  v[7]=1496,
		v[8]=1760, v[9]=2056, v[10]=2382, v[11]=2740, v[12]=3129, v[13]=3549, v[14]=4000;  

	// log coloring (separates high-qdc events better)
	// range: 1477*log(i), scaled for i=15 == 4000
	if (mode == 2)
		v[0]=500,  v[1]=1024, v[2]=1623,  v[3]=2048,  v[4]=2377,  v[5]=2646,  v[6]=2874,  v[7]=3071,
		v[8]=3245, v[9]=3401, v[10]=3542, v[11]=3670, v[12]=3788, v[13]=3898, v[14]=4000;

	// linear coloring (separate qdc events linearly)
	if (mode == 3)
		v[0]=500,  v[1]=750,  v[2]=1000,  v[3]=1250,  v[4]=1500,  v[5]=1750,  v[6]=2000,  v[7]=2250,
		v[8]=2500, v[9]=2750, v[10]=3000, v[11]=3250, v[12]=3500, v[13]=3750, v[14]=4000;

	// set color table
	if (qdc < v[0])                 { TColor *color = gROOT->GetColor(1000);  color->SetRGB(0.00, 0.00, 0.00);  return 1000;}   // black (below threshold)
	if (qdc > v[0] && qdc <=v[1])   { TColor *color = gROOT->GetColor(1001);  color->SetRGB(0.02, 0.00, 0.75);  return 1001;}   // blue side
	if (qdc > v[1] && qdc <=v[2])   { TColor *color = gROOT->GetColor(1002);  color->SetRGB(0.08, 0.00, 0.69);  return 1002;}
	if (qdc > v[2] && qdc <=v[3])   { TColor *color = gROOT->GetColor(1003);  color->SetRGB(0.15, 0.00, 0.64);  return 1003;}
	if (qdc > v[3] && qdc <=v[4])   { TColor *color = gROOT->GetColor(1004);  color->SetRGB(0.21, 0.00, 0.59);  return 1004;}
	if (qdc > v[4] && qdc <=v[5])   { TColor *color = gROOT->GetColor(1005);  color->SetRGB(0.27, 0.01, 0.53);  return 1005;}
	if (qdc > v[5] && qdc <=v[6])   { TColor *color = gROOT->GetColor(1006);  color->SetRGB(0.33, 0.01, 0.48);  return 1006;}
	if (qdc > v[6] && qdc <=v[7])   { TColor *color = gROOT->GetColor(1007);  color->SetRGB(0.40, 0.01, 0.43);  return 1007;}
	if (qdc > v[7] && qdc <=v[8])   { TColor *color = gROOT->GetColor(1008);  color->SetRGB(0.46, 0.01, 0.37);  return 1008;}
	if (qdc > v[8] && qdc <=v[9])   { TColor *color = gROOT->GetColor(1009);  color->SetRGB(0.52, 0.02, 0.32);  return 1009;}
	if (qdc > v[9] && qdc <=v[10])  { TColor *color = gROOT->GetColor(1010);  color->SetRGB(0.58, 0.02, 0.27);  return 1010;}
	if (qdc > v[10] && qdc <=v[11]) { TColor *color = gROOT->GetColor(1011);  color->SetRGB(0.65, 0.02, 0.21);  return 1011;}
	if (qdc > v[11] && qdc <=v[12]) { TColor *color = gROOT->GetColor(1012);  color->SetRGB(0.71, 0.02, 0.16);  return 1012;}
	if (qdc > v[12] && qdc <=v[13]) { TColor *color = gROOT->GetColor(1013);  color->SetRGB(0.77, 0.02, 0.11);  return 1013;}
	if (qdc > v[13] && qdc <=v[14]) { TColor *color = gROOT->GetColor(1014);  color->SetRGB(0.84, 0.02, 0.05);  return 1014;}
	if (qdc > v[14])                { TColor *color = gROOT->GetColor(1015);  color->SetRGB(0.90, 0.03, 0.00);  return 1015;}    // red side
	else return 0;
}


//Panel numbers classified by side on the 0-31 scale
Int_t bottomPanels[12] = {0,1,2,3,4,5,6,7,8,9,10,11};
Int_t topPanels[4] = {17,18,20,21};
Int_t northPanels[4] = {15,16,19,23};
Int_t westPanels[4] = {12,13,14,22};
Int_t eastPanels[4] = {28,29,30,31};
Int_t southPanels[4] = {24,25,26,27};

//mother coordinate system for each panel in x,y,z
Double_t mcscoords[32][3];
// canvas is global
TCanvas *ecan;
//Phi & Theta global to allow for plotting
Double_t phi;
Double_t theta;
Double_t pi = 3.1415926535897;
//---------------------------------------------------------------------------------------------------------------------
void DrawEvent(Int_t qdcVals[], Int_t numberOfPanelsHit, Int_t totalQDC, Int_t runNumber, Int_t eventCount)
{
//--- Definition of a simple geometry
   //gSystem->Load("libGeom");
   TGeoManager *geom = new TGeoManager("Assemblies","Geometry using assemblies");
   
   //--- define some materials
   // (const char* name, Double_t a, Double_t z, Double_t rho, Double_t radlen = 0, Double_t intlen = 0)
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   //TGeoMaterial *matAl = new TGeoMaterial("Al", 26.98,13,2.7);
   //TGeoMaterial *Fe = new TGeoMaterial("Fe",55.845,26,7.87);
   //TGeoMaterial *matSi = new TGeoMaterial("Si", 28.085,14,2.329);
   //TGeoMaterial *Cu = new TGeoMaterial("Cu",63.549,29,8.92);
   
   //--- define some tracking media
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);
   //TGeoMedium *Al = new TGeoMedium("Aluminium",2, matAl);
   
   //--- create the starting point of the volume structure - 
   // half length units. All are in cm.
   Double_t xtop=160;
   Double_t ytop=160;
   Double_t ztop=160;
   TGeoVolume *top = geom->MakeBox("TOP", Vacuum, xtop,ytop,ztop); //gets placed within this box
   geom->SetTopVolume(top);
   geom->SetTopVisible(1);
   
	// declare panels
	TGeoVolume *panel[32];

	// bottom veto panels-------------------------------------------------------------
	// make box for each layer and fill with 6 panels
	// panel 0-5  lower bottom
	// panel 6-11 upper bottom
	// measurements taken from Dr. David Tedeschi on site and
	// Majorana_MuonFluxV21.pdf
	// The southwest side was chosen to be the 'most correct' in terms of offset
	Double_t bottomDetectorWidth = 32.0; //total width of a bottom panel detector
	Double_t totalBottomHeight = 43.18; //upper + lower bottom panel heigth together
	Double_t spaceBetweenDetectorsy = 2.54; //from detector edge to next detector edge in the y axis
	Double_t detectorHeight = 2.54; //from detector edge to next detector edge in the z axis
	
	Double_t xlayer1 = 182.0 / 2; //length of lower bottom detectors in cm
	Double_t xlayer2 = 223.0 / 2; //length of upper bottom detectors in cm
	Double_t ylayer1 = ((bottomDetectorWidth * 6) + (spaceBetweenDetectorsy * 6)) / 2; // = 207.24/2
	Double_t zplace = (203.29) / 2; //where to place the bottom panels on the z axis(hieght of upper panels)
	Double_t upperoffsetx = (zplace - xlayer2 - spaceBetweenDetectorsy/2) + (1.93 / 2); //the offset of the upper bottom panels in x axis
	Double_t loweroffsetx = (zplace - xlayer1 - spaceBetweenDetectorsy/2) + (13.627 / 2); //the offset of the lower bottom panels in x axis
	Double_t upperoffsety = (zplace - ylayer1) + (13.597 / 2); //the offset of the upper bottom panels in y axis
	Double_t loweroffsety = (zplace - ylayer1) - (1.27 / 2); //the offset of the lower bottom panels in the y axis
	Double_t zlayer1 = totalBottomHeight / 4;
	
	// we start with the bottom outside boxes
	TGeoVolume *layerbo1 = geom->MakeBox("layerbo1", Vacuum, xlayer1 + (spaceBetweenDetectorsy / 2),ylayer1,zlayer1); // bottom lower outer boxes
	TGeoVolume *layerbo2 = geom->MakeBox("layerbo2", Vacuum, xlayer2 + (spaceBetweenDetectorsy / 2),ylayer1,zlayer1); // bottom upper outer boxes

	// put 6 panels in layer 1
	Double_t xpanel1 = xlayer1;
	Double_t ypanel1 = ylayer1 / 6.;    
	Double_t zpanel1 = zlayer1;
	Char_t panel_name[50];
	Double_t ypos;

	for (Int_t i=0; i<6;i++){
		sprintf(panel_name,"panel%d",i);
		panel[i]= geom->MakeBox(panel_name, Vacuum, xpanel1,ypanel1,zpanel1);
		ypos = -ylayer1+(2*i+1)*ypanel1;  
		layerbo1->AddNode(panel[i],i, new TGeoTranslation(loweroffsetx,ypos,0));

		sprintf(panel_name,"panel%d",i);
		panel[6+i]= geom->MakeBox(panel_name, Vacuum, xlayer2,ypanel1,zpanel1);
		ypos = -ylayer1+(2*i+1)*ypanel1;  
		layerbo2->AddNode(panel[6+i],i, new TGeoTranslation(upperoffsetx,ypos,0));
	}
  
	top->AddNode(layerbo1,1, new TGeoTranslation(-loweroffsetx,loweroffsety,-zplace-3*zlayer1));  // add bottom layers to mother
	layerbo1->SetVisDaughters(false); // make outside layers invisible
	TGeoRotation *rot1 = new TGeoRotation();  
	rot1->RotateZ(90);  
	top->AddNode(layerbo2,2, new TGeoCombiTrans(-upperoffsety,upperoffsetx,-zplace-zlayer1,rot1)); //add second layer
	layerbo2->SetVisDaughters(false);
	
	//now we do the inside bottom boxes
	TGeoVolume *layerbi1 = geom->MakeBox("layerbi1", Vacuum, xlayer2,ylayer1,zlayer1); // bottom upper inner 
	TGeoVolume *layerbi2 = geom->MakeBox("layerbi2", Vacuum, xlayer1,ylayer1,zlayer1); // bottom lower inner

	// put 6 panels in layer 1
	Double_t xpanel2 = xlayer2;
	Double_t ypanel2 = bottomDetectorWidth/2;
	Double_t zpanel2 = zlayer1-spaceBetweenDetectorsy/2;

	for (Int_t i=0; i<6;i++){
		sprintf(panel_name,"panel%d",i);
		panel[i]= geom->MakeBox(panel_name, Vacuum, xlayer1,ypanel2,detectorHeight/2);
		ypos = -ylayer1+(2*i+1)*ypanel2+(spaceBetweenDetectorsy * i + 1) + loweroffsety;
		layerbi1->AddNode(panel[i],i, new TGeoTranslation(-loweroffsetx,ypos,0));
		mcscoords[i][0] = -loweroffsetx;
		mcscoords[i][1] = ypos;
		mcscoords[i][2] = -zplace-3*zlayer1;

		sprintf(panel_name,"panel%d",i);
		panel[6+i]= geom->MakeBox(panel_name, Vacuum, xlayer2,ypanel2,detectorHeight/2);
		ypos = -ylayer1+(2*i+1)*ypanel2+(spaceBetweenDetectorsy *i + 1) + upperoffsety; 
		layerbi2->AddNode(panel[6+i],i, new TGeoTranslation(upperoffsetx,ypos,0));
		mcscoords[6+i][0] = -ypos;
		mcscoords[6+i][1] = upperoffsetx;
		mcscoords[6+i][2] = -zplace-zlayer1;
	}
	rot1->RotateZ(-90); //rotate back after first rotation  
	top->AddNode(layerbi1,1, new TGeoTranslation(0,0,-zplace-3*zlayer1));  // add bottom layers to mother 
	rot1->RotateZ(90);  //then rotate again
	top->AddNode(layerbi2,2, new TGeoCombiTrans(0,0,-zplace-zlayer1,rot1)); //add second layer
	
	// top veto panels-------------------------------------------------------------
	// panel 17,18 top outer
	// panel 20,21 top inner   
	Double_t topxlayer = 203.29 / 2;
	Double_t topylayer = 203.29 / 2;
	Double_t zlayer = 2.54 / 2;
	
	Double_t xpanel = topxlayer;
	Double_t ypanel = topylayer/6.;    
	Double_t zpanel = zlayer;
	
	TGeoVolume *layert1 = geom->MakeBox("layert1", Vacuum, topxlayer,topylayer,zlayer);   
	TGeoVolume *layert2 = geom->MakeBox("layert2", Vacuum, topxlayer,topylayer,zlayer);   
	ypanel = topylayer/2.;    

	for (Int_t i=0; i<2;i++){
	    sprintf(panel_name,"panel%d",i);
	    panel[17+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	    ypos = -topylayer+(2*i+1)*ypanel;  
	    layert1->AddNode(panel[17+i],i, new TGeoTranslation(0,ypos,0));
	    mcscoords[17+i][0] = 0;
	    mcscoords[17+i][1] = ypos;
	    mcscoords[17+i][2] = topxlayer;

	    sprintf(panel_name,"panel%d",i);
	    panel[20+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	    ypos = -topylayer+(2*i+1)*ypanel;  
	    layert2->AddNode(panel[20+i],i, new TGeoTranslation(0,ypos,0));
	    mcscoords[20+i][0] = ypos;
	    mcscoords[20+i][1] = 0;
	    mcscoords[20+i][2] = topxlayer;
  	}
  	// add top layers to mother
  	top->AddNode(layert1,1, new TGeoTranslation(0,0,topxlayer-zlayer + 4*zlayer));
  	TGeoRotation *rot2 = new TGeoRotation();
  	rot2->RotateZ(-90);  
  	top->AddNode(layert2,2, new TGeoCombiTrans(0,0,topxlayer-(3*zlayer) + 4*zlayer,rot2));

	// north veto panels-------------------------------------------------------------
	// panel 15,16 north outer
	// panel 19,23 north inner
	
	Double_t northxlayer = 203.29 / 2;
	Double_t northylayer = 203.29 / 2;
	Double_t halfHoleLength = 50 / 6; //dont have actual dimensions for yet, currently length of cryovats
	
	TGeoVolume *layern1 = geom->MakeBox("layern1", Vacuum, northxlayer,northylayer,zlayer);   
	TGeoVolume *layern2 = geom->MakeBox("layern2", Vacuum, northxlayer,northylayer,zlayer);   

	xpanel = northxlayer;    
	ypanel = northylayer/2.;    

	// 15
	sprintf(panel_name,"panel%d",15);
	panel[15]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfHoleLength,zpanel);
	ypos = -northylayer+(2*0+1)*ypanel;  
	layern1->AddNode(panel[15],15, new TGeoTranslation(0,ypos,0));
	mcscoords[15][0] = northxlayer;
	mcscoords[15][1] = ypos;
	mcscoords[15][2] = 0;
	
	//16
	sprintf(panel_name,"panel%d",16);
	panel[16]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfHoleLength,zpanel);
	ypos = -northylayer+(2*1+1)*ypanel;  
	layern1->AddNode(panel[16],16, new TGeoTranslation(0,ypos,0));
	mcscoords[16][0] = northxlayer;
	mcscoords[16][1] = ypos;
	mcscoords[16][2] = 0;
	    
	//19
	sprintf(panel_name,"panel%d",19);
	panel[19]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfHoleLength,zpanel);
	ypos = -northylayer+(2*0+1)*ypanel;  
	layern2->AddNode(panel[19],19, new TGeoTranslation(0,ypos,0));
	mcscoords[19][0] = northxlayer;
	mcscoords[19][1] = ypos;
	mcscoords[19][2] = 0;
	//23
	sprintf(panel_name,"panel%d",23);
	panel[23]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfHoleLength,zpanel);
	ypos = -northylayer+(2*1+1)*ypanel;  
	layern2->AddNode(panel[23],23, new TGeoTranslation(0,ypos,0));
	mcscoords[23][0] = northxlayer;
	mcscoords[23][1] = ypos;
	mcscoords[23][2] = 0;

	// add north layers to mother
	TGeoRotation *rot3 = new TGeoRotation();
	rot3->RotateY(90);   
	top->AddNode(layern1,1, new TGeoCombiTrans(northxlayer-zlayer,-halfHoleLength,0,rot3));
	top->AddNode(layern2,2, new TGeoCombiTrans(northxlayer-3*zlayer,halfHoleLength,0,rot3));

	// west veto panels-------------------------------------------------------------
	// panel 12,13 west inner
	// panel 14,22 west outer
	//switch role of x and y to get panel shape
 
	Double_t westxlayer = 203.29 / 2;
	Double_t westylayer = 203.29 / 2;
	 
	TGeoVolume *layerw1 = geom->MakeBox("layerw1", Vacuum, westxlayer,westylayer,zlayer);   
	TGeoVolume *layerw2 = geom->MakeBox("layerw2", Vacuum, westxlayer,westylayer,zlayer);

	xpanel = westxlayer;    
	ypanel = (westylayer-2*zlayer)/2.;

	// 12
	sprintf(panel_name,"panel%d",12);
	panel[12]= geom->MakeBox(panel_name, Vacuum, ypanel+halfHoleLength,xpanel,zpanel);
	ypos = -westylayer+(2*0+1)*ypanel;  
	layerw1->AddNode(panel[12],12, new TGeoTranslation(ypos,0,0));
	mcscoords[12][0] = ypos;
	mcscoords[12][1] = westxlayer;
	mcscoords[12][2] = 0;
	
	//13
	sprintf(panel_name,"panel%d",13);
	panel[13]= geom->MakeBox(panel_name, Vacuum, ypanel-halfHoleLength,xpanel,zpanel);
	ypos = -westylayer+(2*1+1)*ypanel;  
	layerw1->AddNode(panel[13],13, new TGeoTranslation(ypos,0,0));
	mcscoords[13][0] = ypos;
	mcscoords[13][1] = westxlayer;
	mcscoords[13][2] = 0;
	
	//22
	sprintf(panel_name,"panel%d",22);
	panel[22]= geom->MakeBox(panel_name, Vacuum, ypanel-halfHoleLength,xpanel,zpanel);
	ypos = -westylayer+(2*0+1)*ypanel;  
	layerw2->AddNode(panel[22],22, new TGeoTranslation(ypos,0,0));
	mcscoords[22][0] = ypos;
	mcscoords[22][1] = westxlayer;
	mcscoords[22][2] = 0;
	
	//14
	sprintf(panel_name,"panel%d",14);
	panel[14]= geom->MakeBox(panel_name, Vacuum, ypanel+halfHoleLength,xpanel,zpanel);
	ypos = -westylayer+(2*1+1)*ypanel;  
	layerw2->AddNode(panel[14],14, new TGeoTranslation(ypos,0,0));
	mcscoords[14][0] = ypos;
	mcscoords[14][1] = westxlayer;
	mcscoords[14][2] = 0;

	// add west layers to mother
	TGeoRotation *rot4 = new TGeoRotation();
	rot4->RotateX(90);   
	top->AddNode(layerw2,1, new TGeoCombiTrans(-halfHoleLength,westylayer-zlayer,0,rot4));
	top->AddNode(layerw1,2, new TGeoCombiTrans(halfHoleLength,westylayer-3*zlayer,0,rot4));


	///////////////////////// NEW ADDITIONS

	// EAST veto panels-------------------------------------------------------------
	// panel 28,30 EAST inner
	// panel 29,31 EAST outer

	// Z-axis points toward TOP
	// Y-axis points toward EAST
	// X-axis points toward NORTH
	
	Double_t eastxlayer = 203.29 / 2;
	Double_t eastylayer = 203.29 / 2;

	TGeoVolume *layerE1 = geom->MakeBox("layerE1", Vacuum, eastxlayer,eastylayer,zlayer);   
	TGeoVolume *layerE2 = geom->MakeBox("layerE2", Vacuum, eastxlayer,eastylayer,zlayer);   
	xpanel = eastxlayer;    
	ypanel = (eastylayer-2*zlayer)/2.; 

	//28
	sprintf(panel_name,"panel%d",28);
	panel[28]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfHoleLength,zpanel);
	ypos = -eastylayer+(2*0+1)*ypanel;  
	layerE1->AddNode(panel[28],28, new TGeoTranslation(0,ypos,0));
	mcscoords[28][0] = ypos;
	mcscoords[28][1] = -eastylayer;
	mcscoords[28][2] = 0;

	//30
	sprintf(panel_name,"panel%d",30);
	panel[30]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfHoleLength,zpanel);
	ypos = -eastylayer+(2*1+1)*ypanel;  
	layerE1->AddNode(panel[30],30, new TGeoTranslation(0,ypos,0));
	mcscoords[30][0] = ypos;
	mcscoords[30][1] = -eastylayer;
	mcscoords[30][2] = 0;

	//29
	sprintf(panel_name,"panel%d",29);
	panel[29]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfHoleLength,zpanel);
	ypos = -eastylayer+(2*0+1)*ypanel;  
	layerE2->AddNode(panel[29],29, new TGeoTranslation(0,ypos,0));
	mcscoords[29][0] = ypos;
	mcscoords[29][1] = -eastylayer;
	mcscoords[29][2] = 0;

	//31
	sprintf(panel_name,"panel%d",31);
	panel[31]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfHoleLength,zpanel);
	ypos = -eastylayer+(2*1+1)*ypanel;  
	layerE2->AddNode(panel[31],31, new TGeoTranslation(0,ypos,0));
	mcscoords[31][0] = ypos;
	mcscoords[31][1] = -eastylayer;
	mcscoords[31][2] = 0;
	
	//cryovat accessor hole
	TGeoVolume *accessor1 = geom->MakeBox("accessor1", Vacuum, halfHoleLength*2, zlayer*2, halfHoleLength*2);
	top->AddNode(accessor1, 0, new TGeoTranslation(-zlayer*2, -eastylayer+2*zlayer, 0));
	/*
	TGeoCompositeShape *cs1 = new TGeoCompositeShape("cs1","(panel[29] + panel[31])- accessor1");
	TGeoVolume *comp = new TGeoVolume("COMP",cs1);
	top->AddNode(comp,1);
	*/
	// add EAST layers to mother
	TGeoRotation *rot5 = new TGeoRotation();
	rot5->RotateX(90);
	rot5->RotateY(90);   
	top->AddNode(layerE1,1, new TGeoCombiTrans(halfHoleLength,-eastylayer+zlayer,0,rot5));
	top->AddNode(layerE2,2, new TGeoCombiTrans(-halfHoleLength,-eastylayer+3*zlayer,0,rot5));


	// SOUTH veto panels-------------------------------------------------------------
	// panel 24,26 SOUTH inner
	// panel 25,27 SOUTH outer
	
	Double_t southxlayer = 203.29 / 2;
	Double_t southylayer = 203.29 / 2;

	TGeoVolume *layerS1 = geom->MakeBox("layerS1", Vacuum, southxlayer,southylayer,zlayer);   
	TGeoVolume *layerS2 = geom->MakeBox("layerS2", Vacuum, southxlayer,southylayer,zlayer);   

	xpanel = southxlayer;    
	ypanel = (southylayer-4*zlayer)/2.; 

	//24
	sprintf(panel_name,"panel%d",24);
	panel[24]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfHoleLength,zpanel);
	ypos = -southylayer+(2*0+1)*ypanel;  
	layerS1->AddNode(panel[24],24, new TGeoTranslation(0,ypos,0));
	mcscoords[24][0] = -southxlayer;
	mcscoords[24][1] = ypos;
	mcscoords[24][2] = 0;

	//26
	sprintf(panel_name,"panel%d",26);
	panel[26]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfHoleLength,zpanel);
	ypos = -southylayer+(2*1+1)*ypanel;  
	layerS1->AddNode(panel[26],26, new TGeoTranslation(0,ypos,0));
	mcscoords[26][0] = -southxlayer;
	mcscoords[26][1] = ypos;
	mcscoords[26][2] = 0;

	//25
	sprintf(panel_name,"panel%d",25);
	panel[25]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfHoleLength,zpanel);
	ypos = -southylayer+(2*0+1)*ypanel;  
	layerS2->AddNode(panel[25],25, new TGeoTranslation(0,ypos,0));
	mcscoords[25][0] = -southxlayer;
	mcscoords[25][1] = ypos;
	mcscoords[25][2] = 0;

	//27
	sprintf(panel_name,"panel%d",27);
	panel[27]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfHoleLength,zpanel);
	ypos = -southylayer+(2*1+1)*ypanel;  
	layerS2->AddNode(panel[27],27, new TGeoTranslation(0,ypos,0));
	mcscoords[27][0] = -southxlayer;
	mcscoords[27][1] = ypos;
	mcscoords[27][2] = 0;
	
	//cryovat accessor hole
	TGeoVolume *accessor2 = geom->MakeBox("accessor1", Vacuum, zlayer*2, halfHoleLength*2, halfHoleLength*2);
	top->AddNode(accessor2, 0, new TGeoTranslation(-southxlayer+2*zlayer, 0, 0));

	// add SOUTH layers to mother
	TGeoRotation *rot6 = new TGeoRotation();
	rot6->RotateY(90);   
	top->AddNode(layerS1,1, new TGeoCombiTrans(-southxlayer+zlayer,+4*zlayer+halfHoleLength,0,rot6)); // need to move 2x the width of a panel in the -y direction
	top->AddNode(layerS2,2, new TGeoCombiTrans(-southxlayer+3*zlayer,+4*zlayer-halfHoleLength,0,rot6));

	// done adding panels! -------------------------------------------------------------
	
	// create the cryovats--------------------------------------------------------------
	Double_t cryoboxx = 50 / 2;
	Double_t cryoboxy = 90 / 2;
	Double_t cryoboxz = 60 / 2;
	
	TGeoVolume *cryobox = geom->MakeBox("layert1", Vacuum, cryoboxx,cryoboxy,cryoboxz);

	TGeoVolume *cryovat1 = geom->MakeEltu("cryovat", Vacuum, cryoboxy/2, cryoboxy/2, cryoboxz/2);
	cryobox->AddNode(cryovat1,1, new TGeoTranslation(0,cryoboxy/2,0));
	
	TGeoVolume *cryovat2 = geom->MakeEltu("cryovat", Vacuum, cryoboxy/2, cryoboxy/2, cryoboxz/2);
	cryobox->AddNode(cryovat2,2, new TGeoTranslation(0,-cryoboxy/2,0));
	
	
	//cryobox->SetVisDaughters(false);
	top->AddNode(cryobox,1,0);
	
	// draw the tracks--------------------------------------------------------------------
	// 13 = "mu-" in pdg
	TDatabasePDG *pdg = new TDatabasePDG();
	TParticlePDG *muon1 = pdg->GetParticle("mu-");
	// create a new track
	Int_t track_index = geom->AddTrack(1,muon1->PdgCode());
	// get track pointer
	TVirtualGeoTrack *track = gGeoManager->GetTrack(track_index);
	//name track to get line type for particle
	track->SetName(muon1->GetName());


	/*
	//just testing stuff...
	TGeoVolume *hitbox = geom->MakeBox("hitbox", Vacuum, zlayer*2,zlayer*2,zlayer*2);
	TGeoVolume *hole1 = geom->MakeTube("hole1", Vacuum, 10,9,0);
	TGeoVolume *hole2 = geom->MakeTube("hole1", Vacuum, 10,9,0);
	int testpanel = 0;
	int testpanel2 = 11;
	hitbox->AddNode(hole1,1, new TGeoTranslation(mcscoords[testpanel][0],mcscoords[testpanel][1],mcscoords[testpanel][2]));
	hitbox->AddNode(hole2,2, new TGeoTranslation(mcscoords[testpanel2][0],mcscoords[testpanel2][1],mcscoords[testpanel2][2]));
	top->AddNode(hitbox,1,0);
	*/
	
	vector<vector<Double_t>> allHitCoords;
	bool alreadyDone[32] = {false}; //used so the same layer hit isnt added twice
	int numOfPoints = 1;
	int numOfPlanesHit = 0;
	
	//these for loops check if layer partner hit exists, if it does, draws a hole for the layer hit
	for(int i = 0; i < 32; i++) {
		if(qdcVals[i] != 0) {
			for(int j = 0; j < 32; j++) {
				if(qdcVals[j] != 0 && isLayerHit(i, j) == 1 && alreadyDone[j] == false) {
					numOfPlanesHit++;
					vector<Double_t> drawLoc = hitLocation(i, j);
					allHitCoords.push_back(drawLoc);
					
					//cout << mcscoords[i][0] << " " << mcscoords[i][1] << " " << mcscoords[i][2] << endl;
					//cout << mcscoords[j][0] << " " << mcscoords[j][1] << " " << mcscoords[j][2] << endl;
					cout << drawLoc.at(0) << " " << drawLoc.at(1) << " " << drawLoc.at(2) << endl;
			
					TGeoVolume *hole1 = geom->MakeTube("hole1", Vacuum, 9,10,0);
					if(i > 11 && i != 17 && i != 18 && i != 20 && i != 21) { //if i is not top or bottom panels
						TGeoRotation *rot7 = new TGeoRotation();
						rot7->RotateY(90);
						if(i > 27 || i == 12 || i == 13 || i == 14 || i == 22) { //if i is east or west panels
							rot7->RotateY(-90);
							rot7->RotateX(90);
						}
						top->AddNode(hole1,1,new TGeoCombiTrans(drawLoc.at(0), drawLoc.at(1), drawLoc.at(2),rot7));
					}
					else {
						top->AddNode(hole1,1,new TGeoTranslation(drawLoc.at(0), drawLoc.at(1), drawLoc.at(2)));
					}
					alreadyDone[j] = true;
					
					track->AddPoint(drawLoc.at(0), drawLoc.at(1), drawLoc.at(2),i*50);
					numOfPoints++;
				}
			}
		}
		alreadyDone[i] = true;
	}


   //--- close the geometry
   geom->CloseGeometry();


	// draw the event
   geom->SetVisLevel(4);
   geom->SetVisOption(0);
   track->SetLineWidth(5.0);
   track->SetLineStyle(9);
   track->SetLineColor(3);
 
	
	// color particular panels based on qdc value
	for (Int_t k = 0; k<32; k++) {
		panel[k]->SetLineColor(coloring(qdcVals[k],3));
		if (qdcVals[k] > 0) 
			panel[k]->SetLineWidth(3.0);
		else 
			panel[k]->SetLineWidth(1.0);
	}
   
   ecan = new TCanvas("ecan","veto hits",0,0,500,500);
   ecan->cd();
   top->Draw();
   
   
   //draw the tracks for any two sides
	int topSide = 0;
	int bottomSide = 0;
	int northSide = 0;
	int eastSide = 0;
	int southSide = 0;
	int westSide = 0;
   
	for(int i = 0; i < 32; i ++) {
		for(int j = 0; j < 12; j++) {
			if(qdcVals[i] != 0 && bottomPanels[j] == i)
				bottomSide++;
		}
		for(int j = 0; j < 4; j++) {
			if(qdcVals[i] != 0 && topPanels[j] == i)
				topSide++;
			if(qdcVals[i] != 0 && northPanels[j] == i)
				northSide++;
			if(qdcVals[i] != 0 && eastPanels[j] == i)
				eastSide++;
			if(qdcVals[i] != 0 && southPanels[j] == i)
				southSide++;
			if(qdcVals[i] != 0 && westPanels[j] == i)
				westSide++;
		}
   }
   

	if(numOfPlanesHit <= 2 && numberOfPanelsHit == 4 && (bottomSide == 2 || topSide == 2
		|| northSide == 2 || eastSide == 2 || southSide == 2 || westSide == 2)) {
		geom->DrawTracks();
	}
	
   //top->Draw("ogl");
	// add some geometry markers
	TPaveText *pt = new TPaveText(-.9,.9,-.4,.5,"br"); //-.9,-.9,-.4,-.5 for bottom right x1,y1,x2,y2
	//pt->AddText("This is the South Side");
	
	char pavtot[150];
    sprintf(pavtot,"QDC total: %d", totalQDC);
    pt->AddText(pavtot);
        
    char runnum[150];
	sprintf(runnum,"Run number: %d",runNumber);
	pt->AddText(runnum);
		
	char evntcount[150];
	sprintf(evntcount,"Event count: %d",eventCount);
	pt->AddText(evntcount);
	
	int panel1, panel2;
	panel1 = panel2 = 0;
	
	if(numberOfPanelsHit == 2) {	
		//looks for panel hit from the front
		for(int i = 0; i < 32; i++) {
			if(qdcVals[i] != 0) {
				panel1 = i;
				break;
			}
		}
		//continues from where panel1 left off looking for next panel hit
		for(int i = panel1+1; i < 32; i++) {
			if(qdcVals[i] != 0) {
				panel2 = i;
				break;
			}
		}
		
		char nextto[150];
        sprintf(nextto,"Is Next To: %d", isNextTo(panel1, panel2));
		pt->AddText(nextto);
		
		//displays panel number in 1-32 format
		char test1[150];
		sprintf(test1,"panel1: %d",panel1+1);
		pt->AddText(test1);
		
		char test2[150];
		sprintf(test2,"panel2: %d",panel2+1);
		pt->AddText(test2);
	}
	
	//coordinates of each hit location
	for(int i = 0; i < allHitCoords.size(); i++) {
		char coordout[150];
		sprintf(coordout,"Hit Coords: (%.1f,%.1f,%.1f)",allHitCoords[i][0],allHitCoords[i][1],allHitCoords[i][2]);
		pt->AddText(coordout);
	}
	
	//set phi and theta to 361 each time to know if set
	phi = 361;
	theta = 361;
	//some vector info
	if(allHitCoords.size() > 1 && numberOfPanelsHit == 4) {
		TVector3 r1;
		TVector3 r2;
		TVector3 t1;
		r1.SetXYZ(allHitCoords[0][0],allHitCoords[0][1],allHitCoords[0][2]);
		r2.SetXYZ(allHitCoords[1][0],allHitCoords[1][1],allHitCoords[1][2]);
		t1 = r2-r1;
	
		phi = t1.Phi()*(180/pi); //180/pi to get deg rather than rad
		theta = t1.Theta()*(180/pi);
	
		char phiout[150];
		sprintf(phiout,"Phi: %.1f deg",phi); 
		pt->AddText(phiout);
	
		char thetaout[150];
		sprintf(thetaout,"Theta: %.1f deg",theta);
		pt->AddText(thetaout);
	}
	
	pt->Draw();
	
	ecan->Update();

}
//--------------------------------------------------------------

//Matrix to determine which panels are ajacent.
//Input panel number based on the 0-31 scale (panel 1 is panel 0) 
//will return 1 if panels are ajacent.
int isNextTo(int panel1, int panel2) {

	int panelTable[32][32] = {
	{0,1,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,0,1,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,1,0,1,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,1,0,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,1,0,1,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,1,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,1,1,1,1,1,0,1,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,1},
	{1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1},
	{1,1,1,1,1,1,0,1,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1},
	{1,1,1,1,1,1,0,0,1,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0},
	{1,1,1,1,1,1,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,1,0,0},
	{1,1,1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,1,1,0,0},
	{0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,0,0,0,0,0,1,0,1,0,1,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1,0,0,0,1,0,1,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1,0,0,0,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,1,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0},
	{0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,1,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,1,0,0,1,1,1},
	{0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1},
	{0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,1,1,0,1},
	{0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0}
	};

	return panelTable[panel1][panel2];
}

//input two panels on 0-31 scale, returns 1 if same side(both top or both east ect..)
//and diffrent layer(inside & outside)(layer partner), else return 0
int isLayerHit(int panel1, int panel2) {

	int layerTable[32][32] = {
	{0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0}
	};

	return layerTable[panel1][panel2];
}


//Input two panels on 0-31 scale and returns vector with 
//location in mother coordinate system of the hit if on same side and 
//diffrent layer(layer partner), else displays error and returns empty vector.
std::vector<Double_t> hitLocation(int panel1, int panel2) {
	std::vector<Double_t> drawCoords(3);
	
	if(isLayerHit(panel1, panel2) == 1) {
		if(panel1 < panel2) { //if panel1 is the upper panel
			drawCoords.at(0) = mcscoords[panel2][0];
			drawCoords.at(1) = mcscoords[panel1][1];
			drawCoords.at(2) = mcscoords[panel2][2];
			//cout << drawCoords.at(0) << " " << drawCoords.at(1) << " " << drawCoords.at(2) << endl;
		}
		else {
			drawCoords.at(0) = mcscoords[panel1][0];
			drawCoords.at(1) = mcscoords[panel2][1];
			drawCoords.at(2) = mcscoords[panel2][2];
		}
	}
	else {
		cout << "there is an error with hit location" << endl;
	}
		
	return drawCoords;
}

void hitFinder(char* innum1, char* innum2) {
	
	TApplication *myapp=new TApplication("myapp",0,0);
	string input = "";
	int inputnum1 = 0;
	int inputnum2 = 0;
	Int_t runNumber = 0;
	Int_t entry = 0;
	Int_t eventCount = 0;
	Double_t scalerTime = 0.0;
	Int_t qdcVals[32]={0};
	string line = "";
	int lineLength = 0;
	bool found = false;
	bool validNum1 = false;
	bool validNum2 = false;
	
	ifstream VetoHitsFile;
    VetoHitsFile.open(inputFile);
    
	stringstream myStream(innum2);
	if (myStream >> inputnum1) {
		inputnum1 = atoi(innum1);
		validNum1 = true;
	}
	
	if(validNum1 == false) {
		cout << " is not a valid number, please try again" << endl;
		exit(0);
	}
	
	if(validNum1 == true)	 {
		while(getline(VetoHitsFile, line)) {
			lineLength = 0;
			
			stringstream  lineStream(line);
			lineStream >> runNumber >> entry >> eventCount >> scalerTime;

			while(lineStream >> qdcVals[lineLength]) {
				lineLength++;
			}
			if(lineLength != 31) {
				for(int i = lineLength; i < 32; i++) {
					qdcVals[i] = 0;
				}
			}
			
			if(runNumber == inputnum1) {
				found = true;
				break;
			}
		}
		if(VetoHitsFile.eof() && found == false) {
				cout << "Run number not found, please try again" << endl;
				exit(0);
		}
	}
	found = false;
	lineLength = 0;

	stringstream myStream2(innum2);
	if (myStream2 >> inputnum2) {
		inputnum2 = atoi(innum2);
		validNum2 = true;
	}
	if(validNum2 == false) {
		cout << innum2 << " is not a valid number, please try again" << endl;
		exit(0);
	}
	
	if(validNum2 == true && validNum1 == true)	 {
		//go back to the beginning of the file
		VetoHitsFile.clear();
		VetoHitsFile.seekg(0, ios::beg);
		while(getline(VetoHitsFile, line)) {
			lineLength = 0;
			
			stringstream  lineStream(line);
			lineStream >> runNumber >> entry >> eventCount >> scalerTime;
			
			while(lineStream >> qdcVals[lineLength]) {
				lineLength++;
			}
			if(lineLength != 31) {
				for(int i = lineLength; i < 32; i++) {
					qdcVals[i] = 0;
				}
			}
			if(eventCount == inputnum2 && runNumber == inputnum1) {
				DrawEvent(qdcVals, 0, 0, runNumber, eventCount);
				found = true;
				break;
			}
		}
		if(VetoHitsFile.eof() && found == false) {
				cout << "Event count not found, please try again" << endl;
				exit(0);
		}
	}
	
	myapp->Run();
}
 
//--------------------------------------------------------------

//global canvas and graph declarations
TCanvas *graphs = new TCanvas("graphs","Analysis graphs", 900, 0, 900, 600);
TCanvas *QDCcanvas = new TCanvas("QDC","QDC", 900, 0, 900, 600);
TCanvas *projcan = new TCanvas("Projections","Panels hit vs. Multiplicity", 900, 0, 700, 700);
TCanvas *qdctotalcan = new TCanvas("qdctotalcan", "QDC Totals", 900, 0, 1401, 400);
TCanvas *qdctotaln = new TCanvas("qdctotaln", "QDC totals for each n", 1000, 500);
TCanvas *thetaPhican = new TCanvas("thetaPhican","Theta and Phi of particle's tracks", 900, 0, 900, 600);
TCanvas *costhetaPhican = new TCanvas("costhetaPhican","Cos(Theta) and Phi of particle's tracks", 900, 0, 900, 600);

TH1F *graph1 = new TH1F("graph1","Number of panels hit per event", 65, 0, 33);
TH1F *graph2 = new TH1F("graph2","Times each panel hit", 65, 0, 33);
TH2F *graph3 = new TH2F("graph3","QDC per panel for all hits", 65, 0, 4400, 65, 0, 33);
TH2F *graph4 = new TH2F("graph4","Multiplicity vs. QDC total", 100, 0, 32, 100, 0, 20000);
TH2F *graph5 = new TH2F("graph5","QDC per panel for 2 panels hit", 65, 0, 4400, 65, 0, 33);
TH2F *graph6 = new TH2F("graph6","QDC per panel for 4 panels hit", 65, 0, 4400, 65, 0, 33);
TH2F *graph7 = new TH2F("graph7","Panels hit vs. Multiplicity", 65, 0, 33, 65, 0, 33);
TH1F *graph8 = new TH1F("graph8","QDC total for all hits", 100, 0, 18000);
TH1F *graphn2 = new TH1F("graphn2","QDC total for 2 panels hit", 100, 0, 18000);
TH1F *graphn3 = new TH1F("graphn3","QDC total for 3 panels hit", 100, 0, 18000);
TH1F *graphn4 = new TH1F("graphn4","QDC total for 4 panels hit", 100, 0, 18000);
TH1F *graphn5 = new TH1F("graphn5","QDC total for 5 panels hit", 100, 0, 18000);
TH1F *graphn6 = new TH1F("graphn6","QDC total for 6 panels hit", 100, 0, 18000);
TH1F *graphn7 = new TH1F("graphn7","QDC total for 7 panels hit", 100, 0, 18000);
TH1F *graphn8 = new TH1F("graphn8","QDC total for 8 panels hit", 100, 0, 18000);
TH1F *graphn9 = new TH1F("graphn9","QDC total for 9 panels hit", 100, 0, 18000);
TH2F *thetaPhi = new TH2F("thetaPhi","Theta vs. Phi for top and bottom n=4 layer hits", 90, -185, 185, 90, -1, 50);
TH2F *costhetaPhi = new TH2F("costhetaPhi","Cos(Theta) vs. Phi for top and bottom n=4 layer hits", 90, -185, 185, 90, -1.3, 1.3);
	
//fills all plots and prints out wireframe pdfs
void fillPlots(Int_t qdcvals[], Int_t totalQDC, Int_t numberOfPanelsHit, Int_t ievent) {
	
	graph1->Fill(numberOfPanelsHit);
	graph4->Fill(numberOfPanelsHit, totalQDC);
	graph8->Fill(totalQDC);	

	for (Int_t i = 0; i < 32; i++) {
		if(qdcvals[i] != 0) {
			graph2->Fill(i+1);
			graph3->Fill(qdcvals[i], i+1);
			graph7->Fill(numberOfPanelsHit, i+1);
			if(numberOfPanelsHit == 2) {
				graph5->Fill(qdcvals[i], i+1);
			}
			if(numberOfPanelsHit == 4) {
				graph6->Fill(qdcvals[i], i+1);
			}
		}
	}
	//fill qdctotaln canvas and print assembly pdfs
	char pcanname[150]; 
	for(int g = 2; g < 9; g++) {
        if(numberOfPanelsHit == g) {
			switch(g) {
				case 2:
					sprintf(pcanname,"output/2/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn2->Fill(totalQDC);
					break;
                case 3:
					sprintf(pcanname,"output/3/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn3->Fill(totalQDC);
					break;
                case 4:
					sprintf(pcanname,"output/4/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn4->Fill(totalQDC);
					break;
                case 5:
					sprintf(pcanname,"output/5/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn5->Fill(totalQDC);
					break;
                case 6:
					sprintf(pcanname,"output/6/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn6->Fill(totalQDC);
					break;
                case 7:
					sprintf(pcanname,"output/7/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn7->Fill(totalQDC);
					break;
                case 8:
					sprintf(pcanname,"output/8/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn8->Fill(totalQDC);
					break;
                case 9:
					sprintf(pcanname,"output/9/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					graphn9->Fill(totalQDC);
					break;
				default:
					sprintf(pcanname,"output/other/assembly%d.pdf",ievent);
					ecan->Print(pcanname,"pdf");
					break;
			}
		}
	}
	//all above 9 panels hit go to 'other' subdir
	if(numberOfPanelsHit > 9) {
		sprintf(pcanname,"output/other/assembly%d.pdf",ievent);
		ecan->Print(pcanname,"pdf");
	}
	
	bool filled = false;
	
	for(int i = 0; i < 12; i++) {
		for(int j = 0; j < 4; j++) {
			if(phi != 361 && theta != 361 && qdcvals[bottomPanels[i]] != 0 && 
			   qdcvals[topPanels[j]] != 0 && numberOfPanelsHit == 4 && filled == false) {
	
				Double_t costheta = cos(theta);	
				
				thetaPhi->Fill(phi, theta);
				costhetaPhi->Fill(phi, costheta);
				filled = true;
				cout << "filled" << endl;
			}
		}
	}
	
	ecan->cd();
}

void drawPlots() {
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
    graph3->SetYTitle("Panel number");
	graph3->SetXTitle("QDC");
	graph3->SetMarkerStyle(kFullDotSmall);
	//gStyle->SetPalette(kBird);
	graph3->SetStats(0);
    graph3->Draw();
	
	graphs->cd(4);
	graph4->SetXTitle("Multiplicity");
	graph4->SetYTitle("Total QDC");
	graph4->SetMarkerStyle(kFullDotSmall);
	//graph4->GetYaxis()->SetTitleOffset(1.7);
	graph4->SetStats(0);
	graph4->Draw();
    
    graphs->cd(5);
    graph5->SetYTitle("Panel number");
    graph5->SetXTitle("QDC");
    graph5->SetMarkerStyle(kFullDotSmall);
    graph5->SetStats(0);
    graph5->Draw();
    
    graphs->cd(6);
    graph6->SetYTitle("Panel number");
    graph6->SetXTitle("QDC");
    graph6->SetMarkerStyle(kFullDotSmall);
    graph6->SetStats(0);
    graph6->Draw();
	
	//Drawing qdccanvas canvas
	QDCcanvas->Divide(3,2);
	
	QDCcanvas->cd(1);
	graph3->Draw();
	
	TH1D *qdcpanelproj = graph3->ProjectionX();
	QDCcanvas->cd(4);
	qdcpanelproj->SetTitle("QDC distribution for all hits");
	qdcpanelproj->SetYTitle("Count");
	qdcpanelproj->GetYaxis()->SetTitleOffset(1.45);
	qdcpanelproj->Draw("bar");	
	
	QDCcanvas->cd(2);
	graph5->Draw();
	
	TH1D *qdcpanelprojn2 = graph5->ProjectionX();
	QDCcanvas->cd(5);
	qdcpanelprojn2->SetTitle("QDC distribution for 2 panels hit");
	qdcpanelprojn2->SetYTitle("Count");
	qdcpanelprojn2->GetYaxis()->SetTitleOffset(1.45);
	qdcpanelprojn2->Draw("bar");
	
	QDCcanvas->cd(3);
	graph6->Draw();
	
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
	graph7->Draw();
	
	TH1D *countingPanelProj = graph7->ProjectionY();
	projcan->cd(2);
	countingPanelProj->SetTitle("Times each panel hit");
	countingPanelProj->SetYTitle("Count");
	countingPanelProj->Draw("hbar");
	
	TH1D *countingMultProj = graph7->ProjectionX();
	projcan->cd(3);
	countingMultProj->SetTitle("Multiplicity");
	countingMultProj->SetYTitle("Count");
	//countingMultProj->GetYaxis()->SetTitleOffset(1.7);
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
	thetaPhi->SetXTitle("Degrees(Phi)");
	thetaPhi->SetYTitle("Degrees(Theta)");
	thetaPhi->Draw();
	
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
	costhetaPhi->SetXTitle("Degrees(Phi)");
	costhetaPhi->SetYTitle("Cos(Theta)");
	costhetaPhi->Draw();
	
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
}

void printPlots() {
	char graphscan[150];
		sprintf(graphscan,"output/plots/graphs.pdf");
		graphs->Print(graphscan,"pdf");
		
	char qdccan[150];
		sprintf(qdccan,"output/plots/qdcDisplay.pdf");
		QDCcanvas->Print(qdccan,"pdf");
	
	char projcanvas[150];
		sprintf(projcanvas,"output/plots/Projections.pdf");
		projcan->Print(projcanvas,"pdf");
	
	char qdctotcan[150];
		sprintf(qdctotcan,"output/plots/QDCtotal.pdf");
		qdctotalcan->Print(qdctotcan,"pdf");
		
	char qdctotncan[150];
        sprintf(qdctotncan,"output/plots/QDCtotaln.pdf");
        qdctotaln->Print(qdctotncan,"pdf");
        
    char thetaPhicanprint[150];
        sprintf(thetaPhicanprint,"output/plots/thetaPhi.pdf");
        thetaPhican->Print(thetaPhicanprint,"pdf");
        
    char costhetaPhicanprint[150];
        sprintf(costhetaPhicanprint,"output/plots/cosThetaPhi.pdf");
        costhetaPhican->Print(costhetaPhicanprint,"pdf");

}
//--------------------------------------------------------------

void VetoDisplay()
{

	Int_t ievent = 0;

	Int_t runNumber = 0;
	Int_t entry = 0;
	Int_t eventCount = 0;
	Double_t scalerTime = 0.0;
	Int_t qdcVals[32] = {0};
	Int_t numberOfPanelsHit = 0;
	Int_t totalQDC = 0;
	Int_t qdcValSum[32]={0};
	Double_t maxvalsum = 0.0;
	Double_t minvalsum = 0.0;
	Double_t scaleFactor = 0.0;
	string line = "";
	string cell = "";
	Int_t lineLength = 0;
	bool goodPanelNum = false;
	
	//read in vetohit array from file
    ifstream VetoHitsFile;
    VetoHitsFile.open(inputFile);
    //TODO: add fileNotFound error, exit
    while(getline(VetoHitsFile, line)) {
		ievent++;
		//getchar();
		// reset counters and read first part of next line
		numberOfPanelsHit = 0;
		totalQDC = 0;
		lineLength = 0;
		goodPanelNum = false;
		
		stringstream  lineStream(line);
		lineStream >> runNumber >> entry >> eventCount >> scalerTime;
		printf("runNumber:%i  entry:%i  eventCount:%i  scalerTime:%.5f \n",runNumber,entry,eventCount,scalerTime);

		
		while (lineStream.peek() == ' ') { //gets rid of trailing space before qdcvals
			lineStream.get();
		}

		while(getline(lineStream,cell,' ')) {
			try {
			qdcVals[lineLength] = stoi(cell); //parses input values
			}
			catch(exception e) {
				cout << "wrong number of panels. " << lineLength << " given." << endl;
				break;
			}
			lineLength++;
			if(lineLength == 32) { //stops trailing spaces/values after 32 qdcvals
				break;
			}
		}
		/* old way
		while(lineStream >> qdcVals[lineLength]) {
			lineLength++;
		}
		*/
		
		if(lineLength == 32 || lineLength == 24) {
			goodPanelNum = true;
		}
		if(goodPanelNum == true) {
			//if linelength is 24, fill in any missing detectors with 0 qdc
			if(lineLength != 32) { //lineLength 
				for(int i = lineLength; i < 32; i++) {
					qdcVals[i] = 0;
				}
			}
		
			for(int i = 0; i < 32; i++) {
				if(qdcVals[i] != 0) {
					numberOfPanelsHit++;
					totalQDC += qdcVals[i];
					qdcValSum[i] += qdcVals[i];
				}
			}
			/*
			for(int i = 0; i < 32; i++) {
				cout << qdcVals[i] << " ";
			}
			cout << endl;
			*/
			DrawEvent(qdcVals, numberOfPanelsHit, totalQDC, runNumber, eventCount);
			fillPlots(qdcVals, totalQDC, numberOfPanelsHit, ievent);
	}
	}
	
	
	// send qdcValSum as last picture
	minvalsum = qdcValSum[0]; //min equal to first value in sums to allow for comparison
	for(Int_t i = 0; i < 32; i++) {
		if(qdcValSum[i] > maxvalsum) {
			maxvalsum = qdcValSum[i];
		}
		if(qdcValSum[i] < minvalsum) {
			minvalsum = qdcValSum[i];
		}
	}
	//we multiply each sum by a scale factor to get the adjusted sums
	scaleFactor = (4000-500) / (maxvalsum - minvalsum);
	for (Int_t j=0;j<32;j++){
		qdcValSum[j] = qdcValSum[j] * scaleFactor;
		//cout << qdcValSum[j] << endl;
	}
	//a bit buggy as there is no runnumber or eventcount for all sums to draw on pave text
	DrawEvent(qdcValSum, totalQDC, numberOfPanelsHit, runNumber, eventCount); 
	
   	char pcanname[150]; 
   	sprintf(pcanname,"output/other/assembly_sum.pdf");
   	ecan->Print(pcanname,"pdf");

	//lastly, draw and print the plots	
	drawPlots();
	printPlots();

}

int main(int argc, char* argv[]) {
	if(argc == 1) {
		VetoDisplay();
		return 0;
	}
	if(argc == 3) {
		hitFinder(argv[1],argv[2]);
	}
	else {
		cout << "Usage: ./VetoDisplay for a full run of the input file, or "
				<< "./VetoDisplay RunNumber EventCount to see an individual event model" << endl;
	}
}
