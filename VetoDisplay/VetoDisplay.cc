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

//method for coloring the wireframe based on recived QDC per panel
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
int topSide = 0;
int bottomSide = 0;
int northSide = 0;
int eastSide = 0;
int southSide = 0;
int westSide = 0;
int numOfPlanesHit = 0;

//mother coordinate system for each panel in x,y,z, essentially just the center location of each panel
Double_t mcscoords[32][3];
// canvas is global
TCanvas *ecan;
//Phi & Theta global to allow for plotting
Double_t phi;
Double_t theta;
static Double_t pi = 3.1415926535897;
//---------------------------------------------------------------------------------------------------------------------
void DrawEvent(Int_t qdcVals[], Int_t numberOfPanelsHit, Int_t totalQDC, Int_t runNumber, Int_t eventCount)
{
//--- Definition of a simple geometry
   //gSystem->Load("libGeom");
   TGeoManager *geom = new TGeoManager("Assemblies","Geometry using assemblies");
   
   //--- define some materials
   TGeoMaterial *matVacuum = new TGeoMaterial("Vacuum", 0,0,0);
   //--- define some tracking media
   TGeoMedium *Vacuum = new TGeoMedium("Vacuum",1, matVacuum);

   //--- create the starting point of the volume structure - 
   // half length units. All are in cm.
   Double_t xtop=160;
   Double_t ytop=160;
   Double_t ztop=160;
   TGeoVolume *top = geom->MakeBox("TOP", Vacuum, xtop,ytop,ztop); //gets placed within this box
   geom->SetTopVolume(top);
   geom->SetTopVisible(0);
   
	// declare panels
	TGeoVolume *panel[32];
	
	//east side lengths
	Double_t eastxlayerouter = 172.72 / 2;
	Double_t eastxlayerinner = 176.53 / 2;
	Double_t eastylayer = 168.91 / 2;
	Double_t eastHalfOverlapOuter = 15.39 / 2; //half overlap = half the distance that panels on one side would overlap(inside and outside)
	Double_t eastHalfOverlapInner = 10.01 / 2; //i.e. half the length of the cryostat accessor holes on east and south sides
	
	//south side lengths
	Double_t southxlayerouter = 172.72 / 2;
	Double_t southxlayerinner = 176.53 / 2;
	Double_t southylayer = 212.09 / 2;
	Double_t southHalfOverlapOuter = 24.28 / 2;
	Double_t southHalfOverlapInner = 1.12 / 2;
	
	//north side lengths
	Double_t northxlayerouter = 165.1 / 2;
	Double_t northxlayerinner = 170.82 / 2;
	Double_t northylayer = 203.86 / 2;
	Double_t halfAmountOfOverlap = 5.08 / 2;  //used for both north and west side overlaps
	
	//west side lengths
	Double_t westxlayerouter = 165.10 / 2;
	Double_t westxlayerinner = 170.82 / 2;
	Double_t westylayer = 168.91 / 2;
	
	//top lengths
	Double_t topxlayerupper = 168.91 / 2;
	Double_t topylayerupper = 205.74 / 2;
	Double_t topxlayerlower = 212.09 / 2;
	Double_t topylayerlower = 168.91 / 2;
	
	//panel widths
	Double_t zlayer = 2.54 / 2;
	Double_t xpanel = 0.0;
	Double_t ypanel = 0.0;    
	Double_t zpanel = 0.0;
	
	// bottom veto panels-------------------------------------------------------------
	// make box for each layer and fill with 6 panels
	// panel 0-5  lower bottom
	// panel 6-11 upper bottom
	Double_t bottomDetectorWidth = 32.0; //total width of a bottom panel detector
	Double_t totalBottomHeight = 43.18; //upper + lower bottom panel height together
	Double_t spaceBetweenDetectorsy = 2.54; //from detector edge to next detector edge in the y axis
	Double_t detectorHeight = 2.54; //from detector edge to next detector edge in the z axis
	
	Double_t xlayer1 = 182.0 / 2; //length of lower bottom detectors in cm
	Double_t xlayer2 = 223.0 / 2; //length of upper bottom detectors in cm
	Double_t ylayer1 = ((bottomDetectorWidth * 6) + (spaceBetweenDetectorsy * 6)) / 2;
	Double_t zplace = (eastxlayerouter); //where to place the bottom panels on the z axis(hieght of upper panels)
	Double_t upperoffsetx = (southylayer - xlayer2 - spaceBetweenDetectorsy/2) - 1.93; //the offset of the upper bottom panels in x axis
	Double_t loweroffsetx = (eastylayer-2*zlayer - xlayer1 - spaceBetweenDetectorsy/2) + 13.63; //the offset of the lower bottom panels in x axis
	Double_t upperoffsety = (eastylayer-2*zlayer - ylayer1) + 13.96; //the offset of the upper bottom panels in y axis
	Double_t loweroffsety = (southylayer - ylayer1) - 1.27; //the offset of the lower bottom panels in the y axis
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

		sprintf(panel_name,"panel%d",6+i);
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

	// put 6 panels in layer 2
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

		sprintf(panel_name,"panel%d",6+i);
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
	
	TGeoVolume *layert1 = geom->MakeBox("layert1", Vacuum, topxlayerupper,topylayerupper,zlayer);   
	TGeoVolume *layert2 = geom->MakeBox("layert2", Vacuum, topxlayerlower,topylayerlower,zlayer);   

	for (Int_t i=0; i<2;i++){
		xpanel = topxlayerupper;
		ypanel = topylayerupper/2.;    
		zpanel = zlayer;
	
	    sprintf(panel_name,"panel%d",i);
	    panel[17+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	    ypos = -topylayerupper+(2*i+1)*ypanel;  
	    layert1->AddNode(panel[17+i],i, new TGeoTranslation(0,ypos,0));
	    mcscoords[17+i][0] = 0;
	    mcscoords[17+i][1] = ypos;
	    mcscoords[17+i][2] = northxlayerinner+3*zlayer -(westxlayerouter-westxlayerinner) -(eastxlayerouter-westxlayerouter);

		
		xpanel = topxlayerlower;
		ypanel = topylayerlower/2.;
		zpanel = zlayer;

	    sprintf(panel_name,"panel%d",i);
	    panel[20+i]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel,zpanel);
	    ypos = -topylayerlower+(2*i+1)*ypanel;  
	    layert2->AddNode(panel[20+i],i, new TGeoTranslation(0,ypos,0));
	    mcscoords[20+i][0] = ypos;
	    mcscoords[20+i][1] = 0;
	    mcscoords[20+i][2] = northxlayerinner+3*zlayer -(westxlayerouter-westxlayerinner) -(eastxlayerouter-westxlayerouter);
  	}
  	// add top layers to mother
  	top->AddNode(layert1,1, new TGeoTranslation(2*zlayer,0,northxlayerinner+3*zlayer -(westxlayerouter-westxlayerinner) -(eastxlayerouter-westxlayerouter)));
  	TGeoRotation *rot2 = new TGeoRotation();
  	rot2->RotateZ(-90);  
  	top->AddNode(layert2,2, new TGeoCombiTrans(2*zlayer,0,northxlayerinner+zlayer -(westxlayerouter-westxlayerinner) -(eastxlayerouter-westxlayerouter),rot2));

	// EAST veto panels-------------------------------------------------------------
	// panel 28,30 EAST inner
	// panel 29,31 EAST outer

	// Z-axis points toward TOP
	// Y-axis points toward EAST
	// X-axis points toward NORTH
	// when one panel would increase in size, it would have to be translated a direction by half
	// of the size increase for the bottom/sides to persist in the same location, thus the wonky translations

	TGeoVolume *layerE1 = geom->MakeBox("layerE1", Vacuum, eastxlayerouter,eastylayer,zlayer);   
	TGeoVolume *layerE2 = geom->MakeBox("layerE2", Vacuum, eastxlayerinner,eastylayer,zlayer);   
	xpanel = eastxlayerouter;
	ypanel = (eastylayer)/2.; 

	//29
	sprintf(panel_name,"panel%d",29);
	panel[29]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-eastHalfOverlapOuter,zpanel);
	ypos = -eastylayer+(2*0+1)*ypanel;  
	layerE1->AddNode(panel[29],29, new TGeoTranslation(0,ypos,0));
	mcscoords[29][0] = ypos;
	mcscoords[29][1] = -southylayer-zlayer;
	mcscoords[29][2] = 0;

	//31
	sprintf(panel_name,"panel%d",31);
	panel[31]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+eastHalfOverlapOuter,zpanel);
	ypos = -eastylayer+(2*1+1)*ypanel;  
	layerE1->AddNode(panel[31],31, new TGeoTranslation(0,ypos,0));
	mcscoords[31][0] = ypos;
	mcscoords[31][1] = -southylayer-zlayer;
	mcscoords[31][2] = 0;
	
	xpanel = eastxlayerinner;
	
	//28
	sprintf(panel_name,"panel%d",28);
	panel[28]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+eastHalfOverlapInner,zpanel);
	ypos = -eastylayer+(2*0+1)*ypanel;  
	layerE2->AddNode(panel[28],28, new TGeoTranslation(0,ypos,0));
	mcscoords[28][0] = ypos;
	mcscoords[28][1] = -southylayer-zlayer;
	mcscoords[28][2] = 0;

	//30
	sprintf(panel_name,"panel%d",30);
	panel[30]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-eastHalfOverlapInner,zpanel);
	ypos = -eastylayer+(2*1+1)*ypanel;  
	layerE2->AddNode(panel[30],30, new TGeoTranslation(0,ypos,0));
	mcscoords[30][0] = ypos;
	mcscoords[30][1] = -southylayer-zlayer;
	mcscoords[30][2] = 0;
	
	//cryovat accessor hole
	TGeoVolume *accessor1 = geom->MakeBox("accessor1", Vacuum, 25.40 / 2, zlayer*2, 21.59 / 2);
	top->AddNode(accessor1, 0, new TGeoTranslation(0, -southylayer-2*zlayer, 0));

	// add EAST layers to mother
	TGeoRotation *rot5 = new TGeoRotation();
	rot5->RotateX(90);
	rot5->RotateY(90);   
	top->AddNode(layerE1,1, new TGeoCombiTrans(-eastHalfOverlapInner,-southylayer-3*zlayer,0,rot5));
	top->AddNode(layerE2,2, new TGeoCombiTrans(eastHalfOverlapOuter,-southylayer-zlayer,-(eastxlayerouter-eastxlayerinner),rot5));


	// SOUTH veto panels-------------------------------------------------------------
	// panel 24,26 SOUTH inner
	// panel 25,27 SOUTH outer

	TGeoVolume *layerS1 = geom->MakeBox("layerS1", Vacuum, southxlayerouter,southylayer,zlayer);   
	TGeoVolume *layerS2 = geom->MakeBox("layerS2", Vacuum, southxlayerinner,southylayer,zlayer);   

	xpanel = southxlayerouter;    
	ypanel = (southylayer)/2.; 

	//25
	sprintf(panel_name,"panel%d",25);
	panel[25]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+southHalfOverlapOuter,zpanel);
	ypos = -southylayer+(2*0+1)*ypanel;  
	layerS1->AddNode(panel[25],25, new TGeoTranslation(0,ypos,0));
	mcscoords[25][0] = -southxlayerouter;
	mcscoords[25][1] = ypos;
	mcscoords[25][2] = 0;

	//27
	sprintf(panel_name,"panel%d",27);
	panel[27]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-southHalfOverlapOuter,zpanel);
	ypos = -southylayer+(2*1+1)*ypanel;  
	layerS1->AddNode(panel[27],27, new TGeoTranslation(0,ypos,0));
	mcscoords[27][0] = -southxlayerouter;
	mcscoords[27][1] = ypos;
	mcscoords[27][2] = 0;
	
	xpanel = southxlayerinner;
	
	//24
	sprintf(panel_name,"panel%d",24);
	panel[24]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-southHalfOverlapInner,zpanel);
	ypos = -southylayer+(2*0+1)*ypanel;  
	layerS2->AddNode(panel[24],24, new TGeoTranslation(0,ypos,0));
	mcscoords[24][0] = -southxlayerinner;
	mcscoords[24][1] = ypos;
	mcscoords[24][2] = 0;

	//26
	sprintf(panel_name,"panel%d",26);
	panel[26]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+southHalfOverlapInner,zpanel);
	ypos = -southylayer+(2*1+1)*ypanel;  
	layerS2->AddNode(panel[26],26, new TGeoTranslation(0,ypos,0));
	mcscoords[26][0] = -southxlayerinner;
	mcscoords[26][1] = ypos;
	mcscoords[26][2] = 0;
	
	//cryovat accessor hole
	TGeoVolume *accessor2 = geom->MakeBox("accessor1", Vacuum, zlayer*2, 25.40 / 2,  21.59 / 2);
	top->AddNode(accessor2, 0, new TGeoTranslation(-eastylayer, southHalfOverlapOuter, 0));


	// add SOUTH layers to mother
	TGeoRotation *rot6 = new TGeoRotation();
	rot6->RotateY(90);   
	top->AddNode(layerS1,1, new TGeoCombiTrans(-((eastylayer)+zlayer), southHalfOverlapOuter,0,rot6));
	top->AddNode(layerS2,2, new TGeoCombiTrans(-((eastylayer)-zlayer), -southHalfOverlapInner,-(southxlayerouter-southxlayerinner),rot6));
	
	// north veto panels-------------------------------------------------------------
	// panel 15,16 north outer
	// panel 19,23 north inner
	
	TGeoVolume *layern1 = geom->MakeBox("layern1", Vacuum, northxlayerouter,northylayer,zlayer);   
	TGeoVolume *layern2 = geom->MakeBox("layern2", Vacuum, northxlayerinner,northylayer,zlayer);   

	xpanel = northxlayerouter;    
	ypanel = northylayer/2.;    

	// 15
	sprintf(panel_name,"panel%d",15);
	panel[15]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfAmountOfOverlap,zpanel);
	ypos = -northylayer+(2*0+1)*ypanel;  
	layern1->AddNode(panel[15],15, new TGeoTranslation(0,ypos,0));
	mcscoords[15][0] = northxlayerouter;
	mcscoords[15][1] = ypos;
	mcscoords[15][2] = 0;
	
	//16
	sprintf(panel_name,"panel%d",16);
	panel[16]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfAmountOfOverlap,zpanel);
	ypos = -northylayer+(2*1+1)*ypanel;  
	layern1->AddNode(panel[16],16, new TGeoTranslation(0,ypos,0));
	mcscoords[16][0] = northxlayerouter;
	mcscoords[16][1] = ypos;
	mcscoords[16][2] = 0;
	    
	xpanel = northxlayerinner;
	    
	//19
	sprintf(panel_name,"panel%d",19);
	panel[19]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel+halfAmountOfOverlap,zpanel);
	ypos = -northylayer+(2*0+1)*ypanel;  
	layern2->AddNode(panel[19],19, new TGeoTranslation(0,ypos,0));
	mcscoords[19][0] = northxlayerinner;
	mcscoords[19][1] = ypos;
	mcscoords[19][2] = 0;
	//23
	sprintf(panel_name,"panel%d",23);
	panel[23]= geom->MakeBox(panel_name, Vacuum, xpanel,ypanel-halfAmountOfOverlap,zpanel);
	ypos = -northylayer+(2*1+1)*ypanel;  
	layern2->AddNode(panel[23],23, new TGeoTranslation(0,ypos,0));
	mcscoords[23][0] = northxlayerinner;
	mcscoords[23][1] = ypos;
	mcscoords[23][2] = 0;

	// add north layers to mother
	TGeoRotation *rot3 = new TGeoRotation();
	rot3->RotateY(90);   
	top->AddNode(layern1,1, new TGeoCombiTrans((eastylayer)+zlayer,-halfAmountOfOverlap - (southylayer - northylayer),-(eastxlayerouter-northxlayerouter),rot3));
	top->AddNode(layern2,2, new TGeoCombiTrans((eastylayer)-zlayer,halfAmountOfOverlap - (southylayer - northylayer),-(northxlayerouter-northxlayerinner) -(eastxlayerouter-northxlayerouter),rot3));

	// west veto panels-------------------------------------------------------------
	// panel 12,13 west inner
	// panel 14,22 west outer
	//switch role of x and y to get panel shape
	 
	TGeoVolume *layerw1 = geom->MakeBox("layerw1", Vacuum, westxlayerouter,westylayer,zlayer);   
	TGeoVolume *layerw2 = geom->MakeBox("layerw2", Vacuum, westxlayerinner,westylayer,zlayer);

	xpanel = westxlayerouter;    
	ypanel = (westylayer)/2.;

	// 12
	sprintf(panel_name,"panel%d",12);
	panel[12]= geom->MakeBox(panel_name, Vacuum, ypanel+halfAmountOfOverlap,xpanel,zpanel);
	ypos = -westylayer+(2*0+1)*ypanel;  
	layerw1->AddNode(panel[12],12, new TGeoTranslation(ypos,0,0));
	mcscoords[12][0] = ypos;
	mcscoords[12][1] = southylayer-zlayer;
	mcscoords[12][2] = 0;
	
	//13
	sprintf(panel_name,"panel%d",13);
	panel[13]= geom->MakeBox(panel_name, Vacuum, ypanel-halfAmountOfOverlap,xpanel,zpanel);
	ypos = -westylayer+(2*1+1)*ypanel;  
	layerw1->AddNode(panel[13],13, new TGeoTranslation(ypos,0,0));
	mcscoords[13][0] = ypos;
	mcscoords[13][1] = southylayer-zlayer;
	mcscoords[13][2] = 0;
	
	xpanel = westxlayerinner;  
	
	//22
	sprintf(panel_name,"panel%d",22);
	panel[22]= geom->MakeBox(panel_name, Vacuum, ypanel-halfAmountOfOverlap,xpanel,zpanel);
	ypos = -westylayer+(2*0+1)*ypanel;  
	layerw2->AddNode(panel[22],22, new TGeoTranslation(ypos,0,0));
	mcscoords[22][0] = ypos;
	mcscoords[22][1] = southylayer-zlayer;
	mcscoords[22][2] = 0;
	
	//14
	sprintf(panel_name,"panel%d",14);
	panel[14]= geom->MakeBox(panel_name, Vacuum, ypanel+halfAmountOfOverlap,xpanel,zpanel);
	ypos = -westylayer+(2*1+1)*ypanel;  
	layerw2->AddNode(panel[14],14, new TGeoTranslation(ypos,0,0));
	mcscoords[14][0] = ypos;
	mcscoords[14][1] = southylayer-zlayer;
	mcscoords[14][2] = 0;

	// add west layers to mother
	TGeoRotation *rot4 = new TGeoRotation();
	rot4->RotateX(90);   
	top->AddNode(layerw1,1, new TGeoCombiTrans(halfAmountOfOverlap+2*zlayer,southylayer-zlayer,-(eastxlayerouter-westxlayerouter),rot4));
	top->AddNode(layerw2,2, new TGeoCombiTrans(-halfAmountOfOverlap+2*zlayer,southylayer-3*zlayer,-(westxlayerouter-westxlayerinner) -(eastxlayerouter-westxlayerouter),rot4));

	// done adding panels! -------------------------------------------------------------
	
	// create the cryovats--------------------------------------------------------------
	Double_t cryoboxx = 34.34 / 2;
	Double_t cryoboxy = 68.68 / 2;
	Double_t cryoboxz = 44.56 / 2;
	Double_t yoffset = 9.16;
	Double_t halfGapBetween = 3.19;
	
	TGeoVolume *cryobox = geom->MakeBox("cryobox", Vacuum, cryoboxx,cryoboxy,cryoboxz);

	TGeoVolume *cryovat1 = geom->MakeEltu("cryovat1", Vacuum, cryoboxy/2, cryoboxy/2, cryoboxz/2);
	cryobox->AddNode(cryovat1,1, new TGeoTranslation(0,cryoboxy/2 + halfGapBetween - yoffset,0));
	
	TGeoVolume *cryovat2 = geom->MakeEltu("cryovat2", Vacuum, cryoboxy/2, cryoboxy/2, cryoboxz/2);
	cryobox->AddNode(cryovat2,2, new TGeoTranslation(0,-cryoboxy/2 - halfGapBetween - yoffset,0));
	
	top->AddNode(cryobox, 1, 0);
	
	
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

	vector<vector<Double_t>> allHitCoords; //all hit locations stored in 2-D vector
	bool alreadyDone[32] = {false}; //used so the same layer hit isnt added twice
	int numOfPoints = 1;
	
	//these 'for' loops check if layer partner hit exists, if it does, draws a hole for the layer hit
	for(int i = 0; i < 32; i++) {
		if(qdcVals[i] != 0) {
			for(int j = 0; j < 32; j++) {
				if(qdcVals[j] != 0 && isLayerHit(i, j) == 1 && alreadyDone[j] == false) {
					numOfPlanesHit++;
					vector<Double_t> drawLoc = hitLocation(i, j);
					allHitCoords.push_back(drawLoc);

					//if a layer hit exsists, create the "hole"
					TGeoVolume *hole1 = geom->MakeTube("hole1", Vacuum, 9,10,0);
					//orientation of the hole
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
   
   
   //draw the tracks for any two sides ---------------------------------
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

	if( (numOfPlanesHit <= 2 && numberOfPanelsHit == 4) && (bottomSide == 2 || topSide == 2
		|| northSide == 2 || eastSide == 2 || southSide == 2 || westSide == 2) || (numberOfPanelsHit == 5 && topSide == 2 && bottomSide == 2) ) {
		geom->DrawTracks();
	}
	// -----------------------------------------------------------------
	
	// Fill the text box -----------------------------------------------
	TPaveText *pt = new TPaveText(-.9,.9,-.4,.5,"br"); //x1,y1,x2,y2
	//-.9,-.9,-.4,-.5 for bottom left corner coords
	
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
		//looks for first panel hit
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
	
	//set phi and theta to 361 to know if set (0 doesn't work because 0 phi,0 theta is a possible location)
	phi = 361;
	theta = 361;
	//some vector info
	if( (allHitCoords.size() > 1 && numberOfPanelsHit == 4) || (allHitCoords.size() > 1 && numberOfPanelsHit == 5 && topSide == 2 && bottomSide == 2) ) {
		TVector3 r1;
		TVector3 r2;
		TVector3 t1;
		r1.SetXYZ(allHitCoords[0][0],allHitCoords[0][1],allHitCoords[0][2]);
		r2.SetXYZ(allHitCoords[1][0],allHitCoords[1][1],allHitCoords[1][2]); //top panels are always a higher number than the bottom panels
		t1 = r2-r1;
	
		phi = t1.Phi()*(180/pi); //180/pi to get deg rather than rad
		theta = t1.Theta()*(180/pi);
		
		if(phi < 0) {
			phi += 360; //makes negative phi's the correct positive
		}
		
		//check if the track went through the cryostats; -------
		bool passedThroughCryo = false;
		Double_t degrad = pi/180;
		Double_t point[3], dircos[3];
		
		point[0] = allHitCoords[1][0];
		point[1] = allHitCoords[1][1];
		point[2] = allHitCoords[1][2];

		dircos[0] = sin(theta*degrad) * -cos(phi*degrad); //the direction of the track in normallized vector
		dircos[1] = sin(theta*degrad) * -sin(phi*degrad);
		dircos[2] = -cos(theta*degrad);
		
		printf("direction : %f %f %f\n", dircos[0], dircos[1], dircos[2]);

		geom->InitTrack(point,dircos);
		
		while(!geom->IsOutside()) {
			printf("current: %s\n", geom->GetPath());
			string path = geom->GetPath();
			if(path.compare("/TOP_1/cryobox_1") == 0) { 
				passedThroughCryo = true;
			}
			geom->FindNextBoundary();
			cout<<"Distance to Boundary is "<<geom->GetStep()<<" cm."<<endl;
			cout<<"Is step entering: " << geom->IsStepEntering() << endl;;
			geom->Step(); //Steps to the next geometry boundry
		}
		// -------------------------------------------------------------
		
		char passedThroughOut[150];
		sprintf(passedThroughOut,"Passed through Cryostat: %d",passedThroughCryo); 
		pt->AddText(passedThroughOut);
	
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
	{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0},
	{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,1,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,1,1,0,0,1,0,0,1,0,1,0,1,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,1,0},
	{0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,1,0,0,0,0,0,0,1,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,1,1,0,1,1,0,0,0},
	{0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0},
	{0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,1,0,0,1,1,0,0,0,0,1,0,0,1,1,1},
	{0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1},
	{0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,1,0,1,0,1,0,1,0,0,0,0,0,0,1,1,0,1},
	{0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0}
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

//Finds the user's input hit and displays it in 3-D with application.
void hitFinder(char* innum1, char* innum2) {
	
	TApplication *myapp = new TApplication("myapp",0,0);
	string input = "";
	int inputnum1 = 0;
	int inputnum2 = 0;
	Int_t runNumber = 0;
	Int_t entry = 0;
	Int_t eventCount = 0;
	Double_t scalerTime = 0.0;
	Int_t qdcVals[32]={0};
	Int_t numberOfPanelsHit = 0;
	Int_t totalQDC = 0;
	string line = "";
	int lineLength = 0;
	bool found = false;
	bool validNum1 = false;
	bool validNum2 = false;
	
	ifstream VetoHitsFile;
    VetoHitsFile.open(inputFile);
    
	stringstream myStream(innum1);
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
				for(int i = 0; i < 32; i++) {
					if(qdcVals[i] != 0) {
						numberOfPanelsHit++;
						totalQDC += qdcVals[i];
					}
				}
				DrawEvent(qdcVals, numberOfPanelsHit, totalQDC, runNumber, eventCount);
				cout << "Press ctrl+C to exit" << endl;
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
TCanvas *QDCanglecan = new TCanvas("QDCanglecan"," ", 1200, 800);
TCanvas *QDCanglecan2 = new TCanvas("QDCanglecan2"," ", 1200, 800);
TCanvas *QDCanglecan3 = new TCanvas("QDCanglecan3"," ", 1200, 800);
TCanvas *QDCslantcan = new TCanvas("QDCslantcan"," ", 1200, 800);
TCanvas *inThrucan = new TCanvas("inThrucan"," ", 1200, 800);
TCanvas *totalQDCanglecan = new TCanvas("totalQDCanglecan"," ", 1200, 800);
TCanvas *totalQDCanglecan2 = new TCanvas("totalQDCanglecan2"," ", 1200, 800);
TCanvas *totalQDCanglecan3 = new TCanvas("totalQDCanglecan3"," ", 1200, 800);
TCanvas *totalQDCslantcan = new TCanvas("totalQDCslantcan"," ", 1600, 700);

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
TH2F *thetaPhi = new TH2F("thetaPhi","Theta vs. Phi for top and bottom n=4 layer hits", 90, 0, 360, 90, -1, 50);
TH2F *costhetaPhi = new TH2F("costhetaPhi","Cos(Theta) vs. Phi for top and bottom n=4 layer hits", 90, 0, 360, 90, 0, 1.3);
TH2D *QDCangle = new TH2D("QDCangle","QDC vs Theta for top and bottom n=4 layer hits", 180, 0, 45, 180, 0 , 5000);
TH2D *QDCangleTwo = new TH2D("QDCangleTwo","9-13 theta slice", 180, 9, 13, 180, 0 , 5000);
TH2D *QDCangleThree = new TH2D("QDCangleThree","13-17 theta slice", 180, 13, 17, 180, 0 , 5000);
TH2D *QDCangleFour = new TH2D("QDCangleFour","17-23 theta slice", 180, 17, 23, 180, 0 , 5000);
TH2D *QDCangleFive = new TH2D("QDCangleFive","23-32 theta slice", 180, 23, 32, 180, 0 , 5000);
TH2D *QDCangleSix = new TH2D("QDCangleSix","32-45 theta slice", 180, 32, 45, 180, 0 , 5000);
TH2F *QDCslant = new TH2F("QDCslant","QDC vs Slant Depth for top and bottom n=4 layer hits", 100, 1200, 2200, 100, 0, 5000);
TH2F *thetaSlant = new TH2F("thetaSlant","Theta vs Slant Depth for top and bottom n=4 layer hits",100, 1200, 2200, 100, 0, 45);
TH2F *inThruHist = new TH2F("inThruHist","Inches particle passed through detector vs QDC val for top panels only", 180, .95, 1.5, 180, 0, 5000);
TH2D *totalQDCangle = new TH2D("totalQDCangle","Total QDC for top panels vs Theta for top and bottom n=4 layer hits", 180, 0, 90, 180, 0 , 12000);
TH2D *totalQDCangleTwo = new TH2D("totalQDCangleTwo","Total QDC for top panels vs Theta for top and bottom n=4 layer hits 0-9", 180, 0, 9, 180, 0 , 12000);
TH2D *totalQDCangleThree = new TH2D("totalQDCangleThree","Total QDC for top panels vs Theta for top and bottom n=4 layer hits 9-13", 180, 9, 13, 180, 0 , 12000);
TH2D *totalQDCangleFour = new TH2D("totalQDCangleFour","Total QDC for top panels vs Theta for top and bottom n=4 layer hits 13-17", 180, 13, 17, 180, 0 , 12000);
TH2D *totalQDCangleFive = new TH2D("totalQDCangleFive","Total QDC for top panels vs Theta for top and bottom n=4 layer hits 17-23", 180, 17, 23, 180, 0 , 12000);
TH2D *totalQDCangleSix = new TH2D("totalQDCangleSix","Total QDC for top panels vs Theta for top and bottom n=4 layer hits 23-32", 180, 23, 32, 180, 0 , 12000);
TH2D *totalQDCangleSeven = new TH2D("totalQDCangleSeven","Total QDC for top panels vs Theta for top and bottom n=4 layer hits 32-45", 180, 32, 45, 180, 0 , 12000);
TH2F *totalQDCslant = new TH2F("totalQDCslant","Total QDC vs Slant Depth for top and bottom n=4 layer hits", 180, 1200, 2200, 180, 0, 12000);
	
//fills all plots and prints out wireframe pdfs
void fillPlots(Int_t qdcvals[], Int_t totalQDC, Int_t numberOfPanelsHit, Int_t ievent) {
	
	bool totalQDCangleFilled = false;
	bool totalQDCangleFilled2 = false;
	Double_t totalQDCtop = 0.0;
	bool isNotNextToAny = true;
	Int_t whichPanel = 0;
	
	graph1->Fill(numberOfPanelsHit);
	graph4->Fill(numberOfPanelsHit, totalQDC);
	graph8->Fill(totalQDC);	

	//gets total QDC for just the top panels for top and bottom 4 panel events
	for (Int_t i = 0; i < 32; i++) {
		if( qdcvals[i] != 0 && (i == 17 || i == 18 || i == 20 || i == 21) && numOfPlanesHit <= 2 && bottomSide == 2 && topSide == 2 && numberOfPanelsHit == 4) {
			totalQDCtop += qdcvals[i];
		}
	}

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
			
			if(numberOfPanelsHit == 4) {
				if(numOfPlanesHit <= 2 && bottomSide == 2 && topSide == 2) { //top and bottom plane hits
					if(i == 17 || i == 18 || i == 20 || i == 21) { //top panels only
					
					QDCangle->Fill(theta, qdcvals[i]);
					QDCangleTwo->Fill(theta, qdcvals[i]); //just slices of QDCangle, there is almost certainly a better way to do this
					QDCangleThree->Fill(theta, qdcvals[i]);
					QDCangleFour->Fill(theta, qdcvals[i]);
					QDCangleFive->Fill(theta, qdcvals[i]);
					QDCangleSix->Fill(theta, qdcvals[i]);
					
					//TODO:: IMPLEMENT SLANT DEPTH THROUGH ROOT FILE, NOT METHOD
					//Double_t slantDepth = SlantDepth(phi,theta);
					//thetaSlant->Fill(slantDepth, theta);
					//QDCslant->Fill(slantDepth, qdcvals[i]);
					
					Double_t inThru = 1 / cos(theta * (pi/180)); //the track length through the detectors per panel
					inThruHist->Fill(inThru, qdcvals[i]);
				
					if(totalQDCangleFilled == false) {
						totalQDCangle->Fill(theta, totalQDCtop);
						totalQDCangleTwo->Fill(theta, totalQDCtop); //again, just slices of totalQDCangle, there is almost certainly a better way to do this
						totalQDCangleThree->Fill(theta, totalQDCtop);
						totalQDCangleFour->Fill(theta, totalQDCtop);
						totalQDCangleFive->Fill(theta, totalQDCtop);
						totalQDCangleSix->Fill(theta, totalQDCtop);
						totalQDCangleSeven->Fill(theta, totalQDCtop);
						//totalQDCslant->Fill(slantDepth, totalQDCtop);
						totalQDCangleFilled = true;
					}
				
					}
				}
			}
		}
	}
	
	
	/*
	//adding 5 panel hit events to data by excluding a panel based on QDC value or not next to others
	if(numberOfPanelsHit == 5) {
		for(Int_t j = 0; j < 32; j++) {
			if(qdcvals[j] != 0) {
				for(Int_t k = 0; k < 32; k ++) {
					if(qdcvals[k] != 0 && isNextTo(j,k) == 1) {
						isNotNextToAny = false;
					}
				}
				if(isNotNextToAny == true || qdcvals[j] < 500) {
					whichPanel = j;
				}
			}
		}
		if( (isNotNextToAny == true || (qdcvals[whichPanel] < 500 && (whichPanel != 17 || whichPanel != 18 || whichPanel != 20 || whichPanel != 21))) && totalQDCangleFilled2 == false) {
			if(numOfPlanesHit <= 2 && bottomSide == 2 && topSide == 2 && theta != 361) {
				totalQDCangle->Fill(theta, totalQDCtop);
				totalQDCangleTwo->Fill(theta, totalQDCtop);
				totalQDCangleThree->Fill(theta, totalQDCtop);
				totalQDCangleFour->Fill(theta, totalQDCtop);
				totalQDCangleFive->Fill(theta, totalQDCtop);
				totalQDCangleSix->Fill(theta, totalQDCtop);
				totalQDCangleSeven->Fill(theta, totalQDCtop);
				totalQDCangleFilled2 = true;
			}
		}
	}
	*/
	
	bool filled = false;
	//fill the theta and phi graphs
	for(int i = 0; i < 12; i++) {
		for(int j = 0; j < 4; j++) {
			if(phi != 361 && theta != 361 && qdcvals[bottomPanels[i]] != 0 && 
			   qdcvals[topPanels[j]] != 0 && numberOfPanelsHit == 4 && filled == false) {
				
				thetaPhi->Fill(phi, theta);
				Double_t costheta = cos(theta / (180/pi));	
				costhetaPhi->Fill(phi, costheta);
				filled = true;
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
	
	ecan->cd();
}

//the editing and drawing of all the graphs
void drawPlots() {
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
	QDCslantcan->Divide(2,2);
	
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
	//TODO:: GENERATE ROOT FILES FOR GRAPHS, MAKE NEW FILE FOR FITTING
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
	totalQDCanglecan2->Divide(3,2);
	
	totalQDCanglecan2->cd(1);
		totalQDCangleTwo->Draw("colz");
	
	TH1D *totalQDCangleprojy2 = totalQDCangleTwo->ProjectionY("totalqdcangle2",1,180);
	totalQDCanglecan2->cd(4);
		totalQDCangleprojy2->SetTitle("Projection Y");
		totalQDCangleprojy2->SetYTitle("Count");
		totalQDCangleprojy2->Draw("bar");
	
	totalQDCanglecan2->cd(2);
		totalQDCangleThree->Draw("colz");
	
	TH1D *totalQDCangleprojy3 = totalQDCangleThree->ProjectionY("totalqdcangle3",1,180);
	totalQDCanglecan2->cd(5);
		totalQDCangleprojy3->SetTitle("Projection Y");
		totalQDCangleprojy3->SetYTitle("Count");
		totalQDCangleprojy3->Draw("bar");
	
	totalQDCanglecan2->cd(3);
		totalQDCangleFour->Draw("colz");
	
	TH1D *totalQDCangleprojy4 = totalQDCangleFour->ProjectionY("totalqdcangle4",1,180);
	totalQDCanglecan2->cd(6);
		totalQDCangleprojy4->SetTitle("Projection Y");
		totalQDCangleprojy4->SetYTitle("Count");
		totalQDCangleprojy4->Draw("bar");
	
	//totalQDCanglecan3
	totalQDCanglecan3->Divide(3,2);
	
	totalQDCanglecan3->cd(1);
		totalQDCangleFive->Draw("colz");
	
	TH1D *totalQDCangleprojy5 = totalQDCangleFive->ProjectionY("totalqdcangle5",1,180);
	totalQDCanglecan3->cd(4);
		totalQDCangleprojy5->SetTitle("Projection Y");
		totalQDCangleprojy5->SetYTitle("Count");
		totalQDCangleprojy5->Draw("bar");
	
	totalQDCanglecan3->cd(2);
		totalQDCangleSix->Draw("colz");
	
	TH1D *totalQDCangleprojy6 = totalQDCangleSix->ProjectionY("totalqdcangle6",1,180);
	totalQDCanglecan3->cd(5);
		totalQDCangleprojy6->SetTitle("Projection Y");
		totalQDCangleprojy6->SetYTitle("Count");
		totalQDCangleprojy6->Draw("bar");
	
	totalQDCanglecan3->cd(3);
		totalQDCangleSeven->Draw("colz");
	
	TH1D *totalQDCangleprojy7 = totalQDCangleSeven->ProjectionY("totalqdcangle7",1,180);
	totalQDCanglecan3->cd(6);
		totalQDCangleprojy7->SetTitle("Projection Y");
		totalQDCangleprojy7->SetYTitle("Count");
		totalQDCangleprojy7->Draw("bar");
	
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

	char QDCanglecanprint[150];
		sprintf(QDCanglecanprint,"output/plots/QDCangle.pdf");
		QDCanglecan->Print(QDCanglecanprint,"pdf");
		
	char QDCslantcanprint[150];
		sprintf(QDCslantcanprint,"output/plots/QDCslant.pdf");
		QDCslantcan->Print(QDCslantcanprint,"pdf");
		
	char QDCanglecan2print[150];
		sprintf(QDCanglecan2print,"output/plots/QDCangle2.pdf");
		QDCanglecan2->Print(QDCanglecan2print,"pdf");
		
	char QDCanglecan3print[150];
		sprintf(QDCanglecan3print,"output/plots/QDCangle3.pdf");
		QDCanglecan3->Print(QDCanglecan3print,"pdf");
		
	char inThrucanprint[150];
		sprintf(inThrucanprint,"output/plots/inThru.pdf");
		inThrucan->Print(inThrucanprint,"pdf");
		
	char totalQDCanglecanprint[150];
		sprintf(totalQDCanglecanprint,"output/plots/totalQDCangle.pdf");
		totalQDCanglecan->Print(totalQDCanglecanprint,"pdf");
		
	char totalQDCanglecanprint2[150];
		sprintf(totalQDCanglecanprint2,"output/plots/totalQDCangle2.pdf");
		totalQDCanglecan2->Print(totalQDCanglecanprint2,"pdf");
		
	char totalQDCanglecanprint3[150];
		sprintf(totalQDCanglecanprint3,"output/plots/totalQDCangle3.pdf");
		totalQDCanglecan3->Print(totalQDCanglecanprint3,"pdf");
		
	char totalQDCslantcanprint[150];
		sprintf(totalQDCslantcanprint,"output/plots/totalQDCslant.pdf");
		totalQDCslantcan->Print(totalQDCslantcanprint,"pdf");
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
    // MAIN LOOP -------------------------------------------------------
    while(getline(VetoHitsFile, line)) {
		ievent++;
		//getchar();
		// reset counters and read first part of next line
		numberOfPanelsHit = 0;
		totalQDC = 0;
		lineLength = 0;
		goodPanelNum = false;
		topSide = 0;
		bottomSide = 0;
		northSide = 0;
		eastSide = 0;
		southSide = 0;
		westSide = 0;
		numOfPlanesHit = 0;
		
		stringstream  lineStream(line);
		lineStream >> runNumber >> entry >> eventCount >> scalerTime;
		printf("runNumber:%i  entry:%i  eventCount:%i  scalerTime:%.5f \n",runNumber,entry,eventCount,scalerTime);

		
		while (lineStream.peek() == ' ') { //gets rid of trailing space before qdcvals if any
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
			if(lineLength == 32) { //stops trailing spaces/values after 32 qdcvals if any
				break;
			}
		}
	
		if(lineLength == 32 || lineLength == 24) {
			goodPanelNum = true;
		}
		if(goodPanelNum == true) {
			//if linelength is 24, fill in any missing detectors with 0 qdc
			if(lineLength != 32) {
				for(int i = lineLength; i < 32; i++) {
					qdcVals[i] = 0;
				}
			}
			
			//find panels hit, totalqdc per event, and total added qdc per panel
			for(int i = 0; i < 32; i++) {
				if(qdcVals[i] != 0) {
					numberOfPanelsHit++;
					totalQDC += qdcVals[i];
					qdcValSum[i] += qdcVals[i];
				}
			}
		}
	
			
		DrawEvent(qdcVals, numberOfPanelsHit, totalQDC, runNumber, eventCount);
		fillPlots(qdcVals, totalQDC, numberOfPanelsHit, ievent);
		
		// hand picked events ------------------------------------------
		Double_t totalQDCsides = 0.0;
		if( (runNumber == 9954 && eventCount == 293) || (runNumber == 9483 && eventCount == 49) || (runNumber == 10015 && eventCount == 1071)
			|| (runNumber == 9645 && eventCount == 357) || (runNumber == 10124 && eventCount == 465) || (runNumber == 9905 && eventCount == 695)
			|| (runNumber == 10083 && eventCount == 936) ) {
			for(int q = 0; q < 4; q++) {
				totalQDCsides += qdcVals[northPanels[q]];
				totalQDCsides += qdcVals[southPanels[q]];
				totalQDCsides += qdcVals[eastPanels[q]];
				totalQDCsides += qdcVals[westPanels[q]];
			}
		}
		if(totalQDCsides != 0.0) {
			totalQDCangle->Fill(75, totalQDCsides);
		}
		// -------------------------------------------------------------
	}
	// MAIN LOOP END ---------------------------------------------------
	
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
	}
	
	//draw the event on last time (heatmap)
	DrawEvent(qdcValSum, totalQDC, numberOfPanelsHit, runNumber, eventCount); 
	
	//print out the heatmap
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
