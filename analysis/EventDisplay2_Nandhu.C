#include <iostream>
#include <TF1.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLinearFitter.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TRandom.h>
#include <TGButton.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <RQ_OBJECT.h>
#include <TGNumberEntry.h>
#include <TGTextEntry.h>

//Making an event display
class EventDisplay
{
public:
  TList * myGeometryData;
  int global_run_number;
  int Analyse_Secondaries = 1;
  float Theta_min_cut = 2.524;
  float ThetaVerticalCut = 3.02;
  int Photon_min_cut = 0;

  int MaxPMTNo = 50000;
  int MaxPMTHits = 1000;
  float Finger_Edep_Max = 10.0;
  float AnaBar_Edep_Max = 10.0;
  float pedastel_sigma = 2.9;
  int Detector_Offset = 2560;
  int Detector_PMT_Offset = 2500;
  int AnaBar_Offset = 30000;
  int AnaBar_PMT_Offset = 0;

  int Finger_NPhotons_Max = 250;
  int AnaBar_NPhotons_Max = 200;

  const int NUMPADDLE = 14;
  const int NUMBARS = 14;
  const int NUMMODULES = 3;
  const int NUMSIDES = 2;
  const int NUMLAYERS = 2;

  const int NDET = NUMPADDLE * NUMBARS * NUMMODULES * NUMSIDES * NUMLAYERS;

  int NMaxPMT = 14;

  Double_t padLen = 50.2;	//Length of one paddle
  Double_t padWid = 0.54;	//Width of one paddle

  TGMainFrame *display;
  TRootEmbeddedCanvas *screen;
  TGLabel *eventNo;
  int event = 0;
  int run_number = 4000; 
  TFile *f;
  TTree *t;
  TFile *noisef;
  TTree *noise;
  
  Int_t PMT_Nphotons[50000];
  Int_t PMT_Noise_photons[50000];
  float PMT_Time[5000];
  float PMT_Noise_Time[5000];
  
  Int_t Detector_id[50000];
  Int_t Detector_Noise_id[50000];
  Int_t Detector_Nhits;
  TH2Poly *top;
  TH2Poly *bot;
  TH2Poly *topNoise;
  TH2Poly *botNoise;
  TH2Poly *topTemplate;
  TH2Poly *botTemplate;
  TPolyMarker *topCalc;
  TPolyMarker *botCalc;
  TPolyMarker *topCalcNoise;
  TPolyMarker *botCalcNoise;
  TGCheckButton *simulated;
  TGCheckButton *noiseEvents;
  TGCheckButton *time;
  TGNumberEntryField *num;

  int getSide (int fID)
  {
    int iLayer = fID % NUMPADDLE;
    int iBar = ((fID - iLayer) / NUMPADDLE) % NUMBARS;
      return (((fID - iLayer) / NUMPADDLE - iBar) / NUMBARS) % NUMSIDES;
  }

  int getPlane (int fID)
  {
    int iLayer = fID % NUMPADDLE;
    int iBar = ((fID - iLayer) / NUMPADDLE) % NUMBARS;
    int iSide = (((fID - iLayer) / NUMPADDLE - iBar) / NUMBARS) % NUMSIDES;
    int iModule =
      ((((fID - iLayer) / NUMPADDLE -
	 iBar) / NUMBARS) / NUMSIDES) % NUMMODULES;
    return (((((fID - iLayer) / NUMPADDLE -
	       iBar) / NUMBARS) / NUMSIDES) / NUMMODULES) % NUMLAYERS;
  }

  float getXOffsetFromTime (int fID, float time)
  {

    float xoffset;
    int iSide = getSide (fID);
    int iPlane = getPlane (fID);
    TRandom3 *fRand = new TRandom3 (0);

    /* Old calibration - average time
       if (iPlane == 0) {
       double a = -0.00094362480519633;
       double b = -0.03673192847114299;
       double c = 8.609456854512016;
       double xmin = -25.0;
       double ymin = a*xmin*xmin+b*xmin+c-0.003;
       if (time<ymin) {
       xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
       } else {
       xoffset = -25.0+fRand->Uniform(0.0,12.0);
       }
       } else {
       double a = -0.0009339487175907869;
       double b = -0.03579463613239478;
       double c = 9.59488826868313;
       double xmin = -25.0;
       double ymin = a*xmin*xmin+b*xmin+c-0.003;
       if (time<ymin) {
       xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
       } else {
       xoffset = -25.0+fRand->Uniform(0.0,12.0);
       }
       }
     */

    if (iPlane == 0)
      {
	double a = -0.00037499653844674384;
	double b = -0.051184438607694095;
	double c = 3.5929880196450843;
	double discriminant = b * b - 4 * a * (c - time);
	if (discriminant >= 0.0)
	  {
	    xoffset =
	      (-b - sqrt (fabs (b * b - 4 * a * (c - time)))) / (2.0 * a);
	  }
	else
	  {
	    xoffset = -25.0;
	  }
      }
    else
      {
	double a = -0.0004090552401677775;
	double b = -0.050453362706664166;
	double c = 4.537007798185197;
	double discriminant = b * b - 4 * a * (c - time);
	if (discriminant >= 0.0)
	  {
	    xoffset =
	      (-b - sqrt (fabs (b * b - 4 * a * (c - time)))) / (2.0 * a);
	  }
	else
	  {
	    xoffset = -25.0;
	  }
      }

    if (iSide == 1)
      {
	xoffset = -xoffset;
      }

    cout << "==========================" << endl;
    cout << "event = " << event << endl;
    cout << "xoffset = " << xoffset << endl;
    cout << "time = " << time << endl;
    cout << "iPlane = " << iPlane << endl;
    cout << "iSide = " << iSide << endl;
    cout << "fID = " << fID << endl;
    cout << "getSide(fID) = " << getSide (fID) << endl;
    cout << "getPlane(fID) = " << getPlane (fID) << endl;

    return xoffset;

  }

  bool getTrigger (int Detector_Nhits, int *Detector_id)
  {

    bool tophit = false;
    bool bottomhit = false;
    bool fhit = false;
    bool ahit = false;
    bool trigger = false;
    for (int j = 0; j < Detector_Nhits; j++)
      {
	//std::cout << "Detector id = " << Detector_id[j] << std::endl;
	if ((Detector_id[j] == Detector_Offset
	     || Detector_id[j] == Detector_Offset + 1) && !tophit)
	  {
	    tophit = true;
	    //std::cout << "Top hit" << Detector_id[j] << std::endl;
	  }
	if ((Detector_id[j] == Detector_Offset + 2
	     || Detector_id[j] == Detector_Offset + 3) && !bottomhit)
	  {
	    bottomhit = true;
	    //std::cout << "Bottom hit" << Detector_id[j] << std::endl;
	  }
	for (int ibar = 0;
	     ibar < NUMPADDLE * NUMBARS * NUMMODULES * NUMSIDES * NUMLAYERS;
	     ibar++)
	  {
	    if (Detector_id[j] == AnaBar_Offset + ibar)
	      {
		ahit = true;
	      }
	  }
	//if (tophit && bottomhit) {
	if (ahit)
	  {
	    fhit = true;
	    trigger = true;
	  }
      }

    return trigger;
  }

  //Clear all histograms and graphs.
  void Clear(){
    top->ClearBinContents ();
    bot->ClearBinContents ();
    topNoise->ClearBinContents ();
    botNoise->ClearBinContents ();
    topCalc = new TPolyMarker ();
    botCalc = new TPolyMarker ();
    topCalc->SetMarkerStyle (kFullCircle);
    botCalc->SetMarkerStyle (kFullCircle);
    topCalc->SetMarkerColor (kPink);
    botCalc->SetMarkerColor (kPink);
    topCalcNoise = new TPolyMarker ();
    botCalcNoise = new TPolyMarker ();
    topCalcNoise->SetMarkerStyle (kFullCircle);
    botCalcNoise->SetMarkerStyle (kFullCircle);
    topCalcNoise->SetMarkerColor (kBlue);
    botCalcNoise->SetMarkerColor (kBlue);
  }

  void Draw(){
    auto eDisplay = screen->GetCanvas();
    auto legend = new TLegend(0,0,1,1,"");
    legend->SetBorderSize(0);
    //auto legendPad = new TPad("Legend","I'm a TPad",0.43,0.17,0.63,0.37);
    auto legendPad = new TPad("Legend","I'm a TPad",0,0,1,0.05);
    legendPad->Draw();

    //Set Palettes
    //gStyle->SetPalette(kBeach);
    TExec *ex1 = new TExec("ex1"," gStyle->SetPalette(kSolar);");
    TExec *ex2 = new TExec("ex2"," gStyle->SetPalette(kBeach);");
    TExec *ex3 = new TExec("ex3"," gStyle->SetPalette(kBird);");
    eDisplay->cd(1);
    //eDisplay->cd(1)->SetFrameLineColor(0);

    //legend->AddEntry((TObject*)0, "Legend:", "");

    topTemplate->Draw ("COL");
    ex3->Draw();
    topTemplate->Draw ("COL SAME");
    if(simulated->IsDown()){
    	ex1->Draw();
	top->SetLineColor(kOrange);
	top->SetLineWidth(6);
        top->Draw("COL SAME 0");
   	legend->AddEntry(top, "Real Hits", "l");
    }
    if(noiseEvents->IsDown()){
    	ex2->Draw();
	topNoise->SetLineColor(kViolet);
	topNoise->SetLineWidth(6);
    	topNoise->Draw ("SAME COL 0");
   	legend->AddEntry(topNoise, "Noise Hits", "l");
    }
    if(time->IsDown()){
	if(simulated->IsDown()){
    		topCalc->Draw ("SAME *");
   		legend->AddEntry(topCalc, "Calibrated", "p");
	}
	if(noiseEvents->IsDown()){
    		topCalcNoise->Draw ("SAME *");
   		legend->AddEntry(topCalcNoise, "Noise", "p");
	}
    }

    eDisplay->cd(2);
    //eDisplay->cd(2)->SetFrameLineColor(0);

    botTemplate->Draw ("COL");
    ex3->Draw();
    botTemplate->Draw ("COL SAME");
    if(simulated->IsDown()){
    	ex1->Draw();
        bot->Draw("COL SAME 0");
    }
    if(noiseEvents->IsDown()){
    	ex2->Draw();
    	botNoise->Draw ("SAME COL 0");
    }
    if(time->IsDown()){
	if(simulated->IsDown()){
    		botCalc->Draw ("SAME *");
	}
	if(noiseEvents->IsDown()){
    		botCalcNoise->Draw ("SAME *");
	}
    }

    //Remove the stats box
    topTemplate->SetStats (0);
    botTemplate->SetStats (0);
    legendPad->cd();
    legend->SetNColumns(legend->GetNRows());
    legend->Draw();

    //Allow interaction
    eDisplay->Draw ();
    eDisplay->cd();
    eDisplay->Update();
  }

  TPolyMarker* CalculateXPosition(TPolyMarker *calc, int i, int* photons, float* time){
	if (photons[i] > Photon_min_cut){			
		//Only if there are enough photons to justify it
		TVectorD *geo = (TVectorD *) myGeometryData->At (i);
		int xdpos = (*geo)[1] / 10.0;
		int zdpos = (*geo)[3] / 10.0;
		//cout << getXOffsetFromTime(i, PMT_Time[i]);
		//topCalc->SetNextPoint(xdpos + getXOffsetFromTime(i, PMT_Time[i]),zdpos);
		//topCalc->SetNextPoint (xdpos + 25 + getXOffsetFromTime (i, PMT_Time[i]), zdpos);
		calc->SetNextPoint(xdpos + getXOffsetFromTime(i, time[i]), zdpos);
	}
	return calc;
  }

  void Iterate()
  {

    //Show What Event We're looking at.
    std::string title = "Event " + std::to_string (event);
    eventNo->ChangeText(title.c_str());
    eventNo->Resize(eventNo->GetDefaultSize());

    //Reset the histograms so we can fill them again.
    Clear();

    float xdpos;
    float zdpos;

    //Get the event information.
    t->GetEntry (event);

    //Check if the event is valid before making the display.
    /*while (getTrigger (Detector_Nhits, Detector_id) == false)
      {
	event++;
	t->GetEntry (event);
      }*/

    //Fill the bins
    for (int i = 0; i < AnaBar_PMT_Offset + NUMPADDLE * NUMBARS * NUMMODULES * NUMSIDES * NUMLAYERS; i++){
	if (i < 1176)
	  {
	    //Chooses which histograms is filled
	    if(PMT_Nphotons[i] != PMT_Noise_photons[i]){
		topNoise->SetBinContent (i, PMT_Noise_photons[i]);
	   	topCalcNoise = CalculateXPosition(topCalcNoise,i,PMT_Noise_photons, PMT_Noise_Time);
		continue;
	    }

	    top->SetBinContent (i, PMT_Nphotons[i]);
	    //Draw the calculated x position
	    topCalc = CalculateXPosition(topCalc,i,PMT_Nphotons, PMT_Time);
	  }
	else
	  {
	   //Chooses which histograms is filled
	    if(PMT_Nphotons[i] != PMT_Noise_photons[i]){
		botNoise->SetBinContent (i - 1175, PMT_Noise_photons[i]);
		botCalcNoise = CalculateXPosition(botCalcNoise,i,PMT_Noise_photons, PMT_Noise_Time);
		continue;
	    }
	    bot->SetBinContent (i - 1175, PMT_Nphotons[i]);

	    //Draw the calculated x position        
	    botCalc = CalculateXPosition(botCalc,i,PMT_Nphotons, PMT_Time);
	  }
      }
      Draw();
  }

  //Initializes the files
  void FileHandler(){
    auto fileName = "data/AnaBarMC_" + std::to_string (run_number) + ".root";
    auto treeName = "T";

    f = new TFile ((TString) fileName, "READ");
    t = 0;
    f->GetObject (treeName, t);

    //auto noisefileName =
    //  "data/AnaBarMC_" + std::to_string (run_number) + "noise.root";
    //auto noisetreeName = "noise";

    //noisef = new TFile ((TString) noisefileName, "READ");
    //noise = 0;
    //noisef->GetObject (noisetreeName, noise);
    //t->AddFriend (noise); 

    myGeometryData = (TList *) t->GetUserInfo ()->FindObject ("myGeometryData");
  }

  //Initialise the deposition arrays
  void InitializeArrays(){
    
    //Set the arrays to 0
    for ( Int_t i = 0; i < 5000; i++) {
	PMT_Nphotons[i] = 0;;
  	Detector_id[i] = 0;
  	PMT_Time[i] = 0;	
	PMT_Noise_photons[i] = 0;
        PMT_Noise_Time[i] = 0;
    }

    Int_t Detector_Nhits = 0;

    //Get the Photon deposition array
    t->SetBranchAddress ("PMT_Nphotons", &PMT_Nphotons);
    t->SetBranchAddress ("Detector_id", &Detector_id);
    t->SetBranchAddress ("Detector_Nhits", &Detector_Nhits);
    t->SetBranchAddress ("PMT_Time", &PMT_Time);

    t->SetBranchAddress ("PMT_Nphotons_Noise", &PMT_Noise_photons);
    t->SetBranchAddress ("PMT_Time_Noise", &PMT_Noise_Time);
  }

  //Builds the GUI by which we can access the Event Display
  void BuildApplication(){
    display = new TGMainFrame(gClient->GetRoot(),800,800);

    //Show the Event
    TGHorizontalFrame *topbar = new TGHorizontalFrame(display,800,100);
    eventNo = new TGLabel(topbar, "Event 1000");
    eventNo->SetTextJustify(5);
    topbar->AddFrame(eventNo, new TGLayoutHints(kLHintsCenterX,
                                            10,10,10,10));
    display->AddFrame(topbar, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));

    screen = new TRootEmbeddedCanvas("screen",display, 800,800);
    display->AddFrame(screen, new TGLayoutHints(kLHintsExpandX |
                   kLHintsExpandY, 10,10,10,1));

    // Create a horizontal frame widget with check buttons
    TGHorizontalFrame *cframe = new TGHorizontalFrame(display,800,100);
    simulated = new TGCheckButton(cframe,"&Simulated Hits");
    noiseEvents = new TGCheckButton(cframe,"&Noise Hits");
    time = new TGCheckButton(cframe,"&Timing-Distance Calibrated data");
    simulated->SetState(kButtonDown);
    simulated->Connect("Clicked()","EventDisplay",this,"Iterate()");
    noiseEvents->Connect("Clicked()","EventDisplay",this,"Iterate()");
    time->Connect("Clicked()","EventDisplay",this,"Iterate()");
    cframe->AddFrame(simulated, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    cframe->AddFrame(noiseEvents, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    cframe->AddFrame(time, new TGLayoutHints(kLHintsCenterX,5,5,3,4));
    display->AddFrame(cframe, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));
    
    // Create a horizontal frame widget with buttons
    TGHorizontalFrame *hframe = new TGHorizontalFrame(display,800,100);
    TGTextButton *prev = new TGTextButton(hframe,"&prev");
    prev->Connect("Clicked()","EventDisplay",this,"DrawPrev()");
    hframe->AddFrame(prev, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
    num = new TGNumberEntryField (hframe,-1,0, TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive, TGNumberFormat::kNELLimitMax, 0, t->GetEntries());	
    num->Connect("ReturnPressed()","EventDisplay",this,"Goto()");
    hframe->AddFrame(num, new TGLayoutHints(kLHintsCenterX,10,10,3,4));

    TGTextButton *next= new TGTextButton(hframe,"&next");
    next->Connect("Clicked()","EventDisplay",this,"DrawNext()");
    hframe->AddFrame(next, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
    display->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));
    display->SetWindowName("Event Display");

    //Increase the font size
    TGFontPool *pool = gClient->GetFontPool();
    TGFont *font = pool->GetFont("helvetica", 15,  kFontWeightNormal,  kFontSlantRoman);

    FontStruct_t ft = font->GetFontStruct();
    next->SetFont(ft);
    prev->SetFont(ft);
    num->SetFont(ft);
    num->SetHeight(30);

    font = pool->GetFont("helvetica", 13,  kFontWeightNormal,  kFontSlantRoman);
    ft = font->GetFontStruct();
    simulated->SetFont(ft);
    noiseEvents->SetFont(ft);
    time->SetFont(ft);

    font = pool->GetFont("helvetica", 25,  kFontWeightNormal,  kFontSlantRoman);
    ft = font->GetFontStruct();
    eventNo->SetTextFont(ft);
    eventNo->SetWrapLength(500);

    //Display the window
    display->MapSubwindows();
    display->Resize(display->GetDefaultSize());
    display->MapWindow();

    screen->GetCanvas()->Divide(2);
    display->Connect("CloseWindow()", "EventDisplay", this, "DoExit()"); 
    display->DontCallClose();
  }

  //Sets up the geometry of the 2D histogram and Initializes other graphs.
  void BuildGeometry(){

    //Initialize
    top = new TH2Poly ("top", "Top Layer", -80, 100, -200, 200);
    bot = new TH2Poly ("bot", "Bottom Layer", -80, 100, -200, 200);
    topNoise = new TH2Poly ("topNoise", "Top Layer", -80, 100, -200, 200);
    botNoise = new TH2Poly ("botNoise", "Bottom Layer", -80, 100, -200, 200);
    topTemplate = new TH2Poly ("topTemplate", "Top Layer", -80, 100, -200, 200);
    botTemplate = new TH2Poly ("botTemplate", "Bottom Layer", -80, 100, -200, 200);
    topCalc = new TPolyMarker ();
    botCalc = new TPolyMarker ();

    topCalc->SetMarkerStyle (kFullCircle);
    botCalc->SetMarkerStyle (kFullCircle);
    topCalc->SetMarkerColor (kPink);
    botCalc->SetMarkerColor (kPink);

    float xdpos;
    float zdpos;

    //Set up the geometry
    for (Int_t i = 0; i < 1176; i++)
      {
	TVectorD *geo = (TVectorD *) myGeometryData->At (i);
	xdpos = (*geo)[1] / 10.0;
	zdpos = (*geo)[3] / 10.0;
	top->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
	topNoise->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
	topTemplate->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
      }

    for (Int_t i = 1175; i < 2351; i++)
      {
	TVectorD *geo = (TVectorD *) myGeometryData->At (i);
	xdpos = (*geo)[1] / 10.0;
	zdpos = (*geo)[3] / 10.0;
	bot->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
	botNoise->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
	botTemplate->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
      }
  }

  EventDisplay(int run, int value){
    run_number = run;
    event = value;

    FileHandler();
    InitializeArrays();
    BuildApplication();
    BuildGeometry();
 
    Iterate();
  }
  void Goto(){
    	event = (int)num->GetNumber();
	Iterate();
  }
  void DrawPrev(){
	event--;
	Iterate();
  }
  void DrawNext(){
	event++;
	Iterate();
  }
  void DoExit(){
	f->Close ();
    	noisef->Close ();
        display->CloseWindow();
	delete this;
        //gApplication->Terminate();
  }
  
};

void
EventDisplay2_Nandhu (int run_number = 4000, int event = 0)
{
  EventDisplay *viewer = new EventDisplay(run_number,event);
  //viewer.main(run_number,event);
}

