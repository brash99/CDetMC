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

  TGMainFrame *display;
  TRootEmbeddedCanvas *screen;
  TGLabel *eventNo;
  int event = 0;
  int run_number = 4000; 

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

    cout << "========================" << endl;
    cout << "Event" << event << endl;
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


  void main ()
  {
    //global_run_number = run_number;
    std::cout << run_number << std::endl;
    std::string title = "Event " + std::to_string (event);
    eventNo->ChangeText(title.c_str());
    //TString fileName;
    //fileName.Form("data/AnaBarMC_%d.root",run_number);
    auto fileName = "data/AnaBarMC_" + std::to_string (run_number) + ".root";
    auto treeName = "T";

    TFile *f = new TFile ((TString) fileName, "READ");
    TTree *t = 0;
    f->GetObject (treeName, t);

    auto noisefileName =
      "data/AnaBarMC_" + std::to_string (run_number) + "noise.root";
    auto noisetreeName = "noise";

    TFile *noisef = new TFile ((TString) noisefileName, "READ");
    TTree *noise = 0;
    noisef->GetObject (noisetreeName, noise);
    //t->AddFriend (noise); 
    //t->Print();
    myGeometryData = (TList *) t->GetUserInfo ()->FindObject ("myGeometryData");

    //Get the Photon deposition array
    Int_t PMT_Nphotons[50000];
    Int_t Detector_id[50000];
    Float_t PMT_Time[5000];
    Int_t Detector_Nhits;
    //t->SetBranchAddress ("PMT_Nphotons_Noise", &PMT_Nphotons);
    t->SetBranchAddress ("PMT_Nphotons_Noise", &PMT_Nphotons);
    t->SetBranchAddress ("Detector_id", &Detector_id);
    t->SetBranchAddress ("Detector_Nhits", &Detector_Nhits);
    t->SetBranchAddress ("PMT_Time_Noise", &PMT_Time);

    Double_t padLen = 50.2;	//Length of one paddle
    Double_t padWid = 0.54;	//Width of one paddle
    Double_t xdpos;
    Double_t zdpos;

    //std::vector<float> x = v[0].Take<float>("Detector_x");

    //Get the event information.
    t->GetEntry (event);

    //Check if the event is valid before making the display.
    if (getTrigger (Detector_Nhits, Detector_id) == false)
      {
	event++;
	t->GetEntry (event);
      }

    auto top = new TH2Poly ("top", "Top Layer", -80, 100, -200, 200);
    auto bot = new TH2Poly ("bot", "Bottom Layer", -80, 100, -200, 200);
    auto topCalc = new TPolyMarker ();
    auto botCalc = new TPolyMarker ();

    topCalc->SetMarkerStyle (kFullCircle);
    botCalc->SetMarkerStyle (kFullCircle);

    //Set up the geometry
    for (Int_t i = 0; i < 1176; i++)
      {
	TVectorD *geo = (TVectorD *) myGeometryData->At (i);
	xdpos = (*geo)[1] / 10.0;
	zdpos = (*geo)[3] / 10.0;
	top->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
      }

    for (Int_t i = 1175; i < 2351; i++)
      {
	TVectorD *geo = (TVectorD *) myGeometryData->At (i);
	xdpos = (*geo)[1] / 10.0;
	zdpos = (*geo)[3] / 10.0;
	bot->AddBin (xdpos - (padLen / 2), zdpos - (padWid / 2),
		     xdpos + (padLen / 2), zdpos + (padWid / 2));
      }


    //Fill the bins
    for (int i = 0;
	 i <
	 AnaBar_PMT_Offset +
	 NUMPADDLE * NUMBARS * NUMMODULES * NUMSIDES * NUMLAYERS; i++)
      {
	if (i < 1176)
	  {
	    

	    //Draw the calculated x position        
	    if (PMT_Nphotons[i] > Photon_min_cut)
	      {			//Only if there are enough photons to justify it
	        top->SetBinContent (i, PMT_Nphotons[i]);
		TVectorD *geo = (TVectorD *) myGeometryData->At (i);
		xdpos = (*geo)[1] / 10.0;
		zdpos = (*geo)[3] / 10.0;
		//cout << getXOffsetFromTime(i, PMT_Time[i]);
		//topCalc->SetNextPoint(xdpos + getXOffsetFromTime(i, PMT_Time[i]) - getXOffsetFromTime(i, PMT_Time[i]),zdpos);
		cout << "Top" << endl;
		cout << i << endl;
		cout << PMT_Time[i] << endl;
		topCalc->SetNextPoint (xdpos + 
				       getXOffsetFromTime (i, PMT_Time[i]),
				       zdpos);
	      }
	  }
	else
	  {

	    //Draw the calculated x position        
	    if (PMT_Nphotons[i] > Photon_min_cut)
	      {			//Only if there are enough photons to justify it
	        bot->SetBinContent (i - 1175, PMT_Nphotons[i]);
		TVectorD *geo = (TVectorD *) myGeometryData->At (i);
		xdpos = (*geo)[1] / 10.0;
		zdpos = (*geo)[3] / 10.0;
		cout << "Bottom" << endl;
		cout << i << endl;
		cout << PMT_Time[i] << endl;
		//botCalc->SetNextPoint (xdpos + 
		//		       getXOffsetFromTime (i,
		//					   PMT_Time[i]) -
		//		       getXOffsetFromTime (i, PMT_Time[i]),
	        //			       zdpos);
		botCalc->SetNextPoint(xdpos + getXOffsetFromTime(i, PMT_Time[i]),zdpos);
	      }
	  }
      }

    auto eDisplay = screen->GetCanvas();

    //eDisplay->cd();
    //auto eDisplay = new TCanvas("eDisplay",title.c_str(),800,900);
    //eDisplay->Divide(2);
    /*TPad *topside = new TPad ("topside", "Top Layer", 0.0, 0.1, 0.5, 1.0);
    TPad *botside = new TPad ("botside", "Bottom Layer", 0.5, 0.1, 1, 1.0);
    topside->Draw ();
    botside->Draw ();*/
    eDisplay->cd(1);
    //topside->cd ();
    top->Draw ("COL");
    topCalc->Draw ("SAME *");

    eDisplay->cd(2);
    //botside->cd ();
    bot->Draw ("COL");
    botCalc->Draw ("SAME E");

    //Remove the stats box
    top->SetStats (0);
    bot->SetStats (0);

    //Allow interaction
//    eDisplay->Draw ();
    eDisplay->cd();
    eDisplay->Update(); 
    f->Close ();
    noisef->Close ();
    cout << "banana";
  }

  EventDisplay(int run, int value){
    run_number = run;
    event = value;

    display = new TGMainFrame(gClient->GetRoot(),800,800);

    TGHorizontalFrame *topbar = new TGHorizontalFrame(display,800,100);
    eventNo = new TGLabel(topbar, "Title Text");
    topbar->AddFrame(eventNo, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
    display->AddFrame(topbar, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));

    screen = new TRootEmbeddedCanvas("screen",display, 800,800);
    display->AddFrame(screen, new TGLayoutHints(kLHintsExpandX |
                   kLHintsExpandY, 10,10,10,1));
    
    // Create a horizontal frame widget with buttons
    TGHorizontalFrame *hframe = new TGHorizontalFrame(display,800,100);
    TGTextButton *prev = new TGTextButton(hframe,"&prev");
    prev->Connect("Clicked()","EventDisplay",this,"DrawPrev()");
    hframe->AddFrame(prev, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));

    TGTextButton *next= new TGTextButton(hframe,"&next");
    next->Connect("Clicked()","EventDisplay",this,"DrawNext()");
    hframe->AddFrame(next, new TGLayoutHints(kLHintsCenterX,
                                            5,5,3,4));
    display->AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,
                                             2,2,2,2));
    display->SetWindowName("Event Display");
    display->MapSubwindows();
    display->Resize(display->GetDefaultSize());
    display->MapWindow();

    screen->GetCanvas()->Divide(2);
 
    main();
  }
  void DrawPrev(){
	event--;
	main();
  }
  void DrawNext(){
	event++;
	main();
  }
  
};

void
EventDisplay2Noise (int run_number = 4000, int event = 0)
{
  EventDisplay *viewer = new EventDisplay(run_number,event);
  //viewer.main(run_number,event);
}

