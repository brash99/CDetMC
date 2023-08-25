#include <iostream>
#include <TF1.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLinearFitter.h>

//Making an event display

TList* myGeometryData;
int global_run_number;
int Analyse_Secondaries = 1;
float Theta_min_cut = 2.524;
float ThetaVerticalCut = 3.02;
int Photon_min_cut = 75;

static const int MaxPMTNo = 50000;
static const int MaxPMTHits = 1000;
static const float Finger_Edep_Max = 10.0;
static const float AnaBar_Edep_Max = 10.0;
static const float pedastel_sigma = 2.9;
static const Int_t Detector_Offset = 2560;
static const int Detector_PMT_Offset = 2500;
static const int AnaBar_Offset = 30000;
static const int AnaBar_PMT_Offset = 0;

static const int Finger_NPhotons_Max = 250;
static const int AnaBar_NPhotons_Max = 200;

static const int NUMPADDLE = 14;
static const int NUMBARS = 14;
static const int NUMMODULES = 3;
static const int NUMSIDES = 2;
static const int NUMLAYERS = 2;

static const int NDET = NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;

static const int NMaxPMT = 14;

Int_t Detector_Nhits;
Int_t PMT_Nphotons[MaxPMTNo];
Int_t Detector_id[MaxPMTHits];
Float_t PMT_Time[MaxPMTNo];

int getSide(int fID) {
    int iLayer = fID%NUMPADDLE;
    int iBar = ((fID-iLayer)/NUMPADDLE)%NUMBARS;
    return (((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)%NUMSIDES;
}

int getPlane(int fID) {
    int iLayer = fID%NUMPADDLE;
    int iBar = ((fID-iLayer)/NUMPADDLE)%NUMBARS;
    int iSide = (((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)%NUMSIDES;
    int iModule =  ((((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)/NUMSIDES)%NUMMODULES;
    return (((((fID-iLayer)/NUMPADDLE-iBar)/NUMBARS)/NUMSIDES)/NUMMODULES)%NUMLAYERS;
}

float getXOffsetFromTime(int fID, float time) {

    float xoffset;
    int iSide = getSide(fID);
    int iPlane = getPlane(fID);
    TRandom3* fRand = new TRandom3(0);

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
    
    if (iPlane == 0) {
        double a = -0.00037499653844674384;
        double b = -0.051184438607694095;
        double c = 3.5929880196450843;
	double discriminant = b*b - 4*a*(c-time);
	if (discriminant >= 0.0) {
        	xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
	} else {
		xoffset -25.0;
	}
    } else {
        double a = -0.0004090552401677775;
        double b = -0.050453362706664166;
        double c = 4.537007798185197;
	double discriminant = b*b - 4*a*(c-time);
	if (discriminant >= 0.0) {
        	xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
	} else {
		xoffset -25.0;
	}
    }

    if (iSide == 1) {
        xoffset = -xoffset;
    }

    return xoffset;

}

bool getTrigger(int Detector_Nhits, int* Detector_id) {

    bool tophit = false;
    bool bottomhit = false;
    bool fhit = false;
    bool ahit = false;
    bool trigger = false;
    for (int j=0; j<Detector_Nhits; j++) {
        //std::cout << "Detector id = " << Detector_id[j] << std::endl;
        if ((Detector_id[j] == Detector_Offset || Detector_id[j] == Detector_Offset+1) && !tophit) {
            tophit = true;
            //std::cout << "Top hit" << Detector_id[j] << std::endl;
        }
        if ((Detector_id[j] == Detector_Offset+2 || Detector_id[j] == Detector_Offset+3) && !bottomhit) {
            bottomhit = true;
            //std::cout << "Bottom hit" << Detector_id[j] << std::endl;
        }
        for (int ibar=0; ibar<NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; ibar++){
            if (Detector_id[j] == AnaBar_Offset + ibar) {
                ahit = true;
            }
        }
        //if (tophit && bottomhit) {
        if (ahit) {
            fhit = true;
            trigger = true;
        }
    }

    return trigger;
}

void EventDisplay(int run_number = 4000){
	global_run_number = run_number;
	std::cout << run_number << std::endl;
	//TString fileName;
	//fileName.Form("data/AnaBarMC_%d.root",run_number);
	auto fileName = "data/AnaBarMC_"+std::to_string(run_number)+".root";
	auto treeName = "T";

        TFile* f = new TFile((TString)fileName,"READ");
        TTree* t = 0;
	f->GetObject(treeName,t);

	//Get the Photon deposition array
	t->SetBranchAddress("PMT_Nphotons",&PMT_Nphotons);
	t->SetBranchAddress("Detector_id",&Detector_id);
	t->SetBranchAddress("Detector_Nhits",&Detector_Nhits);
	t->SetBranchAddress("PMT_Time",&PMT_Time);

	myGeometryData = (TList*)t->GetUserInfo()->FindObject("myGeometryData");

	Double_t padLen = 50.2; //Length of one paddle
	Double_t padWid = 0.54; //Width of one paddle
	Double_t xdpos;
	Double_t zdpos;

	//std::vector<float> x = v[0].Take<float>("Detector_x");

	for(int event = 0; event < (int)t->GetEntries(); event++){

		//Get the event information.
		t->GetEntry(event);

		//Check if the event is valid before making the display.
		if(getTrigger(Detector_Nhits,Detector_id) == false){
			continue;
		}

                std::cout << "Supposedly good event ... " << getTrigger(Detector_Nhits,Detector_id) << " ... Detector_Nhits = " << Detector_Nhits << " ... Dtector_id = " << Detector_id << std::endl;

	    	auto top = new TH2Poly("top", "Top Layer", -80,100,-200,200);
		auto bot = new TH2Poly("bot", "Bottom Layer", -80,100,-200,200);
		auto topCalc = new TPolyMarker();
		auto botCalc = new TPolyMarker();

		topCalc->SetMarkerStyle(kFullCircle);
		botCalc->SetMarkerStyle(kFullCircle);

		//Set up the geometry
		for(Int_t i = 0; i < 1176; i++){
			TVectorD* geo = (TVectorD*)myGeometryData->At(i);
			xdpos = (*geo)[1]/10.0;
        		zdpos = (*geo)[3]/10.0;
			top->AddBin(xdpos-(padLen/2), zdpos-(padWid/2), xdpos+(padLen/2), zdpos+(padWid/2));
		}

		for(Int_t i = 1175; i < 2351; i++){
			TVectorD* geo = (TVectorD*)myGeometryData->At(i);
			xdpos = (*geo)[1]/10.0;
        		zdpos = (*geo)[3]/10.0;
			bot->AddBin(xdpos-(padLen/2), zdpos-(padWid/2), xdpos+(padLen/2), zdpos+(padWid/2));
		}
			
		//Fill the bins
		for(int i = 0; i < AnaBar_PMT_Offset+NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS; i++){
			if(i < 1176){
				top->SetBinContent(i,PMT_Nphotons[i]);
				
				//Draw the calculated x position	
				if(PMT_Nphotons[i] > Photon_min_cut){ //Only if there are enough photons to justify it
					TVectorD* geo = (TVectorD*)myGeometryData->At(i);
					xdpos = (*geo)[1]/10.0;
        				zdpos = (*geo)[3]/10.0;
					topCalc->SetNextPoint(xdpos + getXOffsetFromTime(i, PMT_Time[i]),zdpos);
                                        std::cout << "Paddle: " << i << " ... x = " << xdpos << " ... offset = " << getXOffsetFromTime(i,PMT_Time[i]) << " ... time = " << PMT_Time[i] << " ... zdpos = " << zdpos << std::endl;
				}
			} else {
				bot->SetBinContent(i-1175,PMT_Nphotons[i]);
				
				//Draw the calculated x position	
				if(PMT_Nphotons[i] > Photon_min_cut){ //Only if there are enough photons to justify it
					TVectorD* geo = (TVectorD*)myGeometryData->At(i);
					xdpos = (*geo)[1]/10.0;
        				zdpos = (*geo)[3]/10.0;
				        botCalc->SetNextPoint(xdpos + getXOffsetFromTime(i, PMT_Time[i]),zdpos);
				}
			}
		}

		std::string title = "Event Display " + std::to_string(event);

		TCanvas *eDisplay = new TCanvas("eDisplay",title.c_str(),800,800);
		eDisplay->Divide(2);

		eDisplay->cd(1);		
		top->Draw("COL");
		topCalc->Draw("SAME E");

		eDisplay->cd(2);
		bot->Draw("COL");	
		botCalc->Draw("SAME E");

		//Remove the stats box
		top->SetStats(0);
		bot->SetStats(0);

		eDisplay->Draw();
		eDisplay->WaitPrimitive();

		//Cleanup
		delete gROOT->FindObject("eDisplay");
		delete gROOT->FindObject("top");
		delete gROOT->FindObject("bot");
	}
}
