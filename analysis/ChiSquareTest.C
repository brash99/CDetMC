#include <iostream>
#include <TF1.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TLinearFitter.h>
#include <TRandom3.h>
#include <TFile.h>
#include <TTree.h>

void print(std::string name, std::vector <int> const &a) {
   std::cout << "The elements of " << name << " are : ";

   for(int i=0; i < a.size(); i++)
   std::cout << a.at(i) << ' ';

   std::cout << '\n';
}


TList* myGeometryData;
int global_run_number;
int Analyse_Secondaries = 1;
float Theta_min_cut = 2.524;
float ThetaVerticalCut = 3.02;
int Photon_min_cut = 75;
int events_to_see = 0;

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

const int NDET = NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS;

int NMaxPMT = 14;

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
		xoffset = -25.0;
	}
    } else {
        double a = -0.0004090552401677775;
        double b = -0.050453362706664166;
        double c = 4.537007798185197;
	double discriminant = b*b - 4*a*(c-time);
	if (discriminant >= 0.0) {
        	xoffset = (-b-sqrt(fabs(b*b-4*a*(c-time))))/(2.0*a);
	} else {
		xoffset = -25.0;
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

//Noise removal algorithm
std::vector<int> RemoveNoise(int event, std::vector<int> layer1, std::vector<int>layer2, float distanceAcceptable = 5.0){
	std::vector<int> output;

        if (event<events_to_see) {
            cout << "layer1.size() = " << layer1.size() << ", layer2.size() = " << layer2.size() << endl;
        }   
                        
        if (event<events_to_see) {
	  for(int i = 0; i < layer1.size(); i++){
		for(int j = 0; j < layer2.size(); j++){
			TVectorD* geo = (TVectorD*)myGeometryData->At(layer1.at(i));
			TVectorD* geo2 = (TVectorD*)myGeometryData->At(layer2.at(j));
			float x1 = (*geo)[1]/10.0;
			float x2 = (*geo2)[1]/10.0;
        		float z1 = (*geo)[3]/10.0;
			float z2 = (*geo2)[3]/10.0;
                        cout << "i = " << i << ", j = " << j << ", x1 = " << x1 << ", x2 = " << x2 << ", z1 = " << z1 << ", z2 = " << z2 << ", dist = " <<  pow(pow((z1-z2),2),0.5) << endl;
                }
          }
        }

	for(int i = 0; i < layer1.size(); i++){
		for(int j = 0; j < layer2.size(); j++){
			TVectorD* geo = (TVectorD*)myGeometryData->At(layer1.at(i));
			TVectorD* geo2 = (TVectorD*)myGeometryData->At(layer2.at(j));
			float x1 = (*geo)[1]/10.0;
			float x2 = (*geo2)[1]/10.0;
        		float z1 = (*geo)[3]/10.0;
			float z2 = (*geo2)[3]/10.0;
			//if ((pow(x1-x2, 2) + pow(z1-z2, 2)) < distanceAcceptable){

			if (pow(pow((z1-z2),2),0.5) < distanceAcceptable){
				if (event < events_to_see) {
					cout << "Pushing back i = " << i << " at x1 = " << x1 << ", z1 = " << z1 << endl;
				}
				output.push_back(layer1.at(i));
				break;
			}
		}
	}
	return output; 
}

void ChiSquareTest(int run_number = 4000){

	global_run_number = run_number;
	std::cout << run_number << std::endl;
	//TString fileName;
	//fileName.Form("data/AnaBarMC_%d.root",run_number);
	auto fileName = "data/AnaBarMC_"+std::to_string(run_number)+".root";
	auto treeName = "T";

        TFile* f = new TFile((TString)fileName,"READ");
        TTree* t = 0;
	f->GetObject(treeName,t);

	//auto noisefileName = "data/AnaBarMC_"+std::to_string(run_number)+"noise.root";
	//auto noisetreeName = "noise";

        //TFile* noisef = new TFile((TString)noisefileName,"READ");
        //TTree* noise = 0;
	//noisef->GetObject(noisetreeName,noise);
	//t->AddFriend(noise);

	//Get the Photon deposition array
	Int_t PMT_Nphotons[50000];
	Int_t PMT_Nphotons_Noise[50000];
	Int_t Detector_id[50000];
	Float_t PMT_Time[5000];
	Float_t PMT_Time_Noise[5000];
	Int_t Detector_Nhits;
	t->SetBranchAddress("PMT_Nphotons",&PMT_Nphotons);
	t->SetBranchAddress("PMT_Nphotons_Noise",&PMT_Nphotons_Noise);
	t->SetBranchAddress("Detector_id",&Detector_id);
	t->SetBranchAddress("Detector_Nhits",&Detector_Nhits);
	t->SetBranchAddress("PMT_Time",&PMT_Time);
	t->SetBranchAddress("PMT_Time_Noise",&PMT_Time_Noise);

	myGeometryData = (TList*)t->GetUserInfo()->FindObject("myGeometryData");

	Double_t padLen = 50.2; //Length of one paddle
	Double_t padWid = 0.54; //Width of one paddle
	Double_t xdpos;
	Double_t zdpos;

	//std::vector<float> x = v[0].Take<float>("Detector_x");

	//Get some data
	int sumTotalNoise = 0;
	int sumUnidentifiedNoise = 0;
	int sumMisidentified = 0;

	for(int event = 0; event < (int)t->GetEntries(); event++){


		//Get the event information.
		t->GetEntry(event);

		//Check if the event is valid before making the display.
		if(getTrigger(Detector_Nhits,Detector_id) == false){
			continue;
		}

		//Initialize
		std::vector<int> paddles; //Triggered paddles from the control
		std::vector<int> noisePaddles; //Input paddles from the noisy data
		std::vector<int> outPaddles; //Cleaned up noise paddles

		std::vector<int> top;
		std::vector<int> bot;

		//Fill with values- only the triggered values
		for(int i=0; i < 2352; i++){
			if(PMT_Nphotons[i] > 0){
				paddles.push_back(i);
				//cout << std::to_string(PMT_Nphotons[i]) + " ";
			}
			if(PMT_Nphotons_Noise[i] > 0){
				noisePaddles.push_back(i);
			}
		}
		//cout << "\n";

		

		//Sort the noisePaddles into top and bottom layers
		for(int i = 0; i < noisePaddles.size();i++){
			if(noisePaddles.at(i) > 1176){
				top.push_back(noisePaddles.at(i));			
			} else{
				bot.push_back(noisePaddles.at(i));
			}			
		}

		//Finish our output values
                if (event<events_to_see) {
                    cout << "top/bottom" << endl;
                }
		top = RemoveNoise(event, top,bot);
                if (event<events_to_see) {
                    cout << "bottom/top" << endl;
                }
		bot = RemoveNoise(event, bot,top);
		outPaddles.insert(outPaddles.end(),top.begin(),top.end());
		outPaddles.insert(outPaddles.end(),bot.begin(),bot.end());

		//Compare our values to what we should have gotten
		int totalNoise = 0;
		int correctNoise = 0;
		int wrongNoise = 0;
		int unIdentified = 0;

		std::vector<int> identifiedNoise;
		//std::vector<int> unIdentifiedNoise;

		//Sort for ease
		if (event < events_to_see) {
			print("paddles",paddles);
			print("noisePaddles",noisePaddles);
			print("outPaddles",outPaddles);
		}
		std::sort(paddles.begin(),paddles.end());
		std::sort(noisePaddles.begin(),noisePaddles.end());
		std::sort(outPaddles.begin(),outPaddles.end()); 

		std::set_difference(noisePaddles.begin(), noisePaddles.end(), outPaddles.begin(), outPaddles.end(), std::inserter(identifiedNoise, std::end(identifiedNoise)));
		//std::set_difference(outPaddles.begin(),outPaddles.end(),paddles.begin(),paddles.end(),unIdentifiedNoise.begin());
		if (event < events_to_see) {
			print("identified noise",identifiedNoise);
		}
		totalNoise = (noisePaddles.size()-paddles.size());	
	
		for(int i = 0; i < identifiedNoise.size(); i++){
			if (event < events_to_see) {
				cout << "PMT_Nphotons[identifiedNoise.at(i)] " << PMT_Nphotons[identifiedNoise.at(i)] << endl;
			}
			if(PMT_Nphotons[identifiedNoise.at(i)] == 0){
				correctNoise++;
			} else {
				wrongNoise++;
			}
		}

		unIdentified = totalNoise - correctNoise;	

                if (event  < events_to_see) {
		cout << "For Event " + std::to_string(event) + " - Total noise: " + std::to_string(totalNoise) + ", Indentified Noise: " + std::to_string(identifiedNoise.size()) + ", Correctly indentified: " + std::to_string(correctNoise) + ", Misindentified: " + std::to_string(wrongNoise) + ", Unindentified: " + std::to_string(unIdentified) + "\n"; 		
		}

		sumTotalNoise+= totalNoise;
		sumUnidentifiedNoise+=unIdentified;
		sumMisidentified+=wrongNoise;

		/*cout << "noisePaddles ";		
		for (int n : noisePaddles)
        		std::cout << n << ' ';
    		std::cout << '\n';

		cout << "outPaddles ";
		for (int n : outPaddles)
        		std::cout << n << ' ';
    		std::cout << '\n'; 

		cout << "paddles ";
		for (int n : paddles)
        		std::cout << n << ' ';
    		std::cout << '\n'; 		

		cout << "identifiedNoise ";
		for (int n : identifiedNoise)
        		std::cout << n << ' ';
    		std::cout << '\n';*/

	}
	cout << "Total Noise for run: " + std::to_string(sumTotalNoise) + ", Unidentified Noise: " + std::to_string(sumUnidentifiedNoise) + ", Misidentified Noise: " + std::to_string(sumMisidentified) + "\n";

	cout << "Misidentification Likelihood: " + std::to_string((float)sumMisidentified/(float)sumTotalNoise*100.0) + "%\nUnidentification Likelihood: " + std::to_string((float)sumUnidentifiedNoise/(float)sumTotalNoise*100.0) + "%\n";
}
