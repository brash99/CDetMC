#include <iostream>
#include <TF1.h>
#include <TVectorD.h>
#include <TMatrixD.h>
#include <TRandom.h>
#include <TLinearFitter.h>
#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TList.h>
#include <TMath.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TStyle.h>
#include <TAxis.h>
#include <TLine.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TChain.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>
#include <TGraph2D.h>
#include <TGraph2DErrors.h>

class noiseMachine{
public:
//Adding Noise
TList* myGeometryData;
int global_run_number;
int Analyse_Secondaries = 1;
float Theta_min_cut = 2.524;
float ThetaVerticalCut = 3.02;
int Photon_min_cut = 75;
int Photon_min_cut_noise = 75;

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

//Adds the noise, used in a separate method to be able be called multiple times
void NoiseMaker(Int_t event, Int_t Photons[], Float_t Time[], int noiseAmount = 1){
	for(int j = 0; j < noiseAmount; j++){
		//Choose a random paddle
		int i = (int)((gRandom->Rndm())*(AnaBar_PMT_Offset+(NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS)));
		Photons[i] = gRandom->Integer(40)+75;
		Time[i] = gRandom->Uniform (3.0, 5.0);
	}
	/*
	if(event<10) {
		cout << "In NoiseMaker" << endl;
		for (int k=0; k< AnaBar_PMT_Offset+(NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS); k++) {
			if (Photons[k] > Photon_min_cut_noise) {
				std::cout << "Event = " << event << "  Photons[" << k << "] = " << Photons[k] << "   Time[" << k << "] = " << Time[k] <<endl; 
			}
		}
	}
	*/
}

void main(int run_number = 4000, int noiseAmount = 1){
	global_run_number = run_number;
	std::cout << run_number << std::endl;
	//TString fileName;
	//fileName.Form("data/AnaBarMC_%d.root",run_number);

	//Open the tree in the file
	auto fileName = "data/AnaBarMC_"+std::to_string(run_number)+".root";
	auto treeName = "T";
        
	TFile* f = new TFile((TString)fileName,"UPDATE");
        TTree* t = 0;
	f->GetObject(treeName,t);

	//Get the Photon deposition array
	Int_t PMT_Nphotons[50000];
	Int_t Detector_id[50000];
	Float_t PMT_Time[5000];
	Int_t Detector_Nhits;
	t->SetBranchAddress("PMT_Nphotons",&PMT_Nphotons);
	t->SetBranchAddress("Detector_id",&Detector_id);
	t->SetBranchAddress("Detector_Nhits",&Detector_Nhits);
	t->SetBranchAddress("PMT_Time",&PMT_Time);

	//Initialize the noise array
	//noise->Branch("PMT_id",     &fPMTNo,     "PMT_id/I" );  
	Int_t Noise_Photons[50000];
	Float_t Noise_Time[5000];
  	auto newBranch1 = t->Branch("PMT_Nphotons_Noise",  &Noise_Photons,  "PMT_Nphotons_Noise[50000]/I" );  
  	auto newBranch2 = t->Branch("PMT_Time_Noise",  &Noise_Time,  "PMT_Time_Noise[5000]/F" ); 


	//Get the values and add the noise
	for(int event = 0; event < (int)t->GetEntries(); event++){
	//for(int event = 0; event < 10; event++){

		//Get the event information.
		t->GetEntry(event);

		/*
		if(event<10) {
			cout << "Before Copy" << endl;
			for (int i=0; i< AnaBar_PMT_Offset+(NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS); i++) {
				if (PMT_Nphotons[i] > Photon_min_cut) {
					std::cout << "Event = " << event << "  PMT_Nphotons[" << i << "] = " << PMT_Nphotons[i] << "   PMT_Time[" << i << "] = " << PMT_Time[i] <<endl; 
				}
			    }
		}
		*/

		//Copy the PMT_Nphotons and PMT_Time arrays to the Noise_Photons and Noise_Time arrays
		std::copy(PMT_Nphotons, PMT_Nphotons+50000, Noise_Photons);
		std::copy(PMT_Time, PMT_Time+5000, Noise_Time);

		/*
		if(event<10) {
			cout << "After Copy" << endl;
			for (int i=0; i< AnaBar_PMT_Offset+(NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS); i++) {
				if (PMT_Nphotons[i] > Photon_min_cut_noise) {
					std::cout << "Event = " << event << "  Noise_Photons[" << i << "] = " << Noise_Photons[i] << "   Noise_Time[" << i << "] = " << Noise_Time[i] <<endl; 
				}
			    }
		}
		*/



		// add noise to the PMT_Nphotons and PMT_Time arrays
		NoiseMaker(event, Noise_Photons,Noise_Time, noiseAmount);
		newBranch1->Fill();
		newBranch2->Fill();


		/*
		if(event<10) {
			cout << "After Adding Noise" << endl;
			for (int i=0; i< AnaBar_PMT_Offset+(NUMPADDLE*NUMBARS*NUMMODULES*NUMSIDES*NUMLAYERS); i++) {
				if (Noise_Photons[i] > Photon_min_cut_noise) {
					std::cout << "Event = " << event << "  Noise_Photons[" << i << "] = " << Noise_Photons[i] << "   Noise_Time[" << i << "] = " << Noise_Time[i] <<endl; 
				}
			    }
		}
		*/
	}
	t->Write("", TObject::kOverwrite);
	f->Close();
}
};

void AddNoise(int run_number=4000, int noiseAmount = 1){
	noiseMachine speakers;
	speakers.main(run_number, noiseAmount);
}

