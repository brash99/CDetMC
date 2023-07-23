#include "TH1.h"
#include "TF1.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TRandom3.h"

#include <iostream>

using namespace std;

static const int MaxHits = 50000;
static const int MaxPMTNo = 20;
static const int MaxPMTHits = 5000;
static const Float_t Finger_Edep_Max = 10.0;
static const Float_t AnaBar_Edep_Max = 10.0;
static const Float_t pedastel_sigma = 3.33;
static const Int_t Detector_Offset = 0;
//static const Int_t Finger_NPhotons_Max = 150;
//static const Int_t AnaBar_NPhotons_Max = 100;
static const Int_t Finger_NPhotons_Max = 250;
static const Int_t AnaBar_NPhotons_Max = 200;

static const Int_t NUMPADDLE = 14;

TTree *tree1;

Float_t Prim_E, Prim_Th, Prim_Ph;
Int_t Prim_pdg;
Int_t Detector_Nhits;
Int_t Detector_pdg[MaxHits];
Int_t Detector_id[MaxHits];
Float_t Detector_x[MaxHits], Detector_y[MaxHits], Detector_z[MaxHits], Detector_t[MaxHits];
Float_t Detector_Ed[MaxHits];  
Int_t PMT_id;
Int_t PMT_Nphotons[MaxPMTNo];
Float_t PMT_Nphotons_Noise[MaxPMTNo];
Int_t PMT_Nphotons_Total;
Float_t PMT_KineticEnergy[MaxPMTNo][MaxPMTHits];

TRandom3* fRand = new TRandom3(-1);



void AnalyseSignalsTrigger(Int_t Analysis_Run_Number = 8881) {

  //-------------------------------------------------------------------
  //Set stuff up for reading
  //-------------------------------------------------------------------
  TString filename;
  filename.Form("data/AnaBarMC_%d.root",Analysis_Run_Number);
  TFile *f1 = new TFile(filename,"READ");
  tree1 = (TTree*)f1->Get("T");

  tree1->SetBranchAddress("Prim_E", &Prim_E);
  tree1->SetBranchAddress("Prim_Th", &Prim_Th);
  tree1->SetBranchAddress("Prim_Ph", &Prim_Ph);
  tree1->SetBranchAddress("Prim_pdg", &Prim_pdg);
  tree1->SetBranchAddress("PMT_id", &PMT_id);
  tree1->SetBranchAddress("PMT_Nphotons", &PMT_Nphotons);
  tree1->SetBranchAddress("PMT_KineticEnergy", &PMT_KineticEnergy);
  tree1->SetBranchAddress("Detector_Nhits", &Detector_Nhits);
  tree1->SetBranchAddress("Detector_pdg", &Detector_pdg);
  tree1->SetBranchAddress("Detector_id", &Detector_id);
  tree1->SetBranchAddress("Detector_x", &Detector_x);
  tree1->SetBranchAddress("Detector_y", &Detector_y);
  tree1->SetBranchAddress("Detector_z", &Detector_z);
  tree1->SetBranchAddress("Detector_t", &Detector_t);
  tree1->SetBranchAddress("Detector_Ed", &Detector_Ed);

  filename.Form("data/AnaBarMC_%d_trigger.root",Analysis_Run_Number);
  TFile *triggerFile = new TFile(filename, "RECREATE");

  TTree *skimTree = tree1->CloneTree(0);


  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 
//  float edeptot[NMaxPMT]; 


//***************** Event loop *****************//

  for (Int_t i = 0; i < nentries; i++) {

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
	    //edeptot[j] = 0.0; 
    	    anabar_hit_paddle[j]=false; 
    }
    tree1->GetEntry(i); 

//    Float_t fMass = 105.70; 
//    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
//    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph);
//    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 

    bool trigger = false;
    bool finger_hit = false;
    bool anabar_hit = false;
    bool anabar_top_hit = false;// unused
    bool anabar_bottom_hit = false; // unused
    int j_finger = 0;
    int j_anabar = 0;
    for (Int_t j=0; j < Detector_Nhits ; j++) {
	//cout << "Detector hit = " << j << " Detector_id[j] = " << Detector_id[j] << endl;
	if (Detector_id[j] == Detector_Offset && !finger_hit) {
		finger_hit = true;
		j_finger = j;
		//cout<<"hit in finger";
	}
	for (Int_t ibar = 1; ibar<15; ibar++){
		if (Detector_id[j+Detector_Offset] == ibar+Detector_Offset) {
		  anabar_hit = true;
		  anabar_hit_paddle[ibar-1]=true;
		  j_anabar = j;
		  //cout << "hit in anabar " << j << endl;
		}
	}
	//if (Detector_id[j] == 14 && !anabar_bottom_hit) {
	//	anabar_bottom_hit = true;
	//	j_anabar = j;
	//}
    }

    //if (finger_hit && anabar_top_hit && anabar_bottom_hit) trigger = true; 
    //if (finger_hit && anabar_top_hit) trigger = true; 
    if (finger_hit && anabar_hit) trigger = true; 

    if (trigger)
	skimTree->Fill();

  }

  triggerFile->cd("T");
  skimTree->Write();


}

