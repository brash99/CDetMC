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
static const Float_t pedastel_sigma = 2.9;
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

// function prototypes to hopefully keep this thing working until I finish reformating it
Double_t langaufun(Double_t *x, Double_t *par);
TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);



void AnalyseSignals(Int_t Analysis_Run_Number = 9996) {

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

// end main method here; put everything else in individual TCanvas methods  //
//**************************************************************************//
}


//****************************************************//
//************* New canvas methods here **************//

TCanvas *plotC1(){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------
  // canvas c1
  TH1F *hFingerX = new TH1F("FingerX","Finger X Position", 100, -120, 120);
  TH1F *hFingerY = new TH1F("FingerY","Finger Y Position", 100, 30, 80);
  TH1F *hFingerZ = new TH1F("FingerZ","Finger Z Position", 100, -140, 60);
  TH1F *hFingerT = new TH1F("FingerT","Finger Time", 100, 0, .4);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 

  for (Int_t i = 0; i < nentries; i++) { 
  //for (Int_t i = 0; i < 3; i++) { 

    //cout << "Event Number: " << i << endl;
    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
    	    anabar_hit_paddle[j]=false; 
    }
    tree1->GetEntry(i); 

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
		//cout<<"hit in finger\n";
	}
	for (Int_t ibar = 1; ibar<15; ibar++){
                //cout << Detector_id[j+Detector_Offset] << " " << (ibar+Detector_Offset) << endl;
		if (Detector_id[j+Detector_Offset] == ibar+Detector_Offset) {
		  anabar_hit = true;
		  anabar_hit_paddle[ibar-1]=true;
		  j_anabar = j;
		  //cout << "hit in anabar" << j << endl;
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

    for (Int_t j=0; j < Detector_Nhits ; j++) {    

	if (trigger) {
                //cout << "Trigger" << endl;
		
		if (Detector_id[j] == Detector_Offset && Detector_pdg[j] == Prim_pdg) {
		        //cout << "Correct id and pdg type" << endl;
                        counter++; // unused
			hFingerX->Fill(Detector_x[j]);
        		hFingerY->Fill(Detector_y[j]);
        		hFingerZ->Fill(Detector_z[j]);
    			hFingerT->Fill(Detector_t[j]);
		}
	}
    }
    //cout << "Counter = " << counter << endl;

  }

  TCanvas *c1 = new TCanvas("c1", "c1", 100,100,500,270);
  c1->Divide(2,2, 0.01, 0.01, 0);

  c1->cd(1);
  hFingerX->Draw();
  c1->cd(2);
  hFingerY->Draw();
  c1->cd(3);
  hFingerZ->Draw();
  c1->cd(4);
  hFingerT->Draw();

  return c1;

}


TCanvas *plotC2 (){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------
  // canvas c2
  TH1F *hPrimE = new TH1F("PrimE","Primary Energy", 100, 0, 25000);
  TH1F *hPrimTh = new TH1F("PrimTh","Primary Theta", 100, 0, TMath::Pi());
  TH1F *hPrimPh = new TH1F("PrimPh","Primary Phi", 100, 0, 2.0*TMath::Pi());
  TH1F *hPrimPdg = new TH1F("PrimPdg","Primary PDG ID", 20, 0, 20);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries(); 

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 

  for (Int_t i = 0; i < nentries; i++) { 

    tree1->GetEntry(i); 

    Float_t fMass;
    if (Prim_pdg == 11) {
	fMass = 0.511;
    } else {
	if (Prim_pdg == 13) {
	    fMass = 105.7;
	} else {
	    if (Prim_pdg == 2212) {
		fMass = 938.28;
	    } else {
		fMass = 939.5;
	    }
        }
    } 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPx        = fMomentum * TMath::Sin(Prim_Th) * TMath::Cos(Prim_Ph); // NOTE: Prim_Th and Prim_Ph are not the theta and phi of the CDet
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph); // coordinate system. They are the traditional spherical theta and 
    Float_t fPz        = fMomentum * TMath::Cos(Prim_Th);			// phi, with theta measured from the +z axis etc. 
    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 
    Float_t fNewPhi;
    if (fPx < 0)
	fNewPhi = TMath::ATan(fPz/fPx) + TMath::Pi();
    else if (fPx > 0 && fPz < 0)
	fNewPhi = TMath::ATan(fPz/fPx) + TMath::TwoPi();
    else
	fNewPhi = TMath::ATan(fPz/fPx);

    hPrimE->Fill(Prim_E);
    hPrimTh->Fill(fNewTheta);
    hPrimPh->Fill(fNewPhi);
    hPrimPdg->Fill(Prim_pdg);

  }

  TCanvas *c2 = new TCanvas("c2", "c2", 600,100,500,270);
  c2->Divide(2,2, 0.01, 0.01, 0);

  c2->cd(1);
  hPrimE->Draw();
  c2->cd(2);
  hPrimTh->Draw();
  c2->cd(3);
  hPrimPh->Draw();
  c2->cd(4);
  hPrimPdg->Draw();

  return c2;

}


TCanvas *plotC3 (){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------
  // canvas c3
  TH1F *hDetectorNhits = new TH1F("DetectorNhits","Detector Number of Hits", 100, 0, 400);
  TH1F *hDetectorPdg = new TH1F("DetectorPdg","Detector PDG ID", 50, -20, 30);
  TH1F *hDetectorID = new TH1F("DetectorID","Detector ID Number", 30, 0, 30);
  TH1F *hPMTID = new TH1F("PMTID","PMT ID Number", 15, 0, 15);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries(); 

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 

  for (Int_t i = 0; i < nentries; i++) { 

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
    	    anabar_hit_paddle[j]=false; 
    }
    tree1->GetEntry(i); 

    hPMTID->Fill(PMT_id);
    hDetectorNhits->Fill(Detector_Nhits);

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

    for (Int_t j=0; j < Detector_Nhits ; j++) {

	if (trigger) {
		
		counter++; // unused
		if (Detector_id[j] == Detector_Offset ) {
        		hDetectorPdg->Fill(Detector_pdg[j]);
        		hDetectorID->Fill(Detector_id[j]);
		}
		if (Detector_id[j] == 1 + Detector_Offset ) {
        		hDetectorPdg->Fill(Detector_pdg[j]);
        		hDetectorID->Fill(Detector_id[j]);
		}
	}
    }

  }

  TCanvas *c3 = new TCanvas("c3", "c3", 1100,100,500,270);
  c3->Divide(2,2, 0.01, 0.01, 0);

  c3->cd(1);
  hDetectorNhits->Draw();
  c3->cd(2);
  gPad->SetLogy();
  hDetectorPdg->Draw();
  c3->cd(3);
  hDetectorID->Draw();
  c3->cd(4);
  hPMTID->Draw();

  return c3;

}


TCanvas *plotC4 (Float_t Theta_min_cut = 0.0, Int_t Analyse_Secondaries = 1){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------

  TH1F *hFingerPMTNphot = new TH1F("FingerPMTNphot","Finger PMT Number of Photons", Finger_NPhotons_Max+10, -10, Finger_NPhotons_Max);
  TH1F *hFingerEd = new TH1F("FingerEd","Finger Energy Deposited", 100, 0.01, Finger_Edep_Max);

  TH1F *hAnaBarPMTNphot[NUMPADDLE];
  TString name, title;
  for(Int_t i = 1; i <= NUMPADDLE; i++){
	name.Form("AnaBarPMTNphotA%d", i);
	title.Form("AnaBar PMT Number of Photons A%d", i);
	hAnaBarPMTNphot[i-1] = new TH1F(name, title, AnaBar_NPhotons_Max*0.9+20, -20, AnaBar_NPhotons_Max*0.9);
  }

  TH1F *hAnaBarEd = new TH1F("AnaBarEd","AnaBar Energy Deposited", 100, 0.01, AnaBar_Edep_Max);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0; // unused

  const int NMaxPMT=14;
  float edeptot[NMaxPMT];

  for (Int_t i = 0; i < nentries; i++) {

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
	    edeptot[j] = 0.0;
    	    anabar_hit_paddle[j]=false;
    }
    float edep0tot = 0.0;
    tree1->GetEntry(i); 

    Float_t fMass;
    if (Prim_pdg == 11) {
	fMass = 0.511;
    } else {
	if (Prim_pdg == 13) {
	    fMass = 105.7;
	} else {
	    if (Prim_pdg == 2212) {
		fMass = 938.28;
	    } else {
		fMass = 939.5;
	    }
        }
    } 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph);
    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 

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

    if (trigger) {
	
	for (Int_t icount = 0;icount<15;icount++){
		PMT_Nphotons_Noise[icount]=PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
	}
	for (Int_t i = 0; i < NUMPADDLE; i++){
		hAnaBarPMTNphot[i]->Fill(PMT_Nphotons_Noise[i]);
	}

        hFingerPMTNphot->Fill(PMT_Nphotons_Noise[14]);
    }

    for (Int_t j=0; j < Detector_Nhits ; j++) {

	if (trigger) {
		
		counter++; // unused
		if (Detector_id[j] == Detector_Offset ) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
			  edep0tot += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
				edep0tot += Detector_Ed[j];
			       }
			} 
		}
		if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
				edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
					edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
		     	       }
			}
		}
	}
    }

    if (trigger) {
    	hFingerEd->Fill(edep0tot); 
    	hAnaBarEd->Fill(edeptot[6]); 
    }

  }

  TCanvas *c4 = new TCanvas("c4", "c4", 1100,1000,500,270);
  c4->Divide(2,2, 0.01, 0.01, 0);
  
  c4->cd(1);
  hFingerEd->Draw();
  c4->cd(2);
  hFingerPMTNphot->Draw();
  c4->cd(3);
  hAnaBarEd->Draw();
  c4->cd(4);
  hAnaBarPMTNphot[0]->Draw();

  return c4;

}


TCanvas *plotC5 (){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------
  // canvas c5
  TH1F *hAnaBarX = new TH1F("AnaBarX","AnaBar X Position", 100, -120, 120);
  TH1F *hAnaBarY = new TH1F("AnaBarY","AnaBar Y Position", 100, -30, 30);
  TH1F *hAnaBarZ = new TH1F("AnaBarZ","AnaBar Z Position", 100, -30, 30);
  TH1F *hAnaBarT = new TH1F("AnaBarT","AnaBar Time", 100, 0, .4);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries(); 

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 

  for (Int_t i = 0; i < nentries; i++) { 

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
    	    anabar_hit_paddle[j]=false; 
    }
    tree1->GetEntry(i); 

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

    for (Int_t j=0; j < Detector_Nhits ; j++) {

	if (trigger) {
		
		counter++; // unused
		if (Detector_id[j] == 1 + Detector_Offset && Detector_pdg[j] == Prim_pdg) {
        		hAnaBarX->Fill(Detector_x[j]);
        		hAnaBarY->Fill(Detector_y[j]);
        		hAnaBarZ->Fill(Detector_z[j]);
    			hAnaBarT->Fill(Detector_t[j]);
		}
	}
    }

  }

  TCanvas *c5 = new TCanvas("c5", "c5", 100,400,500,270);
  c5->Divide(2,2, 0.01, 0.01, 0);

  c5->cd(1);
  hAnaBarX->Draw();
  c5->cd(2);
  hAnaBarY->Draw();
  c5->cd(3);
  hAnaBarZ->Draw();
  c5->cd(4);
  hAnaBarT->Draw();

  return c5;

}


TCanvas *plotC6 (Float_t Theta_min_cut = 0.0, Int_t Analyse_Secondaries = 1){

  // canvas c6
  TH2F *hE1vsE2 = new TH2F("E1vsE2", "AnaBar Edep vs. Finger Edep", 100, 0.01, Finger_Edep_Max, 100, 0.01, AnaBar_Edep_Max);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0; // unused

  const int NMaxPMT=14;
  float edeptot[NMaxPMT];

  for (Int_t i = 0; i < nentries; i++) {

    bool anabar_hit_paddle[NMaxPMT];
    for (Int_t j=0; j<NMaxPMT; j++) {
	    edeptot[j] = 0.0;
    	    anabar_hit_paddle[j]=false; 
    }
    float edep0tot = 0.0; 
    tree1->GetEntry(i); 

    Float_t fMass;
    if (Prim_pdg == 11) {
	fMass = 0.511;
    } else {
	if (Prim_pdg == 13) {
	    fMass = 105.7;
	} else {
	    if (Prim_pdg == 2212) {
		fMass = 938.28;
	    } else {
		fMass = 939.5;
	    }
        }
    } 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph);
    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 

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

    for (Int_t j=0; j < Detector_Nhits ; j++) {

	if (trigger) {
		
		counter++; // unused
		if (Detector_id[j] == Detector_Offset ) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
			  edep0tot += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
				edep0tot += Detector_Ed[j];
			       }
			} 
		}
		if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
				edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
					edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
		     	       }
			}
		}

	}
    }

    if (trigger) {
    	hE1vsE2->Fill(edep0tot,edeptot[6]);
    }

  }

  TCanvas *c6 = new TCanvas("c6", "c6", 600,400,500,270);
  c6->Divide(1,1, 0.01, 0.01, 0);

  c6->cd(1);
  gPad->SetLogz();
  hE1vsE2->Draw("COLZ");

  return c6;

}


TCanvas *plotC7 (Float_t Theta_min_cut = 0.0, Int_t Analyse_Secondaries = 1){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------

  TH2F *hFinger_Edep_vs_Nphot = new TH2F("FingerEdepVsNphot", "Finger Edep vs. Number of Photons", Finger_NPhotons_Max, 0, Finger_NPhotons_Max, 100, 0.01, Finger_Edep_Max);
  TH2F *hNphot0_vs_Nphot1 = new TH2F("AnaBarVsFingerNphot", "AnaBar vs. Finger Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, Finger_NPhotons_Max, 0, Finger_NPhotons_Max);
  TH2F *hAnaBar_Edep_vs_Nphot = new TH2F("AnaBarEdepVsNphot", "AnaBar Edep vs. Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, 100, 0.01, AnaBar_Edep_Max);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries(); 

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 
  float edeptot[NMaxPMT];

  for (Int_t i = 0; i < nentries; i++) { 

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
	    edeptot[j] = 0.0; 
    	    anabar_hit_paddle[j]=false; 
    }
    float edep0tot = 0.0; 
    tree1->GetEntry(i); 

    Float_t fMass;
    if (Prim_pdg == 11) {
	fMass = 0.511;
    } else {
	if (Prim_pdg == 13) {
	    fMass = 105.7;
	} else {
	    if (Prim_pdg == 2212) {
		fMass = 938.28;
	    } else {
		fMass = 939.5;
	    }
        }
    } 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph);
    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 

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
    if (finger_hit && anabar_hit && fNewTheta > 2.524) trigger = true; 
    //if (finger_hit && anabar_hit) trigger = true; 

    if (trigger) {

	PMT_Nphotons_Total=0;
	for (Int_t icount = 0;icount<14;icount++){
		PMT_Nphotons_Total+=PMT_Nphotons[icount];
	}

    	for (Int_t j=0; j < Detector_Nhits ; j++) {
		
		counter++; // unused
		if (Detector_id[j] == Detector_Offset ) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
			  edep0tot += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
				edep0tot += Detector_Ed[j];
			       }
			} 
		}
		if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
				edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
					edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
		     	       }
			}
		}
    	}

    	hFinger_Edep_vs_Nphot->Fill(PMT_Nphotons[14],edep0tot); 
    	hAnaBar_Edep_vs_Nphot->Fill(PMT_Nphotons[6],edeptot[6]); 
    	hNphot0_vs_Nphot1->Fill(PMT_Nphotons_Total,PMT_Nphotons[14]); 
    }

  }

  TCanvas *c7 = new TCanvas("c7", "c7", 400,200,800,500);
  c7->Divide(2,2, 0.01, 0.01, 0);

  TCanvas *c7PE_MeV = new TCanvas("c7PE_MeV", "c7PE_MeV", 350,200,600,600);
  TCanvas *c7Profile = new TCanvas("c7Profile", "c7Profile", 500,200,600,600);

  c7->cd(1);
  hFinger_Edep_vs_Nphot->Draw("COLZ");
  c7->cd(2);
  hAnaBar_Edep_vs_Nphot->Draw("COLZ");
  c7->cd(3);
  hNphot0_vs_Nphot1->Draw("COLZ");
  //c7->cd(4);

  c7Profile->cd();
  TProfile *prof = hAnaBar_Edep_vs_Nphot->ProfileX();

  prof->Fit("pol1");

  c7PE_MeV->cd();
  hAnaBar_Edep_vs_Nphot->Draw("COLZ");

  return c7;

}


TCanvas *plotC8 (Float_t Theta_min_cut = 3.05, Int_t Analyse_Secondaries = 1){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------

  TH1F *hAnaBarPMTNphot[NUMPADDLE];
  TString name, title;
  for(Int_t i = 1; i <= NUMPADDLE; i++){
	name.Form("AnaBarPMTNphotA%d", i);
	title.Form("AnaBar PMT Number of Photons A%d", i);
	hAnaBarPMTNphot[i-1] = new TH1F(name, title, AnaBar_NPhotons_Max*0.9+20, -20, AnaBar_NPhotons_Max*0.9);
  }

  TH2F *hAnaBar_Edep_vs_Nphot = new TH2F("AnaBarEdepVsNphot", "AnaBar Edep vs. Number of Photons", AnaBar_NPhotons_Max, 0, AnaBar_NPhotons_Max, 100, 0.01, AnaBar_Edep_Max);

  TH1F *hAnaBarEd = new TH1F("AnaBarEd","AnaBar Energy Deposited", 100, 0.01, AnaBar_Edep_Max);
  TH1F *hAnaBarEdAll[NUMPADDLE]; // <-- I want to plot this for all 14 paddles to see if there is a trend on where the peak is 
  for (Int_t i = 1; i <= NUMPADDLE; i++){
	name.Form("AnaBarEd%d", i);
	title.Form("AnaBar Energy Deposited A%d", i);
	hAnaBarEdAll[i-1] = new TH1F(name, title, 100, 0.01, AnaBar_Edep_Max);
  }

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries(); 

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 
  float edeptot[NMaxPMT]; 

  for (Int_t i = 0; i < nentries; i++) { 

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
	    edeptot[j] = 0.0; 
    	    anabar_hit_paddle[j]=false; 
    }
    tree1->GetEntry(i); 

    Float_t fMass;
    if (Prim_pdg == 11) {
	fMass = 0.511;
    } else {
	if (Prim_pdg == 13) {
	    fMass = 105.7;
	} else {
	    if (Prim_pdg == 2212) {
		fMass = 938.28;
	    } else {
		fMass = 939.5;
	    }
        }
    } 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph);
    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 

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
    //if (finger_hit && anabar_hit) trigger = true;
    if (finger_hit && anabar_hit && fNewTheta > 2.524) trigger = true; 

    if (trigger) {

	for (Int_t icount = 0;icount < NUMPADDLE;icount++){
		PMT_Nphotons_Noise[icount]=PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
		hAnaBarPMTNphot[icount]->Fill(PMT_Nphotons_Noise[icount]);
	}

    	for (Int_t j=0; j < Detector_Nhits ; j++) {

		counter++; // unused
		if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
				edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
					edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
		     	       }
			}
		}
    	}

	for(Int_t i = 0; i < NUMPADDLE; i ++){
		hAnaBarEdAll[i]->Fill(edeptot[i]);
	}
    	hAnaBarEd->Fill(edeptot[6]);
    	hAnaBar_Edep_vs_Nphot->Fill(PMT_Nphotons[6],edeptot[6]); 
    }


  }

  TCanvas *c8 = new TCanvas("c8", "c8", 500,100,800,500);
  c8->Divide(2,2, 0.01, 0.01, 0);

  TCanvas *cEd = new TCanvas("cEd", "cEd", 600, 100, 800, 500);
  cEd->Divide(4,4);

  TCanvas *cEdOne = new TCanvas("cEdOne", "cEdOne", 450, 100, 800, 500);

  c8->cd(1);
  hAnaBarEd->Draw();
  c8->cd(2);
  hAnaBarPMTNphot[6]->Draw();
  c8->cd(3);
  hAnaBar_Edep_vs_Nphot->Draw("COLZ");

  TF1* function;
  Double_t means[NUMPADDLE], meanErr[NUMPADDLE];

  for(Int_t i = 0; i < NUMPADDLE; i++) {

	if (i+1 == 3)
	    cEdOne->cd();
	else
	    cEd->cd(i+1);
	
	hAnaBarEdAll[i]->Draw();

	//Double_t start = 6.4;
	//Double_t start = 5.8;
	Double_t start = 4.0;

	Double_t par3[3];
	TF1 *gf = new TF1("gf", "gaus", start, AnaBar_Edep_Max);
	hAnaBarEdAll[i]->Fit(gf, "R");
	function = hAnaBarEdAll[i]->GetFunction("gf");
	function->SetLineColor(1);

	means[i] = function->GetParameter(1);
	meanErr[i] = function->GetParError(1);
  }

  for(Int_t i = 0; i < NUMPADDLE; i++){
	cout << "Paddle " << i+1 << ": Mean peak Edep = " << means[i] << " MeV" << endl;
	cout << "    \t Mean Edep error = " << meanErr[i] << " MeV" << endl;
  }

  Double_t sumMeans = 0.0;
  Double_t meanMean;
  Double_t sumMeanErrSqrs = 0.0;
  Double_t sumErr;
  Double_t meanMeanErr;

  for(int i = 0; i < NUMPADDLE; i++){
	sumMeans += means[i];
	sumMeanErrSqrs += meanErr[i]*meanErr[i];
  }
  cout << "Sum of mean error squares = " << sumMeanErrSqrs << endl;

  meanMean = sumMeans/NUMPADDLE;

  sumErr = TMath::Sqrt(sumMeanErrSqrs);
  cout << "Error in sum of means = " << sumErr << endl;
  meanMeanErr = sumErr/NUMPADDLE;

  cout << "Mean peak Edep across all paddles: " << meanMean << " MeV" << endl;
  cout << "Mean peak Edep uncertainty: " << meanMeanErr << " Mev" << endl;


  return c8;

}


TCanvas *plotC9 (Float_t Theta_min_cut = 3.05, Float_t Edep_Threshold = 4.0, Int_t Analyse_Secondaries = 1){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------

  TH1F *hAnaBarPMTNphot[NUMPADDLE];
  TString name, title;
  for(Int_t i = 1; i <= NUMPADDLE; i++){
	name.Form("AnaBarPMTNphotA%d", i);
	title.Form("AnaBar PMT Number of Photons A%d", i);
	hAnaBarPMTNphot[i-1] = new TH1F(name, title, AnaBar_NPhotons_Max+20, -20, AnaBar_NPhotons_Max);
  }

  TH1F *hAnaBarPMTNoiseCutNphot[NUMPADDLE];
  for(Int_t i = 1; i <= NUMPADDLE; i++){
	name.Form("AnaBarPMTNoiseCutNphotA%d", i);
	title.Form("AnaBar PMT Number of Photons A%d", i);
	hAnaBarPMTNoiseCutNphot[i-1] = new TH1F(name, title, AnaBar_NPhotons_Max+20, -20, AnaBar_NPhotons_Max);
	hAnaBarPMTNoiseCutNphot[i-1]->SetLineColor(kRed);
  }

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 
  float edeptot[NMaxPMT]; 

  for (Int_t i = 0; i < nentries; i++) {

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
	    edeptot[j] = 0.0; 
    	    anabar_hit_paddle[j]=false; 
    }
    tree1->GetEntry(i); 

    Float_t fMass;
    if (Prim_pdg == 11) {
	fMass = 0.511;
    } else {
	if (Prim_pdg == 13) {
	    fMass = 105.7;
	} else {
	    if (Prim_pdg == 2212) {
		fMass = 938.28;
	    } else {
		fMass = 939.5;
	    }
        }
    } 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph);
    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 

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
    if (finger_hit && anabar_hit && fNewTheta > 3.05) trigger = true; 
    //if (finger_hit && anabar_hit) trigger = true; 

    if (trigger) {

	for (Int_t icount = 0;icount < NUMPADDLE;icount++){
		PMT_Nphotons_Noise[icount]=PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
		//PMT_Nphotons_Noise[icount]=PMT_Nphotons[icount];
		hAnaBarPMTNphot[icount]->Fill(PMT_Nphotons_Noise[icount]);
	}
    }

    for (Int_t j=0; j < Detector_Nhits ; j++) {

	if (trigger) {
		
		counter++; // unused
		if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
				edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == Prim_pdg && fNewTheta > Theta_min_cut) {
					edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
		     	       }
			}
		}
	}
    }

    if (trigger) {
	for (Int_t i = 0; i < NUMPADDLE; i++){
		if(anabar_hit_paddle[i]&&edeptot[i]>=Edep_Threshold)
			hAnaBarPMTNoiseCutNphot[i]->Fill(PMT_Nphotons_Noise[i]);
	}
    }


  }

  TCanvas *c9 = new TCanvas("c9", "c9", 350,250,800,500);
  c9->Divide(4,4, 0.01, 0.01, 0);

  TCanvas *c9Single = new TCanvas("c9Single", "c9Single", 400, 300, 800, 500);

  printf("Fittingi A1 ...\n"); 
  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  pllo[0]=0.05; pllo[1]=0.5; pllo[2]=1.0; pllo[3]=0.04;
  plhi[0]=10.0; plhi[1]=50.0; plhi[2]=10000.0; plhi[3]=5.0;
  sv[0]=1.8; sv[1]=5.0; sv[2]=1400.0; sv[3]=3.0;
  Double_t chisqr;
  Int_t    ndf;
  Double_t SNRPeak, SNRFWHM;

  for (Int_t i = 0; i < NUMPADDLE; i++){

	if(i == NUMPADDLE - 1)
	    c9Single->cd();
	else
	    c9->cd(i+1);

  	gPad->SetLogy();
  	hAnaBarPMTNphot[i]->Draw();
 	hAnaBarPMTNoiseCutNphot[i]->Draw("SAME");
  	fr[0]=0.7*hAnaBarPMTNphot[i]->GetMean();
  	fr[1]=25.0*hAnaBarPMTNphot[i]->GetMean();
  	TF1 *fitsnr = langaufit(hAnaBarPMTNphot[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
  	langaupro(fp,SNRPeak,SNRFWHM);
  	fitsnr->Draw("SAME");
  }

  return c9;

}


TCanvas *plotC10 (){

  // Histograms
  TH1F *hFingerPMTKE = new TH1F("FingerPMTKE","Photon Wavelength Production Spectrum", 400, 300.0, 700.0);
  TH1F *hAnaBarPMTKEA1 = new TH1F("AnaBarPMTKEA1","Photon Wavelength in WLS at PMT", 400, 300.0, 700.0);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries(); 

  Long64_t counter = 0; // unused

  const int NMaxPMT=14;

  for (Int_t i = 0; i < nentries; i++) { 

    bool anabar_hit_paddle[NMaxPMT];
    for (Int_t j=0; j<NMaxPMT; j++) {
    	    anabar_hit_paddle[j]=false;
    }
    tree1->GetEntry(i); 

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

    if (trigger) {

	for (Int_t iq=1; iq < NUMPADDLE; iq++){
		for (Int_t jq=0; jq<PMT_Nphotons[iq]; jq++){
			//std::cout << "Processing Anabar pmt = " << iq << " hit = " << jq << " Energy = " << PMT_KineticEnergy[iq][jq] << std::endl;
			hAnaBarPMTKEA1->Fill(1240.0/PMT_KineticEnergy[iq][jq]);
		}
	}

	for (Int_t jq=0; jq<PMT_Nphotons[14]; jq++){
		//std::cout << "Processing Finger hit = " << jq << " " << PMT_KineticEnergy[14][jq] << std::endl; 
		hFingerPMTKE->Fill(1240.0/PMT_KineticEnergy[14][jq]);
	}
    }


  }

  TCanvas *c10 = new TCanvas("c10", "c10", 1100,700,500,270);
  c10->Divide(1,2, 0.01, 0.01, 0);

  c10->cd(1);
  hFingerPMTKE->Draw();
  c10->cd(2);
  hAnaBarPMTKEA1->Draw();

  return c10;

}


TCanvas *plotC11 (Float_t Photon_Threshold = 8.0){

  // Histogram
  TH1F *hAnaBarMult = new TH1F("AnaBarMult","Anabar PMT Multiplicity",12,0,12);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 

  for (Int_t i = 0; i < nentries; i++) { 

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
    	    anabar_hit_paddle[j]=false;
    }
    tree1->GetEntry(i); 

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

    if (trigger) {

	for (Int_t icount = 0;icount < NUMPADDLE;icount++){
		PMT_Nphotons_Noise[icount]=PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
	}
	Int_t imult=0;
	for(Int_t icount=0;icount < NUMPADDLE;icount++){
		if(PMT_Nphotons_Noise[icount]>=Photon_Threshold) imult++;
	}
	hAnaBarMult->Fill(imult);
    }


  }

  TCanvas *c11 = new TCanvas("c11", "c11", 100,1000,500,270);
  c11->Divide(1,1, 0.01, 0.01, 0);

  c11->cd(1);
  hAnaBarMult->Draw();

  return c11;

}


TCanvas *plotC12 (){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------
  TH1F *hPrimPx = new TH1F("PrimPx","Primary Px", 100, -2000, 2000);
  TH1F *hPrimPz = new TH1F("PrimPz","Primary Pz", 100, -2000, 2000);
  TH1F *hPrimPy = new TH1F("PrimPy","Primary Py", 100, -25000, 0);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries(); 

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 

  for (Int_t i = 0; i < nentries; i++) { 

    tree1->GetEntry(i); 

    Float_t fMass;
    if (Prim_pdg == 11) {
	fMass = 0.511;
    } else {
	if (Prim_pdg == 13) {
	    fMass = 105.7;
	} else {
	    if (Prim_pdg == 2212) {
		fMass = 938.28;
	    } else {
		fMass = 939.5;
	    }
        }
    } 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPx        = fMomentum * TMath::Sin(Prim_Th) * TMath::Cos(Prim_Ph); // NOTE: Prim_Th and Prim_Ph are not the theta and phi of the CDet
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph); // coordinate system. They are the traditional spherical theta and 
    Float_t fPz        = fMomentum * TMath::Cos(Prim_Th);			// phi, with theta measured from the +z axis etc.

    hPrimPx->Fill(fPx);
    hPrimPy->Fill(fPy);
    hPrimPz->Fill(fPz);

  }

  TCanvas *c12 = new TCanvas("c12", "c12", 600,1000,500,270);
  c12->Divide(2,2, 0.01, 0.01, 0);
  
  c12->cd(1);
  hPrimPx->Draw();
  c12->cd(2);
  hPrimPy->Draw();
  c12->cd(3);
  hPrimPz->Draw();

  return c12;

}



Double_t langaufun(Double_t *x, Double_t *par) {

   //Fit parameters:
   //par[0]=Width (scale) parameter of Landau density
   //par[1]=Most Probable (MP, location) parameter of Landau density
   //par[2]=Total area (integral -inf to inf, normalization constant)
   //par[3]=Width (sigma) of convoluted Gaussian function
   //
   //In the Landau distribution (represented by the CERNLIB approximation), 
   //the maximum is located at x=-0.22278298 with the location parameter=0.
   //This shift is corrected within this function, so that the actual
   //maximum is identical to the MP parameter.

      // Numeric constants
      Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
      Double_t mpshift  = -0.22278298;       // Landau maximum location

      // Control constants
      Double_t np = 100.0;      // number of convolution steps
      Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

      // Variables
      Double_t xx;
      Double_t mpc;
      Double_t fland;
      Double_t sum = 0.0;
      Double_t xlow,xupp;
      Double_t step;
      Double_t i;


      // MP shift correction
      mpc = par[1] - mpshift * par[0]; 

      // Range of convolution integral
      xlow = x[0] - sc * par[3];
      xupp = x[0] + sc * par[3];

      step = (xupp-xlow) / np;

      // Convolution integral of Landau and Gaussian by sum
      for(i=1.0; i<=np/2; i++) {
         xx = xlow + (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);

         xx = xupp - (i-.5) * step;
         fland = TMath::Landau(xx,mpc,par[0]) / par[0];
         sum += fland * TMath::Gaus(x[0],xx,par[3]);
      }

      return (par[2] * step * sum * invsq2pi / par[3]);
}

TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF)
{
   // Once again, here are the Landau * Gaussian parameters:
   //   par[0]=Width (scale) parameter of Landau density
   //   par[1]=Most Probable (MP, location) parameter of Landau density
   //   par[2]=Total area (integral -inf to inf, normalization constant)
   //   par[3]=Width (sigma) of convoluted Gaussian function
   //
   // Variables for langaufit call:
   //   his             histogram to fit
   //   fitrange[2]     lo and hi boundaries of fit range
   //   startvalues[4]  reasonable start values for the fit
   //   parlimitslo[4]  lower parameter limits
   //   parlimitshi[4]  upper parameter limits
   //   fitparams[4]    returns the final fit parameters
   //   fiterrors[4]    returns the final fit errors
   //   ChiSqr          returns the chi square
   //   NDF             returns ndf

   Int_t i;
   Char_t FunName[100];

   sprintf(FunName,"Fitfcn_%s",his->GetName());

   TF1 *ffitold = (TF1*)gROOT->GetListOfFunctions()->FindObject(FunName);
   if (ffitold) delete ffitold;

   TF1 *ffit = new TF1(FunName,langaufun,fitrange[0],fitrange[1],4);
   ffit->SetParameters(startvalues);
   ffit->SetParNames("Width","MP","Area","GSigma");
   
   for (i=0; i<4; i++) {
      ffit->SetParLimits(i, parlimitslo[i], parlimitshi[i]);
   }

   his->Fit(FunName,"RB0");   // fit within specified range, use ParLimits, do not plot

   ffit->GetParameters(fitparams);    // obtain fit parameters
   for (i=0; i<4; i++) {
      fiterrors[i] = ffit->GetParError(i);     // obtain fit parameter errors
   }
   ChiSqr[0] = ffit->GetChisquare();  // obtain chi^2
   NDF[0] = ffit->GetNDF();           // obtain ndf

   return (ffit);              // return fit function

}


Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM) {

   // Seaches for the location (x value) at the maximum of the 
   // Landau-Gaussian convolute and its full width at half-maximum.
   //
   // The search is probably not very efficient, but it's a first try.

   Double_t p,x,fy,fxr,fxl;
   Double_t step;
   Double_t l,lold;
   Int_t i = 0;
   Int_t MAXCALLS = 10000;


   // Search for maximum

   p = params[1] - 0.1 * params[0];
   step = 0.05 * params[0];
   lold = -2.0;
   l    = -1.0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = langaufun(&x,params);
 
      if (l < lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-1);

   maxx = x;

   fy = l/2;


   // Search for right x location of fy

   p = maxx + params[0];
   step = params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;


   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-2);

   fxr = x;


   // Search for left x location of fy

   p = maxx - 0.5 * params[0];
   step = -params[0];
   lold = -2.0;
   l    = -1e300;
   i    = 0;

   while ( (l != lold) && (i < MAXCALLS) ) {
      i++;

      lold = l;
      x = p + step;
      l = TMath::Abs(langaufun(&x,params) - fy);
 
      if (l > lold)
         step = -step/10;
 
      p += step;
   }

   if (i == MAXCALLS)
      return (-3);


   fxl = x;

   FWHM = fxr - fxl;
   return (0);
}


void langaus() {
   // Fill Histogram
   Int_t data[100] = {0,0,0,0,0,0,2,6,11,18,18,55,90,141,255,323,454,563,681,
                    737,821,796,832,720,637,558,519,460,357,291,279,241,212,
                    153,164,139,106,95,91,76,80,80,59,58,51,30,49,23,35,28,23,
                    22,27,27,24,20,16,17,14,20,12,12,13,10,17,7,6,12,6,12,4,
                    9,9,10,3,4,5,2,4,1,5,5,1,7,1,6,3,3,3,4,5,4,4,2,2,7,2,4};
   TH1F *hSNR = new TH1F("snr","Signal-to-noise",400,0,400);

   for (Int_t i=0; i<100; i++) hSNR->Fill(i,data[i]);

}
