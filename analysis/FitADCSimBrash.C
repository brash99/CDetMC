// This is where we will fit the simulated ADC distribution to the real distribution.
// Actually we're going to convert the ADC signals from the real detector into photon numbers and fit the simulated number
// of photons to the data number of photons.

#include "TROOT.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TF1.h"
#include "TH1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLatex.h"
#include "TList.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TString.h"
#include "TStyle.h"
#include "TText.h"
#include "TTree.h"
#include "TGraph.h"
#include "TMath.h"
#include "TRandom3.h"

#include <cstring>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>

using namespace std;

static const int MaxHits = 50000;
static const int MaxPMTNo = 20;
static const int MaxPMTHits = 5000;
static const Float_t Finger_Edep_Max = 10.0;
static const Float_t AnaBar_Edep_Max = 10.0;
static const Float_t pedastel_sigma = 2.71;
static const Int_t Detector_Offset = 0;
static const Int_t Finger_NPhotons_Max = 150;
static const Int_t AnaBar_NPhotons_Max = 200;

static const Int_t NUMPADDLE = 14;
static const Int_t NUMPMT = 14;
static const Int_t NUMPIXEL = 16;
static const Int_t NUMPADDLES = NUMPMT*NUMPIXEL;

static const float adc_charge = 50*1e-15; // 50 fC, corresponding to an ADC chan
static const float e = 1.6e-19; // C, electron charge
static const Int_t xcanvas = 800; // width of canvases
static const Int_t ycanvas = 800; // height of canvases

//static const Int_t PEperMeV = 18; // photo-electrons per MeV
static const Int_t PEperMeV = 17;

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

TTree *T; // tree from root file, raw adc/tdc/hit data
// local arrays to hold data per event from tree
Double_t adc[NUMPADDLES];
Double_t adc_c[NUMPADDLES];
Double_t tdcl[NUMPADDLES];
Double_t tdcl_c[NUMPADDLES];
Double_t tdct[NUMPADDLES];
Double_t tdct_c[NUMPADDLES];
Double_t nahit;
Double_t nthit;
Double_t nhit;

//Latest map set for M1-L
//Int_t pixel1[NUMPMT]={1, 2, 13, 4, 6, 12, 4,1, 1,5,13,1, 1,8};
//Int_t pixel2[NUMPMT]={14, 16,16,6,13,13,16,3,6,9,15,4,16,12};
//Map set for M1-R
Int_t pixel1[NUMPMT]={4, 4, 2, 1, 4, 3, 4,13, 4,13,13,13, 5,13};
Int_t pixel2[NUMPMT]={5, 8,11,13,14,13,16,16,16,16,16,16,16,16};
//Original map of Ralph/Tommy
//Int_t pixel1[NUMPMT]={4, 4, 2, 1, 4, 3, 4,13, 4,13,13,13, 5,13};
//Int_t pixel2[NUMPMT]={5, 8,11,13,14,13,16,16,16,16,16,16,16,16};
//Map to test for missing pixels
//Int_t pixel1[NUMPMT]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};
//Int_t pixel2[NUMPMT]={0,0,0,0,0,0,0,0,0,0,0,0,0,0};

Int_t paddleindex[NUMPADDLES];

Int_t run; // run number used in titles of plots
Int_t n_events_to_analyze; // number of events to analyze ... -1 = all.

// function prototypes for simulation plot fitting
Double_t langaufun(Double_t *x, Double_t *par);
TF1 *langaufit(TH1F *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
Int_t langaupro(Double_t *params, Double_t &maxx, Double_t &FWHM);



void FitADCSimBrash(Int_t Analysis_Run_Number = 88811, Int_t runno = 1550, Int_t n_events=-1) {

  //-------------------------------------------------------------------
  //Set stuff up for reading simulated events
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

  //-------------------------------------------------------------------


  //*************** Set up for reading real data ***************//

  TString file2;
  file2.Form("/home/llorenti/analyzer/replay/rootfiles/scint_%d.root",runno);
  TFile *_file0 = TFile::Open(file2);

  run=runno;

  T = (TTree *)_file0->Get("T"); 
// Setting addresses so we can call them later
  T->SetBranchAddress("C.cdetm1r.adc",&adc);
  T->SetBranchAddress("C.cdetm1r.adc_c",&adc_c);
  T->SetBranchAddress("C.cdetm1r.tdcl",&tdcl);
  T->SetBranchAddress("C.cdetm1r.tdct",&tdct);
  T->SetBranchAddress("C.cdetm1r.tdcl_c",&tdcl_c);
  T->SetBranchAddress("C.cdetm1r.tdct_c",&tdct_c);
  T->SetBranchAddress("C.cdetm1r.nhit",&nhit);
  T->SetBranchAddress("C.cdetm1r.nahit",&nahit);
  T->SetBranchAddress("C.cdetm1r.nthit",&nthit);

  Int_t n_entries = T->GetEntries(); // Checks how many total entries in T
  cout << "Found " << n_entries << " events"<<endl;
  if (n_events==-1){
	n_events_to_analyze = n_entries;
	cout << "Analyzing all events." << endl;
  } else {
	if (n_events < n_entries) {
		n_events_to_analyze = n_events;
		cout << "Analyzing " << n_events << " events." << endl;
  	}else{
		n_events_to_analyze = n_entries;
		cout << "Analyzing " << n_entries << " events." << endl;
	}
  }

  //*************************************************************//

}


//*******************************************************************************//
//****************** Functions for simulation plots, for now ********************//
//*******************************************************************************//

TCanvas *plotC9 (/*Float_t Theta_min_cut = 3.017*/ Float_t Theta_min_cut = 0.0, Float_t Edep_Threshold = 0.0, Int_t Nphot_Neighbor_Cut = 8, Int_t Analyse_Secondaries = 1){

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------

  TH1F *hAnaBarPMTNphot[NUMPADDLE];
  TString name, title;
  for(Int_t i = 1; i <= NUMPADDLE; i++){
	name.Form("AnaBarPMTNphotA%d", i);
	title.Form("AnaBar PMT Number of Photons A%d", i);
	hAnaBarPMTNphot[i-1] = new TH1F(name, title, (AnaBar_NPhotons_Max+20)/4, -20, AnaBar_NPhotons_Max);
  }

  TH1F *hAnaBarPMTNoiseCutNphot[NUMPADDLE];
  for(Int_t i = 1; i <= NUMPADDLE; i++){
	name.Form("AnaBarPMTNoiseCutNphotA%d", i);
	title.Form("AnaBar PMT Number of Photons A%d", i);
	hAnaBarPMTNoiseCutNphot[i-1] = new TH1F(name, title, (AnaBar_NPhotons_Max+20)/4, -20, AnaBar_NPhotons_Max);
	hAnaBarPMTNoiseCutNphot[i-1]->SetLineColor(kRed);
  }

  TH1F *hNewTheta = new TH1F("theta","theta",200,3.14159265/2.0,3.14159265);
  TH2F *hNewThetaPhiCut = new TH2F("thetaphi","thetaphi",200,3.14159265/2.0,3.14159265,200,0,2.0*3.14159265);
  TH1F *hNewThetaCut = new TH1F("theta_cut","theta_cut",200,3.14159265/2.0,3.14159265);

  //-------------------------------------------------------------------
  //Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries = tree1->GetEntries();

  Long64_t counter = 0; // unused

  const int NMaxPMT=14; 
  float edeptot[NMaxPMT]; 

  cout << "Number of Entries = " << nentries << endl;
  int good_triggers = 0;

  for (Int_t i = 0; i < nentries; i++) {

    bool anabar_hit_paddle[NMaxPMT]; 
    for (Int_t j=0; j<NMaxPMT; j++) {
	    edeptot[j] = 0.0; 
    	    anabar_hit_paddle[j]=false; 
    }
    tree1->GetEntry(i); 

    Float_t fMass = 105.70; 
    Float_t fMomentum = sqrt(Prim_E*Prim_E - fMass*fMass); 
    Float_t fPy        = fMomentum * TMath::Sin(Prim_Th) * TMath::Sin(Prim_Ph);
    Float_t fPx        = fMomentum * TMath::Sin(Prim_Th) * TMath::Cos(Prim_Ph);
    Float_t fPz        = fMomentum * TMath::Cos(Prim_Th);
    Float_t fNewTheta = TMath::ACos(fPy/fMomentum); 
    Float_t fNewPhi = TMath::ATan(fPx/fPz);
  
    if (fPx < 0.0){
	fNewPhi = fNewPhi + 3.14159265;
    }
    if (fPx > 0.0 && fPz < 0.0){
	fNewPhi = fNewPhi + 2.0*3.14159265;
    }

  

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

    if (trigger) {
	good_triggers++;

	for (Int_t icount = 0;icount < NUMPADDLE;icount++){
		PMT_Nphotons_Noise[icount]=PMT_Nphotons[icount]+fRand->Gaus(0.0,pedastel_sigma);
		//PMT_Nphotons_Noise[icount]=PMT_Nphotons[icount];
		hAnaBarPMTNphot[icount]->Fill(PMT_Nphotons_Noise[icount]);
	}

    	for (Int_t j=0; j < Detector_Nhits ; j++) {
		
		counter++; // unused
		if (Detector_id[j] > Detector_Offset && Detector_id[j] <= NMaxPMT+Detector_Offset) {
			if (Analyse_Secondaries == 1 && fNewTheta > Theta_min_cut) {
				edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
			}else{ if (Detector_pdg[j] == 13 && fNewTheta > Theta_min_cut) {
					edeptot[Detector_id[j]-1-Detector_Offset] += Detector_Ed[j];
		     	       }
			}
		}
    	}

	hNewTheta->Fill(fNewTheta);
	for (Int_t i = 0; i < NUMPADDLE; i++){
		if(anabar_hit_paddle[i]&&edeptot[i]>=Edep_Threshold && fNewTheta > Theta_min_cut){
		    if(i == 0){
			if(PMT_Nphotons_Noise[i+1] < Nphot_Neighbor_Cut)
			    hAnaBarPMTNoiseCutNphot[i]->Fill(PMT_Nphotons_Noise[i]);
		    }
		    else if(i == NUMPADDLE - 1){
			if(PMT_Nphotons_Noise[i-1] < Nphot_Neighbor_Cut)
			    hAnaBarPMTNoiseCutNphot[i]->Fill(PMT_Nphotons_Noise[i]);
		    }
		    else {
			if(PMT_Nphotons_Noise[i-1] < Nphot_Neighbor_Cut && PMT_Nphotons_Noise[i+1] < Nphot_Neighbor_Cut){
			    hAnaBarPMTNoiseCutNphot[i]->Fill(PMT_Nphotons_Noise[i]);
			    if(PMT_Nphotons_Noise[i] > 5.0*Nphot_Neighbor_Cut) {
				hNewThetaCut->Fill(fNewTheta);
				hNewThetaPhiCut->Fill(fNewTheta,fNewPhi);
			    }
			}
		    }
		}
	}

    }


  }

  cout << "Number of good triggers = " << good_triggers << endl;

  TCanvas *c9 = new TCanvas("c9", "c9", 600,50,800,500);
  c9->Divide(4,4, 0.01, 0.01, 0);

  TCanvas *cSimFitRes = new TCanvas("cSimFitRes", "cSimFitRes", 600,100,800,500);
  TCanvas *cSimFitCounts = new TCanvas("cSimFitCounts", "cSimFitCounts", 600,100,800,500);
  TCanvas *cTheta = new TCanvas("cTheta", "cTheta", 600,100,800,500);
  TCanvas *cThetaPhi = new TCanvas("cThetaPhi", "cThetaPhi", 600,100,800,500);

  printf("Fitting A1 ...\n"); 
  // Setting fit range and start values
  Double_t fr[2];
  Double_t sv[4], pllo[4], plhi[4], fp[4], fpe[4];
  pllo[0]=0.05; pllo[1]=0.5; pllo[2]=1.0; pllo[3]=0.04;
  plhi[0]=10.0; plhi[1]=50.0; plhi[2]=10000.0; plhi[3]=5.0;
  sv[0]=1.8; sv[1]=5.0; sv[2]=1400.0; sv[3]=3.0;
  Double_t chisqr;
  Int_t    ndf;
  Double_t SNRPeak, SNRFWHM;

  TF1 *fcn;
  Double_t constants[NUMPADDLE] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t means[NUMPADDLE] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t sigmas[NUMPADDLE] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t meanErr[NUMPADDLE] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t sigErr[NUMPADDLE] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  Double_t res[NUMPADDLE], resErr[NUMPADDLE];
  Double_t counts[NUMPADDLE], countsErr[NUMPADDLE];
  Double_t upaddle[NUMPADDLE] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14};
  Double_t epaddle[NUMPADDLE] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  Double_t numEntries[NUMPADDLE];

  for (Int_t i = 1; i < NUMPADDLE-1; i++){
	c9->cd(i);
  	gPad->SetLogy();
  	//hAnaBarPMTNphot[i]->Draw();
 	hAnaBarPMTNoiseCutNphot[i]->Draw();

	numEntries[i] = hAnaBarPMTNoiseCutNphot[i]->GetEntries();

	Double_t bc =0.0;
	Double_t bn = 0.0;
	Int_t nbin = AnaBar_NPhotons_Max+20;
	Int_t min = -10;
	Int_t max = AnaBar_NPhotons_Max;

	//TF1 *gfit = new TF1("gfit", "gaus", 0.0, AnaBar_NPhotons_Max*0.8);
	//hAnaBarPMTNoiseCutNphot[i]->Fit(gfit, "R+");

    	for (Int_t j=0.2*nbin; j< nbin ;j++){ // find minimum between pedestal and peak to start the gaussian fit
	  bc = hAnaBarPMTNoiseCutNphot[i]->GetBinContent(j);
	  if( bc == 0.0 || (bc <hAnaBarPMTNoiseCutNphot[i]->GetBinContent(j+1) && bc < hAnaBarPMTNoiseCutNphot[i]->GetBinContent(j+2) && bc <hAnaBarPMTNoiseCutNphot[i]->GetBinContent(j-1) && bc <hAnaBarPMTNoiseCutNphot[i]->GetBinContent(j-2) && bc <hAnaBarPMTNoiseCutNphot[i]->GetBinContent(j-3) && bc <hAnaBarPMTNoiseCutNphot[i]->GetBinContent(j-4))  ){
	    bn = hAnaBarPMTNoiseCutNphot[i]->GetBinCenter(j);
	    break;
	  }
	}       

	Double_t par[3];
	TF1 *g1 = new TF1("g1", "gaus", min, bn);
	hAnaBarPMTNoiseCutNphot[i]->Fit(g1,"R");
	fcn = hAnaBarPMTNoiseCutNphot[i]->GetFunction("g1");
	fcn->SetLineColor(1);

	g1->GetParameters(&par[0]);
	Double_t blow = par[1]+20.0*par[2];
	TF1 *g2 = new TF1("g2","gaus",blow, max);

	hAnaBarPMTNoiseCutNphot[i]->Fit(g2, "R+");

	fcn = hAnaBarPMTNoiseCutNphot[i]->GetFunction("g2");
	fcn->SetLineColor(1);

	means[i] = fcn->GetParameter(1);
	sigmas[i] = fcn->GetParameter(2);
	constants[i] = fcn->GetParameter(0);
	meanErr[i] = fcn->GetParError(1);
	sigErr[i] = fcn->GetParError(2);


//  	fr[0]=0.7*hAnaBarPMTNphot[i]->GetMean();
//  	fr[1]=25.0*hAnaBarPMTNphot[i]->GetMean();
//  	TF1 *fitsnr = langaufit(hAnaBarPMTNphot[i],fr,sv,pllo,plhi,fp,fpe,&chisqr,&ndf);
//  	langaupro(fp,SNRPeak,SNRFWHM);
//  	fitsnr->Draw("SAME");
  }

  for(Int_t i = 1; i < NUMPADDLE-1; i++) {
	
	res[i] = sigmas[i]/means[i];
	resErr[i] = res[i] * TMath::Sqrt( (meanErr[i]*meanErr[i])/(means[i]*means[i]) + (sigErr[i]*sigErr[i])/(sigmas[i]*sigmas[i]) );
	counts[i] = constants[i]*sqrt(2.0*3.14159265)*sigmas[i];
	countsErr[i] = sqrt(counts[i]); 
  }

  res[0]=0.0;
  res[NUMPADDLE-1]=0.0;
  resErr[0]=0.0;
  resErr[NUMPADDLE-1]=0.0;
  counts[0]=0.0;
  counts[NUMPADDLE-1]=0.0;
  countsErr[0]=0.0;
  countsErr[NUMPADDLE-1]=0.0;

  cSimFitRes->cd();

  TGraphErrors *gs = new TGraphErrors(NUMPADDLE, upaddle, res, epaddle, resErr);
  gs->SetMarkerStyle(21);
  gs->SetMarkerColor(4);
  gs->GetXaxis()->SetTitle("Paddle Number");
  gs->GetYaxis()->SetTitle("Sim. Fit Resolution");
  gs->SetTitle("Resolution of Simulated Photoelectron Data");
  gs->Draw("AP");

  cSimFitRes->Update();
  
  cSimFitCounts->cd();

  TGraphErrors *gsc = new TGraphErrors(NUMPADDLE, upaddle, counts, epaddle, countsErr);
  gsc->SetMarkerStyle(21);
  gsc->SetMarkerColor(4);
  gsc->GetXaxis()->SetTitle("Paddle Number");
  gsc->GetYaxis()->SetTitle("Sim. Fit Counts");
  gsc->SetTitle("Counts in Simulated Photoelectron Data");
  gsc->Draw("AP");

  cSimFitCounts->Update();

  cTheta->cd();
  hNewTheta->Draw();
  hNewThetaCut->Draw("SAME");
  cTheta->Update();
  
  cThetaPhi->cd();
  hNewThetaPhiCut->Draw("COLZ");
  cThetaPhi->Update();

  Double_t meanSum = 0.0;
  Double_t avgMean;

  for(int i = 0; i < NUMPADDLE; i++){
	//cout << "Num entries in AnaBarPMTNoiseCutNphot histogram " << i + 1 << ": " << numEntries[i] << endl;
	cout << "Mean PE, paddle " << i + 1 << ": " << means[i] << endl;

	meanSum += means[i];
  }
  avgMean = meanSum/(NUMPADDLE-2);
  cout << "Avg. mean PE: " << avgMean << endl;


  return c9;

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

//*******************************************************************************//
//*******************************************************************************//



//*******************************************************************************//
//******************* Functions for real data plots, for now ********************//
//*******************************************************************************//

void setPaddleIndices(){
	// Setting the location of paddles, taking into account the missing pixels in each PMT
          for(Int_t pmt=0; pmt<NUMPMT;pmt++){
                Int_t ipaddle = (pmt+1)*NUMPADDLE+1;
                for(Int_t pixel=0; pixel<NUMPIXEL;pixel++){
                        Int_t index = pmt*NUMPIXEL+pixel;
                        if (pixel!=pixel1[pmt]-1&&pixel!=pixel2[pmt]-1){
                                ipaddle--;
				paddleindex[ipaddle]=index;
                        }
                }
           }
	   return;
}

Int_t getPaddleIndex(Int_t pmt, Int_t ipaddle){
  //pmt numbers are from 0 to 13 as are ipaddle numbers and the code accounts for removed pixels
	setPaddleIndices();

	Int_t localpaddle = (pmt+1)*NUMPADDLE-ipaddle;

	return paddleindex[localpaddle];
}

Int_t getPaddlePixel(Int_t pmt, Int_t ipaddle){
  //gives the pixel number for the given paddle number of the pmt, again paddle and pmt numbers are from 0 to 13
	setPaddleIndices();

	Int_t localpaddle = (pmt+1)*NUMPADDLE-ipaddle;

	return paddleindex[localpaddle]-(pmt)*NUMPIXEL;
}  

void print_event(Int_t adc_cut=50, Int_t start_event=1, Int_t num_events=1, Int_t tdc_min=750, Int_t tdc_width=300){

        for (Int_t id=start_event;id<start_event+num_events;id++){
          T->GetEntry(id);
	  cout << "Event " << id << endl;
	  for (Int_t pmt=1; pmt<=NUMPMT;pmt++){
	    for (Int_t index=1; index<=NUMPIXEL; index++){
	      Int_t ipixel = (pmt-1)*NUMPIXEL+index-1;
	      if(tdcl[ipixel]>tdc_min&&tdcl[ipixel]<tdc_min+tdc_width){
		cout << "Leading Edge TDC Hit on global pixel " << ipixel << " = " << tdcl[ipixel] << endl;
	      }
	      if(tdct[ipixel]>tdc_min&&tdct[ipixel]<tdc_min+tdc_width){
		cout << "Trailing Edge TDC Hit on global pixel " << ipixel << " = " << tdct[ipixel] << endl;
	      }
	      if(adc_c[ipixel]>adc_cut){
		cout << "ADC Hit on global pixel " << ipixel << " = " << adc_c[ipixel] << endl;
	      }	
	    }
	  }
	}
	
	return;
}

TCanvas *plot_adc_fit(Int_t pmt=7, Int_t tdc_min=850, Int_t tdc_width=100, Int_t adc_neighbor_cut=33, Int_t adc_cut=50, Int_t max=1000){

        TString cut, draw, draw1, title, grtitle;
        title.Form("run_%d_ADC_Fit",run);
        TCanvas *cADCFit = new TCanvas("cADCFit",title,xcanvas,ycanvas);
	title.Form("run_%d_ADC_Mean_Fit",run);
	TCanvas *cADCMeanFit= new TCanvas("cADCMeanFit",title,xcanvas,ycanvas);

	title.Form("run_%d_ADC_to_PE", run);
	TCanvas *cADCToPE = new TCanvas("cADCToPE", title, xcanvas, ycanvas);
	title.Form("run_%d_PE_Fit_Resolution", run);
	TCanvas *cPEFitRes =  new TCanvas("cPEFitRes", title, xcanvas, ycanvas);
	title.Form("run_%d_PE_Total_Counts", run);
	TCanvas *cPECounts =  new TCanvas("cPECounts", title, xcanvas, ycanvas);

	// Plot style stuff
	//TStyle *MyStyle = new TStyle("MyStyle","MyStyle");
	// Set Canvas preferences -> NOTE: this is screwing with the format of the simulation plot =/
	//MyStyle->SetTitleFontSize(0.08);
	//MyStyle->SetTitleX(0.15);
	//MyStyle->SetTitleY(0.99);
	//MyStyle->SetStatW(0.9);
	//MyStyle->SetMarkerStyle(6);
	//gStyle->SetCanvasDefH(xcanvas);
	//gStyle->SetCanvasDefW(ycanvas);
	//gStyle->SetPalette(1);
	//gROOT->SetStyle("MyStyle");


	// Setting up different arrays and objects used to create the two canvases.

	TH1D *htmp[NUMPIXEL];
	TH1I *hPE[NUMPIXEL];
	TF1 *function;
	TF1 *function1;
	Int_t histmax_cut[NUMPIXEL], gaus_cut_plus[NUMPIXEL];
	Float_t gaus_cut_minus[NUMPIXEL];

	Double_t upixel[NUMPIXEL]={1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
	Double_t epixel[NUMPIXEL]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t errors[NUMPIXEL]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	Double_t constants[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t means[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t sigmas[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	Double_t pconstants[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t pmeans[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t psigmas[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	Double_t meanPE[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t sigmaPE[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t constPE[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t meanError[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t sigError[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t constError[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	Double_t resolutions[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t resError[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t countIntegral[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	Double_t ecountIntegral[NUMPIXEL] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

	TString tmpentry;
	//MyStyle->SetStatX(0.9);
        //MyStyle->SetStatY(0.6);
        //MyStyle->SetStatY(0.4);

	setPaddleIndices(); //Setting the geometric paddle locations

	Int_t min=-100; // change max to 1500 for PMT 9
	Int_t nbin=(max-min)/10;

	for(Int_t i = 1; i <= NUMPIXEL; i++)
	  {
	    tmpentry.Form("htmp%d",i);
	    htmp[i-1] = new TH1D(tmpentry,tmpentry,nbin,min,max);
	    htmp[i-1]->SetLineColor(kRed);
	    title.Form("Run %d ADC Fit pmt %d, paddle %d: %d < tdc < %d",run,pmt,i,tdc_min,tdc_min+tdc_width);
	    htmp[i-1]->SetTitle(title);
	  }

	// Filling the histograms using only adc data with good tdc and good adc neighbor cut.

	Int_t nentries=n_events_to_analyze;

	for (Int_t id=1;id<=nentries;id++)
	  {
	    T->GetEntry(id);
	    Int_t ipaddle = (pmt)*NUMPADDLE+1;
	    for (Int_t pixel=0; pixel < NUMPIXEL; pixel++)
	    {
		Int_t index = (pmt-1)*NUMPIXEL+pixel;
		if (pixel!=pixel1[pmt-1]-1 && pixel!=pixel2[pmt-1]-1)
		{
		  ipaddle--;
		  if(tdcl[index] > tdc_min && tdcl[index] < tdc_min+tdc_width)
		  {
		    if (ipaddle < 2) 
		    {
			if (adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut)
			{
			  htmp[pixel]->Fill(adc_c[index]);
			}
		    }
		    else if (ipaddle == NUMPMT*NUMPADDLE) // should this be > or = ?
		    {
			if (adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut)
			{
			    htmp[pixel]->Fill(adc_c[index]);
			}
		    }
		    else
		    {
			if (adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut && adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut)
			{
			  htmp[pixel]->Fill(adc_c[index]);
			}
		    }
		  } // close tdc_width cut loop
		}
	    }
	  }


	// Creating the canvas of adc data with good tdc and fitting the adc with a gaussian or "landau" function.

	cADCFit->Clear();
	cADCMeanFit->Clear();
	cADCToPE->Clear();
	cPEFitRes->Clear();
	cPECounts->Clear();
	cADCFit->Divide(4,4);
	cADCToPE->Divide(4,4);

        Int_t count = 0;
        for (Int_t i=0; i<NUMPIXEL; i++)
	  {
	    if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) 
	      {
         //       means[i]=0.0;
	//	int check = 0;
		cADCFit->cd(count+1);
		gPad->SetLogy();
		cADCFit->Update();
                
		Int_t entries = htmp[i]->GetEntries();
		//float mean = htmp[i]->GetMean(1);
		//float RMS = htmp[i]->GetRMS(1);
                
		htmp[i]->SetStats(0);
		htmp[i]->Draw();

		
		Double_t bc =0.0;
		Double_t bn = 0.0;
		
	    	for (Int_t j=0.15*nbin; j< nbin ;j++){ // find minimum between pedestal and peak to start the gaussian fit
		  bc = htmp[i]->GetBinContent(j);
		  if( bc == 0.0 || (bc <htmp[i]->GetBinContent(j+1) && bc < htmp[i]->GetBinContent(j+2) && bc <htmp[i]->GetBinContent(j-1) && bc <htmp[i]->GetBinContent(j-2) && bc <htmp[i]->GetBinContent(j-3) && bc <htmp[i]->GetBinContent(j-4))  ){
		    bn = htmp[i]->GetBinCenter(j);
		    break;
		  }
		}       

		Double_t par[3];
		TF1 *g1 = new TF1("g1", "gaus", min, bn);
		htmp[i]->Fit(g1,"R");
		function = htmp[i]->GetFunction("g1");
		function->SetLineColor(1);

		pconstants[i] = function->GetParameter(0);
		pmeans[i] = function->GetParameter(1);
		psigmas[i] = function->GetParameter(2);
		
		g1->GetParameters(&par[0]);
		Double_t blow = par[1]+15.0*par[2];
		TF1 *g2 = new TF1("g2","gaus",blow, max);

		htmp[i]->Fit(g2, "R+");

		//htmp[i]->Fit("landau","","", 0, 250);
		function = htmp[i]->GetFunction("g2");
		function->SetLineColor(1);

		constants[i] = function->GetParameter(0);
		means[i] = function->GetParameter(1);
		sigmas[i] = function->GetParameter(2);
              
		count++;
	      
	      }
	  }

	Double_t PEperADC;
	Int_t minPE, maxPE, nbinPE;
	Double_t PEperADCs[NUMPIXEL];
	Int_t minPEs[NUMPIXEL], maxPEs[NUMPIXEL], nbinPEs[NUMPIXEL];

	for(Int_t i = 0; i < NUMPIXEL; i++) {

		PEperADCs[i] = PEperMeV*6.9/means[i];
		//minPEs[i] = 0;
		//maxPEs[i] = TMath::CeilNint(PEperADCs[i]*max);
		//nbinPEs[i] = maxPEs[i];

		minPE = -20;
		maxPE = 200;
		nbinPE = maxPE - minPE;

		tmpentry.Form("hPE%d", i+1);
		//hPE[i] = new TH1I(tmpentry, tmpentry, nbinPEs[i], minPEs[i], maxPEs[i]);
		hPE[i] = new TH1I(tmpentry, tmpentry, nbinPE, minPE, maxPE);
		hPE[i]->SetLineColor(kBlue);
		title.Form("Run %d ADC to PE, pmt %d, paddle/pixel %d", run, pmt, i+1);
		hPE[i]->SetTitle(title);		
	}

	Double_t adcPE;
	Int_t numPE;

	for (Int_t id=1;id<=nentries;id++)// re-loop events to get real histogram of photo-electron number
	  {
	    T->GetEntry(id);
	    Int_t ipaddle = (pmt)*NUMPADDLE+1;
	    for (Int_t pixel=0; pixel < NUMPIXEL; pixel++)
	    {
		Int_t index = (pmt-1)*NUMPIXEL+pixel;
		if (pixel!=pixel1[pmt-1]-1 && pixel!=pixel2[pmt-1]-1)
		{
		  ipaddle--;
		  if(tdcl[index] > tdc_min && tdcl[index] < tdc_min+tdc_width)
		  {
		    if (ipaddle < 2) 
		    {
			if (adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut)
			{
			    adcPE = PEperADCs[pixel]*adc_c[index];
			    if(adcPE - TMath::Floor(adcPE) < 0.5)
				numPE = TMath::FloorNint(adcPE);
			    else
				numPE = TMath::CeilNint(adcPE);
			    hPE[pixel]->Fill(numPE);
			}
		    }
		    else if (ipaddle == NUMPMT*NUMPADDLE) 
		    {
			if (adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut)
			{
			    adcPE = PEperADCs[pixel]*adc_c[index];
			    if(adcPE - TMath::Floor(adcPE) < 0.5)
				numPE = TMath::FloorNint(adcPE);
			    else
				numPE = TMath::CeilNint(adcPE);
			    hPE[pixel]->Fill(numPE);
			}
		    }
		    else
		    {
			if (adc_c[paddleindex[ipaddle-1]] < adc_neighbor_cut && adc_c[paddleindex[ipaddle+1]] < adc_neighbor_cut)
			{
			    adcPE = PEperADCs[pixel]*adc_c[index];
			    if(adcPE - TMath::Floor(adcPE) < 0.5)
				numPE = TMath::FloorNint(adcPE);
			    else
				numPE = TMath::CeilNint(adcPE);
			    hPE[pixel]->Fill(numPE);
			}
		    }
		  } // close tdc_width cut loop
		}
	    }
	  }

	Double_t numEntries[NUMPADDLE];

	count = 0;
	for(Int_t pixel = 0; pixel < NUMPIXEL; pixel++){

	    if (pixel!=pixel1[pmt-1]-1 && pixel!=pixel2[pmt-1]-1)
	    {
		cADCToPE->cd(count+1);
		gPad->SetLogy();
		cADCToPE->Update();

		hPE[pixel]->Draw();

		numEntries[count] = hPE[pixel]->GetEntries();

		Double_t bc =0.0;
		Double_t bn = 0.0;
		
		
	    	//for (Int_t j=0.15*nbinPEs[pixel]; j< nbinPEs[pixel] ;j++){ // find minimum between pedestal and peak to start the gaussian fit
		for (Int_t j=0.45*nbinPE; j< nbinPE ;j++){ // find minimum between pedestal and peak to start the gaussian fit
		  bc = hPE[pixel]->GetBinContent(j);
		  if( bc == 0.0 || (bc <hPE[pixel]->GetBinContent(j+1) && bc < hPE[pixel]->GetBinContent(j+2) && bc <hPE[pixel]->GetBinContent(j-1) && bc <hPE[pixel]->GetBinContent(j-2) && bc <hPE[pixel]->GetBinContent(j-3) && bc <hPE[pixel]->GetBinContent(j-4))  ){
		    bn = hPE[pixel]->GetBinCenter(j);
		    break;
		  }
		}      

		Double_t par[3];
		//TF1 *g1 = new TF1("g1", "gaus", minPEs[pixel], bn);
		TF1 *g1 = new TF1("g1", "gaus", minPE, bn);
		hPE[pixel]->Fit(g1,"R");
		function = hPE[pixel]->GetFunction("g1");
		function->SetLineColor(1);

		g1->GetParameters(&par[0]);
		Double_t blow = par[1]+21.0*par[2];
		//TF1 *g2 = new TF1("g2","gaus",blow, maxPEs[pixel]);
		TF1 *g2 = new TF1("g2","gaus",blow, maxPE);

		hPE[pixel]->Fit(g2, "R+");

		//htmp[i]->Fit("landau","","", 0, 250);
		function = hPE[pixel]->GetFunction("g2");
		function->SetLineColor(1);

		meanPE[pixel] = function->GetParameter(1);
		sigmaPE[pixel] = function->GetParameter(2);
		constPE[pixel] = function->GetParameter(0);
		meanError[pixel] = function->GetParError(1);
		sigError[pixel] = function->GetParError(2);
		constError[pixel] = function->GetParError(0);

		count++;
	    }
	}

	for(Int_t i = 0; i < NUMPIXEL; i++) {

	    if (i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1){
		
		resolutions[i] = sigmaPE[i]/meanPE[i];
		resError[i] = resolutions[i]*TMath::Sqrt( (meanError[i]*meanError[i])/(meanPE[i]*meanPE[i]) + (sigError[i]*sigError[i])/(sigmaPE[i]*sigmaPE[i]) );
		countIntegral[i] = constPE[i]*sqrt(2.0*3.14159265)*sigmaPE[i];
		ecountIntegral[i] = sqrt(countIntegral[i]); 
	    }
	}

	cPEFitRes->cd();

	TGraphErrors *gPE = new TGraphErrors(NUMPIXEL, upixel, resolutions, epixel, resError);
	gPE->SetMarkerStyle(21);
	gPE->SetMarkerColor(4);
	gPE->GetXaxis()->SetTitle("Pixel Number");
	gPE->GetYaxis()->SetTitle("PE Fit Resolution");
	grtitle.Form("run_%d_pmt_%d_PE_fit_res",run,pmt);
	gPE->SetTitle(grtitle);

	gPE->Draw("AP");
	
	cPEFitRes->Update();

	cPECounts->cd();

	TGraphErrors *gPEcounts = new TGraphErrors(NUMPIXEL, upixel, countIntegral, epixel, ecountIntegral);
	gPEcounts->SetMarkerStyle(21);
	gPEcounts->SetMarkerColor(4);
	gPEcounts->GetXaxis()->SetTitle("Pixel Number");
	gPEcounts->GetYaxis()->SetTitle("Total Counts in Peak");
	grtitle.Form("run_%d_pmt_%d_PE_total_counts",run,pmt);
	gPEcounts->SetTitle(grtitle);

	gPEcounts->Draw("AP");

	cPECounts->Update();

	//gPE->GetYaxis()->SetTitleOffset(1.4);

	/*count = 0;
	for(Int_t pixel = 0; pixel < NUMPIXEL; pixel++){

	    if (pixel!=pixel1[pmt-1]-1 && pixel!=pixel2[pmt-1]-1)
	    {
		PEperADC = PEperMeV*8/means[pixel];

		minPE = 0;
		maxPE = TMath::CeilNint(PEperADC*max); 
		nbinPE = maxPE - minPE;

		tmpentry.Form("hPE%d", pixel+1);
		hPE[pixel] = new TH1I(tmpentry, tmpentry, nbinPE, minPE, maxPE);
		hPE[pixel]->SetLineColor(kBlue);
		title.Form("Run %d ADC to PE, pmt %d, paddle/pixel %d", run, pmt, pixel);
		hPE[pixel]->SetTitle(title);

		for(Int_t i = 0; i < nbin; i++){

		    Int_t inc = htmp[pixel]->GetBinContent(i); // cast issue????
		    Double_t adcbin = htmp[pixel]->GetBinCenter(i);

		    if (adcbin >= 0.0){

			Double_t newbin = PEperADC*adcbin;
			Int_t PEbin;

			if(newbin-TMath::Floor(newbin) < 0.5)
			    PEbin = TMath::FloorNint(newbin);
			else
			    PEbin = TMath::CeilNint(newbin);

			hPE[pixel]->Fill(PEbin,inc*1.0);
		    }
		}

		cADCToPE->cd(count+1);
		gPad->SetLogy();
		cADCToPE->Update();

		hPE[pixel]->Draw();

		count++;
	    }

	}*/

	  

	// Printing out means and sigmas to screen with corresponding pixel number.

	Double_t thresh_sum = 0.0;
	Double_t sigma_sum = 0.0;
	Double_t mean_sum = 0.0;
	for (Int_t i=0; i<NUMPIXEL; i++)
	  {
	   if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1)
	    {
	    cout << "Pixel number: " << i+1 << "   \t mean = " << means[i] << "    \t sigma = " << sigmas[i] << endl;
	    cout << "Pixel number: " << i+1 << "   \t ped mean = " << pmeans[i] << "    \t ped sigma = " << psigmas[i] << endl;
	    cout << "Suggested threshold = " << pmeans[i]+3.0*psigmas[i] << endl;
	    thresh_sum += pmeans[i]+3.0*psigmas[i];
	    sigma_sum += psigmas[i];
	    mean_sum += means[i];
	    }
	  }
	Double_t average_thresh = thresh_sum/(NUMPIXEL-2);
	Double_t average_sigma = sigma_sum/(NUMPIXEL-2);
	Double_t average_mean = mean_sum/(NUMPIXEL-2);
	cout << "Suggested Threshold overall = " << average_thresh << endl;
	cout << "Average Sigma = " << average_sigma << endl;
	cout << "Average Mean = " << average_mean << endl;
	cout << "Suggested Calibration = " << PEperMeV*6.9/average_mean << endl;
	cout << "Suggested Smearing = " << average_sigma*PEperMeV*6.9/average_mean << endl;

        title.Form("run_%d_ADC_pmt_%d_tdc_min_%d_max_%d.png",run,pmt,tdc_min,tdc_min+tdc_width);
        cADCFit->Print(title);
        cADCFit->cd(0);

	// End of work on canvas with adc data fitted with a function

	// Creating the canvas of mean adc data with good tdc data.

	cADCMeanFit->cd();

	Double_t mmean[2]={0,0};
	Double_t merror[2]={1,1};
	Double_t mpixel[2]={0,0};
	Double_t mepixel[2]={0,0};
	mpixel[0]=pixel1[pmt-1];
	mpixel[1]=pixel2[pmt-1];

	Double_t yline[2]={100,100};
	Double_t xline[2]={0,17};

	// Calculating errors on the mean adc data and printing out number of entries, means, and errors
	// to screen with corresponding pixel number.

	count = 0;
	for (Int_t i = 0; i < NUMPIXEL; i++)
	{
	    if(i != pixel1[pmt-1]-1 && i != pixel2[pmt-1]-1) 
	    {
		Int_t entries = htmp[i]->GetEntries();
		errors[i] = sigmas[i]/sqrt(entries);
		if(means[i] > 1)
		{
		    cout << "Pixel number: " <<  i+1 << "  \t entries: " << entries << "\t mean: " 
			 << means[i] << "\t error: " << errors[i] << endl;
		}
               	count++;
	    }
	}

	// Creating three different graphs for used pixels, missing pixels, and general formatting.

	TGraphErrors *gr = new TGraphErrors(NUMPIXEL,upixel,means,epixel,errors);
	gr->SetMarkerStyle(21);
	gr->GetXaxis()->SetTitle("Pixel Number");
	gr->GetYaxis()->SetTitle("Mean ADC Fit (Good TDC)");
	gr->GetYaxis()->SetTitleOffset(1.4);
	grtitle.Form("run_%d_pmt_%d_adc_mean",run,pmt);
	gr->SetTitle(grtitle);
	TGraphErrors *gr2 = new TGraphErrors(2,mpixel,mmean,mepixel,merror);
	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(2);
	gr2->SetTitle("");
	TGraph *gr3 = new TGraph(2,xline,yline);
	gr3->SetLineColor(2);
	gr3->SetLineWidth(2);
	gr3->SetLineStyle(2);
	gr3->SetTitle("");

	// Drawing the three different graphs to the adc mean values canvas.

	gr->Draw("AP");
	gr2->Draw("P");
	gr3->Draw("L");
            
        cADCMeanFit->Update();
	cADCMeanFit->cd(0);

	// End of work on canvas with adc mean data.

	// Creating a file to output parameters of function fitted to the ADC data with good TDC.
        TString filename;
        filename.Form("run_%d_pmt_%d_adc_fitresults.txt", run, pmt);
	ofstream outfile(filename);

	outfile << "Run #: " << run << "\t PMT #: " << pmt << " \t  Paramters of fitted ADC data with good TDC. \n";
	outfile << "Pixel # \t Mean Value \t Error \t Sigma \n" << endl;
	for(Int_t i = 0; i < NUMPIXEL; i++)
	  {
	    outfile << i+1 << "\t" << means[i] << "\t" << errors[i] << "\t" << sigmas[i] << endl;
	  }
	cout << "Created output file for parameters." << endl;
	outfile.close();

	for(int i = 0; i < NUMPADDLE; i++){
		cout << "Num entries in hPE histogram " << i+1 << ": " << numEntries[i] << endl;
	}
	
	return cADCMeanFit;
	return cADCFit;
}

//*******************************************************************************//
//*******************************************************************************//
