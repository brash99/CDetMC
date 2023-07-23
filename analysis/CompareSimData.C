#include "TMath.h"

Double_t fitfunction(Double_t *x, Double_t *par) {
	Float_t xx = x[0];
	Float_t binsize = 1000.0/100.0;
	//Double_t f = binsize*par[0]/sqrt(2.0*3.14159265*pow(par[2],2))*exp(-0.5*(xx-par[1])*(xx-par[1])/(par[2]*par[2]));
	Double_t f = par[3]*xx*xx+par[4]*xx+par[5]+binsize*par[0]/sqrt(2.0*3.14159265*pow(par[2],2))*exp(-0.5*(xx-par[1])*(xx-par[1])/(par[2]*par[2]));
	return f;
}


void CompareSimData(Int_t Data_Run_Number = 157, Int_t Simulation_Run_Number = 10) {

  //-------------------------------------------------------------------
  //Set stuff up for reading
  //-------------------------------------------------------------------
  TString simfilename;
  simfilename.Form("data/G4_test%d.root",Simulation_Run_Number);
  TFile *f1 = new TFile(simfilename,"READ");
  TTree *tree1 = (TTree*)f1->Get("T");
  TString datafilename;
  datafilename.Form("data/AnaBar_Processed_%d.root",Data_Run_Number);
  TFile *g1 = new TFile(datafilename,"READ");
  TTree *tree2 = (TTree*)g1->Get("tree1");

  UShort_t Element0, Element1;
  Short_t fRawTDC312_0, fRawTDC314_0;
  Short_t fRawTDC312_1, fRawTDC314_1;
  Short_t fRawTDC312_2, fRawTDC314_2;
  Short_t fRawTDC312_3, fRawTDC314_3;
   
  tree2->SetBranchAddress("Element0", &Element0);
  tree2->SetBranchAddress("Element1", &Element1);
  //tree2->SetBranchAddress("Element2", &Element2);
  //tree2->SetBranchAddress("TDCRef264", &fTDCRef264);
  tree2->SetBranchAddress("RawTDC312_0", &fRawTDC312_0 );
  tree2->SetBranchAddress("RawTDC314_0", &fRawTDC314_0 );
  tree2->SetBranchAddress("RawTDC312_1", &fRawTDC312_1 );
  tree2->SetBranchAddress("RawTDC314_1", &fRawTDC314_1 );
  tree2->SetBranchAddress("RawTDC312_2", &fRawTDC312_2 );
  tree2->SetBranchAddress("RawTDC314_2", &fRawTDC314_2 );
  tree2->SetBranchAddress("RawTDC312_3", &fRawTDC312_3 );
  tree2->SetBranchAddress("RawTDC314_3", &fRawTDC314_3 );
  
  const int MaxHits = 10000;
  Float_t Prim_E, Prim_Th, Prim_Ph;
  Int_t Prim_pdg;
  Int_t Detector_Nhits;
  Int_t Detector_pdg[MaxHits];
  Int_t Detector_id[MaxHits];
  Float_t Detector_x[MaxHits], Detector_y[MaxHits], Detector_z[MaxHits], Detector_t[MaxHits];
  Float_t Detector_Ed[MaxHits];  
  Int_t PMT_id;
  Int_t PMT_Nphotons_Zero;
  Int_t PMT_Nphotons_One;

  tree1->SetBranchAddress("Prim_E", &Prim_E);
  tree1->SetBranchAddress("Prim_Th", &Prim_Th);
  tree1->SetBranchAddress("Prim_Ph", &Prim_Ph);
  tree1->SetBranchAddress("Prim_pdg", &Prim_pdg);
  tree1->SetBranchAddress("PMT_Nphotons_Zero", &PMT_Nphotons_Zero);
  tree1->SetBranchAddress("PMT_Nphotons_One", &PMT_Nphotons_One);
  tree1->SetBranchAddress("PMT_id", &PMT_id);
  tree1->SetBranchAddress("Detector_Nhits", &Detector_Nhits);
  tree1->SetBranchAddress("Detector_pdg", &Detector_pdg);
  tree1->SetBranchAddress("Detector_id", &Detector_id);
  tree1->SetBranchAddress("Detector_x", &Detector_x);
  tree1->SetBranchAddress("Detector_y", &Detector_y);
  tree1->SetBranchAddress("Detector_z", &Detector_z);
  tree1->SetBranchAddress("Detector_t", &Detector_t);
  tree1->SetBranchAddress("Detector_Ed", &Detector_Ed);

  //-------------------------------------------------------------------
  //Create histograms
  //-------------------------------------------------------------------

  TH1F *hTDC0 = new TH1F("TDC0", "Finger TDC 1st Hit", 499, 1, 500);
  TH1F *hTDC2 = new TH1F("TDC2", "AnaBar TDC 1st Hit", 499, 1, 500);
  
  TH1F *hFingerPMTNphot = new TH1F("FingerPMTNphot","Finger PMT Number of Photoelectrons", 100, 0, 200);
  TH1F *hAnaBarPMTNphot = new TH1F("AnaBarPMTNphot","AnaBar PMT Number of Photoelectrons", 100, 0, 200);
  TH1F *hFingerPMTNphotScaled = new TH1F("FingerPMTNphotScaled","Finger PMT Number of Photoelectrons Scaled", 100, 0, 1000);
  TH1F *hAnaBarPMTNphotScaled = new TH1F("AnaBarPMTNphotScaled","AnaBar PMT Number of Photoelectrons Scaled", 100, 0, 1000);
  TH1F *hElement00 = new TH1F("Element00","Finger ADC (Good TDCs)", 100, 0, 1000);
  TH1F *hElement10 = new TH1F("Element10","AnaBar ADC (Good TDCs)", 100, 0, 1000); 

  //Limits and cuts
  //-------------------------------------------------------------------

   Float_t Peds[2] = {240, 195}; // mean + 1sigma of the pedestal gaussian
  
  //-------------------------------------------------------------------
  //Data Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries_data = tree2->GetEntries();

  Double_t counter00_data = 0;
  Double_t counter10_data = 0;
  
  for (Int_t i = 0; i < nentries_data; i++) {
    tree2->GetEntry(i);

    if ((Element0-Peds[0])>0.0 && (Element0-Peds[0]) < 1000.0 &&
	(Element1-Peds[1])>0.0 && (Element1-Peds[1]) < 1000.00 ) {

	hTDC0->Fill( fRawTDC312_0 );
    	hTDC2->Fill( fRawTDC314_0 );
    }
  }
  
  hTDC0->Fit("gaus");
  TF1 *tdc0fit = hTDC0->GetFunction("gaus");
  tdc0_mean = tdc0fit->GetParameter(1);
  tdc0_sigma = tdc0fit->GetParameter(2);
  hTDC2->Fit("gaus");
  TF1 *tdc2fit = hTDC2->GetFunction("gaus");
  tdc2_mean = tdc2fit->GetParameter(1);
  tdc2_sigma = tdc2fit->GetParameter(2);
 
  for (Int_t i = 0; i < nentries_data; i++) {
    tree2->GetEntry(i);

    //Basically same TDC peak for PMT1 and 2.
    if (fRawTDC312_0 > (tdc0_mean-6.0*tdc0_sigma) && fRawTDC312_0 < (tdc0_mean+6.0*tdc0_sigma) ) {
	if (fRawTDC314_0 > tdc2_mean-6.0*tdc2_sigma && fRawTDC314_0 < tdc2_mean+6.0*tdc2_sigma) {
                if ((Element0-Peds[0]) > 0.0 && (Element0-Peds[0]) < 1000.0) {
	        counter00_data++;
		hElement00->Fill(Element0 - Peds[0]);
                }
                if ((Element1-Peds[1]) > 0.0 && (Element1-Peds[1]) < 1000.0) {
	        counter10_data++;
		hElement10->Fill(Element1 - Peds[1]);
                }
	}

    }
  }

  Double_t npoints00 = hElement00->GetEntries();
  Double_t mean00 = hElement00->GetMean();
  Double_t rms00 = hElement00->GetRMS();
  Double_t npoints10 = hElement10->GetEntries();
  Double_t mean10 = hElement10->GetMean();
  Double_t rms10 = hElement10->GetRMS();
  

  //-------------------------------------------------------------------
  // Simulation Event loop
  //-------------------------------------------------------------------
  
  Long64_t nentries_sim = tree1->GetEntries();

  Double_t counter_sim = 0;
  
  for (Int_t i = 0; i < nentries_sim; i++) {
    
    tree1->GetEntry(i);

    bool trigger = false;
    bool finger_hit = false;
    bool anabar_hit = false;
    int j_finger = 0;
    int j_anabar = 0;
    for (Int_t j=0; j < Detector_Nhits ; j++) {
	if (Detector_id[j] == 0 && !finger_hit) {
		finger_hit = true;
		j_finger = j;
	}
	if (Detector_id[j] == 1 && !anabar_hit) {
		anabar_hit = true;
		j_anabar = j;
	}
    }

    if (finger_hit && anabar_hit) trigger = true; 
    if (trigger) {
        counter_sim++;
  	hAnaBarPMTNphot->Fill(PMT_Nphotons_Zero);
        hFingerPMTNphot->Fill(PMT_Nphotons_One);
    }

  }
  
  Double_t simnpoints00 = hFingerPMTNphot->GetEntries();
  Double_t simmean00 = hFingerPMTNphot->GetMean();
  Double_t simrms00 = hFingerPMTNphot->GetRMS();
  Double_t simnpoints10 = hAnaBarPMTNphot->GetEntries();
  Double_t simmean10 = hAnaBarPMTNphot->GetMean();
  Double_t simrms10 = hAnaBarPMTNphot->GetRMS();

  //Double_t finger_nscaling = npoints00/simnpoints00;
  //Double_t anabar_nscaling = npoints10/simnpoints10;
  Double_t finger_nscaling = counter00_data/counter_sim;
  Double_t anabar_nscaling = counter10_data/counter_sim;
  //Double_t finger_xscaling = mean00/simmean00*1.00;
  //Double_t anabar_xscaling = mean10/simmean10*1.00;
  Double_t finger_xscaling = 2.72;
  Double_t anabar_xscaling = 3.09;
  std::cout << "Finger Counter ... nscaling = " << finger_nscaling << "  xscaling = " << finger_xscaling << std::endl;
  std::cout << "AnaBar Counter ... nscaling = " << anabar_nscaling << "  xscaling = " << anabar_xscaling << std::endl;
  
  for (Int_t i = 0; i < nentries_sim; i++) {
    
    tree1->GetEntry(i);

    bool trigger = false;
    bool finger_hit = false;
    bool anabar_hit = false;
    int j_finger = 0;
    int j_anabar = 0;
    for (Int_t j=0; j < Detector_Nhits ; j++) {
	if (Detector_id[j] == 0 && !finger_hit) {
		finger_hit = true;
		j_finger = j;
	}
	if (Detector_id[j] == 1 && !anabar_hit) {
		anabar_hit = true;
		j_anabar = j;
	}
    }

    if (finger_hit && anabar_hit) trigger = true; 
    if (trigger) {
  	hAnaBarPMTNphotScaled->Fill(PMT_Nphotons_Zero*anabar_xscaling,anabar_nscaling);
        hFingerPMTNphotScaled->Fill(PMT_Nphotons_One*finger_xscaling,finger_nscaling);
    }

  }

  
  //-------------------------------------------------------------------
  //Plotting and writing out
  //-------------------------------------------------------------------
  
  
  TCanvas *c4 = new TCanvas("c4", "c4", 100,500,600,400);
  c4->Divide(2,1, 0.01, 0.01, 0);
  
  //c4->cd(1);
  //hFingerPMTNphot->Draw();
  //c4->cd(2);
  //hAnaBarPMTNphot->Draw();
  //c4->cd(3);
  //hElement00->Draw();
  //c4->cd(4);
  //hElement10->Draw();
  c4->cd(1);
  hFingerPMTNphotScaled->SetLineColor(kRed);
  hFingerPMTNphotScaled->Draw();
  hElement00->Draw("SAME");
  c4->cd(2);
  hAnaBarPMTNphotScaled->SetLineColor(kRed);
  hAnaBarPMTNphotScaled->Draw();
  hElement10->Draw("SAME");

}

