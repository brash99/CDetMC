#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"
#include "AnalysisMessenger.hh"

#include "G4UnitsTable.hh"
#include "G4RunManager.hh"
#include "G4Point3D.hh"
#include "G4Transform3D.hh"

#include "TROOT.h"
#include "TApplication.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

//---------------------------------------------------------------------------

AnalysisManager::AnalysisManager()
{
  ZeroArray();

  fOutFileName = TString("output/out_default.root");

  fAnaMessenger = new AnalysisMessenger(this);
}

//---------------------------------------------------------------------------

AnalysisManager::~AnalysisManager()
{
   fROOTtree->Write();
   fROOTfile->Close();
//std::cout<<"Tree has been writen and closed bc analysis manager has been deconstructed"<<std::endl;
}

//---------------------------------------------------------------------------

void AnalysisManager::InitOutput()
{ 

  fROOTfile = new TFile(fOutFileName,"RECREATE","fROOTfile",1);
  fROOTtree = new TTree("T","Output Tree");
  fROOTtree->SetAutoSave();
//std::cout<<"Tree has been created"<<std::endl;
  // Set Primary Branches
  fROOTtree->Branch("Prim_E",      &fPEne,   "Prim_E/F"      );
  fROOTtree->Branch("Prim_X",      &fXvtx,   "Prim_X/F"      );
  fROOTtree->Branch("Prim_Y",      &fYvtx,   "Prim_Y/F"      );
  fROOTtree->Branch("Prim_Z",      &fZvtx,   "Prim_Z/F"      );
  fROOTtree->Branch("Prim_Th",     &fPth,    "Prim_Th/F"     );
  fROOTtree->Branch("Prim_Ph",     &fPph,    "Prim_Ph/F"     );
  fROOTtree->Branch("Prim_pdg",    &fPpdg,   "Prim_pdg/I"    );
  
  // Set PMT Hit Branches
  fROOTtree->Branch("PMT_id",     &fPMTNo,     "PMT_id/I" );  
  fROOTtree->Branch("PMT_Nphotons",  fNphotons,  "PMT_Nphotons[50000]/I" );  
  fROOTtree->Branch("PMT_Time",  fPMTTime,  "PMT_Time[50000]/F" );  
  //fROOTtree->Branch("PMT_KineticEnergy",  fPMTKineticEnergy,  "PMT_KineticEnergy[50000][1000]/F" );  

  // Set Raw Detector Step Hit Branches
  fROOTtree->Branch("Detector_Nhits", &fRAW_Nhits, "Detector_Nhits/I");  
  fROOTtree->Branch("Detector_pdg",   fRAW_pdg,    "Detector_pdg[Detector_Nhits]/I");
  fROOTtree->Branch("Detector_id",    fRAW_id,     "Detector_id[Detector_Nhits]/I");
  fROOTtree->Branch("Detector_x",     fRAW_xpre,   "Detector_x[Detector_Nhits]/F"  );
  fROOTtree->Branch("Detector_y",     fRAW_ypre,   "Detector_y[Detector_Nhits]/F"  );
  fROOTtree->Branch("Detector_z",     fRAW_zpre,   "Detector_z[Detector_Nhits]/F"  );
  fROOTtree->Branch("Detector_t",     fRAW_time,   "Detector_t[Detector_Nhits]/F"  );
  fROOTtree->Branch("Detector_Ed",    fRAW_Edep,   "Detector_Ed[Detector_Nhits]/F" );
//std::cout<<"branches have been created"<<std::endl;
}

//---------------------------------------------------------------------------

void AnalysisManager::ZeroArray()
{
//std::cout<<"Start of ZeroArray"<<std::endl;
  // Primary
  G4ThreeVector zero(0.,0.,0.);
  fPEne   = 9999;
  fXvtx   = 9999;
  fYvtx   = 9999;
  fZvtx   = 9999;
  fPdir   = (zero);
  fPth    = 9999;
  fPph    = 9999;
  fPTime  = 9999;
  fPPDef  = NULL;
  fPpdg   = 9999;

  // PMT
  fPMTNo    = -1;
  for ( Int_t i = 0; i < fMaxPMTNo; i++) {
	fNphotons[i]=0;
	fPMTTime[i]=0;
 	//for (Int_t j = 0; j < fMaxPMTHits; j++) {
	//	fPMTKineticEnergy[i][j] = 0;
	//}	
  }

  // Raw Hits
  fRAW_Nhits  = 0;
  fSteppdef    = NULL;
  fStepp3      = (zero);
  fSteppospre  = (zero);
  fSteppospost = (zero);
  fStepid      = 0;
  fSteptime    = 0;
  fStepedep    = 0;

  for ( Int_t i = 0; i < fMaxhits; i++ ) {
    fRAW_id[i]      = 9999;  
    fRAW_time[i]    = 9999;
    fRAW_Edep[i]    = 9999;
    fRAW_pdg[i]     = 9999;
    fRAW_mass[i]    = 9999;
    fRAW_mom[i]     = 9999;
    fRAW_px[i]      = 9999;
    fRAW_py[i]      = 9999;
    fRAW_pz[i]      = 9999;
    fRAW_xpre[i]    = 9999;
    fRAW_ypre[i]    = 9999;
    fRAW_zpre[i]    = 9999;
    fRAW_xpost[i]   = 9999;
    fRAW_ypost[i]   = 9999;
    fRAW_zpost[i]   = 9999;
    fRAW_Energy[i]  = 9999;
  }
//std::cout<<"End of ZeroArray"<<std::endl;
}

//---------------------------------------------------------------------------

void AnalysisManager::FillArray( Int_t hitn ) 
{
//std::cout<<"The arrays are about to be filled"<<std::endl;
    fRAW_Nhits++;
    fRAW_id[hitn]     = (Int_t)fStepid;
    fRAW_pdg[hitn]    = (Int_t)fSteppdef->GetPDGEncoding();
    fRAW_mass[hitn]   = (Float_t)fSteppdef->GetPDGMass();
    fRAW_time[hitn]   = (Float_t)fSteptime;                                   
    fRAW_mom[hitn]    = (Float_t)fStepp3.mag();                             
    fRAW_px[hitn]     = (Float_t)fStepp3.getX();                             
    fRAW_py[hitn]     = (Float_t)fStepp3.getY();                             
    fRAW_pz[hitn]     = (Float_t)fStepp3.getZ();                             
    fRAW_xpre[hitn]   = (Float_t)fSteppospre.getX();                             
    fRAW_ypre[hitn]   = (Float_t)fSteppospre.getY();                             
    fRAW_zpre[hitn]   = (Float_t)fSteppospre.getZ();                             
    fRAW_xpost[hitn]  = (Float_t)fSteppospost.getX();                             
    fRAW_ypost[hitn]  = (Float_t)fSteppospost.getY();                             
    fRAW_zpost[hitn]  = (Float_t)fSteppospost.getZ();                             
    fRAW_Edep[hitn]   = (Float_t)fStepedep;
    fRAW_Energy[hitn] = TMath::Sqrt( fRAW_mom[hitn]*fRAW_mom[hitn] 
				     + fRAW_mass[hitn]*fRAW_mass[hitn] );

    fRAW_xpre[hitn]   = (fRAW_xpre[hitn] + fRAW_xpost[hitn])/2.;
    fRAW_ypre[hitn]   = (fRAW_ypre[hitn] + fRAW_ypost[hitn])/2.;
    fRAW_zpre[hitn]   = (fRAW_zpre[hitn] + fRAW_zpost[hitn])/2.;
//std::cout<<"The arrays have been filled"<<std::endl;
}

//---------------------------------------------------------------------------

void AnalysisManager::FillTree()
{
//std::cout<<"Tree is about to be filled with primary variables"<<std::endl;
  // Primary Variables
  fPTime  = (Float_t)fPTime;
  fPth    = (Float_t)fPdir.getTheta();                         
  fPph    = (Float_t)fPdir.getPhi();                                                      
  fPEne   = (Float_t)fPEne;                         
  fXvtx   = (Float_t)fXvtx;                         
  fYvtx   = (Float_t)fYvtx;                         
  fZvtx   = (Float_t)fZvtx;                         
  fPpdg   = (Int_t)  fPPDef->GetPDGEncoding();

  fROOTtree->Fill();
  //for (G4int i=0; i<3000;i++){
  //	  std::cout << i << ": " << fNphotons[i] << std::endl;
  //}
  //std::cout<<"Tree has been filled with primary variables"<<std::endl;
}

//---------------------------------------------------------------------------
