#include "TFile.h"
#include "TTree.h"
#include "TRandom3.h"
#include "TDatabasePDG.h"
#include "TPDGCode.h"
#include "TF1.h"
#include "TH1.h"
#include "TBenchmark.h"
#include "TSystem.h"
#include "TMath.h"

#include <iostream>
#include <vector>

// ------------------------------------------------------------------------------------------------

// Functions
void  InitOutput();
void  InitInput();
void  GenerateOneParticle(int fPDGCode, double xvalue, double zvalue);
void  GenerateOneSBSParticle(int iEvent, int runNumber, int nEvents);
void  GenerateOneToyParticle();

// Random number generator
TRandom3*       fRand;

// PDG Database
TDatabasePDG*   fPDG;

// Output for G4
TFile*          fROOTFile;
TTree*          fROOTTree;
TString         fOutFileName;
TString         fInFileName;
Float_t         fVx;
Float_t         fVy;
Float_t         fVz;
Float_t         fPx;
Float_t         fPy;
Float_t         fPz;
Float_t         fP;
Float_t         fM;
Float_t         fE;
Int_t           fPDGCodeTree;

// Input from SBS
TTree* tree1;
int cdet_hit;
vector<int> *sdtrack_idx;
vector<int> *pid;
vector<double> *xpos;
vector<double> *ypos;
vector<double> *zpos;
vector<double> *xmomentum;
vector<double> *ymomentum;
vector<double> *zmomentum;
vector<double> *energy;

// Sampling Functions
TH1*            fMomFlatDist;
TH1*            fMomPowDist;
TH1*            fThetaDist;
TH1*            fPhiDist;
const Float_t   fMomMin  =  0.5;
const Float_t   fMomMax  = 20.0;
const Float_t   fMomMean =  3.5;  
Float_t         fIntRatio;

// ------------------------------------------------------------------------------------------------

void GenParticles( int fPDGCode = 13, int nevents = 100, 
		    int run_number = 2000, double xvalue = 0.0, double zvalue = 0.0) 
{
  
  // Initialise random number generator
  fRand = new TRandom3( run_number );
  
  // Set up PDG Table
  fPDG             = new TDatabasePDG();
  TString pdgtable = gSystem->Getenv( "ROOTSYS" );
  pdgtable.Append( "/etc/pdg_table.txt" );
  fPDG->ReadPDGTable( pdgtable );

  // Initialize input
  TString inname;
  inname.Form("~/CDetMC/macros/gep_12Gev1mil.root");
  fInFileName = inname;
  if (fPDGCode == -1) {
	InitInput();
  }
  
  // Initialise output
  TString fname;
  fname.Form("~/CDetMC/batch/data/AnaBarMC_Gen_%d.root",run_number);
  fOutFileName = fname;
  InitOutput();

  // Initialise sampling functions
  TF1*    momPowFunc  = new TF1("momPowFunc", "x^(-2.7)", fMomMean, fMomMax );
  Float_t meanval     = momPowFunc->Eval( fMomMean );
  Char_t* meanstr     = new Char_t[sizeof("*******")+1];
  sprintf( meanstr, "%1.6f", meanval );
  TF1*    momFlatFunc = new TF1("momFlatFunc", meanstr, fMomMin, fMomMean );
  fIntRatio           = momFlatFunc->Integral( fMomMin, fMomMean ) / 
    ( momPowFunc->Integral( fMomMean, fMomMax ) + momFlatFunc->Integral(fMomMin, fMomMean ) );
  fMomFlatDist   = (TH1*)momFlatFunc->GetHistogram()->Clone("MomFlatDist");
  fMomPowDist    = (TH1*)momPowFunc->GetHistogram()->Clone("MomPowDist");

  TF1* thetaFunc = new TF1("thetaFunc", "cos(x) * cos(x)", TMath::Pi()/2., TMath::Pi()    );
  TF1* phiFunc   = new TF1("phiFunc",   "1",               0.,             TMath::TwoPi() );
  fThetaDist     = (TH1*)thetaFunc->GetHistogram()->Clone("ThetaDist");
  fPhiDist       = (TH1*)phiFunc->GetHistogram()->Clone("PhiDist");

  // Initialise counters
  int   nTotal = 0;  
  TBenchmark* bench  = new TBenchmark();
  bench->Start("Statistics");
  
  // Main event loop
  for( int i = 0; i < nevents; i++ ) 
    {
      nTotal++;
      
      if (fPDGCode == -1 ) {
	std::cout << "SBS" << std::endl;
        GenerateOneSBSParticle(i,run_number,nevents);
      } else {
              if (fPDGCode == -2) {
                GenerateOneToyParticle();
              } else {
                GenerateOneParticle(fPDGCode,xvalue, zvalue);
              }
      }

      fROOTTree->Fill();
      
      if( i % 10 == 0 )
	cout << i << " " << fPDGCode << endl;
    }
  
  // Write output and close file
  fROOTTree->Write();
  fROOTFile->Close();
  
  // Print stats
  bench->Stop("Statistics");
  cout << "\t" <<  nTotal << " total events" << endl;
  bench->Print("Statistics");
  
}

// ------------------------------------------------------------------------------------------------

void InitOutput()
{
  fROOTFile = new TFile(fOutFileName, "RECREATE", "fROOTfile", 1);
  fROOTTree = new TTree("h1", "Generator Output Tree");
  fROOTTree->SetAutoSave();
  const Int_t basket = 64000;
  
  fROOTTree->Branch("X_vtx",   &fVx, "X_vtx/F", basket );
  fROOTTree->Branch("Y_vtx",   &fVy, "Y_vtx/F", basket );
  fROOTTree->Branch("Z_vtx",   &fVz, "Z_vtx/F", basket );
  
  fROOTTree->Branch("Px_p",    &fPx, "Px_p/F",  basket );
  fROOTTree->Branch("Py_p",    &fPy, "Py_p/F",  basket );
  fROOTTree->Branch("Pz_p",    &fPz, "Pz_p/F",  basket );
  fROOTTree->Branch("En_p",    &fE,  "En_p/F",  basket );
  
  fROOTTree->Branch("Mass",    &fM,  "Mass/F",  basket );

  fROOTTree->Branch("PDG", &fPDGCodeTree, "PDG/I",  basket );

}

void InitInput()
{
        TFile *f1 = new TFile(fInFileName,"READ");
        tree1 = (TTree*)f1->Get("T");

        tree1->SetBranchAddress("Earm.CDET_Scint.hit.nhits", &cdet_hit);
        tree1->SetBranchAddress("Earm.CDET_Scint.hit.sdtridx", &sdtrack_idx);
        tree1->SetBranchAddress("SDTrack.PID",&pid);
        tree1->SetBranchAddress("SDTrack.posx",&xpos);
        tree1->SetBranchAddress("SDTrack.posy",&ypos);
        tree1->SetBranchAddress("SDTrack.posz",&zpos);
        tree1->SetBranchAddress("SDTrack.momx",&xmomentum);
        tree1->SetBranchAddress("SDTrack.momy",&ymomentum);
        tree1->SetBranchAddress("SDTrack.momz",&zmomentum);
        tree1->SetBranchAddress("SDTrack.Etot",&energy);

}

// ------------------------------------------------------------------------------------------------

void GenerateOneSBSParticle(int iEvent, int runNumber, int nEvents)
{
	int eventOffset = runNumber%100*nEvents;
	std::cout << eventOffset << " "  << iEvent << std::endl;

        tree1->GetEntry(eventOffset+iEvent);

        double angle = 29.0/180.0*3.14159265;
        double bbdist = 4.50;
        double cdetdist = 4.0735;

        if (cdet_hit>0) {
                fVx =        -(-(*zpos)[(*sdtrack_idx)[0]] * sin(angle) + (*xpos)[(*sdtrack_idx)[0]] * cos(angle))*100;
                fVy =        -((*zpos)[(*sdtrack_idx)[0]] *cos(angle) + (*xpos)[(*sdtrack_idx)[0]] * sin(angle) - cdetdist)*100;
                fVz =         -(*ypos)[(*sdtrack_idx)[0]]*100;
                fPx =   -(-(*zmomentum)[(*sdtrack_idx)[0]] * sin(angle) + (*xmomentum)[(*sdtrack_idx)[0]] * cos(angle))*1000;
                fPy =   -((*zmomentum)[(*sdtrack_idx)[0]] * cos(angle) + (*xmomentum)[(*sdtrack_idx)[0]] * sin(angle))*1000;
                fPz =    -(*ymomentum)[(*sdtrack_idx)[0]]*1000;
                fE =        (*energy)[(*sdtrack_idx)[0]]*1000;
                fPDGCodeTree = (*pid)[(*sdtrack_idx)[0]];
 
               fM = fPDG->GetParticle( fPDGCodeTree )->Mass() * 1000;
  	//	std::cout << fVx << " " << fVy << " " << fVz << std::endl;
  	//	std::cout << fPx << " " << fPy << " " << fPz << std::endl;
  	//	std::cout << fE << " " << fM << " " << fPDGCodeTree << std::endl;
  	//	std::cout << std::endl;

        }

}

void GenerateOneToyParticle()
{

  double xsize = 100.0;
  double ysize = 100.0;
  double mp = 938.2796;
  double ebeam = 11000.0;
  double bbdist = 4.50;
  double angle = 29.0*3.14159265/180.0;

  int module = int(fRand->Uniform(0.0,3.0))+1;
  int topbottom = int(fRand->Uniform(0.0,2.0))+1;

  fVz = bbdist*100.0;

  if (module == 1) {
	  if (topbottom == 1) {
	  	fVx = -xsize/2.0+xsize*fRand->Uniform(0.0,1.0)-12.5;
	  	fVy = -ysize/4.0+ysize/2.0*fRand->Uniform(0.0,1.0)-3.0*ysize/4.0;
	  } else {
	  	fVx = -xsize/2.0+xsize*fRand->Uniform(0.0,1.0)-19.5;
	  	fVy = -ysize/4.0+ysize/2.0*fRand->Uniform(0.0,1.0)-5.0*ysize/4.0;
	  }
  }
  if (module == 2) {
	  fVx = -xsize/2.0+xsize*fRand->Uniform(0.0,1.0)-5.0;
	  fVy = -ysize/2.0+ysize*fRand->Uniform(0.0,1.0);
  }
  if (module == 3) {
	  if (topbottom == 1) {
	  	fVx = -xsize/2.0+xsize*fRand->Uniform(0.0,1.0)-19.5;
	  	fVy = -ysize/4.0+ysize/2.0*fRand->Uniform(0.0,1.0)+5.0*ysize/4.0;
	  } else {
	  	fVx = -xsize/2.0+xsize*fRand->Uniform(0.0,1.0)-12.5;
	  	fVy = -ysize/4.0+ysize/2.0*fRand->Uniform(0.0,1.0)+3.0*ysize/4.0;
	  }
  }

  // Vertex positions of Event 1 in 1000 event g4sbs sample (Angelo), for testing!
  //fVx = -29.956;
  //fVy = -145.003;
  //fVz = 450.0;

  double theta_polar = acos((-fVx*sin(angle)+fVz*cos(angle))/
		  sqrt(fVx*fVx+fVy*fVy+fVz*fVz));

  fE = ebeam*mp/(mp+ebeam*(1.0-cos(theta_polar)));
  fM = fPDG->GetParticle(11)->Mass()*1000;

  fPx = fE*(fVx/sqrt(fVx*fVx+fVy*fVy+fVz*fVz));
  fPy = fE*(fVy/sqrt(fVx*fVx+fVy*fVy+fVz*fVz));
  fPz = fE*(fVz/sqrt(fVx*fVx+fVy*fVy+fVz*fVz));

  // SBS -> CDet Coordinates
  fVx = -fVx;
  double dummy = fVy;
  fVy = -(fVz-470.0);
  fVz = -dummy;
  fPx = -fPx;
  double dummy2 = fPy;
  fPy = -fPz;
  fPz = -dummy2;

  fPDGCodeTree = 11;
                
  //std::cout << module << std::endl;
  //std::cout << fVx << " " << fVy << " " << fVz << std::endl;
  //std::cout << fPx << " " << fPy << " " << fPz << std::endl;
  //std::cout << fE << " " << fM << " " << fPDGCodeTree << std::endl;
  //std::cout << std::endl;

}

void GenerateOneParticle(int fPDGCode, double xvalue, double zvalue)
{

  // Generate vertex position in cm 
  fVx = fRand->Uniform(xvalue-25.0 , xvalue+25.0 );
  fVy = 25.0;
  fVz = fRand->Uniform(zvalue-25.0 , zvalue+25.0 );
  //fVx = fRand->Uniform(xvalue, xvalue+0.00 );
  //fVy = 25.0;
  //fVz = fRand->Uniform(zvalue , zvalue+0.01 );

  // Sample Momentum Distributions (flat from min to mean, p^-2.7 from mean to max)
  if( fRand->Uniform(0.,1) < fIntRatio )
    fP = 1000. * fMomFlatDist->GetRandom();
  else
    fP = 1000. * fMomPowDist->GetRandom();

  // Sample Momentum Distributions (flat from min to mex)
  //fP = 1000.*fRand->Uniform(3.0,3.1);

  // Sample Angular Distributions (cos^2(theta) and flat phi)
  Float_t th = fThetaDist->GetRandom();
  //Float_t th = TMath::Pi()-fRand->Uniform(0.00,0.01);
  Float_t ph = fPhiDist->GetRandom();
  //Float_t th = 3.14159265;
  //Float_t ph = 0.0;
  fPx        = fP * TMath::Sin(th) * TMath::Cos(ph);
  fPz        = fP * TMath::Sin(th) * TMath::Sin(ph);
  fPy        = fP * TMath::Cos(th);
  //fPy        = fP * TMath::Sin(th) * TMath::Sin(ph);
  //fPz        = fP * TMath::Cos(th);
  fM         = fPDG->GetParticle( fPDGCode )->Mass() * 1000;
  fE         = TMath::Sqrt( (fP*fP + fM*fM) );
  fPDGCodeTree = fPDGCode;
  
}

// ------------------------------------------------------------------------------------------------
