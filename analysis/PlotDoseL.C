// -----------------------------------------------------------------------------
//
// This analysis code is written in C++ to be used with CERN ROOT
// It reads in the data produced by the geant4 simulation from the 
// root file (TFile)
// Then it gets the root tree (TTree) in which the data is stored 
// and loops over each step taken by a particle in the phantom.
// It plots the energy deposited as a function of longitudinal distance
// in the phantom for healthy tissue and the tumour.
//
// P423M Lab - David Hamilton
//
// root -b -q PlotDoseL.C++
//
// -----------------------------------------------------------------------------

#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TStyle.h"

#include <iostream>
using namespace std;

void PlotDoseL( const char* finname  = "data/G4_p135.root", 
		const char* partname = "all" )
{

  // -----------------------------------------------------------------------------
  // Open root file, get root tree and associate local and tree variables
  // -----------------------------------------------------------------------------
  
  TFile *fin = new TFile( finname );
  TTree* T   = (TTree*)fin->Get("T");
  
  //Declaration of leaves types
  static const Int_t MaxHits = 10000;
  Float_t         Prim_E;
  Float_t         Prim_Th;
  Float_t         Prim_Ph;
  Int_t           Prim_pdg;
  Int_t           Phantom_Nhits;
  Int_t           Phantom_pdg[MaxHits];
  Int_t           Phantom_id[MaxHits];
  Float_t         Phantom_x[MaxHits];
  Float_t         Phantom_y[MaxHits];
  Float_t         Phantom_z[MaxHits];
  Float_t         Phantom_Ed[MaxHits];

  // Set branch addresses.
  T->SetBranchAddress("Prim_E",&Prim_E);
  T->SetBranchAddress("Prim_Th",&Prim_Th);
  T->SetBranchAddress("Prim_Ph",&Prim_Ph);
  T->SetBranchAddress("Prim_pdg",&Prim_pdg);
  T->SetBranchAddress("Phantom_Nhits",&Phantom_Nhits);
  T->SetBranchAddress("Phantom_pdg",Phantom_pdg);
  T->SetBranchAddress("Phantom_id",Phantom_id);
  T->SetBranchAddress("Phantom_x",Phantom_x);
  T->SetBranchAddress("Phantom_y",Phantom_y);
  T->SetBranchAddress("Phantom_z",Phantom_z);
  T->SetBranchAddress("Phantom_Ed",Phantom_Ed);
  
  Long64_t nentries = T->GetEntries();

  // -----------------------------------------------------------------------------
  // Open root file, get root tree and associate local and tree variables
  // -----------------------------------------------------------------------------

  float voxelsize = 1.0;
  float lmin = -130;
  float lmax = 130;
  float tmin = -10;
  float tmax = 10;
  Int_t nbins = (lmax - lmin)/voxelsize;    
  float edep0tot = 0.0;
  float edep1tot = 0.0;
  float lpos[nbins], edep0[nbins], edep1[nbins];

  for ( Int_t k = 0; k < nbins ; k++ ) {	  
    lpos[k]  = (lmin +  k * voxelsize);
    edep0[k] = 0.0;
    edep1[k] = 0.0;
  }
  
  int partflag = 0;
  if( strcmp(partname, "lepton" ) == 0 ) partflag = 1;
  if( strcmp(partname, "meson" )  == 0 ) partflag = 2;
  if( strcmp(partname, "baryon" ) == 0 ) partflag = 3;
  if( strcmp(partname, "ion" )    == 0 ) partflag = 4;

  // ---------------------------------------------------------------------------

  for ( Long64_t i = 0; i < nentries; i++ ) {
    
    T->GetEntry(i);                       // get the ith detector hit
    if( i % 1000 == 0 )
      cout << i << endl;
    
    for ( Int_t k = 0; k < nbins ; k++ ) {	  
      
      for ( Int_t j = 0; j < Phantom_Nhits ; j++ ) {   // loop over each step in phantom
	
	if( partflag == 1 && abs( Phantom_pdg[j] ) > 100 ) continue; // lepton
	if( partflag == 2 && (abs( Phantom_pdg[j] ) < 100 ||  abs( Phantom_pdg[j] ) > 1000) ) continue; // meson 
	if( partflag == 3 && (abs( Phantom_pdg[j] ) < 1000 ||  abs( Phantom_pdg[j] ) > 10000) ) continue; // baryon
	if( partflag == 4 && abs( Phantom_pdg[j] ) < 10000 ) continue; // ion
	
	if( Prim_Th > 1.57 && (Prim_Ph < -1.57 || Prim_Ph > 1.57) ) {  // top or bottom 
	  if( Phantom_x[j] > tmin && Phantom_x[j] < tmax && Phantom_z[j] > tmin && Phantom_z[j] < tmax ) {    
	    if( Phantom_y[j] > (lmin + k*voxelsize) && Phantom_y[j] < (lmin + (k+1)*voxelsize) ) {
	      if( Phantom_id[j] == 0 ) {
		edep0[k] += Phantom_Ed[j];
		edep0tot += Phantom_Ed[j];
	      }
	      else if( Phantom_id[j] == 1 ) {
		edep1[k] += Phantom_Ed[j];
		edep1tot += Phantom_Ed[j];
	      }
	    }
	  }
	}
	if( Prim_Ph == 0 && (Prim_Th < 1.57 || Prim_Th > 3.1) ) { // front or back
	  if( Phantom_x[j] > tmin && Phantom_x[j] < tmax && Phantom_y[j] > tmin && Phantom_y[j] < tmax ) {    
	    if( Phantom_z[j] > (lmin + k*voxelsize) && Phantom_z[j] < (lmin + (k+1)*voxelsize) ) {
	      if( Phantom_id[j] == 0 ) {
		edep0[k] += Phantom_Ed[j];
		edep0tot += Phantom_Ed[j];
	      }
	      else if( Phantom_id[j] == 1 ) {
		edep1[k] += Phantom_Ed[j];
		edep1tot += Phantom_Ed[j];
	      }
	    }
	  }
	}
      }
    }
  }
  
  // ---------------------------------------------------------------------------

  // Make sure the y-scale is drawn correctly
  edep0[0]       = 0.0;
  edep0[nbins-1] = 0.0;
  
  float nor_max = 0.0;
  float tum_max = 0.0;
  int   nbins_nz = 0;
  for ( Int_t k = 0; k < nbins ; k++ ) {
    edep0[k] = edep0[k]/(float)nentries;
    edep1[k] = edep1[k]/(float)nentries;	  
    if( edep0[k] > nor_max ) 
      nor_max = edep0[k];
    if( edep1[k] > tum_max ) 
      tum_max = edep1[k];
    if( edep0[k] != 0.0 || edep0[k] != 0.0 )
      nbins_nz++;
  }
  
  // ---------------------------------------------------------------------------
  // Create a canvas, make it look pretty and plot the graph
  // ---------------------------------------------------------------------------

  gStyle->SetCanvasColor(0);
  gStyle->SetFrameFillColor(0);

  gStyle->SetPadTopMargin(.05);
  gStyle->SetPadLeftMargin(.15);
  gStyle->SetPadRightMargin(.05);
  gStyle->SetPadBottomMargin(.15);

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetStripDecimals(kFALSE);

  TCanvas* c2 = new TCanvas( "c2", "", 800, 600 );
  c2->SetGrid();
  
  TGraph *g, *g1;
  if( tum_max > nor_max ) {
    g = new TGraph(nbins, lpos, edep1 );
    g->SetMarkerColor( 2 );
    g->SetLineColor( 2 );
  }
  else {
    g = new TGraph(nbins, lpos, edep0 );
    g->SetMarkerColor( 4 );
    g->SetLineColor( 4 );
  }

  g->SetMarkerStyle( 21 );
  g->SetLineWidth( 2 );
  g->GetXaxis()->SetTitleOffset(1.3);
  g->GetYaxis()->SetTitleOffset(1.5);
  g->GetXaxis()->SetTitleFont(42);
  g->GetYaxis()->SetTitleFont(42);
  g->GetXaxis()->SetTitleSize(0.045);
  g->GetYaxis()->SetTitleSize(0.045);
  g->GetXaxis()->SetLabelOffset(0.01);
  g->GetYaxis()->SetLabelOffset(0.01);
  g->GetXaxis()->SetLabelFont(42);
  g->GetYaxis()->SetLabelFont(42);
  g->GetXaxis()->SetLabelSize(0.035);
  g->GetYaxis()->SetLabelSize(0.035);
  g->GetXaxis()->SetNdivisions(505);
  g->GetYaxis()->SetNdivisions(505);
  g->Draw("AL");

  g->GetXaxis()->SetTitle("Longitudinal Distance (mm)");
  g->GetYaxis()->SetTitle("Energy Deposited (MeV)");

  if( tum_max > nor_max ) {
    g1 = new TGraph(nbins, lpos, edep0 );
    g1->SetMarkerColor( 4 );
    g1->SetLineColor( 4 );
  }
  else {
    g1 = new TGraph(nbins, lpos, edep1 );
    g1->SetMarkerColor( 2 );
    g1->SetLineColor( 2 );
  }

  g1->SetMarkerStyle( 20 );
  g1->SetLineWidth( 2 );
  g1->Draw("L");

  // ---------------------------------------------------------------------------
  
  string s1 = (string)finname;
  string outname = s1.substr(0, s1.find(".", 0));
  outname.append("L.C");
  c2->Print(outname.c_str());

  cout << endl << "----------------------------------------------------------------" << endl << endl;
  cout << "Mean voxel energy deposit per incident beam particle = " << (edep0tot+edep1tot)/(float)nentries/(float)nbins_nz << " MeV" << endl;
  cout << "Tumour-to-normal ratio (TNR)                      = " << edep1tot/edep0tot  << endl;
  cout << endl << "----------------------------------------------------------------" << endl << endl;

}
