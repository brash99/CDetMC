#include "G4GeneralParticleSource.hh"
#include "G4Event.hh"
#include "G4String.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleDefinition.hh"
#include "Randomize.hh"

#include "TObjArray.h"
#include "TBranch.h"
#include "TString.h"

#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorMessenger.hh"

using namespace CLHEP;

//---------------------------------------------------------------------------

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction()
{
  fGunMessenger       = new PrimaryGeneratorMessenger(this);

  fGenFile            = NULL;
  fGenTree            = NULL;

  fNevent             = 1;
  fNGenBranches       = 0;
  fDump               = 0;

  fParticleTable      = G4ParticleTable::GetParticleTable();
  fIonTable           = G4IonTable::GetIonTable(); 
  fPDefinition        = NULL;

  fParticleSource     = new G4GeneralParticleSource();
  fParticleGun        = new G4ParticleGun(1);
}

//---------------------------------------------------------------------------

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  if( fGenFile ) {
    fGenFile->Close();
    delete fGenFile;
   }
  delete fParticleGun;
  delete fParticleSource;
  delete fGunMessenger;
}

//---------------------------------------------------------------------------

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  switch(fMode) {

  case EPGA_GPS:

    fParticleSource->GeneratePrimaryVertex(anEvent);

    fVx          = fParticleSource->GetParticlePosition().getX();
    fVy          = fParticleSource->GetParticlePosition().getY();
    fVz          = fParticleSource->GetParticlePosition().getZ();
    fPxp         = fParticleSource->GetParticleMomentumDirection().getX();
    fPyp         = fParticleSource->GetParticleMomentumDirection().getY();
    fPzp         = fParticleSource->GetParticleMomentumDirection().getZ();
    fEp          = fParticleSource->GetParticleEnergy()/GeV;
    fTp          = fParticleSource->GetParticleTime();
    fPDefinition = fParticleSource->GetParticleDefinition();

    break;

  case EPGA_ROOT:
    if(fGenTree) {
      
      fGenTree->GetEvent(fNevent++);

      if( fPDG == 150000 ) {
	fPDefinition = fIonTable->FindIon(6,12,0.0); 
      }
      else if( fPDG == 11 ) {
	fPDefinition = fParticleTable->FindParticle("e-");
      }
      else if( fPDG == 2212 ) {
	fPDefinition = fParticleTable->FindParticle("proton");
      }
      else if( fPDG == 13 ) {
	fPDefinition = fParticleTable->FindParticle("mu-");
      }
      else {
	fPDefinition = fParticleTable->FindParticle("gamma");
      }

      fParticleGun->SetParticlePosition          ( G4ThreeVector( fVx *cm, fVy *cm, fVz *cm) );
      fParticleGun->SetParticleDefinition        ( fPDefinition );
      fParticleGun->SetParticleMomentumDirection ( G4ThreeVector(fPxp, fPyp, fPzp).unit() );
      fParticleGun->SetParticleEnergy            ( (fEp *MeV) );
      fParticleGun->GeneratePrimaryVertex(anEvent);
      
      
      fVx          = fParticleGun->GetParticlePosition().getX()/cm;
      fVy          = fParticleGun->GetParticlePosition().getY()/cm;
      fVz          = fParticleGun->GetParticlePosition().getZ()/cm;
      fPxp         = fParticleGun->GetParticleMomentumDirection().getX();
      fPyp         = fParticleGun->GetParticleMomentumDirection().getY();
      fPzp         = fParticleGun->GetParticleMomentumDirection().getZ();
      fEp          = fParticleGun->GetParticleEnergy()/MeV;
      fTp          = fParticleGun->GetParticleTime();
      fPDefinition = fParticleGun->GetParticleDefinition();

    }
    break;
  default:
    G4cout << "Unknown mode given to PrimiaryGeneratorAction (0 for gps or 1 for root)" << G4endl;
  }			       
}


//---------------------------------------------------------------------------

void PrimaryGeneratorAction::SetUpROOTInput(TString filename)
{
  
  fMode = EPGA_ROOT;
  
  fGenFile = new TFile(filename);
  if(!fGenFile)
    G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find filename" << G4endl;
  
  fGenTree       = dynamic_cast<TTree*>(fGenFile->Get("h1"));
  if(!fGenTree)
    G4cout << "PrimaryGeneratorAction::SetUpRootInput(TString filename) - Didn't find ntuple h1" << G4endl;
    
  fNGenBranches       = fGenTree->GetNbranches();
  TObjArray* objarray = fGenTree->GetListOfBranches();

  for( Int_t i = 0; i < fNGenBranches; i++ ) {
    
    TBranch *branch = dynamic_cast<TBranch*>     (objarray->At(i));
    TString  bname  = TString( const_cast<char*> (branch->GetName()) );

    if( bname == "X_vtx" ) branch->SetAddress( &fVx );
    if( bname == "Y_vtx" ) branch->SetAddress( &fVy );
    if( bname == "Z_vtx" ) branch->SetAddress( &fVz );

    if( bname == "Px_p"  ) branch->SetAddress( &fPxp );
    if( bname == "Py_p"  ) branch->SetAddress( &fPyp );
    if( bname == "Pz_p"  ) branch->SetAddress( &fPzp );
    if( bname == "En_p"  ) branch->SetAddress( &fEp  );

    if( bname == "PDG"  ) branch->SetAddress( &fPDG );

  }
}

//---------------------------------------------------------------------------
