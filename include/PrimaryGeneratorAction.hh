#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh" 
#include "G4ThreeVector.hh"
#include "globals.hh"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TLorentzVector.h"

//---------------------------------------------------------------------------

class PrimaryGeneratorMessenger;
class G4ParticleGun;
class G4GeneralParticleSource;
class G4Event;

enum { EPGA_GPS, EPGA_ROOT};

//---------------------------------------------------------------------------

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
  PrimaryGeneratorAction();    
  ~PrimaryGeneratorAction();

  void GeneratePrimaries (G4Event*);
  void SetUpROOTInput    (TString filename);
  void SetMode(G4int mode){ fMode = mode; }
  void SetDumpBeam(G4int dump){ fDump = dump; }

  G4int GetNEvents()                  { return static_cast<G4int>(fGenTree->GetEntries()); }
  G4int GetMode()                     { return fMode; }
  G4int GetDumpBeam()                 { return fDump; }
  G4ThreeVector GetVertex()           { return G4ThreeVector(fVx,  fVy,  fVz);  }
  G4ThreeVector GetDirection()        { return G4ThreeVector(fPxp, fPyp, fPzp); }
  G4double      GetEnergy()           { return fEp; }
  G4double      GetX()           { return fVx; }
  G4double      GetY()           { return fVy; }
  G4double      GetZ()           { return fVz; }
  G4double      GetTime()             { return fTp; }
  G4ParticleDefinition* GetPrimPDef() { return fPDefinition; }

private:
  PrimaryGeneratorMessenger* fGunMessenger;

  G4ParticleGun*             fParticleGun;    
  G4GeneralParticleSource*   fParticleSource;
  G4ParticleTable*           fParticleTable;
  G4IonTable*                fIonTable; 
  G4ParticleDefinition*      fPDefinition;
  G4int                      fMode;   
  G4int                      fDump;   
  
  TFile*                     fGenFile;    
  TTree*                     fGenTree;  
  Int_t                      fNGenBranches;   
  Int_t                      fNevent;
  Int_t                      fPDG;
  Float_t                    fVx;
  Float_t                    fVy;
  Float_t                    fVz;
  Float_t                    fPxp;
  Float_t                    fPyp;
  Float_t                    fPzp;
  Float_t                    fEp;
  Float_t                    fTp;

};
#endif

//---------------------------------------------------------------------------


