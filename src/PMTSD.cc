#include "PMTSD.hh"
#include "PMTHit.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"

//---------------------------------------------------------------------------

PMTSD::PMTSD(G4String name, G4int Nelements )
  :G4VSensitiveDetector(name)
{
  collectionName.insert(name+G4String("Collection"));
  fCollection = NULL;
  fNelements  = Nelements;
  fNhits      = 0;
  fHCID       = -1;
  fhitID      = new G4int[fNelements];
  fHits       = new G4int[fNelements];
  for(G4int i=0; i<fNelements; i++) fhitID[i] = -1;
  for(G4int i=0; i<fNelements; i++) fHits[i]  = 0; 
}

//---------------------------------------------------------------------------

PMTSD::~PMTSD()
{;}

//---------------------------------------------------------------------------

void PMTSD::Initialize(G4HCofThisEvent*)
{
  fCollection = new PMTHitsCollection(SensitiveDetectorName,collectionName[0]); 

}

//---------------------------------------------------------------------------

G4bool PMTSD::ProcessHits(G4Step*,G4TouchableHistory*)
{
  return false;
}

//---------------------------------------------------------------------------

G4bool PMTSD::ProcessHits_constStep(const G4Step* aStep,
				       G4TouchableHistory* )
{
  //std::cout << "Here I am in PMTSD::ProcessHits_constStep .... " << std::endl;
  if(aStep->GetTrack()->GetDefinition() 
     != G4OpticalPhoton::OpticalPhotonDefinition()) return false;
 
  //std::cout << "Getting PMT number and time .... " << std::endl;
  G4int pmtNumber= aStep->GetPostStepPoint()->GetTouchable()
    ->GetVolume()->GetCopyNo();
  //std::cout << pmtNumber << std::endl;

  G4double pmtTime= aStep->GetPostStepPoint()->GetGlobalTime();
  //std::cout << pmtTime << std::endl;

  
  // Try to get kinetic energy
  G4Track* theTrack = aStep->GetTrack();
  //G4double energy=theTrack->GetKineticEnergy()*1.0E6;
  //std::cout << "PMT No. = " << pmtNumber << ", Photon kinetic energy = " << energy << std::endl;

  //for (G4int iii = 0; iii<10000; iii++){
  //	std::cout << " iii = " << iii << "  fhitID[iii] = " << fhitID[iii] << std::endl;
  //} 
 
  // if this PMT hasn't been hit in this event
  //std::cout << "Accessing fhitID ... " << fhitID[pmtNumber] << " " << energy << std::endl; 
  //if ( fhitID[pmtNumber] == -1 && energy > 0.1) {
  if ( fhitID[pmtNumber] == -1 ) {
    //std::cout << "First PMT hit ... pmtNumber = " << pmtNumber << " energy = " << energy << std::endl; 
    PMTHit* OpHit = new PMTHit;
    OpHit->SetPMTNumber(pmtNumber);
    OpHit->SetPMTTime(pmtTime);
    //OpHit->SetPMTKineticEnergy(0,energy);
    OpHit->IncPhotonCount();

    fhitID[pmtNumber] = fCollection->insert(OpHit) - 1;
    fHits[fNhits++] = pmtNumber;
  }
  //std::cout << "Accessing fhitID through fCollection... " << std::endl; 
  // this is not a new hit
    (*fCollection)[fhitID[pmtNumber]]->IncPhotonCount();
    G4int current_hit_number = (*fCollection)[fhitID[pmtNumber]]->GetPhotonCount();

    G4double current_time = (*fCollection)[fhitID[pmtNumber]]->GetPMTTime();
    G4double new_time = (current_time*(current_hit_number-1)+pmtTime)/current_hit_number;
    (*fCollection)[fhitID[pmtNumber]]->SetPMTTime(new_time);

    //std::cout << "Not a new hit ... pmtNumber = " << pmtNumber << " hitNumber = " << current_hit_number << " old_time = " << current_time << " new_time " << new_time << std::endl; 
    //(*fCollection)[fhitID[pmtNumber]]->SetPMTKineticEnergy(current_hit_number-1,energy);
    //std::cout << "Returning from PMTSD::ProcessHits_constStep" << std::endl; 
  return true;
}
//---------------------------------------------------------------------------

void PMTSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  if(fHCID<0)  fHCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  if(fNhits>0) HCE->AddHitsCollection(fHCID,fCollection);

  for (G4int i=0;i<fNhits;i++) {
    fhitID[fHits[i]] = -1;
    fHits[i]         = 0;
  }
  fNhits = 0;
}

//---------------------------------------------------------------------------

