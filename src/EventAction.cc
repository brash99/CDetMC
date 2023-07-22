#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4ios.hh"
#include "G4SDManager.hh"
#include "G4ThreeVector.hh"

#include "EventAction.hh"
#include "DetectorHit.hh"
#include "PMTHit.hh"
#include "AnalysisManager.hh"
#include "PrimaryGeneratorAction.hh"

//---------------------------------------------------------------------------

EventAction::EventAction( AnalysisManager* ana, PrimaryGeneratorAction* pga )
  :fAnaManager(ana), fPGA(pga), fDetCollID(-1), fPMTCollID(-1)
{;}

//---------------------------------------------------------------------------

EventAction::~EventAction()
{;}

//---------------------------------------------------------------------------

void EventAction::BeginOfEventAction(const G4Event* evt)
{ 
//  std::cout<<"Starting BeginOfEventAction"<<std::endl;
  if( evt->GetEventID() == 0 )
    fAnaManager->InitOutput();
  
  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(fDetCollID < 0 ) 
    fDetCollID = SDman->GetCollectionID("DetSDCollection");
  if(fPMTCollID < 0 ) 
    fPMTCollID = SDman->GetCollectionID("PMTSDCollection");

//  std::cout<<"Ending BeginOfEventAction"<<std::endl;
}

//---------------------------------------------------------------------------

void EventAction::EndOfEventAction(const G4Event* evt)
{
  //std::cout<<"starting EndOfEventAction"<<std::endl;
  G4int event_id   = evt->GetEventID();  
  if ( event_id%10 == 0 ) G4cout <<"Event " << event_id << G4endl;

  //std::cout<<"ZeroArray"<<std::endl;
  fAnaManager->ZeroArray(); 

  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  DetectorHitsCollection *DHC = NULL;
  PMTHitsCollection      *PHC = NULL;
  G4int det_hits = 0;
  G4int pmt_hits = 0;

  //if(HCE) {std::cout<<"HCE is true"<<std::endl;}
  //else{
  //  std::cout<<"HCE is false"<<std::endl;
  //}
  if(HCE) {
    //std::cout<<"HCE is true"<<std::endl;
    if( fDetCollID >= 0 ) DHC = (DetectorHitsCollection*)(HCE->GetHC(fDetCollID));       
    if( fPMTCollID >= 0 ) PHC = (PMTHitsCollection*)(HCE->GetHC(fPMTCollID));
  }

  // Detector
//if( DHC) {std::cout<<"DHC is true"<<std::endl;}
//else{
//  std::cout<<"DHC is false"<<std::endl;
//}
  if(DHC) {
    det_hits = DHC->entries();
    //std::cout << "Number of Detector Hits at end of event = " << det_hits << std::endl;
    //std::cout<<"DHC is true"<<std::endl;
    if( det_hits != 0 ) {
      //std::cout<<"DHC has hit"<<std::endl;
      for( G4int i = 0; i < det_hits; i++) {
	DetectorHit* hit1 = static_cast<DetectorHit*>( DHC->GetHit(i) );
	fAnaManager->SetStepPDef( (G4ParticleDefinition*) hit1->GetPDef() );
	fAnaManager->SetStepP3( (G4ThreeVector) hit1->GetMom() );
	fAnaManager->SetStepPosPre( (G4ThreeVector) hit1->GetPosPre() );
	fAnaManager->SetStepPosPost( (G4ThreeVector) hit1->GetPosPost() );
	fAnaManager->SetStepTime( (G4double) hit1->GetTime() );
	fAnaManager->SetStepEdep( (G4double) hit1->GetEdep() );
	fAnaManager->SetStepID( (G4int) hit1->GetID() );
	fAnaManager->FillArray( i ); 
      //std::cout<<"DHC has hit and FillArray has been called"<<std::endl;
      }
    }
  }

  //std::cout<<"About to check the PMT for hits"<<std::endl;
  //if( pmt_hits != 0 ) {std::cout<<"pmt_hits != 0"<<std::endl;}
  //else{
  // std::cout<<"pmt_hits == 0"<<std::endl;
  //}
  // PMT
  if(PHC) {
    pmt_hits = PHC->entries();
    //std::cout << "Number of PMT Hits at end of event = " << pmt_hits << std::endl;
    if( pmt_hits != 0 ) {
      for( G4int j = 0; j < pmt_hits; j++) {
        //std::cout << "Point 1" << std::endl;
	PMTHit* hit2 = static_cast<PMTHit*>( PHC->GetHit(j) );
        //std::cout << "Point 2" << std::endl;
	//fAnaManager->SetPMTKE( hit2 );
        //std::cout << "Point 3" << std::endl;
	fAnaManager->SetPhotonCount( (G4int) hit2->GetPMTNumber(), (G4int) hit2->GetPhotonCount() );
	fAnaManager->SetPhotonTime( (G4int) hit2->GetPMTNumber(), (G4double) hit2->GetPMTTime() );
        //std::cout << j << " " << hit2->GetPMTNumber() << " "  << hit2->GetPhotonCount() << " "  << hit2->GetPMTTime() << std::endl; 
        //std::cout << "Point 4" << std::endl;
	fAnaManager->SetPMTNumber( (G4int) hit2->GetPMTNumber() );
        //std::cout << "Point 5" << std::endl;
        //std::cout<<"PMT numbers have been set"<<std::endl;
	//if (hit2->GetPMTNumber() < 2504) {
	//	std::cout << "hit " << j << " pmt number " << hit2->GetPMTNumber() << " number of photons = " << hit2->GetPhotonCount() << std::endl;
	//}
      }
    }
  }
  //std::cout<<"About to check the Primary info stuffs"<<std::endl;
  // Primary
  //if( det_hits != 0 ) {std::cout<<"det_hits != 0"<<std::endl;}
  //else{
  //  std::cout<<"det_hits == 0"<<std::endl;
  //}
  if( det_hits != 0 ) {
    fAnaManager->SetPrimaryDirection ( (G4ThreeVector)fPGA->GetDirection() );
    fAnaManager->SetPrimaryEnergy    ( (G4double)fPGA->GetEnergy() );
    fAnaManager->SetPrimaryX    ( (G4double)fPGA->GetX() );
    fAnaManager->SetPrimaryY    ( (G4double)fPGA->GetY() );
    fAnaManager->SetPrimaryZ    ( (G4double)fPGA->GetZ() );
    fAnaManager->SetPrimaryTime      ( (G4double)fPGA->GetTime() );
    fAnaManager->SetPrimaryPDef      ( (G4ParticleDefinition*)fPGA->GetPrimPDef() );
    
    fAnaManager->FillTree(); 
    //std::cout<<"Primay things have been set and tree has been filled"<<std::endl;
  }
  //std::cout<<"End of EndOfEventAction"<<std::endl;
}

//---------------------------------------------------------------------------
