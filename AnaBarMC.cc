#include <stdexcept>
#include <iostream>
#include "globals.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "AnalysisManager.hh"
#include "G4TrajectoryDrawByParticleID.hh"

//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif

//---------------------------------------------------------------------------

int main(int argc, char** argv)
{

  G4RunManager* runManager = new G4RunManager;
  PhysicsList*  phys       = new PhysicsList();
  runManager->SetUserInitialization(phys);
  
  DetectorConstruction*   detCon     = new DetectorConstruction();
  AnalysisManager*        anaManager = new AnalysisManager(detCon);
  runManager->SetUserInitialization(detCon);
  
  PrimaryGeneratorAction* pga        = new PrimaryGeneratorAction();
  runManager->SetUserAction(pga);
  SteppingAction*         step       = new SteppingAction();
  runManager->SetUserAction(step);
  EventAction*            event      = new EventAction( anaManager, pga );
  runManager->SetUserAction(event);

  G4UImanager * UI         = G4UImanager::GetUIpointer();
  G4VisManager* visManager = 0;

  if (argc == 1)   // Define UI session for interactive mode.
    {
      G4cout << "Interactive mode!!!!!!!!!!!!!" << G4endl;
//#ifdef G4VIS_USE
      visManager = new G4VisExecutive;
      visManager->Initialize();

//#endif
      G4UIsession * session = 0;
#ifdef G4UI_USE_TCSH
      session = new G4UIterminal(new G4UItcsh);
#else
      session = new G4UIterminal();
#endif
      session->SessionStart();
      delete session;
    }
  else           // Batch mode
    {
        G4cout << "Batch mode!!!!!!!!!!!!!" << G4endl;
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
      if( pga->GetMode()==EPGA_ROOT ) {
	G4String commandr = "/run/beamOn ";
	G4int nev = pga->GetNEvents();
	char snev[50];
	snprintf(snev, 50, "%d", nev);
	UI->ApplyCommand(commandr+snev);
      }
    }

  if(visManager) delete visManager;
  delete anaManager;
  delete runManager;

  return 0;
}

//---------------------------------------------------------------------------
