#include <stdexcept>
#include <iostream>
#include "globals.hh"

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UIExecutive.hh"
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

  G4UIExecutive* ui = nullptr;
  if (argc == 1) {
    ui = new G4UIExecutive(argc, argv);
  }

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

  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();

  G4UImanager* UIManager = G4UImanager::GetUIpointer();

  if (ui) {
    UIManager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UIManager->ApplyCommand(command+fileName);
  }

  delete visManager;
  delete anaManager;
  delete runManager;
  return 0;

}

//---------------------------------------------------------------------------
