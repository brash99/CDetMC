#include <stdexcept>
#include <iostream>
#include "globals.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "AnalysisManager.hh"


#include "G4RunManagerFactory.hh"
#include "G4String.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"

//---------------------------------------------------------------------------

int main(int argc, char** argv)
{
    // detect interactive mode (if no arguments) and define UI session
    G4UIExecutive* ui = nullptr;
    if(argc == 1)
        ui = new G4UIExecutive(argc, argv);

    auto runManager = G4RunManagerFactory::CreateRunManager();

    PhysicsList*  phys       = new PhysicsList();
    runManager->SetUserInitialization(phys);
  
    DetectorConstruction*   detCon     = new DetectorConstruction();
    AnalysisManager*        anaManager = new AnalysisManager(detCon);
    runManager->SetUserInitialization(detCon);


    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();

    G4UImanager* UImanager = G4UImanager::GetUIpointer();

    if(ui)
    {
        // interactive mode
        UImanager->ApplyCommand("/control/execute macros/vis.mac");
        ui->SessionStart();
        delete ui;
    }
    else
    {
        // batch mode
        G4String command  = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }

    if (visManager) delete visManager;
    delete anaManager;
    delete runManager;

    return 0;
}

//---------------------------------------------------------------------------
