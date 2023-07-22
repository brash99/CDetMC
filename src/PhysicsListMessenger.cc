#include "PhysicsListMessenger.hh"

#include "PhysicsList.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

//---------------------------------------------------------------------------

PhysicsListMessenger::PhysicsListMessenger(PhysicsList* pPhys)
  :fPPhysicsList(pPhys)
{   
  fPhysDir = new G4UIdirectory("/AnaBarMC/physics/");
  fPhysDir->SetGuidance("physics control.");

  fGammaCutCmd = new G4UIcmdWithADoubleAndUnit("/AnaBarMC/physics/setGCut",this);  
  fGammaCutCmd->SetGuidance("Set gamma cut.");
  fGammaCutCmd->SetParameterName("Gcut",false);
  fGammaCutCmd->SetUnitCategory("Length");
  fGammaCutCmd->SetRange("Gcut>0.0");
  fGammaCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fElectCutCmd = new G4UIcmdWithADoubleAndUnit("/AnaBarMC/physics/setECut",this);  
  fElectCutCmd->SetGuidance("Set electron cut.");
  fElectCutCmd->SetParameterName("Ecut",false);
  fElectCutCmd->SetUnitCategory("Length");
  fElectCutCmd->SetRange("Ecut>0.0");
  fElectCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fProtoCutCmd = new G4UIcmdWithADoubleAndUnit("/AnaBarMC/physics/setPCut",this);  
  fProtoCutCmd->SetGuidance("Set positron cut.");
  fProtoCutCmd->SetParameterName("Pcut",false);
  fProtoCutCmd->SetUnitCategory("Length");
  fProtoCutCmd->SetRange("Pcut>0.0");
  fProtoCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);  

  fAllCutCmd = new G4UIcmdWithADoubleAndUnit("/AnaBarMC/physics/setCuts",this);
  fAllCutCmd->SetGuidance("Set cut for all.");
  fAllCutCmd->SetParameterName("cut",false);
  fAllCutCmd->SetUnitCategory("Length");
  fAllCutCmd->SetRange("cut>0.0");
  fAllCutCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fPListCmd = new G4UIcmdWithAString("/AnaBarMC/physics/addPhysics",this);
  fPListCmd->SetGuidance("Add physics list.");
  fPListCmd->SetParameterName("PList",false);
  fPListCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fOpticalCmd = new G4UIcmdWithAnInteger("/AnaBarMC/physics/optical",this);
  fOpticalCmd->SetGuidance("Set whether optical processes should be on");
  fOpticalCmd->SetParameterName("Optical",false);
  fOpticalCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  
  fHadronicCmd = new G4UIcmdWithAnInteger("/AnaBarMC/physics/hadronic",this);
  fHadronicCmd->SetGuidance("Set whether hadronic processes should be on");
  fHadronicCmd->SetParameterName("Hadronic",false);
  fHadronicCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

}

//---------------------------------------------------------------------------

PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fPhysDir;
  delete fGammaCutCmd;
  delete fElectCutCmd;
  delete fProtoCutCmd;
  delete fAllCutCmd;
  delete fPListCmd;
  delete fOpticalCmd;
  delete fHadronicCmd;
}

//---------------------------------------------------------------------------

void PhysicsListMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if( command == fGammaCutCmd )
   { fPPhysicsList->SetCutForGamma(fGammaCutCmd->GetNewDoubleValue(newValue));}

  if( command == fElectCutCmd )
   { fPPhysicsList->SetCutForElectron(fElectCutCmd->GetNewDoubleValue(newValue));}

  if( command == fProtoCutCmd )
   { fPPhysicsList->SetCutForPositron(fProtoCutCmd->GetNewDoubleValue(newValue));}

  if( command == fAllCutCmd )
    {
      G4double cut = fAllCutCmd->GetNewDoubleValue(newValue);
      fPPhysicsList->SetCutForGamma(cut);
      fPPhysicsList->SetCutForElectron(cut);
      fPPhysicsList->SetCutForPositron(cut);
    }

  if( command == fPListCmd )
   { fPPhysicsList->AddPhysicsList(newValue);}

  if( command == fOpticalCmd )
    { fPPhysicsList->SetOpticalProcesses(fOpticalCmd->GetNewIntValue(newValue)); }

  if( command == fHadronicCmd )
    { fPPhysicsList->SetHadronicProcesses(fHadronicCmd->GetNewIntValue(newValue)); }
}

//---------------------------------------------------------------------------

