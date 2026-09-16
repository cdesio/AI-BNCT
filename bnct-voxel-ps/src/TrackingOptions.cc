#include "TrackingOptions.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIdirectory.hh"

TrackingOptions::TrackingOptions()
  : fKillSecondaries(false),
    fRecordSecondaries(true),
    fDirectory(new G4UIdirectory("/tracking/")),
    fKillSecondariesCommand(nullptr),
    fRecordSecondariesCommand(nullptr)
{
  fDirectory->SetGuidance("Tracking/scoring controls.");

  fKillSecondariesCommand = new G4UIcmdWithABool("/tracking/killSecondaries", this);
  fKillSecondariesCommand->SetGuidance("Kill secondary tracks before transport.");
  fKillSecondariesCommand->SetParameterName("killSecondaries", false);
  fKillSecondariesCommand->SetDefaultValue(false);
  fKillSecondariesCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fRecordSecondariesCommand = new G4UIcmdWithABool("/tracking/recordSecondaries", this);
  fRecordSecondariesCommand->SetGuidance("Record secondary steps/phase-space entries if they are transported.");
  fRecordSecondariesCommand->SetParameterName("recordSecondaries", false);
  fRecordSecondariesCommand->SetDefaultValue(true);
  fRecordSecondariesCommand->AvailableForStates(G4State_PreInit, G4State_Idle);
}

TrackingOptions::~TrackingOptions()
{
  delete fRecordSecondariesCommand;
  delete fKillSecondariesCommand;
  delete fDirectory;
}

void TrackingOptions::SetNewValue(G4UIcommand* command, G4String value)
{
  if (command == fKillSecondariesCommand) {
    fKillSecondaries = fKillSecondariesCommand->GetNewBoolValue(value);
  } else if (command == fRecordSecondariesCommand) {
    fRecordSecondaries = fRecordSecondariesCommand->GetNewBoolValue(value);
  }
}
