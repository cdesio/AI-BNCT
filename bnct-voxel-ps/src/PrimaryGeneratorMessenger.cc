#include "PrimaryGeneratorMessenger.hh"

#include "PrimaryGeneratorAction.hh"

#include "G4UIcmdWith3Vector.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIdirectory.hh"

PrimaryGeneratorMessenger::PrimaryGeneratorMessenger(PrimaryGeneratorAction* generator)
  : fGenerator(generator),
    fDirectory(new G4UIdirectory("/primary/")),
    fParticleCommand(nullptr),
    fEnergyCommand(nullptr),
    fPositionCommand(nullptr),
    fDirectionCommand(nullptr)
{
  fDirectory->SetGuidance("Primary particle controls.");

  fParticleCommand = new G4UIcmdWithAString("/primary/particle", this);
  fParticleCommand->SetGuidance("Set primary particle: alpha, lithium+++, lithium++, lithium+, lithium.");
  fParticleCommand->SetParameterName("particle", false);
  fParticleCommand->SetDefaultValue("alpha");
  fParticleCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEnergyCommand = new G4UIcmdWithADoubleAndUnit("/primary/energy", this);
  fEnergyCommand->SetGuidance("Set primary kinetic energy.");
  fEnergyCommand->SetParameterName("energy", false);
  fEnergyCommand->SetRange("energy>=0.");
  fEnergyCommand->SetDefaultValue(1.47);
  fEnergyCommand->SetDefaultUnit("MeV");
  fEnergyCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fPositionCommand = new G4UIcmdWith3VectorAndUnit("/primary/position", this);
  fPositionCommand->SetGuidance("Set primary position.");
  fPositionCommand->SetParameterName("x", "y", "z", false);
  fPositionCommand->SetDefaultUnit("nm");
  fPositionCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fDirectionCommand = new G4UIcmdWith3Vector("/primary/direction", this);
  fDirectionCommand->SetGuidance("Set primary direction.");
  fDirectionCommand->SetParameterName("dx", "dy", "dz", false);
  fDirectionCommand->AvailableForStates(G4State_PreInit, G4State_Idle);
}

PrimaryGeneratorMessenger::~PrimaryGeneratorMessenger()
{
  delete fDirectionCommand;
  delete fPositionCommand;
  delete fEnergyCommand;
  delete fParticleCommand;
  delete fDirectory;
}

void PrimaryGeneratorMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  if (command == fParticleCommand) {
    fGenerator->SetParticleName(value);
  } else if (command == fEnergyCommand) {
    fGenerator->SetEnergy(fEnergyCommand->GetNewDoubleValue(value));
  } else if (command == fPositionCommand) {
    fGenerator->SetPosition(fPositionCommand->GetNew3VectorValue(value));
  } else if (command == fDirectionCommand) {
    fGenerator->SetDirection(fDirectionCommand->GetNew3VectorValue(value));
  }
}
