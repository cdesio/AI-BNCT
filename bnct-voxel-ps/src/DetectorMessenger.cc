#include "DetectorMessenger.hh"

#include "DetectorConstruction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIparameter.hh"

#include <sstream>

DetectorMessenger::DetectorMessenger(DetectorConstruction* detector)
  : fDetector(detector),
    fDirectory(new G4UIdirectory("/det/")),
    fVoxelSizeCommand(nullptr),
    fMaxStepCommand(nullptr),
    fCheckOverlapsCommand(nullptr),
    fGridCommand(nullptr)
{
  fDirectory->SetGuidance("Detector and voxel grid controls.");

  fVoxelSizeCommand = new G4UIcmdWithADoubleAndUnit("/det/setVoxelSize", this);
  fVoxelSizeCommand->SetGuidance("Set cubic voxel size.");
  fVoxelSizeCommand->SetParameterName("size", false);
  fVoxelSizeCommand->SetRange("size>0.");
  fVoxelSizeCommand->SetDefaultValue(300.0);
  fVoxelSizeCommand->SetDefaultUnit("nm");
  fVoxelSizeCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fMaxStepCommand = new G4UIcmdWithADoubleAndUnit("/det/setMaxStep", this);
  fMaxStepCommand->SetGuidance("Set user max step in voxels. Use a large value to avoid extra geometric step limiting.");
  fMaxStepCommand->SetParameterName("maxStep", false);
  fMaxStepCommand->SetRange("maxStep>0.");
  fMaxStepCommand->SetDefaultValue(30.0);
  fMaxStepCommand->SetDefaultUnit("nm");
  fMaxStepCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fCheckOverlapsCommand = new G4UIcmdWithABool("/det/checkOverlaps", this);
  fCheckOverlapsCommand->SetGuidance("Enable or disable Geant4 geometry overlap checks.");
  fCheckOverlapsCommand->SetParameterName("checkOverlaps", false);
  fCheckOverlapsCommand->SetDefaultValue(false);
  fCheckOverlapsCommand->AvailableForStates(G4State_PreInit, G4State_Idle);

  fGridCommand = new G4UIcommand("/det/setGrid", this);
  fGridCommand->SetGuidance("Set grid dimensions: nx ny nz.");
  auto* nx = new G4UIparameter("nx", 'i', false);
  nx->SetParameterRange("nx>0");
  fGridCommand->SetParameter(nx);
  auto* ny = new G4UIparameter("ny", 'i', false);
  ny->SetParameterRange("ny>0");
  fGridCommand->SetParameter(ny);
  auto* nz = new G4UIparameter("nz", 'i', false);
  nz->SetParameterRange("nz>0");
  fGridCommand->SetParameter(nz);
  fGridCommand->AvailableForStates(G4State_PreInit, G4State_Idle);
}

DetectorMessenger::~DetectorMessenger()
{
  delete fGridCommand;
  delete fCheckOverlapsCommand;
  delete fMaxStepCommand;
  delete fVoxelSizeCommand;
  delete fDirectory;
}

void DetectorMessenger::SetNewValue(G4UIcommand* command, G4String value)
{
  if (command == fVoxelSizeCommand) {
    fDetector->SetVoxelSize(fVoxelSizeCommand->GetNewDoubleValue(value));
    return;
  }

  if (command == fMaxStepCommand) {
    fDetector->SetMaxStep(fMaxStepCommand->GetNewDoubleValue(value));
    return;
  }

  if (command == fCheckOverlapsCommand) {
    fDetector->SetCheckOverlaps(fCheckOverlapsCommand->GetNewBoolValue(value));
    return;
  }

  if (command == fGridCommand) {
    std::istringstream input(value);
    G4int nx = 0;
    G4int ny = 0;
    G4int nz = 0;
    input >> nx >> ny >> nz;
    fDetector->SetGrid(nx, ny, nz);
  }
}
