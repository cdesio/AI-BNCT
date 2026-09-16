#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "G4Run.hh"

RunAction::RunAction(const RunConfig& config)
  : G4UserRunAction(), fConfig(config)
{}

void RunAction::BeginOfRunAction(const G4Run*)
{
  auto* analysis = G4AnalysisManager::Instance();
  analysis->SetDefaultFileType("root");
  analysis->SetVerboseLevel(0);
  analysis->SetNtupleDirectoryName("Ntuples");
  analysis->OpenFile(fConfig.outputStem);
  CreateNtuples();
}

void RunAction::EndOfRunAction(const G4Run*)
{
  auto* analysis = G4AnalysisManager::Instance();
  analysis->Write();
  analysis->CloseFile();
  analysis->Clear();
}

void RunAction::CreateNtuples()
{
  auto* analysis = G4AnalysisManager::Instance();
  analysis->SetFirstNtupleId(0);

  analysis->CreateNtuple("primary", "Primary event setup");
  analysis->CreateNtupleIColumn(0, "EventID");
  analysis->CreateNtupleIColumn(0, "SeedID");
  analysis->CreateNtupleSColumn(0, "ParticleType");
  analysis->CreateNtupleDColumn(0, "InitialEnergy_eV");
  analysis->CreateNtupleDColumn(0, "InitialPosX_nm");
  analysis->CreateNtupleDColumn(0, "InitialPosY_nm");
  analysis->CreateNtupleDColumn(0, "InitialPosZ_nm");
  analysis->CreateNtupleDColumn(0, "InitialDirX");
  analysis->CreateNtupleDColumn(0, "InitialDirY");
  analysis->CreateNtupleDColumn(0, "InitialDirZ");
  analysis->CreateNtupleDColumn(0, "VoxelSize_nm");
  analysis->CreateNtupleIColumn(0, "Nx");
  analysis->CreateNtupleIColumn(0, "Ny");
  analysis->CreateNtupleIColumn(0, "Nz");
  analysis->FinishNtuple(0);

  analysis->CreateNtuple("steps", "Voxel step tracking");
  analysis->CreateNtupleIColumn(1, "EventID");
  analysis->CreateNtupleIColumn(1, "SeedID");
  analysis->CreateNtupleIColumn(1, "TrackID");
  analysis->CreateNtupleIColumn(1, "ParentID");
  analysis->CreateNtupleIColumn(1, "StepNumber");
  analysis->CreateNtupleSColumn(1, "ParticleType");
  analysis->CreateNtupleIColumn(1, "PDGCode");
  analysis->CreateNtupleDColumn(1, "KEPre_eV");
  analysis->CreateNtupleDColumn(1, "KEPost_eV");
  analysis->CreateNtupleDColumn(1, "Edep_eV");
  analysis->CreateNtupleDColumn(1, "StepLength_nm");
  analysis->CreateNtupleDColumn(1, "PrePosX_nm");
  analysis->CreateNtupleDColumn(1, "PrePosY_nm");
  analysis->CreateNtupleDColumn(1, "PrePosZ_nm");
  analysis->CreateNtupleDColumn(1, "PostPosX_nm");
  analysis->CreateNtupleDColumn(1, "PostPosY_nm");
  analysis->CreateNtupleDColumn(1, "PostPosZ_nm");
  analysis->CreateNtupleIColumn(1, "VoxelID");
  analysis->CreateNtupleIColumn(1, "VoxelX");
  analysis->CreateNtupleIColumn(1, "VoxelY");
  analysis->CreateNtupleIColumn(1, "VoxelZ");
  analysis->CreateNtupleSColumn(1, "ProcessName");
  analysis->CreateNtupleSColumn(1, "VolumeName");
  analysis->FinishNtuple(1);

  analysis->CreateNtuple("phase_space", "Voxel entry phase space");
  analysis->CreateNtupleIColumn(2, "EventID");
  analysis->CreateNtupleIColumn(2, "SeedID");
  analysis->CreateNtupleIColumn(2, "TrackID");
  analysis->CreateNtupleIColumn(2, "ParentID");
  analysis->CreateNtupleSColumn(2, "ParticleType");
  analysis->CreateNtupleIColumn(2, "PDGCode");
  analysis->CreateNtupleDColumn(2, "PositionX_nm");
  analysis->CreateNtupleDColumn(2, "PositionY_nm");
  analysis->CreateNtupleDColumn(2, "PositionZ_nm");
  analysis->CreateNtupleDColumn(2, "LocalPositionX_nm");
  analysis->CreateNtupleDColumn(2, "LocalPositionY_nm");
  analysis->CreateNtupleDColumn(2, "LocalPositionZ_nm");
  analysis->CreateNtupleDColumn(2, "DirectionX");
  analysis->CreateNtupleDColumn(2, "DirectionY");
  analysis->CreateNtupleDColumn(2, "DirectionZ");
  analysis->CreateNtupleDColumn(2, "KineticEnergy_eV");
  analysis->CreateNtupleDColumn(2, "Time_ns");
  analysis->CreateNtupleIColumn(2, "VoxelID");
  analysis->CreateNtupleIColumn(2, "VoxelX");
  analysis->CreateNtupleIColumn(2, "VoxelY");
  analysis->CreateNtupleIColumn(2, "VoxelZ");
  analysis->FinishNtuple(2);
}
