#include "EventAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "VoxelGrid.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4SystemOfUnits.hh"

EventAction::EventAction(const RunConfig& config,
                         const DetectorConstruction* detector,
                         const PrimaryGeneratorAction* primary)
  : G4UserEventAction(),
    fConfig(config),
    fDetector(detector),
    fPrimary(primary)
{}

void EventAction::BeginOfEventAction(const G4Event* event)
{
  auto* analysis = G4AnalysisManager::Instance();
  const auto& grid = fDetector->GetGrid();
  const auto& pos = fPrimary->GetPosition();
  const auto& dir = fPrimary->GetDirection();

  analysis->FillNtupleIColumn(0, 0, event->GetEventID());
  analysis->FillNtupleIColumn(0, 1, fConfig.seed);
  analysis->FillNtupleSColumn(0, 2, fPrimary->GetParticleName());
  analysis->FillNtupleDColumn(0, 3, fPrimary->GetEnergy() / eV);
  analysis->FillNtupleDColumn(0, 4, pos.x() / nm);
  analysis->FillNtupleDColumn(0, 5, pos.y() / nm);
  analysis->FillNtupleDColumn(0, 6, pos.z() / nm);
  analysis->FillNtupleDColumn(0, 7, dir.x());
  analysis->FillNtupleDColumn(0, 8, dir.y());
  analysis->FillNtupleDColumn(0, 9, dir.z());
  analysis->FillNtupleDColumn(0, 10, grid.GetVoxelSize() / nm);
  analysis->FillNtupleIColumn(0, 11, grid.GetNx());
  analysis->FillNtupleIColumn(0, 12, grid.GetNy());
  analysis->FillNtupleIColumn(0, 13, grid.GetNz());
  analysis->AddNtupleRow(0);
}
