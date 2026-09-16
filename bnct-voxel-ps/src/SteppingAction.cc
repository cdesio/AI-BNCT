#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "TrackingOptions.hh"
#include "VoxelGrid.hh"

#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4SystemOfUnits.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VProcess.hh"

#include <array>
#include <cmath>

SteppingAction::SteppingAction(const RunConfig& config,
                               const DetectorConstruction* detector,
                               const TrackingOptions* trackingOptions)
  : G4UserSteppingAction(),
    fConfig(config),
    fDetector(detector),
    fTrackingOptions(trackingOptions),
    fCurrentEventID(-1)
{
  if (!fConfig.binaryOutput.empty()) {
    fBinaryPhaseSpace.open(fConfig.binaryOutput, std::ios::out | std::ios::binary);
    if (!fBinaryPhaseSpace.is_open()) {
      G4ExceptionDescription message;
      message << "Could not open binary phase-space output: "
              << fConfig.binaryOutput;
      G4Exception("SteppingAction::SteppingAction", "BNCTPS001",
                  FatalException, message);
    }
  }
}

SteppingAction::~SteppingAction()
{
  if (fBinaryPhaseSpace.is_open()) {
    fBinaryPhaseSpace.close();
  }
}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!step || !step->GetTrack()) return;
  if (!fTrackingOptions->RecordSecondaries() && step->GetTrack()->GetParentID() != 0) {
    return;
  }

  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  if (!pre || !post) return;

  const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  const G4int eventID = event ? event->GetEventID() : -1;
  if (eventID != fCurrentEventID) {
    fCurrentEventID = eventID;
    fSeenVoxelEntries.clear();
  }

  if (IsVoxelStep(pre)) {
    const G4int copyNo = pre->GetTouchableHandle()->GetCopyNumber();
    FillStepRow(step, copyNo);
  }

  if (IsVoxelEntry(step)) {
    const G4int copyNo = post->GetTouchableHandle()->GetCopyNumber();
    if (MarkVoxelEntrySeen(eventID, step->GetTrack()->GetTrackID(), copyNo)) {
      FillPhaseSpaceRow(step, post, copyNo);
    }
  }
}

bool SteppingAction::IsVoxelStep(const G4StepPoint* point) const
{
  if (!point || !point->GetPhysicalVolume()) return false;
  return point->GetPhysicalVolume()->GetName() == "voxel";
}

bool SteppingAction::IsVoxelEntry(const G4Step* step) const
{
  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  if (!pre || !post) return false;
  if (!IsVoxelStep(post)) return false;

  if (!IsVoxelStep(pre)) return true;

  return pre->GetTouchableHandle()->GetCopyNumber()
         != post->GetTouchableHandle()->GetCopyNumber();
}

bool SteppingAction::MarkVoxelEntrySeen(G4int eventID, G4int trackID, G4int copyNo)
{
  if (eventID != fCurrentEventID) {
    fCurrentEventID = eventID;
    fSeenVoxelEntries.clear();
  }

  return fSeenVoxelEntries.insert({trackID, copyNo}).second;
}

G4int SteppingAction::GetDNAReplayParticleID(const G4ParticleDefinition* particle) const
{
  if (!particle) return -1;

  const auto& name = particle->GetParticleName();
  const auto pdg = particle->GetPDGEncoding();
  if (name == "e-" || pdg == 11) return 1;
  if (name == "gamma" || pdg == 22) return 2;
  if (name == "alpha" || name == "alpha+" || name == "helium" || pdg == 1000020040) return 3;
  if (name == "e+" || pdg == -11) return 11;
  if (name == "lithium+++" || name == "lithium++" || name == "lithium+"
      || name == "lithium" || pdg == 1000030070) {
    return 53;
  }

  return -1;
}

void SteppingAction::FillStepRow(const G4Step* step, G4int copyNo) const
{
  auto* analysis = G4AnalysisManager::Instance();
  const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  const auto* track = step->GetTrack();
  const auto* pre = step->GetPreStepPoint();
  const auto* post = step->GetPostStepPoint();
  const auto* particle = track->GetParticleDefinition();
  const auto& grid = fDetector->GetGrid();

  G4int ix = 0;
  G4int iy = 0;
  G4int iz = 0;
  grid.GetIndices(copyNo, ix, iy, iz);

  G4String processName = "none";
  if (post->GetProcessDefinedStep()) {
    processName = post->GetProcessDefinedStep()->GetProcessName();
  }

  analysis->FillNtupleIColumn(1, 0, event->GetEventID());
  analysis->FillNtupleIColumn(1, 1, fConfig.seed);
  analysis->FillNtupleIColumn(1, 2, track->GetTrackID());
  analysis->FillNtupleIColumn(1, 3, track->GetParentID());
  analysis->FillNtupleIColumn(1, 4, track->GetCurrentStepNumber());
  analysis->FillNtupleSColumn(1, 5, particle->GetParticleName());
  analysis->FillNtupleIColumn(1, 6, particle->GetPDGEncoding());
  analysis->FillNtupleDColumn(1, 7, pre->GetKineticEnergy() / eV);
  analysis->FillNtupleDColumn(1, 8, post->GetKineticEnergy() / eV);
  analysis->FillNtupleDColumn(1, 9, step->GetTotalEnergyDeposit() / eV);
  analysis->FillNtupleDColumn(1, 10, step->GetStepLength() / nm);
  analysis->FillNtupleDColumn(1, 11, pre->GetPosition().x() / nm);
  analysis->FillNtupleDColumn(1, 12, pre->GetPosition().y() / nm);
  analysis->FillNtupleDColumn(1, 13, pre->GetPosition().z() / nm);
  analysis->FillNtupleDColumn(1, 14, post->GetPosition().x() / nm);
  analysis->FillNtupleDColumn(1, 15, post->GetPosition().y() / nm);
  analysis->FillNtupleDColumn(1, 16, post->GetPosition().z() / nm);
  analysis->FillNtupleIColumn(1, 17, copyNo);
  analysis->FillNtupleIColumn(1, 18, ix);
  analysis->FillNtupleIColumn(1, 19, iy);
  analysis->FillNtupleIColumn(1, 20, iz);
  analysis->FillNtupleSColumn(1, 21, processName);
  analysis->FillNtupleSColumn(1, 22, pre->GetPhysicalVolume()->GetName());
  analysis->AddNtupleRow(1);
}

void SteppingAction::FillPhaseSpaceRow(const G4Step* step, const G4StepPoint* point,
                                       G4int copyNo) const
{
  auto* analysis = G4AnalysisManager::Instance();
  const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  const auto* track = step->GetTrack();
  const auto* particle = track->GetParticleDefinition();
  const auto& grid = fDetector->GetGrid();
  const auto position = point->GetPosition();
  const auto direction = point->GetMomentumDirection();
  auto local = position;
  const auto touchable = point->GetTouchableHandle();
  if (touchable && touchable->GetHistory()) {
    local = touchable->GetHistory()->GetTopTransform().TransformPoint(position);
  }

  G4int ix = 0;
  G4int iy = 0;
  G4int iz = 0;
  grid.GetIndices(copyNo, ix, iy, iz);

  analysis->FillNtupleIColumn(2, 0, event->GetEventID());
  analysis->FillNtupleIColumn(2, 1, fConfig.seed);
  analysis->FillNtupleIColumn(2, 2, track->GetTrackID());
  analysis->FillNtupleIColumn(2, 3, track->GetParentID());
  analysis->FillNtupleSColumn(2, 4, particle->GetParticleName());
  analysis->FillNtupleIColumn(2, 5, particle->GetPDGEncoding());
  analysis->FillNtupleDColumn(2, 6, position.x() / nm);
  analysis->FillNtupleDColumn(2, 7, position.y() / nm);
  analysis->FillNtupleDColumn(2, 8, position.z() / nm);
  analysis->FillNtupleDColumn(2, 9, local.x() / nm);
  analysis->FillNtupleDColumn(2, 10, local.y() / nm);
  analysis->FillNtupleDColumn(2, 11, local.z() / nm);
  analysis->FillNtupleDColumn(2, 12, direction.x());
  analysis->FillNtupleDColumn(2, 13, direction.y());
  analysis->FillNtupleDColumn(2, 14, direction.z());
  analysis->FillNtupleDColumn(2, 15, point->GetKineticEnergy() / eV);
  analysis->FillNtupleDColumn(2, 16, point->GetGlobalTime() / ns);
  analysis->FillNtupleIColumn(2, 17, copyNo);
  analysis->FillNtupleIColumn(2, 18, ix);
  analysis->FillNtupleIColumn(2, 19, iy);
  analysis->FillNtupleIColumn(2, 20, iz);
  analysis->AddNtupleRow(2);

  WriteBinaryPhaseSpaceRow(step, point, copyNo, local);
}

void SteppingAction::WriteBinaryPhaseSpaceRow(const G4Step* step,
                                              const G4StepPoint* point,
                                              G4int copyNo,
                                              const G4ThreeVector& local) const
{
  if (!fBinaryPhaseSpace.is_open()) return;

  const auto* event = G4EventManager::GetEventManager()->GetConstCurrentEvent();
  const auto* track = step->GetTrack();
  const auto* particle = track->GetParticleDefinition();
  const G4int particleID = GetDNAReplayParticleID(particle);
  if (particleID < 0) return;

  const auto position = point->GetPosition();
  const auto direction = point->GetMomentumDirection();
  std::array<G4double, 18> record = {
    local.x(),
    local.y(),
    local.z(),
    direction.x(),
    direction.y(),
    direction.z(),
    point->GetKineticEnergy(),
    static_cast<G4double>(event ? event->GetEventID() : -1),
    static_cast<G4double>(particleID),
    static_cast<G4double>(copyNo),
    point->GetGlobalTime(),
    static_cast<G4double>(particleID),
    position.x(),
    position.y(),
    position.z(),
    static_cast<G4double>(fConfig.seed),
    static_cast<G4double>(track->GetTrackID()),
    static_cast<G4double>(track->GetParentID())
  };

  fBinaryPhaseSpace.write(reinterpret_cast<const char*>(record.data()),
                          static_cast<std::streamsize>(record.size() * sizeof(G4double)));
}
