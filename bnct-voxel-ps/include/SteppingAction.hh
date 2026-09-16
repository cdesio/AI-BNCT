#ifndef BNCT_STEPPING_ACTION_HH
#define BNCT_STEPPING_ACTION_HH

#include "G4ThreeVector.hh"
#include "G4UserSteppingAction.hh"
#include "RunConfig.hh"

#include <fstream>
#include <set>
#include <utility>

class DetectorConstruction;
class G4ParticleDefinition;
class G4Step;
class G4StepPoint;
class G4Track;
class TrackingOptions;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(const RunConfig& config,
                 const DetectorConstruction* detector,
                 const TrackingOptions* trackingOptions);
  ~SteppingAction() override;

  void UserSteppingAction(const G4Step* step) override;

private:
  bool IsVoxelStep(const G4StepPoint* point) const;
  bool IsVoxelEntry(const G4Step* step) const;
  bool MarkVoxelEntrySeen(G4int eventID, G4int trackID, G4int copyNo);
  G4int GetDNAReplayParticleID(const G4ParticleDefinition* particle) const;
  void FillStepRow(const G4Step* step, G4int copyNo) const;
  void FillPhaseSpaceRow(const G4Step* step, const G4StepPoint* point, G4int copyNo) const;
  void WriteBinaryPhaseSpaceRow(const G4Step* step, const G4StepPoint* point,
                                G4int copyNo, const G4ThreeVector& local) const;

  const RunConfig& fConfig;
  const DetectorConstruction* fDetector;
  const TrackingOptions* fTrackingOptions;
  mutable std::ofstream fBinaryPhaseSpace;
  G4int fCurrentEventID;
  std::set<std::pair<G4int, G4int>> fSeenVoxelEntries;
};

#endif
