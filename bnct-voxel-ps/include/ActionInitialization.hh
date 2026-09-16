#ifndef BNCT_ACTION_INITIALIZATION_HH
#define BNCT_ACTION_INITIALIZATION_HH

#include "G4VUserActionInitialization.hh"
#include "RunConfig.hh"

class DetectorConstruction;
class TrackingOptions;

class ActionInitialization : public G4VUserActionInitialization {
public:
  ActionInitialization(const RunConfig& config,
                       const DetectorConstruction* detector,
                       const TrackingOptions* trackingOptions);
  ~ActionInitialization() override = default;

  void BuildForMaster() const override;
  void Build() const override;

private:
  const RunConfig& fConfig;
  const DetectorConstruction* fDetector;
  const TrackingOptions* fTrackingOptions;
};

#endif
