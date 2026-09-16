#ifndef BNCT_STACKING_ACTION_HH
#define BNCT_STACKING_ACTION_HH

#include "G4UserStackingAction.hh"

class TrackingOptions;

class StackingAction : public G4UserStackingAction {
public:
  explicit StackingAction(const TrackingOptions* options);
  ~StackingAction() override = default;

  G4ClassificationOfNewTrack ClassifyNewTrack(const G4Track* track) override;

private:
  const TrackingOptions* fOptions;
};

#endif
