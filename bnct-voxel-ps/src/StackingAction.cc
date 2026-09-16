#include "StackingAction.hh"

#include "TrackingOptions.hh"

#include "G4Track.hh"

StackingAction::StackingAction(const TrackingOptions* options)
  : G4UserStackingAction(), fOptions(options)
{}

G4ClassificationOfNewTrack StackingAction::ClassifyNewTrack(const G4Track* track)
{
  if (fOptions->KillSecondaries() && track && track->GetParentID() != 0) {
    return fKill;
  }
  return fUrgent;
}
