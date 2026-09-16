#include "ActionInitialization.hh"

#include "DetectorConstruction.hh"
#include "EventAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "StackingAction.hh"
#include "SteppingAction.hh"
#include "TrackingOptions.hh"

ActionInitialization::ActionInitialization(const RunConfig& config,
                                           const DetectorConstruction* detector,
                                           const TrackingOptions* trackingOptions)
  : G4VUserActionInitialization(),
    fConfig(config),
    fDetector(detector),
    fTrackingOptions(trackingOptions)
{}

void ActionInitialization::BuildForMaster() const
{
  SetUserAction(new RunAction(fConfig));
}

void ActionInitialization::Build() const
{
  auto* primary = new PrimaryGeneratorAction();
  SetUserAction(primary);
  SetUserAction(new RunAction(fConfig));
  SetUserAction(new EventAction(fConfig, fDetector, primary));
  SetUserAction(new SteppingAction(fConfig, fDetector, fTrackingOptions));
  SetUserAction(new StackingAction(fTrackingOptions));
}
