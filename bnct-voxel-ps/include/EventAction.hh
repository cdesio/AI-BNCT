#ifndef BNCT_EVENT_ACTION_HH
#define BNCT_EVENT_ACTION_HH

#include "G4UserEventAction.hh"
#include "RunConfig.hh"

class DetectorConstruction;
class PrimaryGeneratorAction;

class EventAction : public G4UserEventAction {
public:
  EventAction(const RunConfig& config,
              const DetectorConstruction* detector,
              const PrimaryGeneratorAction* primary);
  ~EventAction() override = default;

  void BeginOfEventAction(const G4Event* event) override;
  void EndOfEventAction(const G4Event*) override {}

private:
  const RunConfig& fConfig;
  const DetectorConstruction* fDetector;
  const PrimaryGeneratorAction* fPrimary;
};

#endif
