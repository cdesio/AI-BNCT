#ifndef BNCT_RUN_ACTION_HH
#define BNCT_RUN_ACTION_HH

#include "G4UserRunAction.hh"
#include "RunConfig.hh"

class RunAction : public G4UserRunAction {
public:
  explicit RunAction(const RunConfig& config);
  ~RunAction() override = default;

  void BeginOfRunAction(const G4Run*) override;
  void EndOfRunAction(const G4Run*) override;

private:
  void CreateNtuples();

  const RunConfig& fConfig;
};

#endif
