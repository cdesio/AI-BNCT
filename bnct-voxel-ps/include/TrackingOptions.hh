#ifndef BNCT_TRACKING_OPTIONS_HH
#define BNCT_TRACKING_OPTIONS_HH

#include "G4UImessenger.hh"

class G4UIcmdWithABool;
class G4UIdirectory;

class TrackingOptions : public G4UImessenger {
public:
  TrackingOptions();
  ~TrackingOptions() override;

  void SetNewValue(G4UIcommand* command, G4String value) override;

  bool KillSecondaries() const { return fKillSecondaries; }
  bool RecordSecondaries() const { return fRecordSecondaries; }

private:
  bool fKillSecondaries;
  bool fRecordSecondaries;
  G4UIdirectory* fDirectory;
  G4UIcmdWithABool* fKillSecondariesCommand;
  G4UIcmdWithABool* fRecordSecondariesCommand;
};

#endif
