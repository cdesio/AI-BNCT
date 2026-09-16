#ifndef BNCT_PRIMARY_GENERATOR_MESSENGER_HH
#define BNCT_PRIMARY_GENERATOR_MESSENGER_HH

#include "G4UImessenger.hh"

class G4UIcmdWith3Vector;
class G4UIcmdWith3VectorAndUnit;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithAString;
class G4UIdirectory;
class PrimaryGeneratorAction;

class PrimaryGeneratorMessenger : public G4UImessenger {
public:
  explicit PrimaryGeneratorMessenger(PrimaryGeneratorAction* generator);
  ~PrimaryGeneratorMessenger() override;

  void SetNewValue(G4UIcommand* command, G4String value) override;

private:
  PrimaryGeneratorAction* fGenerator;
  G4UIdirectory* fDirectory;
  G4UIcmdWithAString* fParticleCommand;
  G4UIcmdWithADoubleAndUnit* fEnergyCommand;
  G4UIcmdWith3VectorAndUnit* fPositionCommand;
  G4UIcmdWith3Vector* fDirectionCommand;
};

#endif
