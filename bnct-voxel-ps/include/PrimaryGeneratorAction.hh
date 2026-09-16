#ifndef BNCT_PRIMARY_GENERATOR_ACTION_HH
#define BNCT_PRIMARY_GENERATOR_ACTION_HH

#include "G4ThreeVector.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"

class G4Event;
class G4ParticleGun;
class PrimaryGeneratorMessenger;

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction();
  ~PrimaryGeneratorAction() override;

  void GeneratePrimaries(G4Event* event) override;

  void SetParticleName(const G4String& name);
  void SetEnergy(G4double energy);
  void SetPosition(const G4ThreeVector& position);
  void SetDirection(const G4ThreeVector& direction);

  const G4String& GetParticleName() const { return fParticleName; }
  G4double GetEnergy() const { return fEnergy; }
  const G4ThreeVector& GetPosition() const { return fPosition; }
  const G4ThreeVector& GetDirection() const { return fDirection; }

private:
  void ConfigureGun();

  G4ParticleGun* fGun;
  PrimaryGeneratorMessenger* fMessenger;
  G4String fParticleName;
  G4double fEnergy;
  G4ThreeVector fPosition;
  G4ThreeVector fDirection;
};

#endif
