#include "PrimaryGeneratorAction.hh"

#include "PrimaryGeneratorMessenger.hh"

#include "G4Alpha.hh"
#include "G4Event.hh"
#include "G4Exception.hh"
#include "G4IonTable.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(),
    fGun(new G4ParticleGun(1)),
    fMessenger(new PrimaryGeneratorMessenger(this)),
    fParticleName("alpha"),
    fEnergy(1.47 * MeV),
    fPosition(0.0, 0.0, 0.0),
    fDirection(0.0, 0.0, 1.0)
{
  ConfigureGun();
}

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fMessenger;
  delete fGun;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  ConfigureGun();
  fGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::SetParticleName(const G4String& name)
{
  fParticleName = name;
}

void PrimaryGeneratorAction::SetEnergy(G4double energy)
{
  fEnergy = energy;
}

void PrimaryGeneratorAction::SetPosition(const G4ThreeVector& position)
{
  fPosition = position;
}

void PrimaryGeneratorAction::SetDirection(const G4ThreeVector& direction)
{
  fDirection = direction.unit();
}

void PrimaryGeneratorAction::ConfigureGun()
{
  G4ParticleDefinition* particle = nullptr;
  if (fParticleName == "alpha") {
    particle = G4Alpha::AlphaDefinition();
  } else if (fParticleName.find("lithium") != G4String::npos) {
    particle = G4IonTable::GetIonTable()->GetIon(3, 7, 0.0);
  } else {
    particle = G4ParticleTable::GetParticleTable()->FindParticle(fParticleName);
  }

  if (!particle) {
    G4ExceptionDescription msg;
    msg << "Unknown primary particle: " << fParticleName;
    G4Exception("PrimaryGeneratorAction::ConfigureGun", "UnknownParticle",
                FatalException, msg);
  }

  fGun->SetParticleDefinition(particle);
  if (fParticleName == "lithium+++") {
    fGun->SetParticleCharge(3.0 * eplus);
  } else if (fParticleName == "lithium++") {
    fGun->SetParticleCharge(2.0 * eplus);
  } else if (fParticleName == "lithium+") {
    fGun->SetParticleCharge(1.0 * eplus);
  } else if (fParticleName == "lithium") {
    fGun->SetParticleCharge(0.0);
  }
  fGun->SetParticleEnergy(fEnergy);
  fGun->SetParticlePosition(fPosition);
  fGun->SetParticleMomentumDirection(fDirection);
}
