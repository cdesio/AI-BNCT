#include "PhysicsList.hh"

#include "G4Alpha.hh"
#include "G4Electron.hh"
#include "G4EmParameters.hh"
#include "G4Gamma.hh"
#include "G4GenericIon.hh"
#include "G4IonConstructor.hh"
#include "G4Positron.hh"
#include "G4PhysicsListHelper.hh"
#include "G4ProcessManager.hh"
#include "G4ProductionCutsTable.hh"
#include "G4StepLimiter.hh"
#include "G4SystemOfUnits.hh"
#include "G4hMultipleScattering.hh"
#include "G4ionIonisation.hh"
#include "G4NuclearStopping.hh"

PhysicsList::PhysicsList(G4bool transportOnly)
  : fTransportOnly(transportOnly)
{
  SetDefaultCutValue(1.0 * micrometer);
  SetVerboseLevel(0);
  UseCoupledTransportation();

  auto* emParameters = G4EmParameters::Instance();
  emParameters->SetVerbose(0);
  emParameters->SetWorkerVerbose(0);
  emParameters->SetMinEnergy(1 * keV);
  emParameters->SetMaxEnergy(20 * MeV);
  emParameters->SetNumberOfBinsPerDecade(5);
  emParameters->SetFluo(false);
  emParameters->SetAuger(false);
  emParameters->SetPixe(false);
  emParameters->SetDeexcitationIgnoreCut(false);

  G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(1 * keV, 20 * MeV);
}

void PhysicsList::ConstructParticle()
{
  G4Gamma::GammaDefinition();
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
  G4IonConstructor().ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
  AddTransportation();

  if (fTransportOnly) {
    return;
  }

  auto* helper = G4PhysicsListHelper::GetPhysicsListHelper();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while ((*particleIterator)()) {
    auto* particle = particleIterator->value();
    const auto& name = particle->GetParticleName();
    if (name == "alpha" || name == "GenericIon") {
      helper->RegisterProcess(new G4hMultipleScattering(), particle);
      helper->RegisterProcess(new G4ionIonisation(), particle);
      helper->RegisterProcess(new G4NuclearStopping(), particle);
      if (particle->GetProcessManager()) {
        particle->GetProcessManager()->AddDiscreteProcess(new G4StepLimiter());
      }
    }
  }
}
