#ifndef BNCT_PHYSICS_LIST_HH
#define BNCT_PHYSICS_LIST_HH

#include "G4VModularPhysicsList.hh"

class PhysicsList : public G4VModularPhysicsList {
public:
  explicit PhysicsList(G4bool transportOnly = false);
  ~PhysicsList() override = default;

  void ConstructParticle() override;
  void ConstructProcess() override;

private:
  G4bool fTransportOnly;
};

#endif
