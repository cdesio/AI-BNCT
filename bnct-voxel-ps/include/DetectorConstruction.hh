#ifndef BNCT_DETECTOR_CONSTRUCTION_HH
#define BNCT_DETECTOR_CONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "VoxelGrid.hh"

class DetectorMessenger;
class G4LogicalVolume;
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction();
  ~DetectorConstruction() override;

  G4VPhysicalVolume* Construct() override;

  void SetVoxelSize(G4double value);
  void SetGrid(G4int nx, G4int ny, G4int nz);
  void SetMaxStep(G4double value);
  void SetCheckOverlaps(G4bool value);

  const VoxelGrid& GetGrid() const { return fGrid; }
  G4double GetMaxStep() const { return fMaxStep; }
  G4bool GetCheckOverlaps() const { return fCheckOverlaps; }

private:
  VoxelGrid fGrid;
  G4double fMaxStep;
  G4bool fCheckOverlaps;
  DetectorMessenger* fMessenger;
  G4LogicalVolume* fVoxelLogical;
};

#endif
