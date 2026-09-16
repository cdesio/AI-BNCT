#include "VoxelParameterisation.hh"

#include "VoxelGrid.hh"

#include "G4Box.hh"
#include "G4VPhysicalVolume.hh"

VoxelParameterisation::VoxelParameterisation(const VoxelGrid* grid)
  : fGrid(grid)
{}

void VoxelParameterisation::ComputeTransformation(const G4int copyNo,
                                                  G4VPhysicalVolume* physVol) const
{
  physVol->SetTranslation(fGrid->GetLocalCenter(copyNo));
  physVol->SetRotation(nullptr);
}

void VoxelParameterisation::ComputeDimensions(G4Box& box, const G4int,
                                              const G4VPhysicalVolume*) const
{
  const G4double half = 0.5 * fGrid->GetVoxelSize();
  box.SetXHalfLength(half);
  box.SetYHalfLength(half);
  box.SetZHalfLength(half);
}
