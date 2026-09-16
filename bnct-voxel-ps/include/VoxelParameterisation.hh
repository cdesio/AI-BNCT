#ifndef BNCT_VOXEL_PARAMETERISATION_HH
#define BNCT_VOXEL_PARAMETERISATION_HH

#include "G4VPVParameterisation.hh"

class G4Box;
class G4VPhysicalVolume;
class VoxelGrid;

class VoxelParameterisation : public G4VPVParameterisation {
public:
  explicit VoxelParameterisation(const VoxelGrid* grid);
  ~VoxelParameterisation() override = default;

  void ComputeTransformation(const G4int copyNo,
                             G4VPhysicalVolume* physVol) const override;
  void ComputeDimensions(G4Box& box, const G4int copyNo,
                         const G4VPhysicalVolume* physVol) const override;

private:
  const VoxelGrid* fGrid;
};

#endif
