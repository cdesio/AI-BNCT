#ifndef BNCT_VOXEL_GRID_HH
#define BNCT_VOXEL_GRID_HH

#include "G4ThreeVector.hh"
#include "globals.hh"

class VoxelGrid {
public:
  VoxelGrid();

  void SetVoxelSize(G4double value);
  void SetGrid(G4int nx, G4int ny, G4int nz);

  G4double GetVoxelSize() const { return fVoxelSize; }
  G4int GetNx() const { return fNx; }
  G4int GetNy() const { return fNy; }
  G4int GetNz() const { return fNz; }
  G4int GetNumberOfVoxels() const { return fNx * fNy * fNz; }

  G4double GetSizeX() const { return fNx * fVoxelSize; }
  G4double GetSizeY() const { return fNy * fVoxelSize; }
  G4double GetSizeZ() const { return fNz * fVoxelSize; }

  G4int GetVoxelID(G4int ix, G4int iy, G4int iz) const;
  void GetIndices(G4int copyNo, G4int& ix, G4int& iy, G4int& iz) const;
  G4ThreeVector GetLocalCenter(G4int copyNo) const;

private:
  G4double fVoxelSize;
  G4int fNx;
  G4int fNy;
  G4int fNz;
};

#endif
