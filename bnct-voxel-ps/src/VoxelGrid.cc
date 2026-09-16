#include "VoxelGrid.hh"

#include "G4SystemOfUnits.hh"
#include "G4Exception.hh"

VoxelGrid::VoxelGrid()
  : fVoxelSize(300.0 * nm), fNx(4), fNy(4), fNz(40)
{}

void VoxelGrid::SetVoxelSize(G4double value)
{
  if (value <= 0.0) {
    G4Exception("VoxelGrid::SetVoxelSize", "InvalidVoxelSize", FatalException,
                "Voxel size must be positive.");
  }
  fVoxelSize = value;
}

void VoxelGrid::SetGrid(G4int nx, G4int ny, G4int nz)
{
  if (nx <= 0 || ny <= 0 || nz <= 0) {
    G4Exception("VoxelGrid::SetGrid", "InvalidGrid", FatalException,
                "All grid dimensions must be positive.");
  }
  fNx = nx;
  fNy = ny;
  fNz = nz;
}

G4int VoxelGrid::GetVoxelID(G4int ix, G4int iy, G4int iz) const
{
  return ix + fNx * (iy + fNy * iz);
}

void VoxelGrid::GetIndices(G4int copyNo, G4int& ix, G4int& iy, G4int& iz) const
{
  ix = copyNo % fNx;
  const G4int tmp = copyNo / fNx;
  iy = tmp % fNy;
  iz = tmp / fNy;
}

G4ThreeVector VoxelGrid::GetLocalCenter(G4int copyNo) const
{
  G4int ix = 0;
  G4int iy = 0;
  G4int iz = 0;
  GetIndices(copyNo, ix, iy, iz);

  const G4double x = -0.5 * GetSizeX() + (ix + 0.5) * fVoxelSize;
  const G4double y = -0.5 * GetSizeY() + (iy + 0.5) * fVoxelSize;
  const G4double z = -0.5 * GetSizeZ() + (iz + 0.5) * fVoxelSize;
  return G4ThreeVector(x, y, z);
}
