/// \file SBPoint.cc
/// \brief Implementation of the SBPoint class

#include "SBPoint_.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

SBPoint::SBPoint(unsigned int pId, int64_t pPos, int64_t strand, int64_t source):
fId(pId),
fPosition(pPos),
fpCluster(0),
fSBsource(source)
{
  fStrand = strand;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SBPoint::~SBPoint()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator != (const SBPoint& pPt) const
{
  return pPt.fId != fId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator == (const SBPoint& pPt) const
{
  return pPt.fId == fId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator < (const SBPoint& pPt) const
{
  return pPt.fId < fId;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

bool SBPoint::operator > (const SBPoint& pPt) const
{
  return pPt.fId > fId;
}