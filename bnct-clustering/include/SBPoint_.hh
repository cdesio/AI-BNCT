//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// Based on an  example provided by the Geant4-DNA collaboration (clustering)
// Any report or published results obtained using the Geant4-DNA software
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
//
// Origninal example authors: Henri Payno and Yann Perrot
//
// Calculation of direct damage moved to python script
// Additions to the example include: adding indirect damage (python script) 
// changing to bp index, adding DSB cluster size, separating damage into direct, 
// indirect, mixed, hybrid and total and DSB distance.
//
//
/// \file SBPoint.hh
/// \brief Definition of the SBPoint class

#ifndef SB_POINT_HH
#define SB_POINT_HH

#include <assert.h>
#include <stdint.h>

#include <cinttypes>

class ClusterSBPoints;
/// \brief define an object to represent a DNA strand break (SB)
class SBPoint
{
public:
  /**
   * SB point constructor
   * @param pId SB unique id
   * @param pPos sugar molecule position (index)
   * @param strand DNA strand of damage (0,1)
   * @param source damage source, 1 - direct, 2 - indirect
   */
  SBPoint(unsigned int pId, int64_t pPos, int64_t strand, int64_t source);

  /**
   * SB point destructor
   */
  ~SBPoint();


  /**
   * Get SB id 
   * @return SB id
   */
  int64_t GetID() const
  {
    return fId;
  }

  /**
   * Get SB positon 
   * @return SB position
   */
  int64_t GetPosition() const {
    return fPosition;
  }

  /**
   * Get cluster which SB is associated with 
   * @return SB cluster
   */
  ClusterSBPoints* GetCluster() const {
    return fpCluster;
  }

  /**
   * Get stand of SB  
   * @return SB strand
   */
  int64_t GetTouchedStrand() const {
    return fStrand;
  }

  /**
   * Set cluster which SB belongs to
   * @param pCluster cluster for the SB
   */
  void SetCluster(ClusterSBPoints* pCluster)
  {
    assert(pCluster); fpCluster = pCluster;
  }

  /**
   * Check if SB already belongs to a cluster
   * @return bool
   */
  bool HasCluster() const {
    return fpCluster != 0;
  }

  /**
   * Remove clsuter which SB belongs to
   */
  void CleanCluster() {
    fpCluster = 0;
  }

  /**
   * Get source of SB
   * @return SB source, 1 - direct , 2 - indirect
   */
  int64_t GetStrandSource(){return fSBsource;}

  bool operator != (const SBPoint& ) const;
  bool operator == (const SBPoint& ) const;
  bool operator < (const SBPoint& ) const;
  bool operator > (const SBPoint& ) const;

private:

  uint64_t fId;             //SB ID
  int64_t fPosition;      //Position (copy number)
  ClusterSBPoints* fpCluster;   // Associated cluster
  int64_t fStrand;                // Strand
  int64_t fSBsource{0}; // 0 unknown, 1 direct, 2 indirect


};

#endif // SB_POINT_HH