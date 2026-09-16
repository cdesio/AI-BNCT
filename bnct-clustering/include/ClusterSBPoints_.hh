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
/// \file ClusterSBPoints.hh
/// \brief Definition of the ClusterSBPoints class

#ifndef CLUSTER_SB_POINT_HH
#define CLUSTER_SB_POINT_HH

#include "SBPoint_.hh"

#include <set>
#include <vector>

/// \brief define an object of a cluster of SB Points
class ClusterSBPoints
{
public:
  /**
   * SB cluster constructor
   * @param pPoints SB points in the cluster
   * @param pContinuous whether DNA geometry is continuous
   */
  ClusterSBPoints(std::set<SBPoint *> pPoints, bool fContinuous);

  /**
   * SB cluster destructor
   */
  ~ClusterSBPoints();

  /**
   * check if cluster is a DSB
   * @return bool, is DSB
   */
  bool IsDSB() const
  {
    return fIsDoubleSB;
  }

  /**
   * check if cluster is a SSB
   * @return bool, is SSB
   */
  bool IsSSB() const
  {
    return !IsDSB();
  }

  /**
   * get number of points in cluster
   * @return cluster size
   */
  uint64_t GetSize() const
  {
    return fpRegisteredSBPoints.size();
  }

  /**
   * add a SB to the cluster
   */
  void AddSBPoint(SBPoint *pSBPoint);

  /**
   * Get centre of cluster
   * @return index of the centre of the cluster
   */
  int64_t GetBarycenter() const;

  /**
   * SB is within minimum distance of cluster centre
   * @param SBPoint SB to check
   * @param pMinBP distance to check within
   * @return bool
   */
  bool HasInBarycenter(const SBPoint *, double pMinBP);

  /**
   * Merge SB cluster with another cluster
   * @param pCluster cluster to merge with
   */
  void MergeWith(ClusterSBPoints *pCluster);

  /**
   * source of damage causing cluster
   * @return source of damage 0 unknown, 1 direct, 2 indirect, 3 hybrid, 4 mixed
   */
  int64_t GetClusterSource() { return fSource; }

  /**
   * Get strand breaks within the cluster
   * @return set of strand breaks
   */
  std::set<SBPoint *> GetRegistredSBPoints() const
  {
    return fpRegisteredSBPoints;
  }

  /**
   * Remove points from cluster
   */
  void Clear();

private:
  /**
   * Check if cluster is DSB after adding SB
   */
  void UpdateDoubleStrand();

  /**
   * Update cluster source after SB is added
   */
  void UpdateStrandSource();

  bool fIsDoubleSB;   // is a double strand break
  int64_t fSource{0}; // 0 unknown, 1 direct, 2 indirect, 3 hybrid, 4 mixed

  std::set<SBPoint *> fpRegisteredSBPoints; // set of SB in the cluster
  bool fContinuous{1}; // whether DNA structure is continuous 
};

#endif // CLUSTER_SB_POINT_HH