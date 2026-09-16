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
/// \file ClusteringAlgorithm.hh
/// \brief Definition of the ClusteringAlgorithm class


#ifndef ClusteringAlgorithm_H
#define ClusteringAlgorithm_H 1

#include "ClusterSBPoints_.hh"
#include "SBPoint_.hh"

#include <map>

class ClusteringAlgorithm
{
public:
  /**
   * Clustering algorithm constructor
   * @param pEps distance between damages to be in a cluster
   * @param pContinuous whether DNA geometry is continuous
   */
  ClusteringAlgorithm(int64_t  pEps, bool pContinuous);

  /**
   * Clustering algorithm destructor
   */
  ~ClusteringAlgorithm();


  /**
   * Register a DNA damage (copy, strand, source)
   * @param copy copy number of the DNA damage
   * @param strand strand of the damage (0 or 1)
   * @param source source of the damage (1 - direct, 2 - indirect)
   */
  void RegisterDamage(int64_t copy, int64_t strand, int64_t source);

  /**
   * Run clustering algorithm 
   * @return map of cluster size distributiuon
   */
  std::map<int64_t ,int64_t > RunClustering();

  /**
   * Clean all data structures
   */
  void  Purge();

  /**
   * Return the number of SSB of a given source
   * 0 - both, 1 - direct, 2 - indirect
   * @return number of SSB
   */
  int64_t  GetSSB(int64_t SBsource) const;

  /**
   * Return the number of complex SSB of a given source
   * 0 - both, 1 - direct, 2 - indirect
   * @return number of cSSB
   */
  int64_t  GetComplexSSB(int64_t SBsource) const;

  /**
   * Return the number of DSB of a given source
   * 0 - both, 1 - direct, 2 - indirect
   * @return number of DSB
   */
  int64_t  GetDSB(int64_t SBsource) const;

  /**
   * Return the number of individual breaks of a given source
   * 0 - both, 1 - direct, 2 - indirect
   * @return number of SB
   */
  int64_t  GetTotalSB(int64_t SBsource) const;

  /**
   * Calculate cluster size distribution 
   * @return map representing cluster size distribution 
   * first int : cluster size (1 = SSB)
   * second int : counts
   */
  std::map<int64_t ,int64_t > GetClusterSizeDistribution();

  /**
   * Calculate cluster size distribution for DSB only based on source
   * 0 - all, 1 - direct, 2 - indirect, 3 - hybrid, 4 - mixed
   * @return map representing cluster size distribution 
   * first int : cluster size (min 2 as DSB)
   * second int : counts
   */
  std::map<int64_t ,int64_t > GetDSBClusterSizeDistribution(int64_t SBsource);

  /**
   * Calculate distance between DSB clusters based on source
   * 0 - all, 1 - direct, 2 - indirect, 3 - hybrid, 4 - mixed
   * @return vector of distances
   */
  std::vector<int64_t> GetDSBClusterDistanceDistribution(int64_t SBsource);


private:
  /**
   * Check if a SB point can be merged to a cluster, and do it
   * @return bool
   */
  bool FindCluster(SBPoint* pPt);

  /**
   * Check if two points can be merged onto the same cluster
   * @return bool
   */
  bool AreOnTheSameCluster(int64_t ,int64_t ,int64_t , bool);

  /**
   * Merge clusters
   */
  void MergeClusters();

  /**
   * Add SSB to cluster
   */
  void IncludeUnassociatedPoints();

  int64_t  fEps;         // distance to merge SBPoints
  std::vector<SBPoint*> fpSetOfPoints;  // Data structure containing all SB points
  std::vector<ClusterSBPoints*> fpClusters;  // Data structure containing all clusters
  std::vector<ClusterSBPoints*> fpClustersDSB;  // Data structure containing all DSB clusters
  unsigned int fNextSBPointID;  // ID of the next SB point
  bool fContinuous{1}; // whether DNA structure is continuous

};

#endif
