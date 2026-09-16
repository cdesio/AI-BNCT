#include <vector>
#include "ClusteringAlgorithm.hh"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>
#include <stdint.h>
#include <stdlib.h>
#include <iostream>
namespace py = pybind11;

using namespace std;

/**
 * Pybind module to calculate the clusters of DNA damage from the locations of individual strand breaks, due to direct, indirect and both effects combined
 * 
 * @param numEvt a vector containining the event number of all strand breaks
 * @param eventsListDirect a vector containining the event number of each direct strand break
 * @param copyListDirect a vector containining the copy number of each direct strand break
 * @param strandListDirect a vector containining the strand number (0 or 1) of each direct strand break
 * @param eventsListIndirect a vector containining the event number of each direct strand break
 * @param copyListIndirect a vector containining the copy number of each direct strand break
 * @param strandListIndirect a vector containining the strand number (0 or 1) of each direct strand break
 * @param fContinuous whether geometry is continuous or discontinuous
 * @return vector containing [SBresults - number of each cluster type per event, clusterFrequency - DSB cluster size (1-10) per event, clusterSizes - distance between DSB clusters per event]
 */
std::vector<vector<vector<int64_t>>> clustering(std::vector<int64_t> numEvt, std::vector<int64_t> eventsListDirect, std::vector<int64_t> copyListDirect, std::vector<int64_t> strandListDirect, std::vector<int64_t> eventsListIndirect, std::vector<int64_t> copyListIndirect, std::vector<int64_t> strandListIndirect, bool fContinuous)
{
    std::vector<vector<int64_t>> SBresults; // vector to store results - number of each cluster type per event
    std::vector<vector<int64_t>> clusterFrequency; // vector to store results - DSB cluster size (1-10) per event
    std::vector<vector<int64_t>> clusterSizes; // vector to store results - distance between DSB clusters per event

    // start clustering object for direct, indirect, and total DNA strand breaks
    ClusteringAlgorithm *fpClusteringDirect = new ClusteringAlgorithm(10, fContinuous);  
    ClusteringAlgorithm *fpClusteringIndirect = new ClusteringAlgorithm(10, fContinuous); 
    ClusteringAlgorithm *fpClusteringTotal = new ClusteringAlgorithm(10, fContinuous);    

    // Perform clustering per event
    for (auto e: numEvt)
    {
        // temporary vector to store the number of each cluster type for event e
        // 0 - event number
        // 1 - TotalSBdirect - all strand breaks from direct effect
        // 2 - SSBdirect - SSB from direct effect
        // 3 - cSSBdirect - complext SSB from direct effect
        // 4 - DSBdirect - DSB from direct effect
        // 5 - TotalSBindirect - all strand breaks from indirect effect
        // 6 - SSBindirect - SSB from indirect effect
        // 7 - cSSBindirect - complex SSB from indirect effect
        // 8 - DSBindirect - DSB from indirect effect
        // 9 - TotalSBtotal - all strand breaks from both effects
        // 10 - SSBtotal - SSB from both effects
        // 11 - cSSBtotal - complex SSB from both effects
        // 12 - DSBtotal - DSB from both effects

        std::vector<int64_t> tempResults(18, 0);

        // temporary vector to store the number of DSB clusters of a given size per event
        // 0: event number
        // 1 - 10: direct DSB with cluster size (1-10)
        // 11 - 20: indirect DSB with cluster size (1-10) 
        // 21 - 30: hybrid DSB with cluster size (1-10) 
        // 31 - 40: mixed DSB with cluster size (1-10) 
        // 41 - 50: all DSB with cluster size (1-10) 
        // 
        std::vector<int64_t> clusterFrequencyTemp(51, 0); 

        bool eventFound =false;

        // add damages to clustering object from direct effects
        for (int64_t i = 0; i < eventsListDirect.size(); ++i)
        {
            if (e == eventsListDirect[i])
            {
                fpClusteringDirect->RegisterDamage(copyListDirect[i], strandListDirect[i], 1);
                fpClusteringTotal->RegisterDamage(copyListDirect[i], strandListDirect[i], 1); 
                eventFound=true;
            }
        }

        // add damages to clustering object from indirect effects
        for (int64_t i = 0; i < eventsListIndirect.size(); ++i)
        {
            if (e == eventsListIndirect[i])
            {
                fpClusteringIndirect->RegisterDamage(copyListIndirect[i], strandListIndirect[i], 2);
                fpClusteringTotal->RegisterDamage(copyListIndirect[i], strandListIndirect[i], 2); 
                eventFound=true;
            }
        }

        if (!eventFound) continue; // if no DNA strand breaks in event e, continue

        // run clustering algorithm to obtain cluster sizes, direct damage
        std::map<int64_t, int64_t> sizeDistributionDirect = fpClusteringDirect->RunClustering();
        tempResults[0] = e;
        tempResults[1] = fpClusteringDirect->GetTotalSB(1);
        tempResults[2] = fpClusteringDirect->GetSSB(1);
        tempResults[3] = fpClusteringDirect->GetComplexSSB(1);
        tempResults[4] = fpClusteringDirect->GetDSB(1);

        fpClusteringDirect->Purge();

        sizeDistributionDirect.clear();

        // run clustering algorithm to obtain cluster sizes, indirect damage
        std::map<int64_t, int64_t> sizeDistributionIndirect = fpClusteringIndirect->RunClustering();

        tempResults[5] = fpClusteringIndirect->GetTotalSB(2);
        tempResults[6] = fpClusteringIndirect->GetSSB(2);
        tempResults[7] = fpClusteringIndirect->GetComplexSSB(2);
        tempResults[8] = fpClusteringIndirect->GetDSB(2);

        fpClusteringIndirect->Purge();

        sizeDistributionIndirect.clear();

        // run clustering algorithm to obtain cluster sizes, both damage types
        std::map<int64_t, int64_t> sizeDistribution = fpClusteringTotal->RunClustering();
        tempResults[9] = fpClusteringTotal->GetTotalSB(0);
        tempResults[10] = fpClusteringTotal->GetSSB(0);
        tempResults[11] = fpClusteringTotal->GetComplexSSB(0);
        tempResults[12] = fpClusteringTotal->GetDSB(0);

        // Complex DSBs contain more than two strand breaks in one DSB cluster.
        for (int64_t source = 0; source <= 4; ++source)
        {
            const auto sizes = fpClusteringTotal->GetDSBClusterSizeDistribution(source);
            for (const auto& [size, count] : sizes)
            {
                if (size > 2)
                {
                    const int64_t offset = source == 0 ? 4 : source - 1;
                    tempResults[13 + offset] += count;
                }
            }
        }


        // Calculate number of DSB which cluster size 0-10, direct DSB
        clusterFrequencyTemp[0] = e;
        int c = 1;
        std::map<int64_t, int64_t> sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(1); //direct

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        // Calculate number of DSB which cluster size 0-10, indirect DSB
        c = 11;
        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(2); // indirect

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        // Calculate number of DSB which cluster size 0-10, hybrid DSB
        c = 21;
        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(3); //hybrid

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        // Calculate number of DSB which cluster size 0-10, mixed DSB
        c = 31;

        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(4); //mixed

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());
        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        // Calculate number of DSB which cluster size 0-10, all DSB.
        c = 41;

        sizeDistributionDSB = fpClusteringTotal->GetDSBClusterSizeDistribution(0);

        fpClusteringTotal->GetDSBClusterDistanceDistribution(0);

        while ( (!sizeDistributionDSB.empty())&&(sizeDistributionDSB.begin()->first <=10 ))
        {
            clusterFrequencyTemp[c+sizeDistributionDSB.begin()->first-1] = sizeDistributionDSB.begin()->second; //-1 zero indexing
            sizeDistributionDSB.erase(sizeDistributionDSB.begin());

        } 
        
        sizeDistribution.clear();
        sizeDistributionDSB.clear();

        // Calculate distance between DSB clusters (in base pairs)
        std::vector<int64_t> sizeDistanceDSB = fpClusteringTotal->GetDSBClusterDistanceDistribution(0);

        fpClusteringTotal->Purge();
        
        SBresults.push_back(tempResults);
        clusterFrequency.push_back(clusterFrequencyTemp);
        clusterSizes.push_back(sizeDistanceDSB);

    }
    std::vector<vector<vector<int64_t>>> results{SBresults, clusterFrequency, clusterSizes};

    delete fpClusteringDirect;
    delete fpClusteringIndirect;
    delete fpClusteringTotal;

    return results;
}

PYBIND11_MODULE(clustering, m)
{
    m.doc() = "clustering"; // optional module docstring

    m.def("clustering", &clustering, "A function that performs clustering");
}
