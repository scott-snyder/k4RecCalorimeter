#ifndef RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H
#define RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H

// std
#include <cstdint>
#include <map>
#include <sys/types.h>
#include <utility>
#include <vector>
#include <optional>

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GaudiKernel/ToolHandle.h"

// Key4HEP
#include "RecCaloCommon/ICaloReadNeighboursMap.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/INoiseConstTool.h"
#include "k4FWCore/DataHandle.h"
#include "k4Interface/IGeoSvc.h"

// EDM4HEP
namespace edm4hep {
class CalorimeterHit;
class CalorimeterHitCollection;
class ClusterCollection;
} // namespace edm4hep

// DD4HEP
namespace dd4hep {
namespace DDSegmentation {
  class Segmentation;
  class BitFieldCoder;
} // namespace DDSegmentation
} // namespace dd4hep

/** @class CaloTopoClusterFCCee k4RecCalorimeter/RecFCCeeCalorimeter/src/components/CaloTopoClusterFCCee.h
 *
 *  Algorithm building the topological clusters for the energy reconstruction, following ATLAS note
 *  ATL-LARG-PUB-2008-002.
 *  1. Finds the seeds for the cluster, looking for cells exceeding the signal/noise ratio that is given by "seedSigma".
 *  2. The vector of seeds are sorted energy.
 *  3. Adds the neighbouring cells to the cluster in case their signal/noise ratio is larger than the "neighbourSigma".
 *  4. The found and added neighbours function as next seeds and their neighbours are added until no more cells exceed
 * the threshold.
 *  5. In the last step the neighbours that did not exceed the threshold the first time are tested on
 * "lastNeighbourSigma". In case that a neighbour is found that has already been assigned to another cluster, both
 * clusters are merged and assigned to the "older" clusterID, this is the one originating from a higher seed energy. The
 * iteration over neighburing cellIDs is continued.
 *
 *  @author Coralie Neubueser
 *  @author Giovanni Marchiori - algorithm rewritten for significant speed-up
 */

class CaloTopoClusterFCCee : public Gaudi::Algorithm {
public:
  using CellID = dd4hep::DDSegmentation::CellID;

  CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc);

  /**
   *
   */
  virtual StatusCode initialize() override;

  virtual StatusCode execute(const EventContext&) const override;

private:
  /// internal cell representation used for clustering, to avoid cloning EDM objects repeatedly
  struct FastCell {
    CellID cellID;
    float energy;
    float x;
    float y;
    float z;
    uint8_t type; // 0=unused,1=seed,2=neighbour,3=lastNeighbour
    float SoverN;
    unsigned icoll;
    unsigned ihit;
  };
  using FastCluster = std::vector<FastCell>;
  using FastClusterMap = std::map<uint32_t, FastCluster>; // TODO make it unordered or a vector

  /**
   * @brief Map from cell ids to cell state:
   *          0 = no cell
   *          > 0: cell index + 1
   *          < 0: -cluster index - 1
   *
   * This is implemented one of two ways.  If there is an ICaloIndexer
   * available and we're dealing with more than 1% of the total cells,
   * then we store the states using a std::vector (full representation).
   * Otherwise, we use a std::unordered_map (sparse representation).
   */
  class CellsMap
  {
  public:
    /**
     * @brief Constructor.
     * @param allCells All cells being clustered.
     * @param indexer ICaloIndexer object for the detectors being handled,
     *                or nullptr if one is not available.
     */
    CellsMap (const std::vector<FastCell>& allCells,
              const k4::recCalo::ICaloIndexer* indexer);


    /**
     * @brief Look up a cell state given a cell id.
     * Returns a reference to 0 if the cell wasn't in the input set.
     * This should not be overwritten!
     */
    int32_t& find (CellID cellid)
    {
      if (m_indexer) {
        unsigned ndx = m_indexer->index (cellid);
        // Neighbour tool may return invalid cells...
        if (ndx == k4::recCalo::ICaloIndexer::INVALID)
          return m_zero;
        return m_cellVec.at(ndx);
      }
      else {
        auto it = m_cellMap.find (cellid);
        if (it == m_cellMap.end()) {
          return m_zero;
        }
        return it->second;
      }
    }

    /// Map of cellid->state used in the sparse representation.
    std::unordered_map<CellID, int32_t> m_cellMap;

    /// Dummy.  In the sparse representation, we return a reference to this
    /// for cells that are not present.
    int32_t m_zero = 0;

    /// Indexer object.  If this is non-null, we're using the full
    /// representation.
    const k4::recCalo::ICaloIndexer* m_indexer = nullptr;

    /// Vector of cell states for the full representation.
    std::vector<int32_t> m_cellVec;
  };

 
  /** Build clusters from the found seeds.
   * First the function initialises a cluster in the preClusterCollection for the seed cells,
   * then it calls the CaloTopoClusterFCCee::searchForNeighbours function to retrieve the vector of next cellIDs to add
   * and loop over to find neighbours. The iteration of search for neighbours is continued until no more neihgbours are
   * found. Then a last round of adding neighbouring cells to the cluster is run where the parameter lastNeighbourSigma
   * is applied.
   *   @param[in] seedCells, collection of seeding cells (vector of fastcells)
   *   @param[in] allCells, collection of all cells (map cellID -> fastcell)
   *   @param[in] clusters, collection of clusters to be filled by the algorithm (map of clusterID -> FastCluster)
   */
  StatusCode buildClusters(const std::vector<FastCell>& seedCells,
                           const std::vector<FastCell>& allCells,
                           CellsMap& allCellsMap,
                           FastClusterMap& clusters) const;
  /** Search for neighbours and add them to cluster collection
   *   @param[in] cellID, the cell ID for which to find the neighbours
   *   @param[in] clusterID, the current cluster ID
   *   @param[in] nSigma, the signal/noise ratio to be exceeded by the neighbouring cell to be added to cluster
   *   @param[in] allCells, map of all cells (CellID -> FastCell)
   *   @param[in] usedCells, map of used cells (CellID -> cluster ID)
   *   @param[in] clusters, map that is filled with clusterID pointing to the associated cells, in a pair of
   *              cluster index and cell collection
   *   @param[in] clusterMembers, map (cluster ID -> set of CellIDs) that is filled by the algorithm to keep track of
   * clustered cells
   *   @param[in] allowClusterMerge, bool to allow for clusters to be merged
   *   return vector of cellID of found neighbours
   */
  std::vector<CellID> searchForNeighbours(const CellID cellID,
                                          uint& clusterID,
                                          int nSigma,
                                          const std::vector<FastCell>& allCells,
                                          CellsMap& allCellsMap,
                                          FastClusterMap& clusters,
                                          std::unordered_map<uint32_t, std::unordered_set<CellID>>& clusterMembers,
                                          bool allowClusterMerge) const;

private:
  /// List of input cell collections
  Gaudi::Property<std::vector<std::string>> m_cellCollections{
      this, "cells", {}, "Names of CalorimeterHit collections to read"};
  /// the vector of input k4FWCore::DataHandles for the input cell collections
  mutable std::vector<k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> > m_cellCollectionHandles;
  // Cluster collection (output)
  mutable k4FWCore::DataHandle<edm4hep::ClusterCollection> m_clusterCollection{"clusters", Gaudi::DataHandle::Writer,
                                                                               this};
  // Cluster cells in collection (output)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_clusterCellsCollection{
      "clusterCells", Gaudi::DataHandle::Writer, this};

  Gaudi::Property<std::vector<int>> m_caloIDs{this, "calorimeterIDs", {}, "Corresponding list of calorimeter IDs"};

  /// Handle for the cells noise tool
  ToolHandle<k4::recCalo::INoiseConstTool> m_noiseTool
    {this, "noiseTool", "TopoCaloNoisyCells", "Handle for the cells noise tool"};
  /// Handle for neighbours tool
  ToolHandle<k4::recCalo::ICaloReadNeighboursMap> m_neighboursTool
    {this, "neigboursTool", "TopoCaloNeighbours", "Handle for tool to retrieve cell neighbours"};
  // flag to use a pre-calculated neighbor map
  Gaudi::Property<bool> m_useNeighborMap{this, "useNeighborMap", true, "use pre-calculated neighbor map"};
  // use GeoSvc when the neighbor map is not present
  SmartIF<IGeoSvc> m_geoSvc;
  // name of the readout: only needed if useNeighborMap is set to false
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "",
                                             "name of the readout (needed if useNeighborMap=false)"};
  // pointer to the segmentation object
  dd4hep::DDSegmentation::Segmentation* m_segmentation = nullptr;

  /// Seed threshold in sigma
  Gaudi::Property<int> m_seedSigma{this, "seedSigma", 4, "number of sigma in noise threshold"};
  /// Neighbour threshold in sigma
  Gaudi::Property<int> m_neighbourSigma{this, "neighbourSigma", 2, "number of sigma in noise threshold"};
  /// Last neighbour threshold in sigma
  Gaudi::Property<int> m_lastNeighbourSigma{this, "lastNeighbourSigma", 0, "number of sigma in noise threshold"};
  /// Cluster energy threshold
  Gaudi::Property<float> m_minClusterEnergy{this, "minClusterEnergy", 0., "minimum cluster energy"};

  /// System encoding string
  Gaudi::Property<std::string> m_systemEncoding{this, "systemEncoding", "system:4", "System encoding string"};

  /// Flag if a new output cell collection of clustered cells should be created
  Gaudi::Property<bool> m_createClusterCellCollection{this, "createClusterCellCollection", false};
  /// General decoder to encode the calorimeter sub-system to determine which
  /// positions tool to use
  std::optional<dd4hep::DDSegmentation::BitFieldCoder> m_decoder;
  int m_indexSystem;

  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc
    { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };
};
#endif /* RECFCCEECALORIMETER_CALOTOPOCLUSTERFCCEE_H */
