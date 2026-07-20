#include "CaloTopoClusterFCCee.h"

// std
#include <algorithm>
#include <memory>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

#include "k4FWCore/MetadataUtils.h"
#include "k4FWCore/GaudiChecks.h"

#include "RecCaloCommon/phihelper.h"

// EDM4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/ClusterCollection.h"
#include "edm4hep/Constants.h"

// DD4hep
#include "DD4hep/Readout.h"

DECLARE_COMPONENT(CaloTopoClusterFCCee)

CaloTopoClusterFCCee::CaloTopoClusterFCCee(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc) {
  declareProperty("clusters", m_clusterCollection, "Handle for calo clusters (output collection)");
  declareProperty("clusterCells", m_clusterCellsCollection, "Handle for clusters (output collection)");
}

StatusCode CaloTopoClusterFCCee::initialize() {

  K4_GAUDI_CHECK( Algorithm::initialize() );

  K4_GAUDI_CHECK( m_indexerSvc.retrieve() );

  // create handles for input cell collections
  for (const std::string& col : m_cellCollections) {
    debug() << "Creating handle for input cell (CalorimeterHit) collection : " << col << endmsg;
    try {
      m_cellCollectionHandles.emplace_back(col, Gaudi::DataHandle::Reader, this);;
    } catch (...) {
      error() << "Error creating handle for input collection: " << col << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // use pre-calculated neighbor map i.e. TTree to retrieve neighbors
  if (m_useNeighborMap) {
    // retrieve cells neighbours tool
    if (!m_neighboursTool.retrieve()) {
      error() << "Unable to retrieve the cells neighbours tool!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  // use DDSegmentation to retrieve neighbors
  if (!m_useNeighborMap) {
    m_geoSvc = service("GeoSvc");

    if (!m_geoSvc) {
      error() << "Unable to locate Geometry Service. "
              << "Make sure you have GeoSvc in the configuration." << endmsg;
      return StatusCode::FAILURE;
    }

    if (m_geoSvc->getDetector()->readouts().find(m_readoutName) == m_geoSvc->getDetector()->readouts().end()) {
      error() << "Readout <<" << m_readoutName << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }

    // get segmentation
    m_segmentation = m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation();
  }

  // retrieve cells noise tool
  if (!m_noiseTool.retrieve()) {
    error() << "Unable to retrieve the cells noise tool!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  // setup system decoder
  m_decoder.emplace (m_systemEncoding);
  m_indexSystem = m_decoder->index("system");

  // initialise the list of metadata for the clusters
  std::vector<std::string> shapeParameterNames = {"dR_over_E"};
  k4FWCore::putCollectionParameter(m_clusterCollection.objKey(), edm4hep::labels::ShapeParameterNames,
                                   shapeParameterNames, this);

  if (m_createClusterCellCollection) {
    std::vector<int> IDs;
    for (int ID : m_caloIDs) {
      IDs.push_back(ID);
    }

    std::vector<std::string> colls;
    for (const std::string& coll : m_cellCollections) {
      colls.push_back(coll);
    }

    if (IDs.size() == colls.size()) {
      k4FWCore::putCollectionParameter(m_clusterCollection.objKey(), "inputSystemIDs", IDs, this);
      k4FWCore::putCollectionParameter(m_clusterCollection.objKey(), "inputCellCollections", colls, this);
    } else {
      warning() << "Sizes of input cell and systemID collections of tower tool are different, no metadata written"
                << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}


CaloTopoClusterFCCee::CellsMap::CellsMap (const std::vector<FastCell>& allCells,
                                          const k4::recCalo::ICaloIndexer* indexer)
{
  size_t ncells = allCells.size();
  if (indexer && ncells > 0.01 * indexer->cellIDs().size()) {
    m_indexer = indexer;
    m_cellVec.resize (indexer->cellIDs().size());
    for (size_t icell = 0; icell < ncells; ++icell) {
      unsigned ndx = indexer->index (allCells[icell].cellID);
      m_cellVec.at(ndx) = icell+1;
    }
  }
  else {
    for (size_t icell = 0; icell < ncells; ++icell) {
      m_cellMap.emplace (allCells[icell].cellID, icell+1);
    }
  }
}


StatusCode CaloTopoClusterFCCee::execute(const EventContext&) const {

  // create output collections
  edm4hep::ClusterCollection* outClusters = m_clusterCollection.createAndPut();
  edm4hep::CalorimeterHitCollection* outClusterCells = nullptr;
  if (m_createClusterCellCollection) {
    outClusterCells = m_clusterCellsCollection.createAndPut();
  }

  // get input collection with calorimeter cells and build cell cache and flat cell map
  std::vector<FastCell> allCells;
  allCells.reserve(2000000);
  std::vector<int> caloIDs;

  std::vector<const edm4hep::CalorimeterHitCollection*> colls;
  for (size_t icoll = 0; auto& hdl : m_cellCollectionHandles) {
    const edm4hep::CalorimeterHitCollection* coll = hdl.get();
    colls.push_back (coll);
    for (size_t ihit = 0; const edm4hep::CalorimeterHit& hit : *coll) {
      int caloID = m_decoder->get(hit.getCellID(), m_indexSystem);
      if (std::ranges::find (caloIDs, caloID) == caloIDs.end()) {
        caloIDs.push_back (caloID);
      }

      // create fast flat cell
      CellID cID = hit.getCellID();
      float energy = hit.getEnergy();
      auto pos = hit.getPosition();
      auto [rms, offset] = m_noiseTool->getNoisePerCell(cID);
      float sovern = (rms > 0.) ? (std::fabs(energy - offset) / rms) : 999999.;
      allCells.emplace_back (cID, energy, (float)pos.x, (float)pos.y, (float)pos.z, 0, sovern, static_cast<unsigned>(icoll), static_cast<unsigned>(ihit));
      ++ihit;
    }
    ++icoll;
  }

  // skip event if no cells to cluster
  if (allCells.empty()) {
    debug() << "No active cells, skipping event..." << endmsg;
    return StatusCode::SUCCESS;
  }
  debug() << "Number of active cells                               : " << allCells.size() << endmsg;

  // Try to find an indexer object.  Ok if null --- we'll fall back
  // to using an unordered_map.
  const k4::recCalo::ICaloIndexer* indexer = m_indexerSvc->indexer (caloIDs, true);

  CellsMap allCellsMap (allCells, indexer);

  // find seeds (cells with S/N > seedSigma)
  // and sort by energy in reversed order
  std::vector<FastCell> seedCellsVec;
  seedCellsVec.reserve(allCells.size() / 10);
  for (const FastCell& cell : allCells) {
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "cellID   = " << cell.cellID << endmsg;
    if (cell.SoverN > m_seedSigma) {
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Found seed" << endmsg;
      seedCellsVec.push_back(cell);
    }
  }
  std::sort(seedCellsVec.begin(), seedCellsVec.end(),
            [](const FastCell& a, const FastCell& b) { return a.energy > b.energy; });

  debug() << "Number of seeds found                                : " << seedCellsVec.size() << endmsg;

  // build clusters (find neighbouring cells)
  // cluster maps clusterID to a FastCluster (vector of FastCells)
  debug() << "Building clusters" << endmsg;

  FastClusterMap clusters;
  StatusCode sc = buildClusters(seedCellsVec, allCells, allCellsMap, clusters);

  if (sc.isFailure()) {
    error() << "Unable to build the clusters!" << endmsg;
    return StatusCode::FAILURE;
  }

  // keep only clusters with sufficient energy and build EDM output clusters
  debug() << "Building EDM clusters from " << clusters.size() << " clusters" << endmsg;

  double checkTotEnergy = 0.;
  double checkTotEnergyAboveThreshold = 0.;
  int clusterWithMixedCells = 0;

  for (const auto& [clusterId, cluster] : clusters) {

    double clusterEnergy = 0.;
    std::unordered_map<int, int> system;
    system.reserve(4);

    for (const FastCell& fastcell : cluster) {

      clusterEnergy += fastcell.energy;
      unsigned systemId = m_decoder->get(fastcell.cellID, m_indexSystem);
      system[int(systemId)]++;
    }

    checkTotEnergy += clusterEnergy;
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "Cluster energy: " << clusterEnergy << endmsg;
    if (clusterEnergy < m_minClusterEnergy) {
      continue;
    }

    // build cluster
    debug() << "Building cluster with ID: " << clusterId << endmsg;
    edm4hep::MutableCluster outCluster;

    // set cluster energy
    outCluster.setEnergy(clusterEnergy);
    checkTotEnergyAboveThreshold += clusterEnergy;

    // loop over the cells attached to the cluster to calculate cluster barycenter and attach cells to cluster
    double clusterPosX = 0.;
    double clusterPosY = 0.;
    double clusterPosZ = 0.;

    double sumCellPhi = 0.;
    double sumCellTheta = 0.;
    double phi0 = 0;

    double deltaR = 0.;

    std::vector<double> cellPhi;
    std::vector<double> cellTheta;
    std::vector<double> cellEnergy;

    cellPhi.reserve(cluster.size());
    cellTheta.reserve(cluster.size());
    cellEnergy.reserve(cluster.size());

    using k4::recCalo::deltaPhi, k4::recCalo::wrapToPi;

    for (const FastCell& fastcell : cluster) {

      const edm4hep::CalorimeterHit& cell =
        colls.at(fastcell.icoll)->at(fastcell.ihit);

      double energy = fastcell.energy;

      clusterPosX += fastcell.x * energy;
      clusterPosY += fastcell.y * energy;
      clusterPosZ += fastcell.z * energy;

      double phi = std::atan2(fastcell.y, fastcell.x);
      double theta = std::atan2(std::hypot(fastcell.x, fastcell.y), fastcell.z);

      cellPhi.push_back(phi);
      cellTheta.push_back(theta);
      cellEnergy.push_back(energy);

      if (sumCellPhi == 0)
        phi0 = phi;
      sumCellPhi += deltaPhi(phi, phi0) * energy;
      sumCellTheta += theta * energy;

      // attach cell
      if (m_createClusterCellCollection) {

        edm4hep::MutableCalorimeterHit newcell = cell.clone();
        newcell.setType(fastcell.type);
        outClusterCells->push_back(newcell);
        outCluster.addToHits(newcell);

      } else {
        outCluster.addToHits(cell);
      }
    }

    // calculate cluster barycenter
    if (clusterEnergy > 0.0) {
      outCluster.setPosition(
          edm4hep::Vector3f(clusterPosX / clusterEnergy, clusterPosY / clusterEnergy, clusterPosZ / clusterEnergy));

      sumCellPhi = wrapToPi(sumCellPhi / clusterEnergy + phi0);
      sumCellTheta /= clusterEnergy;

      for (size_t i = 0; i < cellEnergy.size(); ++i) {
        deltaR +=
          std::hypot(cellTheta[i] - sumCellTheta, deltaPhi(cellPhi[i], sumCellPhi)) * cellEnergy[i];
      }
      outCluster.addToShapeParameters(deltaR / clusterEnergy);
    } else {
      outCluster.setPosition(edm4hep::Vector3f(0., 0., 0.));
      outCluster.addToShapeParameters(0.0);
    }

    outClusters->push_back(outCluster);

    if (system.size() > 1)
      clusterWithMixedCells++;
  }

  if (msgLevel() <= MSG::DEBUG) {
    debug() << "Number of clusters:                                 " << outClusters->size() << endmsg;
    debug() << "Number of clusters with cells in multiple systems:  " << clusterWithMixedCells << endmsg;
    debug() << "Total energy of clusters:                           " << checkTotEnergy << endmsg;
    debug() << "Total energy of clusters above threshold:           " << checkTotEnergyAboveThreshold << endmsg;
    if (m_createClusterCellCollection) {
      debug() << "Leftover cells :                                    " << allCells.size() - outClusterCells->size()
              << endmsg;
    }
  }

  return StatusCode::SUCCESS;
}

StatusCode CaloTopoClusterFCCee::buildClusters(const std::vector<FastCell>& seedCells,
                                               const std::vector<FastCell>& allCells,
                                               CellsMap& allCellsMap,
                                               FastClusterMap& clusters) const {

  if (msgLevel() <= MSG::VERBOSE)
    verbose() << "Initial number of seeds to loop over: " << seedCells.size() << endmsg;

  std::unordered_map<uint32_t, std::unordered_set<CellID>> clusterMembers;
  clusterMembers.reserve(seedCells.size());

  // loop over every seeds in calo to build a cluster (or merge with another cluster if appropriate)
  uint32_t seedCounter = 0;
  for (const FastCell& seedCell : seedCells) {
    seedCounter++;
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "Looking at seed: " << seedCounter << endmsg;

    CellID seedId = seedCell.cellID;
    int32_t& cellState = allCellsMap.find(seedId);
    if (cellState < 0) {
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Seed already assigned to another cluster" << endmsg;
      continue;
    }

    uint32_t clusterId = seedCounter;

    // seed insertion (type = 1)
    FastCluster& cluster = clusters[clusterId];
    cluster.reserve(128);
    cluster.push_back(seedCell);
    cluster.back().type = 1;
    cellState = -clusterId-1;
    clusterMembers[clusterId].insert(seedId);

    std::vector<std::vector<CellID>> nextNeighbours(100);
    nextNeighbours[0] =
      searchForNeighbours(seedId, clusterId, m_neighbourSigma, allCells, allCellsMap, clusters, clusterMembers, true);
    if (msgLevel() <= MSG::VERBOSE)
      verbose() << "Found " << nextNeighbours[0].size() << " neighbours.." << endmsg;

    // first loop over seeds neighbours
    int it = 0;
    while (nextNeighbours[it].size() > 0) {
      it++;
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "it: " << it << endmsg;
      nextNeighbours.emplace_back(std::vector<uint64_t>{});
      for (CellID& id : nextNeighbours[it - 1]) {
        if (id == 0) {
          error() << "Building of cluster is stopped due to missing cell ID "
                     "in neighbours map!"
                  << endmsg;
          return StatusCode::FAILURE;
        }
        if (msgLevel() <= MSG::VERBOSE)
          verbose() << "Next neighbours assigned to cluster ID: " << clusterId << endmsg;
        std::vector<CellID> additionalNeighbours =
          searchForNeighbours(id, clusterId, m_neighbourSigma, allCells, allCellsMap, clusters, clusterMembers, true);
        nextNeighbours[it].insert(nextNeighbours[it].end(), additionalNeighbours.begin(), additionalNeighbours.end());
      }
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Found " << nextNeighbours[it].size() << " more neighbours.." << endmsg;
    }

    // last try with different condition on neighbours
    if (nextNeighbours[it].size() == 0) {
      // loop over all clustered cells
      FastCluster& aCluster = clusters[clusterId];
      for (size_t i = 0; i < aCluster.size(); ++i) {
        const FastCell& cell = aCluster[i];
        if (cell.type <= 2) {
          CellID cID = cell.cellID;
          if (msgLevel() <= MSG::VERBOSE)
            verbose() << "Add neighbours of " << cID << " in last round with thr = " << m_lastNeighbourSigma.value()
                      << " x sigma." << endmsg;
          std::vector<CellID> lastNeighbours = searchForNeighbours(cID, clusterId, m_lastNeighbourSigma, allCells, allCellsMap, clusters,
                                                                     clusterMembers, false);
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}

auto CaloTopoClusterFCCee::searchForNeighbours(
    const CellID cellID,
    uint& clusterID,
    int nSigma,
    const std::vector<FastCell>& allCells,
    CellsMap& allCellsMap,
    FastClusterMap& clusters,
    std::unordered_map<uint32_t, std::unordered_set<CellID>>& clusterMembers,
    bool allowClusterMerge) const -> std::vector<CellID>
{
  std::vector<CellID> additionalNeighbours;

  // retrieve neighbours
  std::vector<CellID> neighboursVec;
  if (m_useNeighborMap) {

    neighboursVec = m_neighboursTool->neighbours(cellID);

  } else {

    std::set<dd4hep::DDSegmentation::CellID> outputNeighbors;
    m_segmentation->neighbours(cellID, outputNeighbors);
    neighboursVec.assign(outputNeighbors.begin(), outputNeighbors.end());
  }

  if (neighboursVec.empty()) {
    error() << "No neighbours for cellID " << cellID << endmsg;
    return {0};
  }

  // loop over neighbours
  if (msgLevel() <= MSG::VERBOSE)
    verbose() << "For cluster: " << clusterID << " , cell " << cellID << endmsg;
  for (const CellID neighbourID : neighboursVec) {

    int32_t& cellState = allCellsMap.find(neighbourID);

    // CASE 1: unused cell -> candidate addition
    if (cellState > 0) {
      if (msgLevel() <= MSG::VERBOSE)
        verbose() << "Found neighbour with CellID: " << neighbourID << endmsg;

      const FastCell& hit = allCells[cellState-1];
      bool addNeighbour = (hit.SoverN > nSigma) || (nSigma == 0);

      if (addNeighbour) {
        if (msgLevel() <= MSG::VERBOSE)
          verbose() << "Neighbour kept, hit = " << hit.cellID << endmsg;
        int cellType = (nSigma == m_lastNeighbourSigma) ? 3 : 2;
        FastCluster& cluster = clusters[clusterID];
        cluster.push_back(hit);
        cluster.back().type = cellType;
        cellState = -clusterID-1;
        clusterMembers[clusterID].insert(neighbourID);
        additionalNeighbours.emplace_back(neighbourID);
      } else {
        if (msgLevel() <= MSG::VERBOSE)
          verbose() << "Neighbour NOT kept, hit = " << hit.cellID << endmsg;
      }
    }

    // CASE 2: already used -> possible merge
    else if (cellState < 0 && (-cellState-1) != static_cast<int>(clusterID) && allowClusterMerge) {

      uint32_t targetCluster = -cellState-1;

      FastCluster& src = clusters[clusterID];
      FastCluster& dst = clusters[targetCluster];
      std::unordered_set<CellID>& dstMembers = clusterMembers[targetCluster];

      if (msgLevel() <= MSG::VERBOSE) {
        verbose() << "Neighbour " << neighbourID << " was found in cluster " << targetCluster << ", cluster "
                  << clusterID << " will be merged!" << endmsg;
        verbose() << "Assigning all cells ( " << clusters[clusterID].size() << " ) to Cluster " << targetCluster
                  << " with ( " << clusters[targetCluster].size() << " ). " << endmsg;
      }

      // merge all cells
      for (const FastCell& c : src) {

        allCellsMap.find(c.cellID) = -targetCluster - 1;

        if (dstMembers.insert(c.cellID).second) {
          dst.push_back(c);
        }
      }

      clusters.erase(clusterID);
      clusterMembers.erase(clusterID);
      // changed clusterId -> if more neighbours are found, correct assignment
      clusterID = targetCluster;
      // found neighbour for next search
      additionalNeighbours.emplace_back(neighbourID);
      // end loop to ensure correct cluster assignment
      break;
    }
  }

  return additionalNeighbours;
}
