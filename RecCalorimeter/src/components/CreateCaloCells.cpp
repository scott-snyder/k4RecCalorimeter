#include "CreateCaloCells.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

// DD4hep
#include "DD4hep/DetType.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"
#include "TGeoManager.h"

// edm4hep
#include "edm4hep/CalorimeterHit.h"

#include <stdexcept>
#include <algorithm>


DECLARE_COMPONENT(CreateCaloCells)


/**
 * @brief Constructor.
 * @param full If true, we initialize for all cells.
 * @param cellids List of all possible cell ids.  Only used in full mode.
 * @param indexer Helper to convert from cell id to the indix within
 *                @c cellids.  Only used in full mode.
 */
CreateCaloCells::CaloCells::CaloCells (bool full,
                                       std::span<const uint64_t> cellids,
                                       const k4::recCalo::ICaloIndexer* indexer)
  : m_mode (full ? FULL : SPARSE),
    m_indexer (indexer)
{
  if (full) {
    // Full mode: initialize @c m_cells with all cell ids, and resize
    // @c m_ihits for call cells.
    m_cells.reserve (cellids.size());
    for (uint64_t id : cellids) m_cells.emplace_back (id, 0);
    m_ihits.resize (cellids.size(), INVALID_IHIT);
  }
  else {
    // Sparse mode.  Don't fill in anything yet, but reserve some space
    // for cells.
    m_cells.reserve (2000);
  }
}


/**
 * @brief Return the index within @c m_cells for a given cell id.
 * @param cellid The cell id to find.
 *
 * May return INVALID_ICELL if the cell id doesn't exist in the container.
 * Will throw if we have a filtered full container.
 */
auto CreateCaloCells::CaloCells::indexByID (uint64_t cellid) -> index_t
{
  if (m_mode == FULL) {
    return m_indexer->index (cellid);
  }
  else if (m_mode == SPARSE) {
    return m_indices.try_emplace (cellid, std::make_pair (INVALID_ICELL, INVALID_IHIT)).first->second.first;
  }
  throw std::runtime_error ("indexByID used on filtered cells");
}


/**
 * @brief Return a @c CaloCell wrapper for a given cell id, adding
 *        the cell if needed.
 * @param cellid The cell id to find.
 *
 * May add a new cell.
 * Will throw if we have a filtered full container.
 */
auto CreateCaloCells::CaloCells::cellByID (uint64_t cellid) -> CaloCell
{
  if (m_mode == FULL) {
    index_t icell =  m_indexer->index (cellid);
    if (icell == INVALID_ICELL) [[unlikely]] {
      throw std::out_of_range ("cellByID");
    }
    return CaloCell (m_cells[icell].first, m_cells[icell].second,
                     m_ihits[icell]);
  }

  else if (m_mode == SPARSE) {
    auto& p = m_indices.try_emplace (cellid, std::make_pair (INVALID_ICELL, INVALID_IHIT) ).first->second;
    if (p.first == INVALID_ICELL) {
      // Add a new cell.
      p.first = m_cells.size();
      m_cells.emplace_back (cellid, 0);
      p.second = INVALID_IHIT;
    }
    index_t icell = p.first;
    return CaloCell (m_cells[icell].first, m_cells[icell].second, p.second);
  }

  std::abort();
}


/**
 * @brief Return a @c CaloCell wrapper for a given either the
 *        cell id or the index within @c m_cells.
 * @param icell Index within @c m_cells, or @c INVALID_CELL.
 * @param cellid The cell id to find.
 *
 * Uses @c icell for a @c FULL container, otherwise @c cellid.
 */
inline
auto CreateCaloCells::CaloCells::cellByIndexOrID (index_t icell, uint64_t cellid) -> CaloCell
{
  if (m_mode == FULL)
    return CaloCell (m_cells[icell].first, m_cells[icell].second,
                     m_ihits[icell]);
  return cellByID (cellid);
}


/**
 * @brief Return hit index for a given cell in @c m_cells.
 * @param icell Index of the cell in @c m_cells.
 *
 * Returns the corresponding hit index, or @c INVALID_IHIT.
 */
size_t CreateCaloCells::CaloCells::ihitByIndex (index_t icell) const
{
  if (m_mode == FULL) {
    return m_ihits[icell];
  }

  else if (m_mode == SPARSE) {
    auto it = m_indices.find (m_cells[icell].first);
    if (it == m_indices.end()) return INVALID_IHIT;
    return it->second.second;
  }

  // Filtered
  index_t jcell = m_indexer->index (m_cells[icell].first);
  if (jcell == INVALID_ICELL) return INVALID_IHIT;
  return m_ihits[jcell];
}


/**
 * @brief Sort the sells in order of increasing cell id.
 */
void CreateCaloCells::CaloCells::sort()
{
  if (m_mode != SPARSE) return;

  std::ranges::sort (m_cells);

  // Then reset indices if we've 
  index_t ncells = m_cells.size();
  for (index_t icell = 0; icell < ncells; ++icell) {
    size_t cellID = m_cells[icell].first;
    m_indices[cellID].first = icell;
  }
}


/**
 * @brief Called after filtering.
 * @param orig_size The container size before filtering.
 */
void CreateCaloCells::CaloCells::afterFilter (size_t orig_size)
{
  // If we've filtered a full container, change start to FILTERED.
  if (m_mode == FULL && m_cells.size() != orig_size)
    m_mode = FILTERED;
}

                                       

CreateCaloCells::CreateCaloCells(const std::string& name, ISvcLocator* svcLoc)
    : Gaudi::Algorithm(name, svcLoc), m_geoSvc("GeoSvc", name) {
  declareProperty("hits", m_hits, "Hits from which to create cells (input)");
  declareProperty("cells", m_cells, "The created calorimeter cells (output)");
  declareProperty("links", m_links, "The links between hits and cells (output)");

  declareProperty("calibTool", m_calibTool, "Handle for tool to calibrate Geant4 energy to EM scale tool");
  declareProperty("noiseTool", m_noiseTool, "Handle for the calorimeter cells noise tool");
  declareProperty("geometryTool", m_geoTool, "Handle for the geometry tool");
}

StatusCode CreateCaloCells::initialize() {
  K4RECCALORIMETER_CHECK( Gaudi::Algorithm::initialize() );

  info() << "CreateCaloCells initialized" << endmsg;
  info() << "do calibration : " << m_doCellCalibration << endmsg;
  info() << "add cell noise      : " << m_addCellNoise << endmsg;
  info() << "remove cells below threshold : " << m_filterCellNoise << endmsg;
  info() << "add position information to the cell : " << m_addPosition << endmsg;
  info() << "emulate crosstalk : " << m_addCrosstalk << endmsg;


  // Initialization of tools
  K4RECCALORIMETER_CHECK( m_indexerSvc.retrieve() );

  // Cell crosstalk tool
  if (m_addCrosstalk) {
    K4RECCALORIMETER_CHECK( m_crosstalksTool.retrieve() );
  }
  // Calibrate Geant4 energy to EM scale tool
  if (m_doCellCalibration) {
    K4RECCALORIMETER_CHECK( m_calibTool.retrieve() );
  }
  // Cell noise tool
  if (m_addCellNoise || m_filterCellNoise) {
    K4RECCALORIMETER_CHECK( m_noiseTool.retrieve() );
    // Geometry settings
    K4RECCALORIMETER_CHECK( m_geoTool.retrieve() );
  }

  if (m_addCellNoise) {
    // Construct cell indices.
    {
      // Use the geoTool to get a collection of all CellIDs.
      m_cellIDs = m_geoTool->cellIDs();
      {
        std::string fname = name() + ".cellids";
        FILE* f = fopen(fname.c_str(), "w");
        for (size_t id : m_cellIDs)
          fprintf (f, "%012lx\n", id);
        fclose (f);
      }
    }
  }
  if (m_addPosition) {
    dd4hep::VolumeManager vman_glob = m_geoSvc->getDetector()->volumeManager();
    int id = m_geoTool->id();
    if (id >= 0)
      m_volman = vman_glob.subdetector (id);
    else
      m_volman = vman_glob;
  }

  if (m_cellPos.isEnabled()) {
    K4RECCALORIMETER_CHECK( m_cellPos.retrieve() );
  }

  // Copy over the CellIDEncoding string from the input collection to the output collection
  auto hitsEncoding = m_hitsCellIDEncoding.get_optional();
  K4RECCALORIMETER_CHECK( hitsEncoding.has_value() );
  m_cellsCellIDEncoding.put(hitsEncoding.value());

  m_decoder = dd4hep::DDSegmentation::BitFieldCoder(hitsEncoding.value());
  m_systemIndex = m_decoder.index ("system");
  m_layerIndex = m_decoder.index ("layer");

  // Find the calorimeter type word for each system ID.
  findCaloTypes();

  if (m_links.objKey().empty()) {
    m_links.updateKey (m_cells.objKey() + "SimCaloHitLinks");
  }

  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCells::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
  const edm4hep::SimCalorimeterHitCollection* hits = m_hits.get();
  debug() << "Input Hit collection size: " << hits->size() << endmsg;

  // Find calorimeter type.
  const k4::recCalo::ICaloIndexer* indexer = nullptr;
  int calotype = 0;
  if (!hits->empty()) {
    uint64_t cellid = hits->begin()->getCellID();
    unsigned detid = m_decoder.get (cellid, m_systemIndex);
    if (detid < m_caloTypes.size()) calotype = m_caloTypes[detid];
    if (calotype == 0) {
      error() << "detector id " << detid << " is not a calorimeter" << endmsg;
    }
    indexer = m_indexerSvc->indexer (detid, !m_addCellNoise);
    if (m_addCellNoise) {
      if (!indexer) {
        error() << "Cannot find indexer for detid " << detid << endmsg;
        return StatusCode::FAILURE;
      }
    }
  }


  // 0. Clear all cells
  CaloCells cells (m_addCellNoise, m_cellIDs, indexer);

  // 1. Merge energy deposits into cells
  // If running with noise map already was prepared. Otherwise it is being
  // created below
  for (size_t ihit = 0; const auto& hit : *hits) {
    // If we have an indexer, check that the cell ID is valid; skip if not.
    // This can happen, for example, if we've made the cryostat active.
    uint64_t cellid = hit.getCellID();
    index_t cellIndex = CaloCells::INVALID_ICELL;
    if (indexer) {
      cellIndex = indexer->index (cellid);
    }
    if (cellIndex != CaloCells::INVALID_ICELL || !indexer) {
      // Add the new cell if needed and accumulate energy.
      CaloCell cell = cells.cellByIndexOrID (cellIndex, cellid);
      cell.energy += hit.getEnergy();
      if (cell.ihit == CaloCells::INVALID_IHIT) cell.ihit = ihit;
    }
    ++ihit;
  }
  debug() << "Number of calorimeter cells after merging of hits: " << cells.size() << endmsg;

  // 2. Emulate cross-talk (if asked)
  if (m_addCrosstalk) {
    addCrosstalk (cells);
  }

  // 3. Calibrate simulation energy to EM scale
  if (m_doCellCalibration) {
    m_calibTool->calibrate(cells.m_cells);
  }

  // 4. Add noise to all cells
  if (m_addCellNoise) {
    m_noiseTool->addRandomCellNoise(cells.m_cells);
  }

  if (m_discritN >= 0) {
    auto discrit = [&] (float e) -> float
      {
        if (m_discritN == 0) return e;
        if (e < m_discritMin) return m_discritMin;
        if (e > m_discritMax) return m_discritMax;
        float range = m_discritMax - m_discritMin;
        float x = (e - m_discritMin) / range;
        return static_cast<int>(x*m_discritN + 0.5) / static_cast<float>(m_discritN) * range + m_discritMin;
      };
    auto calib = cells.m_cells;
    for (auto& p : calib) p.second = 1;
    if (m_doCellCalibration) m_calibTool->calibrate(calib);
    for (size_t i = 0; i < cells.m_cells.size(); ++i) {
      float c = calib[i].second;
      auto&p = cells.m_cells[i];
      p.second = c * discrit(p.second/c);
    }
  }

  // 5. Filter cells
  if (m_filterCellNoise) {
    size_t orig_size = cells.size();
    m_noiseTool->filterCellNoise(cells.m_cells);
    cells.afterFilter (orig_size);
  }

  // 6. Copy information to CaloHitCollection
  edm4hep::CalorimeterHitCollection* edmCellsCollection = new edm4hep::CalorimeterHitCollection();

  // Make sure cells are sorted.  (A no-op if we have all cells; they'll
  // already be sorted in that case.)
  cells.sort();

  // This is actually kind of expensive --- hoist it out of the loop.
  bool haveCellPos = m_cellPos.isEnabled();

  for (index_t icell = 0; icell < cells.size(); ++icell) {
    double energy = cells.m_cells[icell].second;
    if (energy == 0 && !m_addCellNoise) continue;
    auto newCell = edmCellsCollection->create();
    newCell.setEnergy(energy);
    uint64_t cellid = cells.m_cells[icell].first;
    newCell.setCellID(cellid);

    static constexpr double inv_mm = 1 / dd4hep::mm;

    size_t ihit = cells.ihitByIndex (icell);

    if (haveCellPos && (m_addPosition || ihit == CaloCells::INVALID_IHIT)) {
      // Recalculate cell position given cell using positioning tool.
      // We have the tool, and either the position isn't available
      // in the input, or we were requested to recalculate it.
      dd4hep::Position pos = m_cellPos->xyzPosition(cellid) * inv_mm;
      newCell.setPosition({static_cast<float>(pos.x()),
          static_cast<float>(pos.y()),
          static_cast<float>(pos.z())});
    }

    else if (m_addPosition) {
      // We were requested to recalculate the positions, but the
      // positioning tool is not available.
      // Calculate a position based on the volume center.
      uint64_t volid = cellid;
      if (const dd4hep::DDSegmentation::Segmentation* seg = m_geoTool->segmentation()) {
        volid = seg->volumeID (volid);
      }
      auto detelement = m_volman.lookupDetElement(volid);
      double inLocal[] = {0, 0, 0};
      const auto outGlobal = detelement.nominal().localToWorld(inLocal);
      edm4hep::Vector3f position = edm4hep::Vector3f(outGlobal.X() * inv_mm,
                                                     outGlobal.Y() * inv_mm,
                                                     outGlobal.Z()* inv_mm);
      newCell.setPosition(position);
    }

    else if (ihit != CaloCells::INVALID_IHIT) {
      // This cell was present in the input.  Take the position from there.
      newCell.setPosition((*hits)[ihit].getPosition());
    }

    // Otherwise, the position will be left set to 0.

    int layer = m_decoder.get(cellid, m_layerIndex);
    newCell.setType (calotype + 10000 * layer);
  }

  // create hits<->cell links
  edm4hep::CaloHitSimCaloHitLinkCollection* edmCellHitLinksCollection = new edm4hep::CaloHitSimCaloHitLinkCollection();
  for (const auto& hit : *hits) {
    index_t icell = cells.indexByID (hit.getCellID());
    if (icell != CaloCells::INVALID_ICELL) {
      // create Sim<->Reco hit associations
      auto link = edmCellHitLinksCollection->create();
      link.setFrom((*edmCellsCollection)[icell]);
      link.setTo(hit);
    }
  }

  // push the CaloHitCollection to event store
  m_cells.put(edmCellsCollection);
  m_links.put(edmCellHitLinksCollection);

  debug() << "Output Cell collection size: " << edmCellsCollection->size() << endmsg;

  return StatusCode::SUCCESS;
}


/// Build m_caloTypes, giving calorimeter type per system ID.
void CreateCaloCells::findCaloTypes()
{
  for (const auto& p : m_geoSvc->getDetector()->detectors()) {
    dd4hep::DetElement det (p.second);
    dd4hep::DetType detType(det.typeFlag());
    int id = det.id();
    if (detType.is(dd4hep::DetType::CALORIMETER)) {
      int calotype = 0;
      int caloid = 0;
      int layout = 0;
      if (detType.is(dd4hep::DetType::ELECTROMAGNETIC)) {
        calotype = 0;
        caloid = 1;
      } else if (detType.is(dd4hep::DetType::HADRONIC)) {
        calotype = 1;
        caloid = 2;
      } else if (detType.is(dd4hep::DetType::MUON)) {
        calotype = 2;
        caloid = 3;
      } else {
        warning() << "Detector type for calorimeter " << id << " is neither ELECTROMAGNETIC, HADRONIC nor MUON" << endmsg;
      }
      if (detType.is(dd4hep::DetType::BARREL)) {
        layout = 1;
      } else if (detType.is(dd4hep::DetType::ENDCAP)) {
        layout = 2;
      } else {
        warning() << "Detector type for calorimeter " << id << " is neither BARREL nor ENDCAP" << endmsg;
      }

      if (static_cast<int>(m_caloTypes.size()) <= id) {
        m_caloTypes.resize (id+1);
      }
      m_caloTypes[id] = calotype + 10 * caloid + 1000 * layout;
    }
  }
}


void CreateCaloCells::addCrosstalk (CaloCells& cells) const
{
  // Derive the cross-talk contributions without affecting yet the nominal energy
  // (one has to emulate crosstalk based on cells free from any cross-talk contributions)
  CaloCells::CellData_t cells_orig = cells.m_cells;
  // yet the nominal energy
  // loop over cells with nominal energies
  for (size_t jcell = 0; jcell < cells_orig.size(); ++jcell) {
    double this_energy = cells_orig[jcell].second;
    if (this_energy == 0) continue;
    uint64_t this_cellId = cells.m_cells[jcell].first;
    auto vec_neighbours = m_crosstalksTool->getNeighbours(this_cellId); // a span of neighbour IDs
    auto vec_crosstalks = m_crosstalksTool->getCrosstalks(this_cellId); // a span of crosstalk coefficients
    // loop over crosstalk neighbours of the cell under study
    for (unsigned int i_cell = 0; i_cell < vec_neighbours.size(); i_cell++) {
      // signal transfer = energy deposit brought by EM shower hits * crosstalk coefficient
      double signal_transfer = this_energy * vec_crosstalks[i_cell];
      // for the cell under study, record the signal transfer that will be subtracted from its final cell energy
      cells.m_cells[jcell].second -= signal_transfer;
      // for the crosstalk neighbour, record the signal transfer that will be added to its final cell energy
      CaloCell other_cell = cells.cellByID (vec_neighbours[i_cell]);
      other_cell.energy += signal_transfer;
    }
  }
}
