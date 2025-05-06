#pragma GCC optimize "-O0"
#include "CreateCaloCells.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/k4_check.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"
#include "TGeoManager.h"

// edm4hep
#include "edm4hep/CalorimeterHit.h"

DECLARE_COMPONENT(CreateCaloCells)

CreateCaloCells::CreateCaloCells(const std::string& name, ISvcLocator* svcLoc) :
Gaudi::Algorithm(name, svcLoc), m_geoSvc("GeoSvc", name) {
  declareProperty("hits", m_hits, "Hits from which to create cells (input)");
  declareProperty("cells", m_cells, "The created calorimeter cells (output)");
  declareProperty("links", m_links, "The links between hits and cells (output)");

  declareProperty("calibTool", m_calibTool, "Handle for tool to calibrate Geant4 energy to EM scale tool");
  declareProperty("noiseTool", m_noiseTool, "Handle for the calorimeter cells noise tool");
  declareProperty("geometryTool", m_geoTool, "Handle for the geometry tool");
}

StatusCode CreateCaloCells::initialize() {
  K4_CHECK( Gaudi::Algorithm::initialize() );

  info() << "CreateCaloCells initialized" << endmsg;
  info() << "do calibration : " << m_doCellCalibration << endmsg;
  info() << "add cell noise      : " << m_addCellNoise << endmsg;
  info() << "remove cells below threshold : " << m_filterCellNoise << endmsg;
  info() << "add position information to the cell : " << m_addPosition << endmsg;
  info() << "emulate crosstalk : " << m_addCrosstalk << endmsg;

  // Initialization of tools
  // Cell crosstalk tool
  if (m_addCrosstalk) {
    K4_CHECK( m_crosstalksTool.retrieve() );
  }
  // Calibrate Geant4 energy to EM scale tool
  if (m_doCellCalibration) {
    K4_CHECK( m_calibTool.retrieve() );
  }
  // Cell noise tool
  if (m_addCellNoise || m_filterCellNoise) {
    K4_CHECK( m_noiseTool.retrieve() );
    // Geometry settings
    K4_CHECK( m_geoTool.retrieve() );
    // Prepare map of all existing cells in calorimeter to add noise to all
    K4_CHECK( m_geoTool->prepareEmptyCells(m_cellsMap) );
    verbose() << "Initialised empty cell map with size " << m_cellsMap.size() << endmsg;
    // noise filtering erases cells from the cell map after each event, so we need
    // to backup the empty cell map for later reuse
    if (m_addCellNoise && m_filterCellNoise) {
      m_emptyCellsMap = m_cellsMap;
    }

    // Construct cell indices.
    {
      // Use the geoTool to get a collection of all CellIDs.
      std::unordered_map<uint64_t, double> cellsMap;
      K4_CHECK( m_geoTool->prepareEmptyCells(cellsMap) );

      // Now make a sorted list of them.
      auto r = cellsMap | std::views::transform ([](auto x){return x.first;});
      std::vector<size_t> cellIDs (std::ranges::begin(r), std::ranges::end(r));
      std::ranges::sort (cellIDs);

      // And index them.
      for (size_t i = 0; size_t id : cellIDs)
        m_cellsIndexMap[id] = i++;
    }
  }
  if (m_addPosition){
    m_volman = m_geoSvc->getDetector()->volumeManager();
  }

  if (m_cellPos.isEnabled()) {
    K4_CHECK( m_cellPos.retrieve() );
  }

  // Copy over the CellIDEncoding string from the input collection to the output collection
  auto hitsEncoding = m_hitsCellIDEncoding.get_optional();
  K4_CHECK( hitsEncoding.has_value() );
  m_cellsCellIDEncoding.put(hitsEncoding.value());

  if (m_links.objKey().empty()) {
    m_links.updateKey (m_cells.objKey() + "SimCaloHitLinks");
  }

  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCells::execute(const EventContext&) const {
  // Get the input collection with Geant4 hits
  const edm4hep::SimCalorimeterHitCollection* hits = m_hits.get();
  debug() << "Input Hit collection size: " << hits->size() << endmsg;

  CellsInfo cells (m_cellsIndexMap.empty() ? 2048 :  m_cellsIndexMap.size());
  CellsIndex cellsIndex (m_cellsIndexMap);

  // 0. Clear all cells
  if (m_addCellNoise) {
    // if cells are not filtered, the map has same size in each event, equal to the total number
    // of cells in the calorimeter, so we can just reset the values to 0
    // if cells are filtered, during each event they are removed from the cellsMap, so one has to
    // restore the initial map of all empty cells
    if (!m_filterCellNoise)
      std::for_each(m_cellsMap.begin(), m_cellsMap.end(), [](std::pair<const uint64_t, double>& p) { p.second = 0; });
    else
      m_cellsMap = m_emptyCellsMap;
  } else {
    m_cellsMap.clear();
  }


  // 1. Merge energy deposits into cells
  // If running with noise map already was prepared. Otherwise it is being
  // created below
  for (size_t ihit = 0; const auto& hit : *hits) {
    verbose() << "CellID : " << hit.getCellID() << endmsg;
    m_cellsMap[hit.getCellID()] += hit.getEnergy();
    size_t& icell = cellsIndex.index (hit.getCellID());
    if (icell == INVALID) {
      icell = cells.add (hit.getCellID(), hit.getEnergy(), ihit);
    }
    else {
      cells.energy(icell) += hit.getEnergy();
    }
    ++ihit;
  }
  debug() << "Number of calorimeter cells after merging of hits: " << m_cellsMap.size() << endmsg;

  // 2. Emulate cross-talk (if asked)
  if(m_addCrosstalk) {
    // Derive the cross-talk contributions without affecting yet the nominal energy
    // (one has to emulate crosstalk based on cells free from any cross-talk contributions)
    m_CrosstalkCellsMap.clear(); // this is a temporary map to hold energy exchange due to cross-talk, without affecting yet the nominal energy
    // loop over cells with nominal energies
    for (const auto& this_cell : m_cellsMap) {
      uint64_t this_cellId = this_cell.first;
      auto vec_neighbours = m_crosstalksTool->getNeighbours(this_cellId); // a vector of neighbour IDs
      auto vec_crosstalks = m_crosstalksTool->getCrosstalks(this_cellId); // a vector of crosstalk coefficients
      // loop over crosstalk neighbours of the cell under study
      for (unsigned int i_cell=0; i_cell<vec_neighbours.size(); i_cell++) {
        // signal transfer = energy deposit brought by EM shower hits * crosstalk coefficient
        double signal_transfer = this_cell.second * vec_crosstalks[i_cell];
        // for the cell under study, record the signal transfer that will be subtracted from its final cell energy
        m_CrosstalkCellsMap[this_cellId] -= signal_transfer;
        // for the crosstalk neighbour, record the signal transfer that will be added to its final cell energy
        m_CrosstalkCellsMap[vec_neighbours[i_cell]] += signal_transfer;
      }
    }

    // apply the cross-talk contributions on the nominal cell-energy map
    for (const auto& this_cell : m_CrosstalkCellsMap) {
      m_cellsMap[this_cell.first] += this_cell.second;
      size_t& icell = cellsIndex.index (this_cell.first);
      if (icell == INVALID) {
        icell = cells.add (this_cell.first, this_cell.second);
      }
      else {
        cells.energy(icell) += this_cell.second;
      }
    }
    
  }

  // 3. Calibrate simulation energy to EM scale
  if (m_doCellCalibration) {
    m_calibTool->calibrate(m_cellsMap);
    m_calibTool->calibrate(cells.m_cells);
  }

  // 4. Add noise to all cells
  if (m_addCellNoise) {
    m_noiseTool->addRandomCellNoise(m_cellsMap);
    m_noiseTool->addRandomCellNoise(cells.m_cells);
  }

  // 5. Filter cells
  if (m_filterCellNoise) {
    m_noiseTool->filterCellNoise(m_cellsMap);
    m_noiseTool->filterCellNoise(cells.m_cells);
  }

  // 6. Copy information to CaloHitCollection
  edm4hep::CalorimeterHitCollection* edmCellsCollection = new edm4hep::CalorimeterHitCollection();

  if (cells.size() != m_cellsMap.size()) std::abort();

  for (size_t icell = 0; icell < cells.size(); ++icell) {
    double energy = cells.energy(icell);
    if (energy == 0 && !m_addCellNoise) continue;
    auto newCell = edmCellsCollection->create();
    newCell.setEnergy(energy);
    uint64_t cellid = cells.cellID(icell);
    newCell.setCellID(cellid);

    {
      auto it = m_cellsMap.find (cellid);
      if (it == m_cellsMap.end()) std::abort();
      if (it->second != energy) std::abort();
    }

    static constexpr double inv_mm = 1 / dd4hep::mm;

    bool hasInputHit = cells.hasInputHit(icell);

    if (m_cellPos.isEnabled() && (m_addPosition || !hasInputHit)) {
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
      auto detelement = m_volman.lookupDetElement(cellid);
      const auto& transformMatrix = detelement.nominal().worldTransformation();
      double outGlobal[3];
      double inLocal[] = {0, 0, 0};
      transformMatrix.LocalToMaster(inLocal, outGlobal);
      edm4hep::Vector3f position = edm4hep::Vector3f(outGlobal[0] * inv_mm,
                                                     outGlobal[1] * inv_mm,
                                                     outGlobal[2] * inv_mm);
      newCell.setPosition(position);
    }

    else if (hasInputHit) {
      // This cell was present in the input.  Take the position from there.
      newCell.setPosition((*hits)[cells.ihit(icell)].getPosition());
    }

    // Otherwise, the position will be left set to 0.
  }
#if 0
  for (const auto& cell : m_cellsMap) {
    if (m_addCellNoise || cell.second != 0) {
      auto newCell = edmCellsCollection->create();
      newCell.setEnergy(cell.second);
      uint64_t cellid = cell.first;
      newCell.setCellID(cellid);
      // recalc with tool: have tool && (addPosition || (!addPosition && cell not in input)
      // recalc w/o tool: no tool && addPosition
      // take from input: no tool && !addPosition && cell in input
      // 0: no tool && !addPosition && cell not in input

      // If addPosition: recalcuate, either using or not using the tool
      // If not addPosition and cell in input: take pos from input
      // If not addPosition and cell not in input:
      //   if tool available recalc from tool, else 0
      if (m_addPosition){
        if (m_cellPos.isEnabled()) {
          static constexpr double inv_mm = 1 / dd4hep::mm;
          dd4hep::Position pos = m_cellPos->xyzPosition(cellid) * inv_mm;
          newCell.setPosition({static_cast<float>(pos.x()),
                               static_cast<float>(pos.y()),
                               static_cast<float>(pos.z())});
        }
        else {
          auto detelement = m_volman.lookupDetElement(cellid);
          const auto& transformMatrix = detelement.nominal().worldTransformation();
          double outGlobal[3];
          double inLocal[] = {0, 0, 0};
          transformMatrix.LocalToMaster(inLocal, outGlobal);
          edm4hep::Vector3f position = edm4hep::Vector3f(outGlobal[0] / dd4hep::mm, outGlobal[1] / dd4hep::mm, outGlobal[2] / dd4hep::mm);
          newCell.setPosition(position);
        }
      }
    }
  }
#endif

  // XXX Avoid N^2!
  // create hits<->cell links
  edm4hep::CaloHitSimCaloHitLinkCollection* edmCellHitLinksCollection = new edm4hep::CaloHitSimCaloHitLinkCollection();
  for (const auto& cell : *edmCellsCollection) {
    auto cellID = cell.getCellID();
    for (const auto& hit : *hits) {
      auto hitID = hit.getCellID();
      if (hitID == cellID) {
        // create Sim<->Reco hit associations
        auto link = edmCellHitLinksCollection->create();
        link.setFrom(cell);
        link.setTo(hit);
      }
    }
  }

  // push the CaloHitCollection to event store
  m_cells.put(edmCellsCollection);
  m_links.put(edmCellHitLinksCollection);

  debug() << "Output Cell collection size: " << edmCellsCollection->size() << endmsg;

  return StatusCode::SUCCESS;
}

StatusCode CreateCaloCells::finalize() { return Gaudi::Algorithm::finalize(); }
