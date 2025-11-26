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
  K4_GAUDI_CHECK( Gaudi::Algorithm::initialize() );

  info() << "CreateCaloCells initialized" << endmsg;
  info() << "do calibration : " << m_doCellCalibration << endmsg;
  info() << "add cell noise      : " << m_addCellNoise << endmsg;
  info() << "remove cells below threshold : " << m_filterCellNoise << endmsg;
  info() << "add position information to the cell : " << m_addPosition << endmsg;
  info() << "emulate crosstalk : " << m_addCrosstalk << endmsg;

  // Initialization of tools
  // Cell crosstalk tool
  if (m_addCrosstalk) {
    K4_GAUDI_CHECK( m_crosstalksTool.retrieve() );
  }
  // Calibrate Geant4 energy to EM scale tool
  if (m_doCellCalibration) {
    K4_GAUDI_CHECK( m_calibTool.retrieve() );
  }
  // Cell noise tool
  if (m_addCellNoise || m_filterCellNoise) {
    K4_GAUDI_CHECK( m_noiseTool.retrieve() );
    // Geometry settings
    K4_GAUDI_CHECK( m_geoTool.retrieve() );

    // Construct cell indices.
    {
      // Use the geoTool to get a collection of all CellIDs.
      std::vector<size_t> cellIDs = m_geoTool->cellIDs();

      // And index them.
      for (size_t i = 0; size_t id : cellIDs)
        m_cellsIndexMap[id] = i++;
    }
  }
  if (m_addPosition) {
    dd4hep::VolumeManager vman_glob = m_geoSvc->getDetector()->volumeManager();
    m_volman = vman_glob.subdetector (m_geoTool->id());
  }

  if (m_cellPos.isEnabled()) {
    K4_GAUDI_CHECK( m_cellPos.retrieve() );
  }

  // Copy over the CellIDEncoding string from the input collection to the output collection
  auto hitsEncoding = m_hitsCellIDEncoding.get_optional();
  K4_GAUDI_CHECK( hitsEncoding.has_value() );
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
    cells.m_cells.resize (m_cellsIndexMap.size());
    for (const auto& p : m_cellsIndexMap) {
      cells.m_cells.at(p.second).first = p.first;
      cellsIndex.index(p.first) = p.second;
    }
  }

  // 1. Merge energy deposits into cells
  // If running with noise map already was prepared. Otherwise it is being
  // created below
  for (size_t ihit = 0; const auto& hit : *hits) {
    verbose() << "CellID : " << hit.getCellID() << endmsg;
    size_t& icell = cellsIndex.index (hit.getCellID(), ihit);
    if (icell == INVALID) {
      icell = cells.add (hit.getCellID(), hit.getEnergy());
    }
    else {
      cells.energy(icell) += hit.getEnergy();
    }
    ++ihit;
  }
  debug() << "Number of calorimeter cells after merging of hits: " << cells.size() << endmsg;

  // 2. Emulate cross-talk (if asked)
  if (m_addCrosstalk) {
    // Derive the cross-talk contributions without affecting yet the nominal energy
    // (one has to emulate crosstalk based on cells free from any cross-talk contributions)
    CellsInfo cells_orig = cells;
                                 // yet the nominal energy
    // loop over cells with nominal energies
    for (size_t jcell = 0; jcell < cells_orig.size(); ++jcell) {
      uint64_t this_cellId = cells_orig.cellID(jcell);
      double this_energy = cells_orig.energy(jcell);
      auto vec_neighbours = m_crosstalksTool->getNeighbours(this_cellId); // a vector of neighbour IDs
      auto vec_crosstalks = m_crosstalksTool->getCrosstalks(this_cellId); // a vector of crosstalk coefficients
      // loop over crosstalk neighbours of the cell under study
      for (unsigned int i_cell = 0; i_cell < vec_neighbours.size(); i_cell++) {
        // signal transfer = energy deposit brought by EM shower hits * crosstalk coefficient
        double signal_transfer = this_energy * vec_crosstalks[i_cell];
        // for the cell under study, record the signal transfer that will be subtracted from its final cell energy
        cells.energy(jcell) -= signal_transfer;
        // for the crosstalk neighbour, record the signal transfer that will be added to its final cell energy
        size_t& kcell = cellsIndex.index (vec_neighbours[i_cell]);
        if (kcell == INVALID) {
          kcell = cells.add (vec_neighbours[i_cell], signal_transfer);
        }
        else {
          cells.energy(kcell) += signal_transfer;
        }
      }
    }
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
    m_noiseTool->filterCellNoise(cells.m_cells);
  }

  // 6. Copy information to CaloHitCollection
  edm4hep::CalorimeterHitCollection* edmCellsCollection = new edm4hep::CalorimeterHitCollection();

  cellsIndex.sort (cells);

  for (size_t icell = 0; icell < cells.size(); ++icell) {
    double energy = cells.energy(icell);
    if (energy == 0 && !m_addCellNoise) continue;
    auto newCell = edmCellsCollection->create();
    newCell.setEnergy(energy);
    uint64_t cellid = cells.cellID(icell);
    newCell.setCellID(cellid);

    static constexpr double inv_mm = 1 / dd4hep::mm;

    size_t ihit = cellsIndex.ihit(cellid);

    if (m_cellPos.isEnabled() && (m_addPosition || ihit == INVALID)) {
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

    else if (ihit != INVALID) {
      // This cell was present in the input.  Take the position from there.
      newCell.setPosition((*hits)[ihit].getPosition());
    }

    // Otherwise, the position will be left set to 0.
  }

  // create hits<->cell links
  edm4hep::CaloHitSimCaloHitLinkCollection* edmCellHitLinksCollection = new edm4hep::CaloHitSimCaloHitLinkCollection();
  for (const auto& hit : *hits) {
    auto hitID = hit.getCellID();
    size_t icell = cellsIndex.index(hitID);
    // create Sim<->Reco hit associations
    auto link = edmCellHitLinksCollection->create();
    link.setFrom((*edmCellsCollection)[icell]);
    link.setTo(hit);
  }

  // push the CaloHitCollection to event store
  m_cells.put(edmCellsCollection);
  m_links.put(edmCellHitLinksCollection);

  debug() << "Output Cell collection size: " << edmCellsCollection->size() << endmsg;

  return StatusCode::SUCCESS;
}
