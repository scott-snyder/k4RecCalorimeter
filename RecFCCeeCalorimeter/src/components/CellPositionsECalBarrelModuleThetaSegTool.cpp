#include "CellPositionsECalBarrelModuleThetaSegTool.h"

// EDM
#include "edm4hep/CalorimeterHitCollection.h"
#include "RecCaloCommon/GaudiChecks.h"

#include <cmath>

DECLARE_COMPONENT(CellPositionsECalBarrelModuleThetaSegTool)

StatusCode CellPositionsECalBarrelModuleThetaSegTool::initialize() {
  K4_GAUDI_CHECK( AlgTool::initialize() );
  K4_GAUDI_CHECK( m_geoSvc.retrieve() );
  K4_GAUDI_CHECK( m_indexerSvc.retrieve() );
  K4_GAUDI_CHECK( m_constantsSvc.retrieve() );

  // get segmentation
  dd4hep::Segmentation segmentation = m_geoSvc->getDetector()->readout(m_readoutName).segmentation();
  m_segmentation = dynamic_cast<dd4hep::DDSegmentation::FCCSWGridModuleThetaMerged_k4geo*>(
      segmentation.segmentation());
  if (m_segmentation == nullptr) {
    error() << "There is no module-theta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  debug() << "Found merged module-theta segmentation" << endmsg;
  for (int iLayer = 0; iLayer < m_segmentation->nLayers(); iLayer++) {
    info() << "Layer : " << iLayer << " theta merge : " << m_segmentation->mergedThetaCells(iLayer)
           << " module merge : " << m_segmentation->mergedModules(iLayer) << endmsg;
    if (m_segmentation->mergedThetaCells(iLayer) < 1) {
      error() << "Number of cells merged along theta should be >= 1!!!!" << endmsg;
    }
    if (m_segmentation->mergedModules(iLayer) < 1) {
      error() << "Number of modules merged should be >= 1!!!!" << endmsg;
    }
  }

  int detID = segmentation.detector()->id;
  m_indexer = m_indexerSvc->indexer (detID);

  std::string dataKey = m_readoutName.value() + "-cellPositions";
  const PositionData* data = m_constantsSvc->getObj<PositionData> (dataKey);
  if (!data) {
    dd4hep::VolumeManager volman_glob = m_geoSvc->getDetector()->volumeManager();
    dd4hep::VolumeManager volman = volman_glob.subdetector (detID);

    std::span<const uint64_t> ids = m_indexer->cellIDs();

    PositionData positions;
    positions.resize (ids.size());
    for (uint64_t id : ids) {
      unsigned index = m_indexer->index (id);
      dd4hep::DDSegmentation::CellID volumeId = m_segmentation->volumeID(id);
      dd4hep::VolumeManagerContext* vc = volman.lookupContext(volumeId);
      dd4hep::DDSegmentation::Vector3D inSeg = m_segmentation->position(id);
      positions.at(index) = vc->localToWorld(dd4hep::Position(inSeg));
    }
    K4_GAUDI_CHECK( m_constantsSvc->putObj (dataKey, std::move (positions)) );
    data = m_constantsSvc->getObj<PositionData> (dataKey);
    K4_GAUDI_CHECK( data != nullptr );
  }
  m_positions = *data;
  return StatusCode::SUCCESS;
}

void CellPositionsECalBarrelModuleThetaSegTool::getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                                                             edm4hep::CalorimeterHitCollection& outputColl) const {

  debug() << "Input collection size : " << aCells.size() << endmsg;

  // Loop through input cell collection, call xyzPosition method for each cell
  // and assign position to cloned hit to be saved in outputColl
  for (const auto& cell : aCells) {
    auto outSeg = CellPositionsECalBarrelModuleThetaSegTool::xyzPosition(cell.getCellID());
    auto edmPos = edm4hep::Vector3f();
    edmPos.x = outSeg.x() / dd4hep::mm;
    edmPos.y = outSeg.y() / dd4hep::mm;
    edmPos.z = outSeg.z() / dd4hep::mm;

    auto positionedHit = cell.clone();
    positionedHit.setPosition(edmPos);
    outputColl.push_back(positionedHit);

    // Debug information about cell position
    debug() << "Cell energy (GeV) : " << positionedHit.getEnergy() << "\tcellID " << positionedHit.getCellID()
            << endmsg;
    debug() << "Position of cell (mm) : \t" << outSeg.x() / dd4hep::mm << "\t" << outSeg.y() / dd4hep::mm << "\t"
            << outSeg.z() / dd4hep::mm << "\n"
            << endmsg;
  }
  debug() << "Output positions collection size: " << outputColl.size() << endmsg;
}

dd4hep::Position CellPositionsECalBarrelModuleThetaSegTool::xyzPosition(const uint64_t& aCellId) const {

  // find position of volume corresponding to first of group of merged cells
  unsigned index = m_indexer->index (aCellId);
  if (index >= m_positions.size()) throw std::out_of_range ("CellPositionsECalBarrelModuleThetaSegTool::xyzPosition");
  return m_positions[index];
}

int CellPositionsECalBarrelModuleThetaSegTool::layerId(const uint64_t& aCellId) const {
  return m_segmentation->layer(aCellId);
}
