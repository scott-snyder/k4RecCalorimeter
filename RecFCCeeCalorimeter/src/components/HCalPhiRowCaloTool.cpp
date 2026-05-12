/**
 * @file RecFCCeeCalorimeter/src/components/HCalPhiRowCaloTool.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date May, 2026
 * @brief Calorimeter tool for Allegro HCal (with row indexing).
 */

#include "HCalPhiRowCaloTool.h"
#include "DD4hep/Detector.h"
#include "RecCaloCommon/IDMapIndexer.h"
#include "detectorSegmentations/FCCSWHCalPhiRow_k4geo.h"
#include "k4FWCore/GaudiChecks.h"


DECLARE_COMPONENT(HCalPhiRowCaloTool)

/** Gaudi initialize method.
 */
StatusCode HCalPhiRowCaloTool::initialize() {
  // Look up the detector ID.
  std::string detidname = "DetID_";
  if (readoutName() == "HCalBarrelReadoutPhiRow")
    detidname += "HCAL_Barrel";
  else
    detidname += "HCAL_Endcap";

  K4_GAUDI_CHECK(geoSvc().retrieve());
  m_id = geoSvc()->getDetector()->constant<int>(detidname);
  if (m_id < 0) {
    error() << "Can't find detector ID " << detidname << endmsg;
    return StatusCode::FAILURE;
  }

  K4_GAUDI_CHECK(CalorimeterToolBase::initialize());
  return StatusCode::SUCCESS;
}

/** Return the subdetector ID.
 */
int HCalPhiRowCaloTool::id() const {  return m_id; }

/** Fill vector with all existing cells for this geometry.
 */
StatusCode HCalPhiRowCaloTool::collectCells(std::vector<uint64_t>& cells) const {
  bool is_barrel = (readoutName() == "HCalBarrelReadoutPhiRow");
  cells.reserve(is_barrel ? 1031680 : 1164800);
  const auto* seg =
      dynamic_cast<const dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*> (readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWHCalPhiRow_k4geo "
               "segmentation."
            << endmsg;
    return StatusCode::FAILURE;
  }

  const dd4hep::DDSegmentation::BitFieldCoder& decoder = *readout().idSpec().decoder();

  dd4hep::DDSegmentation::CellID cID = 0;
  decoder.set(cID, "system", m_id);

  size_t layer_id = decoder.index(seg->fieldNameLayer());
  size_t row_id = decoder.index(seg->fieldNameRow());
  size_t phi_id = decoder.index(seg->fieldNamePhi());
  size_t pseudolayer_id = is_barrel ? 0 : decoder.index(seg->fieldNamePseudoLayer());

  int numLayers = 0;
  for (int n : seg->numLayers())
    numLayers += n;
  for (int layer = 0; layer < numLayers; ++layer) {
    decoder.set(cID, layer_id, layer);
    for (int row : seg->cellIndexes(layer)) {
      decoder.set(cID, row_id, row);
      if (!is_barrel)
        decoder.set(cID, pseudolayer_id, seg->definePseudoLayer(cID));
      for (int phi = 0; phi < seg->phiBins(); ++phi) {
        decoder.set(cID, phi_id, phi);
        cells.push_back(cID);
      }
    }
  }

  return StatusCode::SUCCESS;
}


/** Return a new indexer object for this subdetector.
 */
std::unique_ptr<k4::recCalo::ICaloIndexer> HCalPhiRowCaloTool::indexer() const {
  const auto* seg =
      dynamic_cast<const dd4hep::DDSegmentation::FCCSWHCalPhiRow_k4geo*> (readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWHCalPhiTow_k4geo "
               "segmentation."
            << endmsg;
    return nullptr;
  }

  bool is_barrel = (readoutName() == "HCalBarrelReadoutPhiRow");

  using Indexer_t = k4::recCalo::IDMapIndexer<3>;
  dd4hep::IDDescriptor idSpec = readout().idSpec();
  // Sorry for the ugly formatting, but clang-format insists.
  std::vector<Indexer_t::FieldDesc_t> fields{Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameLayer())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNameRow())),
                                             Indexer_t::IDMap_t::makeDesc(*idSpec.field(seg->fieldNamePhi())) };

  const dd4hep::BitFieldElement& sysField = *idSpec.field("system");
  if (sysField.offset() != 0) {
    throw std::runtime_error("Bad system field offset; must be zero");
  }

  std::vector<Indexer_t::FieldDesc_t> ignoredFields;
  if (!is_barrel) {
    ignoredFields.push_back (Indexer_t::IDMap_t::makeDesc (*idSpec.field(seg->fieldNamePseudoLayer())));
  }

  // Again, clang-format complains if this is written legibly...
  return std::make_unique<Indexer_t> (this->id(), sysField.width(), fields, cellIDs(), is_barrel ? 4200000 : 1000000,
                                      ignoredFields);
}
