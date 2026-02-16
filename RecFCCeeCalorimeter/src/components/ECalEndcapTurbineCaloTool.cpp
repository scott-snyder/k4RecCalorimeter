/**
 * @file RecFCCeeCalorimeter/src/components/ECalEndcapTurbineCaloTool.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Calorimeter tool for Allegro ECal endcap.
 */


#include "ECalEndcapTurbineCaloTool.h"
#include "RecCaloCommon/IDMapIndexer.h"
#include "detectorSegmentations/FCCSWEndcapTurbine_k4geo.h"
#include "DD4hep/detail/DetectorInterna.h"


DECLARE_COMPONENT(ECalEndcapTurbineCaloTool)


ECalEndcapTurbineCaloTool::WheelParams::WheelParams
  (const dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo& seg, int iwheel)
{
  nmodules = seg.nModules (iwheel);
  nrho = seg.numCellsRho (iwheel);
  nz = seg.numCellsZ (iwheel);
}


/** Fill vector with all existing cells for this geometry.
 */
StatusCode ECalEndcapTurbineCaloTool::collectCells(std::vector<uint64_t>& cells) const
{
  cells.reserve (1203200);

  const auto* seg =
    dynamic_cast<const dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*> (readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWEndcapTurbine_k4geo "
               "segmentation."
            << endmsg;
    return StatusCode::FAILURE;
  }

  const dd4hep::DDSegmentation::BitFieldCoder& decoder =
    *readout().idSpec().decoder();
  int id = readout().segmentation().detector()->id;

  dd4hep::DDSegmentation::CellID cID = 0;
  decoder.set(cID, "system", id);

  size_t side_id = decoder.index (seg->fieldNameSide());
  size_t wheel_id = decoder.index (seg->fieldNameWheel());
  size_t module_id = decoder.index (seg->fieldNameModule());
  size_t rho_id = decoder.index (seg->fieldNameRho());
  size_t z_id = decoder.index (seg->fieldNameZ());
  size_t layer_id = decoder.index (seg->fieldNameLayer());

  int nwheels = 3;
  std::vector<WheelParams> p;
  for (int iwheel = 0; iwheel < nwheels; ++iwheel) {
    p.emplace_back (*seg, iwheel);
  }

  int sides[2] = {-1, 1};
  for (int iside : sides) {
    decoder.set (cID, side_id, iside);
    for (int iwheel = 0; iwheel < nwheels; ++iwheel) {
      const WheelParams& w = p[iwheel];
      decoder.set (cID, wheel_id, iwheel);
      for (int imodule = 0; imodule < w.nmodules; ++imodule) {
        decoder.set (cID, module_id, imodule);
        for (int irho = 0; irho < w.nrho; ++irho) {
          decoder.set (cID, rho_id, irho);
          for (int iz = 0; iz < w.nz; ++iz) {
            decoder.set (cID, z_id, iz);
            decoder.set (cID, layer_id, seg->expLayer (iwheel, irho, iz));
            cells.push_back (cID);
          }
        }
      }
    }
  }

  return StatusCode::SUCCESS;
}


/** Return a new indexer object for this subdetector.
 */
std::unique_ptr<k4::recCalo::ICaloIndexer>
ECalEndcapTurbineCaloTool::indexer() const
{
  const auto* seg =
    dynamic_cast<const dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo*> (readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWECalEndcapTurbine_k4geo "
               "segmentation."
            << endmsg;
    return nullptr;
  }

  using Indexer_t = k4::recCalo::IDMapIndexer<4>;
  dd4hep::IDDescriptor idSpec = readout().idSpec();

  // Combine side+wheel together into a single field.
  const dd4hep::BitFieldElement* bfe_side = idSpec.field ("side");
  const dd4hep::BitFieldElement* bfe_wheel = idSpec.field ("wheel");
  if (bfe_side->offset() + bfe_side->width() != bfe_wheel->offset()) {
    error() << "side and wheel fields not adjacent" << endmsg;
    return nullptr;
  }
  Indexer_t::FieldDesc_t field_sidewheel
    (bfe_side->offset(), bfe_side->width() + bfe_wheel->width());

  std::vector<Indexer_t::FieldDesc_t> fields
    { field_sidewheel,
      Indexer_t::IDMap_t::makeDesc (*idSpec.field(seg->fieldNameRho())),
      Indexer_t::IDMap_t::makeDesc (*idSpec.field(seg->fieldNameZ())),
      Indexer_t::IDMap_t::makeDesc (*idSpec.field(seg->fieldNameModule()))
    };
  std::vector<Indexer_t::FieldDesc_t> ignoredFields
    { Indexer_t::IDMap_t::makeDesc (*idSpec.field(seg->fieldNameLayer()))
    };

  const dd4hep::BitFieldElement& sysField = *idSpec.field("system");
  if (sysField.offset() != 0) {
    throw std::runtime_error ("Bad system field offset; must be zero");
  }

  return std::make_unique<Indexer_t> (this->id(),
                                      sysField.width(),
                                      fields, cellIDs(),
                                      6.2 * 1024 * 1024,
                                      ignoredFields);
}
