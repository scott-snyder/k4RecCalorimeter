#include "HCalPhiThetaCaloTool.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "DD4hep/detail/DetectorInterna.h"
#include <algorithm>
#include <functional>


DECLARE_COMPONENT(HCalPhiThetaCaloTool)


/** Fill vector with all existing cells for this geometry.
 */
StatusCode HCalPhiThetaCaloTool::collectCells(std::vector<uint64_t>& cells) const
{
  const auto* seg =
    dynamic_cast<const dd4hep::DDSegmentation::FCCSWHCalPhiTheta_k4geo*> (readout().segmentation().segmentation());
  if (!seg) {
    error() << "Unable to cast segmentation pointer!!!! Tool only applicable to FCCSWHCalPhiTheta_k4geo "
               "segmentation."
            << endmsg;
    return StatusCode::FAILURE;
  }

  const dd4hep::DDSegmentation::BitFieldCoder& decoder =
    *readout().idSpec().decoder();
  int id = readout().segmentation().detector()->id;

  dd4hep::DDSegmentation::CellID cID = 0;
  decoder.set(cID, "system", id);

  size_t layer_id = decoder.index (seg->fieldNameLayer());
  size_t theta_id = decoder.index (seg->fieldNameTheta());
  size_t phi_id = decoder.index (seg->fieldNamePhi());

  int numLayers = 0;
  for (int n : seg->numLayers()) numLayers += n;
  for (int layer = 0; layer < numLayers; ++layer) {
    const std::vector<int>& thetaBins = seg->thetaBins (layer);
    decoder.set(cID, layer_id, layer);
    for (int theta : thetaBins) {
      decoder.set(cID, theta_id, theta);
      for (int phi = 0; phi < seg->phiBins(); ++phi) {
        decoder.set(cID, phi_id, phi);
        cells.push_back(cID);
      }
    }
  }

  return StatusCode::SUCCESS;
}
