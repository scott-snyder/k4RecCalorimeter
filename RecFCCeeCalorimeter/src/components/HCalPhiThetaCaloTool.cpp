#include "HCalPhiThetaCaloTool.h"
#include "detectorSegmentations/FCCSWHCalPhiTheta_k4geo.h"
#include "DD4hep/detail/DetectorInterna.h"
#include <algorithm>
#include <functional>

// k4FWCore
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/k4_check.h"

DECLARE_COMPONENT(HCalPhiThetaCaloTool)

StatusCode HCalPhiThetaCaloTool::initialize() {
  K4_CHECK( AlgTool::initialize() );
  K4_CHECK( CalorimeterToolBase::getReadout (m_readoutName) );
  return StatusCode::SUCCESS;
}


StatusCode HCalPhiThetaCaloTool::collectCells(std::function<void(uint64_t)> cellFunc) const
{
  std::cout << "hhh0 " << name() << "\n";
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

  int numLayers = std::ranges::fold_left (seg->numLayers(), 0,
                                          std::plus<int>());
  std::cout << "hhh1 " << numLayers << " " << seg->numLayers().size() << "\n";
  for (int layer = 0; layer < numLayers; ++layer) {
    const std::vector<int>& thetaBins = seg->thetaBins (layer);
    std::cout << "hhh " << id << " " << layer << " "
              << thetaBins.size() << " [";
    for (int th : thetaBins) std::cout << th << ", ";
    std::cout << "]\n";
    decoder.set(cID, layer_id, layer);
    for (int theta : thetaBins) {
      decoder.set(cID, theta_id, theta);
      for (int phi = 0; phi < seg->phiBins(); ++phi) {
        decoder.set(cID, phi_id, phi);
        cellFunc(cID);
      }
    }
  }

  return StatusCode::SUCCESS;
}
