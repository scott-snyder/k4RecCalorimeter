#include "CalibrateCaloHitsTool.h"
#include "k4FWCore/GaudiChecks.h"

DECLARE_COMPONENT(CalibrateCaloHitsTool)

StatusCode CalibrateCaloHitsTool::initialize() {
  K4_GAUDI_CHECK( AlgTool::initialize() );

  info() << "Calibration constant: 1/sampling fraction=" << m_invSamplingFraction << endmsg;
  return StatusCode::SUCCESS;
}

void CalibrateCaloHitsTool::calibrate(std::unordered_map<uint64_t, double>& aHits) const {
  // Loop through energy deposits, multiply energy to get cell energy at electromagnetic scale
  for (auto& p : aHits) {
    p.second *= m_invSamplingFraction;
  }
}

void CalibrateCaloHitsTool::calibrate(std::vector<std::pair<uint64_t, double> >& aHits) const {
  // Loop through energy deposits, multiply energy to get cell energy at electromagnetic scale
  for (auto& p : aHits) {
    p.second *= m_invSamplingFraction;
  }
}

