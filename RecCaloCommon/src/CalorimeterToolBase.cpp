#include "RecCaloCommon/CalorimeterToolBase.h"
#include "k4Interface/IGeoSvc.h"
#include "k4FWCore/k4_check.h"
#include "DD4hep/Detector.h"
#include <algorithm>


StatusCode CalorimeterToolBase::initialize()
{
  K4_CHECK( AlgTool::initialize() );
  K4_CHECK( m_geoSvc.retrieve() );
  return StatusCode::SUCCESS;
}


StatusCode CalorimeterToolBase::getReadout (const std::string& readoutName)
{
  if (!readoutName.empty()) {
    // Check if readouts exist
    info() << "Readout: " << readoutName << endmsg;
    const dd4hep::Detector* det = geoSvc().getDetector();
    auto it = det->readouts().find(readoutName);
    if (it == det->readouts().end()) {
      error() << "Readout <<" << readoutName << ">> does not exist." << endmsg;
      return StatusCode::FAILURE;
    }
    m_readout = it->second;
  }
  return StatusCode::SUCCESS;
}


const IGeoSvc& CalorimeterToolBase::geoSvc() const
{
  return *m_geoSvc;
}


StatusCode CalorimeterToolBase::prepareEmptyCells(std::unordered_map<uint64_t, double>& aCells) const
{
  auto pushCell = [&aCells] (uint64_t cellId) { aCells.emplace (cellId, 0); };
  return collectCells(pushCell);
}


std::vector<uint64_t> CalorimeterToolBase::cellIDs() const
{
  std::vector<uint64_t> cells;
  auto pushCell = [&cells] (uint64_t cellId) { cells.push_back (cellId); };
  if (collectCells(pushCell).isSuccess()) {
    std::ranges::sort(cells);
    const auto ret = std::ranges::unique(cells);
    cells.erase(ret.begin(), ret.end());
  }
  else {
    cells.clear();
  }
  return cells;
}


