#include "NoiseCaloCellsVsThetaFromFileTool.h"
#include "k4FWCore/k4_check.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"

DECLARE_COMPONENT(NoiseCaloCellsVsThetaFromFileTool)

NoiseCaloCellsVsThetaFromFileTool::NoiseCaloCellsVsThetaFromFileTool(const std::string& type,
                                                                     const std::string& name,
                                                                     const IInterface* parent)
: base_class (type, name, parent)
{
  // Override some property defaults from the base class.
  m_readoutName = "ECalBarrelThetaModuleMerged";
}

StatusCode NoiseCaloCellsVsThetaFromFileTool::initialize() {
  m_elecNoiseHistoName = m_elecNoiseRMSHistoName;

  K4_GAUDI_CHECK( m_cellPositionsTool.retrieve() );

  K4_GAUDI_CHECK( ReadNoiseFromFileTool::initialize() );

  K4_GAUDI_CHECK( m_randSvc = service<IRndmGenSvc> ("RndmGenSvc", false) );
  K4_GAUDI_CHECK( m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)) );

  debug() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;

  return StatusCode::SUCCESS;
}

template <class C>
void NoiseCaloCellsVsThetaFromFileTool::addRandomCellNoiseT(C& aCells) const {
  for (auto& p : aCells) {
    p.second += getNoiseOffsetPerCell(p.first);
    p.second += (getNoiseRMSPerCell(p.first) * m_gauss.shoot());
  }
}

void NoiseCaloCellsVsThetaFromFileTool::addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  addRandomCellNoiseT(aCells);
}

void NoiseCaloCellsVsThetaFromFileTool::addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  addRandomCellNoiseT (aCells);
}

template <typename C>
void NoiseCaloCellsVsThetaFromFileTool::filterCellNoiseT(C& aCells) const {
  // Erase a cell if it has energy below a threshold from the vector
  if (m_useAbsInFilter) {
    std::erase_if(aCells,
                  [&](auto& p) { return std::abs(p.second-getNoiseOffsetPerCell(p.first)) < m_filterThreshold * getNoiseRMSPerCell(p.first); });
  }
  else {
    std::erase_if(aCells,
                  [&](auto& p) { return p.second < getNoiseOffsetPerCell(p.first) + m_filterThreshold * getNoiseRMSPerCell(p.first); });
  }
}

void NoiseCaloCellsVsThetaFromFileTool::filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  filterCellNoiseT (aCells);
}

void NoiseCaloCellsVsThetaFromFileTool::filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  filterCellNoiseT (aCells);
}


StatusCode NoiseCaloCellsVsThetaFromFileTool::initBinning (NoiseData& data,
                                                           const k4::recCalo::ICaloIndexer& indexer) const
{
  data.m_bins.resize (indexer.cellIDs().size());

  TH1F* h_rms = &data.m_histoElecNoiseRMS.at(0);
  TH1F* h_offset = m_setNoiseOffset ? &data.m_histoElecNoiseOffset.at(0) : nullptr;
  for (uint64_t id : indexer.cellIDs()) {
    unsigned ndx = indexer.index (id);
    double cellTheta = m_cellPositionsTool->xyzPosition(id).Theta();
    int ibinRMS = h_rms->FindFixBin(cellTheta);
    int ibinOffset = h_offset ? h_offset->FindFixBin(cellTheta) : 0;
    data.m_bins.at(ndx) = std::make_pair(ibinRMS, ibinOffset);
  }

  return StatusCode::SUCCESS;
}
