#include "NoiseCaloCellsFromFileBaseTool.h"
#include "k4FWCore/GaudiChecks.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"

// ROOT
#include "TFile.h"
#include "TSystem.h"

StatusCode NoiseCaloCellsFromFileBaseTool::initialize() {
  m_elecNoiseHistoName = m_elecNoiseRMSHistoName;

  K4_GAUDI_CHECK(base_class::initialize());
  K4_GAUDI_CHECK(m_randSvc = service<IRndmGenSvc>("RndmGenSvc", false));
  K4_GAUDI_CHECK(m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)));

  debug() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;

  return StatusCode::SUCCESS;
}

template <class C>
void NoiseCaloCellsFromFileBaseTool::addRandomCellNoiseT(C& aCells, CLHEP::RandGauss& r) const {
  for (auto& p : aCells) {
    p.second += getNoiseOffsetPerCell(p.first);
    p.second += (getNoiseRMSPerCell(p.first) * r.shoot());
  }
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::unordered_map<CellID, double>& aCells) const {
  CLHEP::Ranlux64Engine e;
  initEvent(e);
  CLHEP::RandGauss r(e);

  using p_t = std::pair<uint64_t, double>;
  std::vector<p_t> cells (aCells.begin(), aCells.end());
  std::ranges::sort (cells, [](const p_t& a, const p_t& b) { return a.first < b.first; });
  addRandomCellNoiseT(cells, r);
  for (const p_t& p : cells) aCells[p.first] = p.second;
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::vector<std::pair<CellID, double>>& aCells) const {
  CLHEP::Ranlux64Engine e;
  initEvent(e);
  CLHEP::RandGauss r(e);
  addRandomCellNoiseT(aCells, r);
}

template <typename C>
void NoiseCaloCellsFromFileBaseTool::filterCellNoiseT(C& aCells) const {
  // Erase a cell if it has energy below a threshold from the vector
  if (m_useAbsInFilter) {
    std::erase_if(aCells, [&](auto& p) {
      return std::abs(p.second - getNoiseOffsetPerCell(p.first)) < m_filterThreshold * getNoiseRMSPerCell(p.first);
    });
  } else {
    std::erase_if(aCells, [&](auto& p) {
      return p.second < getNoiseOffsetPerCell(p.first) + m_filterThreshold * getNoiseRMSPerCell(p.first);
    });
  }
}

void NoiseCaloCellsFromFileBaseTool::filterCellNoise(std::unordered_map<CellID, double>& aCells) const {
  filterCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileBaseTool::filterCellNoise(std::vector<std::pair<CellID, double>>& aCells) const {
  filterCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileBaseTool::initEvent(CLHEP::Ranlux64Engine& e) const
{
  const edm4hep::EventHeaderCollection* ehs = m_header.get();
  edm4hep::EventHeader eh = ehs->at(0);
  // FIXME: Use all bits of event number.
  long seeds[] = {static_cast<long>(eh.getRunNumber()),
                  static_cast<long>(eh.getEventNumber()),
                  0x76439862,
                  0};
  if (seeds[0] == 0) seeds[0] = 0x7fffffff;
  if (seeds[1] == 0) seeds[1] = 0x7fffffff;
  e.setSeeds (seeds);
}

