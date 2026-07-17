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

static FILE* flog = 0;
template <class C>
void NoiseCaloCellsFromFileBaseTool::addRandomCellNoiseT(C& aCells) const {
  if (!flog) {
    flog = fopen ("f.log", "w");
  }
  fprintf (flog, "addRandomCellNoise\n");
  for (auto& p : aCells) {
    double orig = p.second;
    double offs = getNoiseOffsetPerCell(p.first);
    double rms = getNoiseRMSPerCell(p.first);
    double noise = (rms * m_gauss.shoot());
    p.second += offs;
    p.second += noise;
    fprintf (flog, "%lu %f %f %f %f %f\n",
             p.first, orig, offs, rms, noise, p.second);
  }
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::unordered_map<CellID, double>& aCells) const {
  addRandomCellNoiseT(aCells);
}

void NoiseCaloCellsFromFileBaseTool::addRandomCellNoise(std::vector<std::pair<CellID, double>>& aCells) const {
  addRandomCellNoiseT(aCells);
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

