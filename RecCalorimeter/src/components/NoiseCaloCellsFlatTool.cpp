#include "NoiseCaloCellsFlatTool.h"
#include "k4FWCore/k4_check.h"
#include <GaudiKernel/StatusCode.h>

DECLARE_COMPONENT(NoiseCaloCellsFlatTool)

StatusCode NoiseCaloCellsFlatTool::initialize() {
  K4_GAUDI_CHECK( AlgTool::initialize() );
  K4_GAUDI_CHECK( m_randSvc = service<IRndmGenSvc> ("RndmGenSvc", true) );
  K4_GAUDI_CHECK( m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)) );

  info() << "RMS of the cell noise: " << m_cellNoiseRMS * 1.e3 << " MeV" << endmsg;
  info() << "Offset of the cell noise: " << m_cellNoiseOffset * 1.e3 << " MeV" << endmsg;
  info() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;
  return StatusCode::SUCCESS;
}

template <typename C>
void NoiseCaloCellsFlatTool::addRandomCellNoiseT (C& aCells) const
{
  for (auto& p : aCells) {
    p.second += m_cellNoiseOffset + (m_gauss.shoot() * m_cellNoiseRMS);
  }
}

void NoiseCaloCellsFlatTool::addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  using p_t = std::pair<uint64_t, double>;
  std::vector<p_t> cells (aCells.begin(), aCells.end());
  std::ranges::sort (cells, [](const p_t& a, const p_t& b) { return a.first < b.first; });
  addRandomCellNoiseT (aCells);
  for (const p_t& p : cells) aCells[p.first] = p.second;
}

void NoiseCaloCellsFlatTool::addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  addRandomCellNoiseT (aCells);
}

template <typename C>
void NoiseCaloCellsFlatTool::filterCellNoiseT (C& aCells) const
{
  // Erase a cell if it has energy below a threshold
  double threshold = m_cellNoiseOffset + m_filterThreshold * m_cellNoiseRMS;
  std::erase_if (aCells, [threshold] (auto& p) { return p.second < threshold; });
}

void NoiseCaloCellsFlatTool::filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const
{
  filterCellNoiseT (aCells);
}

void NoiseCaloCellsFlatTool::filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const
{
  filterCellNoiseT (aCells);
}
