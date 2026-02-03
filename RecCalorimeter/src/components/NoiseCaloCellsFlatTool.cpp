#include "NoiseCaloCellsFlatTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"
#include <GaudiKernel/StatusCode.h>

DECLARE_COMPONENT(NoiseCaloCellsFlatTool)

StatusCode NoiseCaloCellsFlatTool::initialize() {
  K4RECCALORIMETER_CHECK(AlgTool::initialize());
  K4RECCALORIMETER_CHECK(m_randSvc = service<IRndmGenSvc>("RndmGenSvc", true));
  K4RECCALORIMETER_CHECK(m_gauss.initialize(m_randSvc, Rndm::Gauss(0., 1.)));

  info() << "RMS of the cell noise: " << m_cellNoiseRMS * 1.e3 << " MeV" << endmsg;
  info() << "Offset of the cell noise: " << m_cellNoiseOffset * 1.e3 << " MeV" << endmsg;
  info() << "Filter noise threshold: " << m_filterThreshold << "*sigma" << endmsg;
  return StatusCode::SUCCESS;
}

void NoiseCaloCellsFlatTool::initEvent(CLHEP::Ranlux64Engine& e) const
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

template <typename C>
void NoiseCaloCellsFlatTool::addRandomCellNoiseT (C& aCells, CLHEP::RandGauss& r) const
{
  std::map<uint64_t, typename C::iterator> m;
  for (auto i = aCells.begin(); i != aCells.end(); ++i) m[i->first] = i;
  for (auto& p : m) {
    p.second->second += m_cellNoiseOffset;
    p.second->second += r.fire() * m_cellNoiseRMS;
  }
#if 0
  for (auto& p : aCells) {
    p.second += m_cellNoiseOffset + (m_gauss.shoot() * m_cellNoiseRMS);
  }
#endif
}

void NoiseCaloCellsFlatTool::addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  CLHEP::Ranlux64Engine e;
  initEvent(e);
  CLHEP::RandGauss r(e);
  using p_t = std::pair<uint64_t, double>;
  std::vector<p_t> cells (aCells.begin(), aCells.end());
  std::ranges::sort (cells, [](const p_t& a, const p_t& b) { return a.first < b.first; });
  addRandomCellNoiseT (aCells, r);
  for (const p_t& p : cells) aCells[p.first] = p.second;
}

void NoiseCaloCellsFlatTool::addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const {
  CLHEP::Ranlux64Engine e;
  initEvent(e);
  CLHEP::RandGauss r(e);
  addRandomCellNoiseT(aCells, r);
}

template <typename C>
void NoiseCaloCellsFlatTool::filterCellNoiseT(C& aCells) const {
  // Erase a cell if it has energy below a threshold
  double threshold = m_cellNoiseOffset + m_filterThreshold * m_cellNoiseRMS;
  std::erase_if(aCells, [threshold](auto& p) { return p.second < threshold; });
}

void NoiseCaloCellsFlatTool::filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const {
  filterCellNoiseT(aCells);
}

void NoiseCaloCellsFlatTool::filterCellNoise(std::vector<std::pair<uint64_t, double>>& aCells) const {
  filterCellNoiseT(aCells);
}
