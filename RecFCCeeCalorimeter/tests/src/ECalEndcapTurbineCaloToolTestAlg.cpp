/**
 * @file RecFCCeeCalorimeter/tests/src/ECalEndcapTurbineCaloToolTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Test for ECalEndcapTurbineCaloTool
 */

#undef NDEBUG
#include "k4Interface/ICalorimeterTool.h"
#include "k4FWCore/k4_check.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include <span>
#include <cassert>


namespace k4::recCalo {


class ECalEndcapTurbineCaloToolTestAlg
  : public Algorithm
{
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  ToolHandle<ICalorimeterTool> m_tool
  { this, "ECalEndcapTurbineTool", "ECalEndcapTurbineCaloTool", "" };
};


DECLARE_COMPONENT(k4::recCalo::ECalEndcapTurbineCaloToolTestAlg);


StatusCode ECalEndcapTurbineCaloToolTestAlg::initialize()
{
  K4_GAUDI_CHECK( m_tool.retrieve() );

  K4_GAUDI_CHECK (m_tool->readoutName() == "ECalEndcapTurbine");
  K4_GAUDI_CHECK (m_tool->id() == 5);

  std::span<const uint64_t> ids = m_tool->cellIDs();
  size_t ncells = ids.size();
  K4_GAUDI_CHECK (ncells == 1203200);

  std::unique_ptr<k4::recCalo::ICaloIndexer> indexer = m_tool->indexer();
  K4_GAUDI_CHECK (indexer != nullptr);
  K4_GAUDI_CHECK( indexer->detIDBits() == 4 );

  K4_GAUDI_CHECK (indexer->cellIDs().size() == ncells);
  for (size_t i = 0; i < ncells; i++) {
    K4_GAUDI_CHECK (ids[i] == indexer->cellIDs()[i]);
    K4_GAUDI_CHECK (indexer->index(ids[i]) == i);
  }

  return StatusCode::SUCCESS;
}


StatusCode ECalEndcapTurbineCaloToolTestAlg::execute()
{
  return StatusCode::SUCCESS;
}


} // namespace k4::recCalo
