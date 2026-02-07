/**
 * @file RecFCCeeCalorimeter/tests/src/HCalPhiThetaCaloToolTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Test for HCalPhiThetaCaloTool
 */

#undef NDEBUG
#include "RecCaloCommon/ICalorimeterTool.h"
#include "k4FWCore/GaudiChecks.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include <span>
#include <cassert>


namespace k4::recCalo {


class HCalPhiThetaCaloToolTestAlg
  : public Algorithm
{
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  StatusCode testTool (const k4::recCalo::ICalorimeterTool& tool, size_t exp_ncells) const;

  ToolHandle<k4::recCalo::ICalorimeterTool> m_barrelTool
  { this, "HCalBarrelTool", "HCalPhiThetaCaloTool", "" };
  ToolHandle<k4::recCalo::ICalorimeterTool> m_endcapTool
  { this, "HCalEndcapTool", "HCalPhiThetaCaloTool", "" };
};


DECLARE_COMPONENT(k4::recCalo::HCalPhiThetaCaloToolTestAlg);


StatusCode HCalPhiThetaCaloToolTestAlg::initialize()
{
  K4_GAUDI_CHECK( m_barrelTool.retrieve() );
  K4_GAUDI_CHECK( m_endcapTool.retrieve() );

  K4_GAUDI_CHECK (m_barrelTool->readoutName() == "HCalBarrelReadout");
  K4_GAUDI_CHECK (m_barrelTool->id() == 8);
  K4_GAUDI_CHECK (m_endcapTool->readoutName() == "HCalEndcapReadout");
  K4_GAUDI_CHECK (m_endcapTool->id() == 9);

  K4_GAUDI_CHECK( testTool (*m_barrelTool, 210944) );
  K4_GAUDI_CHECK( testTool (*m_endcapTool, 80896) );

  return StatusCode::SUCCESS;
}


StatusCode HCalPhiThetaCaloToolTestAlg::testTool (const k4::recCalo::ICalorimeterTool& tool,
                                                  size_t exp_ncells) const
{
  std::span<const uint64_t> ids = tool.cellIDs();
  size_t ncells = ids.size();
  K4_GAUDI_CHECK (ncells == exp_ncells);

  std::unique_ptr<ICaloIndexer> indexer = tool.indexer();
  K4_GAUDI_CHECK (indexer != nullptr);
  K4_GAUDI_CHECK( indexer->detIDBits() == 4 );

  K4_GAUDI_CHECK (indexer->cellIDs().size() == ncells);
  for (size_t i = 0; i < ncells; i++) {
    K4_GAUDI_CHECK (ids[i] == indexer->cellIDs()[i]);
    K4_GAUDI_CHECK (indexer->index(ids[i]) == i);
  }

  return StatusCode::SUCCESS;
}


StatusCode HCalPhiThetaCaloToolTestAlg::execute()
{
  return StatusCode::SUCCESS;
}


} // namespace k4::recCalo
