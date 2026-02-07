/**
 * @file RecFCCeeCalorimeter/tests/src/HCalPhiThetaCaloToolTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Test for HCalPhiThetaCaloTool
 */

#undef NDEBUG
#include "k4Interface/ICalorimeterTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"
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
  StatusCode testTool (const ICalorimeterTool& tool, size_t exp_ncells) const;

  ToolHandle<ICalorimeterTool> m_barrelTool
  { this, "HCalBarrelTool", "HCalPhiThetaCaloTool", "" };
  ToolHandle<ICalorimeterTool> m_endcapTool
  { this, "HCalEndcapTool", "HCalPhiThetaCaloTool", "" };
};


DECLARE_COMPONENT(k4::recCalo::HCalPhiThetaCaloToolTestAlg);


StatusCode HCalPhiThetaCaloToolTestAlg::initialize()
{
  K4RECCALORIMETER_CHECK( m_barrelTool.retrieve() );
  K4RECCALORIMETER_CHECK( m_endcapTool.retrieve() );

  K4RECCALORIMETER_CHECK (m_barrelTool->readoutName() == "HCalBarrelReadout");
  K4RECCALORIMETER_CHECK (m_barrelTool->id() == 8);
  K4RECCALORIMETER_CHECK (m_endcapTool->readoutName() == "HCalEndcapReadout");
  K4RECCALORIMETER_CHECK (m_endcapTool->id() == 9);

  K4RECCALORIMETER_CHECK( testTool (*m_barrelTool, 210944) );
  K4RECCALORIMETER_CHECK( testTool (*m_endcapTool, 80896) );

  return StatusCode::SUCCESS;
}


StatusCode HCalPhiThetaCaloToolTestAlg::testTool (const ICalorimeterTool& tool,
                                                  size_t exp_ncells) const
{
  std::span<const uint64_t> ids = tool.cellIDs();
  size_t ncells = ids.size();
  K4RECCALORIMETER_CHECK (ncells == exp_ncells);

  std::unique_ptr<k4::recCalo::ICaloIndexer> indexer = tool.indexer();
  K4RECCALORIMETER_CHECK (indexer != nullptr);
  K4RECCALORIMETER_CHECK( indexer->detIDBits() == 4 );

  K4RECCALORIMETER_CHECK (indexer->cellIDs().size() == ncells);
  for (size_t i = 0; i < ncells; i++) {
    K4RECCALORIMETER_CHECK (ids[i] == indexer->cellIDs()[i]);
    K4RECCALORIMETER_CHECK (indexer->index(ids[i]) == i);
  }

  return StatusCode::SUCCESS;
}


StatusCode HCalPhiThetaCaloToolTestAlg::execute()
{
  return StatusCode::SUCCESS;
}


} // namespace k4::recCalo
