/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
/**
 * @file RecFCCeeCalorimeter/tests/src/TubeLayerModuleThetaCaloToolTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Test for TubeLayerModuleThetaCaloToolTest
 */

#undef NDEBUG
#include "k4Interface/ICalorimeterTool.h"
#include "k4FWCore/k4_check.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"
#include <span>
#include <cassert>


namespace k4::recCalo {


class TubeLayerModuleThetaCaloToolTestAlg
  : public Algorithm
{
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  ToolHandle<ICalorimeterTool> m_tool
  { this, "CalorimeterTool", "TubeLayerModuleThetaCaloTool", "" };
};


DECLARE_COMPONENT(k4::recCalo::TubeLayerModuleThetaCaloToolTestAlg);


StatusCode TubeLayerModuleThetaCaloToolTestAlg::initialize()
{
  K4_GAUDI_CHECK( m_tool.retrieve() );
  assert (m_tool->readoutName() == "ECalBarrelModuleThetaMerged");
  assert (m_tool->id() == 4);

  std::span<const uint64_t> ids = m_tool->cellIDs();
  size_t ncells = ids.size();
  assert (ncells == 2041344);

  std::unique_ptr<k4::recCalo::ICaloIndexer> indexer = m_tool->indexer();
  assert (indexer != nullptr);

  assert (indexer->cellIDs().size() == ncells);
  for (size_t i = 0; i < ncells; i++) {
    assert (ids[i] == indexer->cellIDs()[i]);
    assert (indexer->index(ids[i]) == i);
  }

  return StatusCode::SUCCESS;
}


StatusCode TubeLayerModuleThetaCaloToolTestAlg::execute()
{
  return StatusCode::SUCCESS;
}


} // namespace k4::recCalo

