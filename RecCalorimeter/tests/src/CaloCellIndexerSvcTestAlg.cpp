/**
 * @file RecCalorimeter/tests/src/CaloCellIndexerSvcTestAlg.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Test for CaloCellIndexerSvc
 */


#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "k4FWCore/k4_check.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/ServiceHandle.h"
#include <span>


namespace k4::recCalo {


class CaloCellIndexerSvcTestAlg
  : public Algorithm
{
public:
  using Algorithm::Algorithm;

  virtual StatusCode initialize() override;
  virtual StatusCode execute() override;

private:
  ServiceHandle<ICaloCellIndexerSvc> m_svc
  { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };

  Gaudi::Property<int> m_detID
  { this, "DetID", 4, "" };
};


DECLARE_COMPONENT(k4::recCalo::CaloCellIndexerSvcTestAlg);


StatusCode CaloCellIndexerSvcTestAlg::initialize()
{
  K4_GAUDI_CHECK( m_svc.retrieve() );

  K4_GAUDI_CHECK( m_svc->indexer (999) == nullptr );

  const ICaloIndexer* indexer = m_svc->indexer (m_detID);
  K4_GAUDI_CHECK( indexer != nullptr );
  K4_GAUDI_CHECK( indexer->detIDs().size() == 1 );
  K4_GAUDI_CHECK( indexer->detIDs()[0] == m_detID );

  std::span<const uint64_t> ids = indexer->cellIDs();
  for (size_t i = 0; i < ids.size(); ++i) {
    K4_GAUDI_CHECK( indexer->index(ids[i]) == i );
  }

  return StatusCode::SUCCESS;
}


StatusCode CaloCellIndexerSvcTestAlg::execute()
{
  return StatusCode::SUCCESS;
}


} // namespace k4::recCalo

