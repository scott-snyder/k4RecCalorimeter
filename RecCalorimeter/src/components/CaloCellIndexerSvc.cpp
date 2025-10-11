/**
 * @file CaloCellIndexerSvc.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for calorimeter cell indexers.
 */


#include "CaloCellIndexerSvc.h"
#include "k4FWCore/k4_check.h"


DECLARE_COMPONENT(k4::recCalo::CaloCellIndexerSvc);


namespace k4::recCalo {


/**
 * @brief Gaudi initialize method.
 */
StatusCode CaloCellIndexerSvc::initialize()
{
  K4_GAUDI_CHECK( Service::initialize() );
  K4_GAUDI_CHECK( m_geoTools.retrieve() );

  // Make indexers for all tools that support it.
  for (ToolHandle<ICalorimeterTool>& tool : m_geoTools) {
    std::unique_ptr<k4::recCalo::ICaloIndexer> indexer = tool->indexer();
    if (indexer) {
      int detid = tool->id();
      if (m_indexers.size() <= static_cast<size_t>(detid)) {
        m_indexers.resize (detid+1);
      }
      if (m_indexers[detid]) {
        error() << "Duplicated detector ID " << detid << endmsg;
        return StatusCode::FAILURE;
      }
      m_indexers[detid] = std::move(indexer);
    }
  }

  return StatusCode::SUCCESS;
}


/**
 * @brief Return indexer for a given subdetector.
 * @param detID Subdetector ID for the desired indexer.
 *
 * Returns a pointer to the indexer or nullptr if there isn't one defined.
 */
const k4::recCalo::ICaloIndexer* CaloCellIndexerSvc::indexer (int detID) const
{
  if (static_cast<size_t>(detID) < m_indexers.size()) {
    return m_indexers[detID].get();
  }
  return nullptr;
}


} // namespace k4::recCalo

