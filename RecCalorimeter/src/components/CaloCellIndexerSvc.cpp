/**
 * @file CaloCellIndexerSvc.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for calorimeter cell indexers.
 */


#include "CaloCellIndexerSvc.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"
#include "RecCaloCommon/MultiIndexer.h"


DECLARE_COMPONENT(k4::recCalo::CaloCellIndexerSvc);


namespace k4::recCalo {


/**
 * @brief Gaudi initialize method.
 */
StatusCode CaloCellIndexerSvc::initialize()
{
  K4RECCALORIMETER_CHECK( Service::initialize() );
  K4RECCALORIMETER_CHECK( m_constantsSvc.retrieve() );
  K4RECCALORIMETER_CHECK( m_geoTools.retrieve() );

  std::lock_guard lock (m_mutex);

  // Make indexers for all tools that support it.
  for (ToolHandle<ICalorimeterTool>& tool : m_geoTools) {
    std::unique_ptr<ICaloIndexer> indexer = tool->indexer();
    if (indexer) {
      int detid = tool->id();
      if (detid > MAX_DETID) {
        error() << "Out-of-range detector ID " << detid << endmsg;
        return StatusCode::FAILURE;
      }
      mask_t mask = static_cast<mask_t>(1) << detid;
      if (m_indexers.find(mask) != m_indexers.end()) {
        error() << "Duplicate detector ID " << detid << endmsg;
        return StatusCode::FAILURE;
      }
      m_indexers[mask] = std::move(indexer);
    }
  }

  return StatusCode::SUCCESS;
}


/**
 * @brief Return indexer for a given subdetector.
 * @param detID Subdetector ID for the desired indexer.
 * @param quiet If true, don't print an error if we don't find an indexer.
 *
 * Returns a pointer to the indexer or nullptr if there isn't one defined.
 */
const ICaloIndexer* CaloCellIndexerSvc::indexer (int detID,
                                                 bool quiet /*= false*/) const
{
  std::lock_guard lock (m_mutex);

  if (detID > MAX_DETID) {
    error() << "Out-of-range detector ID " << detID << endmsg;
    return nullptr;
  }
  mask_t mask = static_cast<mask_t>(1) << detID;
  auto it = m_indexers.find (mask);
  if (it != m_indexers.end()) {
    return it->second.get();
  }
  if (!quiet) {
    error() << "Cannot find indexer for detID " << detID << endmsg;
  }
  return nullptr;
}


/**
 * @brief Return indexer for a given set of subdetectors.
 */
const ICaloIndexer* CaloCellIndexerSvc::indexer (std::span<const int> detIDs,
                                                 bool quiet /*= false*/)
{
  // Degenerate cases
  if (detIDs.empty()) return nullptr;
  if (detIDs.size() == 1) {
    return indexer (detIDs[0], quiet);
  }

  std::lock_guard lock (m_mutex);

  // Build detID mask.
  mask_t mask = 0;
  for (int id : detIDs) {
    if (id > MAX_DETID) {
      error() << "Out-of-range detector ID " << id << endmsg;
      return nullptr;
    }
    mask |= (static_cast<mask_t> (1) << id);
  }

  // Do we already have this combination?
  auto it = m_indexers.find (mask);
  if (it != m_indexers.end()) {
    return it->second.get();
  }

  // Make a new MultiIndexer.
  std::vector<const ICaloIndexer*> indexers;
  size_t detIDBits = 0;
  for (int id : detIDs) {
    const ICaloIndexer* indexer = this->indexer (id, quiet);
    if (!indexer) return nullptr;
    indexers.push_back (indexer);
    if (detIDBits == 0)
      detIDBits = indexer->detIDBits();
    else if (detIDBits != indexer->detIDBits()) {
      error() << "Inconsistent detIDBits" << endmsg;
      return nullptr;
    }
  }

  auto mi = std::make_unique<MultiIndexer> (detIDBits,
                                            indexers,
                                            *m_constantsSvc);

  return m_indexers.emplace (mask, std::move(mi)).first->second.get();
}


} // namespace k4::recCalo

