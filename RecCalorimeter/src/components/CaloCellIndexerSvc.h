// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCalorimeter/src/components/CaloCellIndexerSvc.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for calorimeter cell indexers.
 */


#ifndef RECCALORIMETER_CALOCELLINDEXERSVC_H
#define RECCALORIMETER_CALOCELLINDEXERSVC_H


#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "k4Interface/ICalorimeterTool.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ToolHandle.h"
#include "GaudiKernel/ServiceHandle.h"
#include <mutex>
#include <climits>


namespace k4::recCalo {


/**
 * @brief Holder for calorimeter cell indexers.
 *
 * This holds the indexer objects for different subdetectors.
 * We put it in a service to avoid duplicating them across different tools.
 *
 * One might think of having this be part of ICaloCellConstantsSvc,
 * but we run into initialization loops in that case.
 */
class CaloCellIndexerSvc
  : public extends<Service, ICaloCellIndexerSvc>
{
public:
  using base_class::base_class;


  /**
   * @brief Gaudi initialize method.
   */
  virtual StatusCode initialize() override;


  /**
   * @brief Return indexer for a given subdetector.
   * @param detID Subdetector ID for the desired indexer.
   * @param quiet If true, don't print an error if we don't find an indexer.
   *
   * Returns a pointer to the indexer or nullptr if there isn't one defined.
   */
  const ICaloIndexer* indexer (int detID, bool quiet = false) const override;


  /**
   * @brief Return indexer for a given set of subdetectors.
   * @param detIDs Subdetector IDs to index.
   * @param quiet If true, don't print an error if we don't find an indexer.
   *
   * Returns a pointer to the indexer or nullptr if there isn't one defined.
   */
  virtual const ICaloIndexer* indexer (std::span<const int> detIDs,
                                       bool quiet = false) override;


private:
  ToolHandleArray<ICalorimeterTool> m_geoTools
  { this, "GeoTools", {} };

  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
  { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc" };

  /// Map from bit mask of detector IDs to indexer objects.
  using mask_t = uint64_t;
  constexpr static int MAX_DETID = sizeof(mask_t) * CHAR_BIT;
  std::map<mask_t, std::unique_ptr<ICaloIndexer> > m_indexers;

  /// Guard access to the map.
  mutable std::recursive_mutex m_mutex;
};


} // namespace k4::recCalo


#endif // not RECCALORIMETER_CALOCELLINDEXERSVC_H
