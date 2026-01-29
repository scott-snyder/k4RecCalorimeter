// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCalorimeter/src/components/CaloCellIndexerSvc.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for calorimeter cell indexers.
 */


#ifndef RECCALORIMETER_CALOCELLINDEXERSVC_H
#define RECCALORIMETER_CALOCELLINDEXERSVC_H


#include "k4Interface/ICalorimeterTool.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "GaudiKernel/Service.h"
#include "GaudiKernel/ToolHandle.h"


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
   *
   * Returns a pointer to the indexer or nullptr if there isn't one defined.
   */
  const k4::recCalo::ICaloIndexer* indexer (int detID) const override;


private:
  ToolHandleArray<ICalorimeterTool> m_geoTools
  { this, "GeoTools", {} };

  /// Indexer objects.
  std::vector<std::unique_ptr<ICaloIndexer> > m_indexers;
};


} // namespace k4::recCalo


#endif // not RECCALORIMETER_CALOCELLINDEXERSVC_H
