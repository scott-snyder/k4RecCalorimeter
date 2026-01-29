// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/ICaloCellIndexerSvc.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Holder for calorimeter cell indexers.
 */


#ifndef RECCALOCOMMON_ICALOCELLINDEXERSVC_H
#define RECCALOCOMMON_ICALOCELLINDEXERSVC_H


#include "k4Interface/ICaloIndexer.h"
#include "GaudiKernel/IInterface.h"


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
class ICaloCellIndexerSvc
  : virtual public IInterface
{
public:
  DeclareInterfaceID (ICaloCellIndexerSvc, 1, 0);


  /**
   * @brief Return indexer for a given subdetector.
   * @param detID Subdetector ID for the desired indexer.
   *
   * Returns a pointer to the indexer or nullptr if there isn't one defined.
   */
  virtual const k4::recCalo::ICaloIndexer* indexer (int detID) const = 0;
};


} // namespace k4::recCalo



#endif // not RECCALOCOMMON_ICALOCELLINDEXERSVC_H
