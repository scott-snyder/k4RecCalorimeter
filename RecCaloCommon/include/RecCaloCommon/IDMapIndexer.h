// This file's extension implies that it's C, but it's really -*- C++ -*-.
/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
/**
 * @file RecCaloCommon/IDMapIndexer.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Implementation of ICaloIndexer using IDMap.
 */


#ifndef RECCALOCOMMON_IDMAPINDEXER_H
#define RECCALOCOMMON_IDMAPINDEXER_H


#include "RecCaloCommon/IDMap.h"
#include "k4Interface/ICaloIndexer.h"


namespace k4::recCalo {


/**
 * @brief Implementation of ICaloIndexer using IDMap.
 *
 * This is an implementation of @c ICaloIndexer for a specific subdetector
 * based on @c IDMap with a fixed number of fields.
 */
template <unsigned NFIELDS>
class IDMapIndexer : public ICaloIndexer
{
public:
  /// Type of an index.
  using index_t = ICaloIndexer::index_t;

  /// Flag indicating an invalid index.
  static constexpr index_t INVALID = static_cast<index_t>(-1);

  /// Type of the mapping.
  using IDMap_t = IDMapN<index_t, NFIELDS>;

  /// Type describing a field; used to initialize the mapping.
  using FieldDesc_t = typename IDMap_t::FieldDesc_t;


  /**
   * @brief Constructor.
   * @param fields Set of fields to use from the identifiers.
   *               Must be sufficient to make identifiers unique.
   *               For best results, should be listed in order of increasing
   *               bit width.
   *               The size must be exactly @c NFIELDS.
   * @param ids Set of all identifiers to be indexed.
   *            Should be sorted for best results.
   * param sizeHint If non-zero, this is an estimate of the total size,
   *                in bytes, required by this mapping.  This will be
   *                used to reserve an appropriate size for the data vector.
   */
  IDMapIndexer (std::span<const FieldDesc_t> fields,
                std::span<const uint64_t> ids,
                size_t sizeHint = 0);


  /**
   * @brief Return the index of an identifier.
   * @param id The identifier to look for.
   *
   * Returns the index of @c id in @c cellIDs(), or @c INVALID.
   */
  virtual index_t index (uint64_t id) const override final;


  /**
   * @brief Return the set of all identifiers that we index.
   */
  virtual std::span<const uint64_t> cellIDs() const override final;

  
private:
  /// The mapping.
  IDMap_t m_map;

  /// The set of cells that we index.
  std::span<const uint64_t> m_cellIDs;
};


template <unsigned NFIELDS>
IDMapIndexer<NFIELDS>::IDMapIndexer (std::span<const FieldDesc_t> fields,
                                     std::span<const uint64_t> ids,
                                     size_t sizeHint /*= 0*/)
  : m_map (fields, INVALID,ids,
           [](size_t i) { return i; },
           sizeHint),
    m_cellIDs (ids)
{
}


/**
 * @brief Return the index of an identifier.
 */
template <unsigned NFIELDS>
inline
auto IDMapIndexer<NFIELDS>::index (uint64_t id) const -> index_t
{
  return m_map.lookup (id);
}


/**
 * @brief Return the set of all identifiers that we index.
 */
template <unsigned NFIELDS>
inline
std::span<const uint64_t> IDMapIndexer<NFIELDS>::cellIDs() const
{
  return m_cellIDs;
}


} // namespace k4::recCalo

#endif // not RECCALOCOMMON_IDMAPINDEXER_H
