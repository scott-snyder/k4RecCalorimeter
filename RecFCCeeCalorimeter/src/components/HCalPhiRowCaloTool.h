// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecFCCeeCalorimeter/src/components/HCalPhiRowCaloTool.h
 * @author scott snyder <snyder@bnl.gov>
 * @date May, 2026
 * @brief Calorimeter tool for Allegro HCal (with row indexing).
 */


#ifndef RECFCCEECALORIMETER_HCALPHIROWCALOTOOL_H
#define RECFCCEECALORIMETER_HCALPHIROWCALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class HCalPhiRowCaloTool
 *
 *  Manage IDs for HCal (with row indexing).
 */
class HCalPhiRowCaloTool : public CalorimeterToolBase
{
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~HCalPhiRowCaloTool() = default;

  /** Gaudi initialize method.
   */
  virtual StatusCode initialize() override final;

  /** Return the subdetector ID.
   */
  virtual int id() const final override;

  /** Return a new indexer object for this subdetector.
   */
  virtual std::unique_ptr<k4::recCalo::ICaloIndexer> indexer() const override final;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const override final;

private:
  /// Detector ID.
  int m_id = -1;
};

#endif /* RECFCCEECALORIMETER_HCALPHIROWCALOTOOL_H */
