#ifndef RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H
#define RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class HCalPhiThetaCaloTool
 *
 *  Generate all cell IDs for HCal.
 */
class HCalPhiThetaCaloTool : public CalorimeterToolBase
{
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~HCalPhiThetaCaloTool() = default;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const override final;
};

#endif /* RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H */
