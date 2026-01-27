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
  HCalPhiThetaCaloTool (const std::string& type,
                        const std::string& name,
                        const IInterface* parent);
  //using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~HCalPhiThetaCaloTool() = default;

  virtual unsigned index(uint64_t cellID) const;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const override final;
};

#endif /* RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H */
