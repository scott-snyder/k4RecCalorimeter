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

  virtual StatusCode initialize() override final;


protected:
  virtual StatusCode collectCells(std::function<void(uint64_t)> cellFunc) const override final;

private:
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", ""};
};

#endif /* RECFCCEECALORIMETER_HCALPHITHETACALOTOOL_H */
