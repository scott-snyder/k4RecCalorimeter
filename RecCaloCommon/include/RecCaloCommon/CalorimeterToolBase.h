// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecCaloCommon/CalorimeterToolBase.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Aug, 2025
 * @brief Common base for ICalorimeterTool implementations.
 */


#ifndef RECCALOCOMMON_CALORIMETERTOOLBASE_H
#define RECCALOCOMMON_CALORIMETERTOOLBASE_H


#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ServiceHandle.h"
#include "k4Interface/ICalorimeterTool.h"
#include "DD4hep/Readout.h"


class IGeoSvc;



/** @class CalorimeterToolBase RecCaloCommon/include/CalorimeterToolBase.h
 *
 * This factors out the implementations of prepareEmptyCells/cellIDs
 * in terms of a common method collectCells().
 */
class CalorimeterToolBase : public extends<AlgTool, ICalorimeterTool>
{
public:
  using base_class::base_class;

  virtual StatusCode initialize() override;

  virtual StatusCode prepareEmptyCells(std::unordered_map<uint64_t, double>& aCells) const override final;

  virtual std::vector<uint64_t> cellIDs() const override final;


protected:
  StatusCode getReadout (const std::string& readoutName);

  const dd4hep::Readout readout() const { return m_readout; }

  virtual StatusCode collectCells(std::function<void(uint64_t)> cellFunc) const = 0;


private:
  const IGeoSvc& geoSvc() const;

  ServiceHandle<IGeoSvc> m_geoSvc { this, "GeoSvc", "GeoSvc" };

  dd4hep::Readout m_readout;
};


#endif // not RECCALOCOMMON_CALORIMETERTOOLBASE_H
