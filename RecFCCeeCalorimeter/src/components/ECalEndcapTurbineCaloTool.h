// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file RecFCCeeCalorimeter/src/components/ECalEndcapTurbineCaloTool.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Calorimeter tool for Allegro ECal endcap.
 */


#ifndef RECFCCEECALORIMETER_ECALENDCAPTURBINECALOTOOL_H
#define RECFCCEECALORIMETER_ECALENDCAPTURBINECALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"


namespace dd4hep::DDSegmentation {
class FCCSWEndcapTurbine_k4geo;
}


/** @class ECalEndcapTurbineCaloTool
 *
 *  Manage IDs for HCal.
 */
class ECalEndcapTurbineCaloTool : public CalorimeterToolBase
{
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~ECalEndcapTurbineCaloTool() = default;


  /** Return a new indexer object for this subdetector.
   */
  virtual std::unique_ptr<ICaloIndexer> indexer() const override final;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const override final;


private:
  struct WheelParams
  {
    WheelParams (const dd4hep::DDSegmentation::FCCSWEndcapTurbine_k4geo& seg,
                 int iwheel);
    int nmodules;
    int nrho;
    int nz;
  };
};


#endif // not RECFCCEECALORIMETER_ECALENDCAPTURBINECALOTOOL_H
