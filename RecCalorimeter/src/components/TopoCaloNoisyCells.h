#ifndef RECCALORIMETER_TOPOCALONOISYCELLS_H
#define RECCALORIMETER_TOPOCALONOISYCELLS_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// k4FWCore
#include "k4Interface/INoiseConstTool.h"
#include "k4Interface/ICalorimeterTool.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include <utility>

class IGeoSvc;
class TFile;

/** @class TopoCaloNoisyCells Reconstruction/RecCalorimeter/src/components/TopoCaloNoisyCells.h
 *TopoCaloNoisyCells.h
 *
 *  Tool that reads a ROOT file containing the TTree with branchs "cellId", "noiseLevel", and "noiseOffset".
 *  This tool reads the tree, creates a map, and allows a lookup of noise level and mean noise of a cell, by its cellID.
 *
 *  @author Coralie Neubueser
 */

class TopoCaloNoisyCells : public extends<AlgTool, INoiseConstTool> {
public:
  using base_class::base_class;
  virtual ~TopoCaloNoisyCells() = default;
  /** Read a root file and the stored TTree of cellIDs to noise values.
   * return StatusCode
   */
  virtual StatusCode initialize() override final;

  /** Expected noise per cell in terms of sigma of Gaussian distibution.
   *   @param[in] aCellId of the cell of interest.
   *   return double.
   */
  virtual double getNoiseRMSPerCell(uint64_t aCellId) const override final;

  /** Expected noise per cell in terms of mean of distibution.
   *   @param[in] aCellId of the cell of interest.
   *   return double.
   */
  virtual double getNoiseOffsetPerCell(uint64_t aCellId) const override final;

  /** Expected noise per cell.
   *   @param[in] aCellId of the cell of interest.
   *   return [rms, offset]
   */
  virtual std::pair<double, double>
  getNoisePerCell(uint64_t aCellId) const override final;

private:
  /// Name
  Gaudi::Property<std::string> m_fileName{this, "fileName",
                                          "/afs/cern.ch/user/c/cneubuse/public/FCChh/cellNoise_map_segHcal.root"};
  ToolHandle<ICalorimeterTool> m_geoTool{this, "geometryTool", ""};
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
  { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "" };

  // rms, offset
  using NoiseData = std::vector<std::pair<double, double> >;
  const NoiseData* m_data = nullptr;

  NoiseData readData (TFile& inFile) const;
};

#endif /* RECCALORIMETER_TOPOCALONOISYCELLS_H */
