#ifndef RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILETOOL_H
#define RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/IRndmGenSvc.h"
#include "GaudiKernel/RndmGenerators.h"

// k4geo
// #include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"

// k4FWCore
#include "k4Interface/ICellPositionsTool.h"
#include "k4Interface/INoiseCaloCellsTool.h"
#include "k4Interface/INoiseConstTool.h"
#include "k4FWCore/k4_check.h"
#include "edm4hep/EventHeaderCollection.h"
#include "k4FWCore/DataHandle.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "CLHEP/Random/RandGauss.h"

#include "ReadNoiseFromFileTool.h"

class IGeoSvc;

// Root
class TH1F;

/** @class NoiseCaloCellsVsThetaFromFileTool
 *
 *  Tool for calorimeter noise
 *  Access noise constants from TH1F histogram (noise vs. theta)
 *  createRandomCellNoise: Create random CaloHits (gaussian distribution) for the vector of cells
 *  filterCellNoise: remove cells with energy below threshold*sigma from the vector of cells
 * The tool needs a cell positioning tool to translate cellID to cell theta.
 * In alternative, the tool could be rewritten to use a specific segmentation class for the cellID->theta
 * translation, but it would be coupled to a specific readout.
 * To avoid all these calls to the positioning tool, one could either
 * - save directly the noise histograms as histos of noise vs thetaID
 * - or, keep histos of noise vs theta, but change the interfaces and the tool to accept
 *   cells rather than cellIDs as input. One would then get theta from the cells.
 *
 *  @author Giovanni Marchiori
 *  @date   2024-07
 *
 */

class NoiseCaloCellsVsThetaFromFileTool : public extends<ReadNoiseFromFileTool, INoiseCaloCellsTool> {
public:
  NoiseCaloCellsVsThetaFromFileTool (const std::string& type,
                                     const std::string& name,
                                     const IInterface* parent);

  virtual ~NoiseCaloCellsVsThetaFromFileTool() = default;
  virtual StatusCode initialize() override final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::unordered_map<uint64_t, double>& aCells) const override final;

  /** @brief Create random CaloHits (gaussian distribution) for the vector of cells (aCells).
   * Vector of cells must contain all cells in the calorimeter with their cellIDs.
   */
  virtual void addRandomCellNoise(std::vector<std::pair<uint64_t, double> >& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::unordered_map<uint64_t, double>& aCells) const override final;

  /** @brief Remove cells with energy below threshold*sigma from the vector of cells
   */
  virtual void filterCellNoise(std::vector<std::pair<uint64_t, double> >& aCells)    const override final;


protected:
  virtual StatusCode initBinning (NoiseData& data,
                                  const ICaloIndexer& indexer) const override;


private:
  template <typename C>
  void addRandomCellNoiseT (C& aCells, CLHEP::RandGauss& r) const;
  template <typename C>
  void filterCellNoiseT (C& aCells) const;

  /// Handle for tool to get cell positions
  ToolHandle<ICellPositionsTool> m_cellPositionsTool{this, "cellPositionsTool", "CellPositionsDummyTool", "Handle for tool to retrieve cell positions"};


  Gaudi::Property<std::string> m_elecNoiseRMSHistoName{this, "elecNoiseRMSHistoName", "h_elecNoise_layer",
                                                       "Name of electronics noise RMS histogram"};

  /// Energy threshold (cells with Ecell < filterThreshold*m_cellNoise removed)
  Gaudi::Property<double> m_filterThreshold{
      this, "filterNoiseThreshold", 3,
      "Energy threshold (cells with Ecell < offset + filterThreshold*m_cellNoise removed)"};
  /// Change the cell filter condition to remove only cells with abs(Ecell) < offset + filterThreshold*m_cellNoise
  /// removed) This avoids to keep only 'one side'  of the noise fluctuations and prevents biasing cluster energy
  /// towards higher energies
  Gaudi::Property<bool> m_useAbsInFilter{
      this, "useAbsInFilter", false,
      "If set, cell filtering condition becomes: drop cell if abs(Ecell-offset) < filterThreshold*m_cellNoise"};

  /// Random Number Service
  SmartIF<IRndmGenSvc> m_randSvc;
  /// Gaussian random number generator used for the generation of random noise hits
  Rndm::Numbers m_gauss;

  mutable k4FWCore::DataHandle<edm4hep::EventHeaderCollection> m_header{"EventHeader", Gaudi::DataHandle::Reader, this};

  void initEvent (CLHEP::Ranlux64Engine& e) const;
};

#endif /* RECFCCEECALORIMETER_NOISECALOCELLSVSTHETAFROMFILETOOL_H */
