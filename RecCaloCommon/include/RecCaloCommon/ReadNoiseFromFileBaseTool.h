// This file's extension implies that it's C, but it's really -*- C++ -*-.
/*
 * Copyright (C) 2002-2026 CERN for the benefit of the ATLAS collaboration.
 */
/**
 * @file RecCaloCommon/ReadNoiseFromFileBaseTool.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Jul, 2026
 * @brief 
 */


#ifndef RECCALOCOMMON_READNOISEFROMFILEBASETOOL_H
#define RECCALOCOMMON_READNOISEFROMFILEBASETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolHandle.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Interfaces
#include "RecCaloCommon/INoiseConstTool.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"

// k4geo
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"

#include <memory>

class IGeoSvc;

// Root
class TH1;

/** @class ReadNoiseFromFileTool
 *
 *  Tool to read the stored noise constant per cell in the calorimeters
 *  Access noise constants from TH1F histogram (noise vs. |eta|)
 *
 *  @author Jana Faltova, Coralie Neubueser
 *  @date   2018-01
 *
 */

class ReadNoiseFromFileBaseTool : public extends<AlgTool, k4::recCalo::INoiseConstTool> {
public:
  using base_class::base_class;
  virtual ~ReadNoiseFromFileBaseTool() = default;

  virtual StatusCode initialize() override;

  /// Find the appropriate noise constant from the histogram
  virtual double getNoiseRMSPerCell(CellID aCellID) const override final;
  virtual double getNoiseOffsetPerCell(CellID aCellID) const override final;
  virtual std::pair<double, double> getNoisePerCell(CellID aCellID) const override final;

protected:
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalHitsPhiEta", "Name of the detector readout"};
  /// Noise offset, if false, mean is set to 0
  Gaudi::Property<bool> m_setNoiseOffset{this, "setNoiseOffset", true, "Set a noise offset per cell"};
  /// Name of electronics noise histogram
  Gaudi::Property<std::string> m_elecNoiseHistoName{this, "elecNoiseHistoName", "h_elecNoise_layer",
                                                    "Name of electronics noise histogram"};
  /// Handle to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};

  struct NoiseData {
    NoiseData() = default;
    NoiseData(const NoiseData&);
    NoiseData(NoiseData&&) = default;
    NoiseData& operator=(const NoiseData&) = delete;
    NoiseData& operator=(NoiseData&&) = default;
    /// Histograms with pileup RMS (index in array - radial layer)
    std::vector<std::unique_ptr<TH1> > m_histoPileupNoiseRMS;
    /// Histograms with electronics noise RMS (index in array - radial layer)
    std::vector<std::unique_ptr<TH1> > m_histoElecNoiseRMS;

    /// Histograms with pileup offset (index in array - radial layer)
    std::vector<std::unique_ptr<TH1> > m_histoPileupOffset;
    /// Histograms with electronics noise offset (index in array - radial layer)
    std::vector<std::unique_ptr<TH1> > m_histoElecNoiseOffset;

    // rms, offset
    std::vector<std::pair<unsigned, unsigned> > m_bins;
  };

  virtual StatusCode initBinning (NoiseData& data,
                                  const k4::recCalo::ICaloIndexer& indexer) const = 0;

  // Decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  int m_index_activeField = -1;


private:
  /// Open file and read noise histograms in the memory
  StatusCode initNoiseFromFile(NoiseData& data) const;

  double getNoiseRMSPerCell(int ibin, unsigned cellLayer) const;
  double getNoiseOffsetPerCell(int ibin, unsigned cellLayer) const;

  /// Add pileup contribution to the electronics noise? (only if read from file)
  Gaudi::Property<bool> m_addPileup{this, "addPileup", true,
                                    "Add pileup contribution to the electronics noise? (only if read from file)"};

  /// Name of the file with noise constants
  Gaudi::Property<std::string> m_noiseFileName{this, "noiseFileName", "", "Name of the file with noise constants"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer",
                                                 "Name of active layers for sampling calorimeter"};
  /// Name of pileup histogram
  Gaudi::Property<std::string> m_pileupHistoName{this, "pileupHistoName", "h_pileup_layer", "Name of pileup histogram"};
  /// Name of electronics noise offset histogram
  Gaudi::Property<std::string> m_elecNoiseOffsetHistoName{this, "elecNoiseOffsetHistoName", "h_mean_pileup_layer",
                                                          "Name of electronics noise offset histogram"};
  /// Name of pileup offset histogram
  Gaudi::Property<std::string> m_pileupOffsetHistoName{this, "pileupOffsetHistoName", "h_pileup_layer",
                                                       "Name of pileup offset histogram"};

  /// For config compatibility.
  Gaudi::Property<uint> m_numRadialLayers{this, "numRadialLayers", 0, "Number of radial layers.  Deprecated: use numHistograms instead."};
  /// Number of radial layers/wheels
  Gaudi::Property<uint> m_numHistograms{this, "numHistograms", 3, "Number of histograms"};
 
  /// Factor to apply to the noise values to get them in GeV if e.g. they were produced in MeV
  Gaudi::Property<float> m_scaleFactor{this, "scaleFactor", 1, "Factor to apply to the noise values"};

  const NoiseData* m_data = nullptr;

  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc
  { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
  { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "" };
  const k4::recCalo::ICaloIndexer* m_indexer = nullptr;
};


#endif // not RECCALOCOMMON_READNOISEFROMFILEBASETOOL_H
