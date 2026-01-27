#ifndef RECCALORIMETER_READNOISEFROMFILETOOL_H
#define RECCALORIMETER_READNOISEFROMFILETOOL_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"
#include "GaudiKernel/ToolHandle.h"

// k4FWCore
#include "k4FWCore/DataHandle.h"

// Interfaces
#include "RecCaloCommon/INoiseConstTool.h"
#include "RecCaloCommon/ICellPositionsTool.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"

// k4geo
#include "detectorSegmentations/FCCSWGridPhiEta_k4geo.h"

class IGeoSvc;

// Root
class TH1F;

/** @class ReadNoiseFromFileTool
 *
 *  Tool to read the stored noise constant per cell in the calorimeters
 *  Access noise constants from TH1F histogram (noise vs. |eta|)
 *
 *  @author Jana Faltova, Coralie Neubueser
 *  @date   2018-01
 *
 */

class ReadNoiseFromFileTool : public extends<AlgTool, k4::recCalo::INoiseConstTool> {
public:
  using base_class::base_class;
  virtual ~ReadNoiseFromFileTool() = default;

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

  struct NoiseData {
    /// Histograms with pileup RMS (index in array - radial layer)
    std::vector<TH1F> m_histoPileupNoiseRMS;
    /// Histograms with electronics noise RMS (index in array - radial layer)
    std::vector<TH1F> m_histoElecNoiseRMS;

    /// Histograms with pileup offset (index in array - radial layer)
    std::vector<TH1F> m_histoPileupOffset;
    /// Histograms with electronics noise offset (index in array - radial layer)
    std::vector<TH1F> m_histoElecNoiseOffset;

    // rms, offset
    std::vector<std::pair<unsigned, unsigned> > m_bins;
  };

  virtual StatusCode initBinning (NoiseData& data,
                                  const k4::recCalo::ICaloIndexer& indexer) const;


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

  /// Number of radial layers
  Gaudi::Property<uint> m_numRadialLayers{this, "numRadialLayers", 3, "Number of radial layers"};

  /// Factor to apply to the noise values to get them in GeV if e.g. they were produced in MeV
  Gaudi::Property<float> m_scaleFactor{this, "scaleFactor", 1, "Factor to apply to the noise values"};

  const NoiseData* m_data = nullptr;

  /// Handle to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc"};
  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc
  { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
  { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "" };
  // Decoder
  dd4hep::DDSegmentation::BitFieldCoder* m_decoder;
  const k4::recCalo::ICaloIndexer* m_indexer = nullptr;
  int m_index_activeField = -1;
};

#endif /* RECCALORIMETER_READNOISEFROMFILETOOL_H */
