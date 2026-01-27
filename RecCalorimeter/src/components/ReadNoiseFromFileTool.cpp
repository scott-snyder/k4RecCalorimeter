#include "ReadNoiseFromFileTool.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

// k4geo
#include "detectorCommon/DetUtils_k4geo.h"

// k4FWCore
#include "k4Interface/IGeoSvc.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Readout.h"
#include "DDSegmentation/Segmentation.h"

// ROOT
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TSystem.h"

DECLARE_COMPONENT(ReadNoiseFromFileTool)

StatusCode ReadNoiseFromFileTool::initialize() {
  K4RECCALORIMETER_CHECK( AlgTool::initialize() );

  K4RECCALORIMETER_CHECK( m_geoSvc.retrieve() );
  K4RECCALORIMETER_CHECK( m_indexerSvc.retrieve() );
  K4RECCALORIMETER_CHECK( m_constantsSvc.retrieve() );

  // Take readout bitfield decoder from GeoSvc
  dd4hep::Readout readout = m_geoSvc->getDetector()->readout(m_readoutName);
  m_decoder = readout.idSpec().decoder();

  m_index_activeField = m_decoder->index(m_activeFieldName);

  int detID = readout.segmentation().detector()->id;
  m_indexer = m_indexerSvc->indexer (detID);
  K4RECCALORIMETER_CHECK( m_indexer != nullptr );

  // get noise constants.
  m_data = m_constantsSvc->getObj<NoiseData> (m_noiseFileName);
  if (!m_data) {
    NoiseData noise;
    K4RECCALORIMETER_CHECK( ReadNoiseFromFileTool::initNoiseFromFile(noise) );
    K4RECCALORIMETER_CHECK( m_constantsSvc->putObj (m_noiseFileName, std::move (noise)) );
    m_data = m_constantsSvc->getObj<NoiseData> (m_noiseFileName);
    K4RECCALORIMETER_CHECK( m_data != nullptr );
  }

  return StatusCode::SUCCESS;
}

StatusCode ReadNoiseFromFileTool::initNoiseFromFile(NoiseData& data) const
{
  // Check if file exists
  if (m_noiseFileName.empty()) {
    error() << "Name of the file with the noise values not provided!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (gSystem->AccessPathName(m_noiseFileName.value().c_str())) {
    error() << "Provided file with the noise values not found!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  }
  std::unique_ptr<TFile> noiseFile(TFile::Open(m_noiseFileName.value().c_str(), "READ"));
  if (noiseFile->IsZombie()) {
    error() << "Unable to open the file with the noise values!" << endmsg;
    error() << "File path: " << m_noiseFileName.value() << endmsg;
    return StatusCode::FAILURE;
  } else {
    info() << "Using the following file with the noise values: " << m_noiseFileName.value() << endmsg;
  }

  std::string elecNoiseLayerHistoName, pileupLayerHistoName;
  std::string elecNoiseOffsetLayerHistoName, pileupOffsetLayerHistoName;
  // Read the histograms with electronics noise and pileup from the file
  for (unsigned i = 0; i < m_numRadialLayers; i++) {
    elecNoiseLayerHistoName = m_elecNoiseHistoName + std::to_string(i + 1);
    debug() << "Getting histogram with a name " << elecNoiseLayerHistoName << endmsg;
    data.m_histoElecNoiseRMS.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(elecNoiseLayerHistoName.c_str())));
    if (data.m_histoElecNoiseRMS.at(i).GetNbinsX() < 1) {
      error() << "Histogram  " << elecNoiseLayerHistoName
              << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (m_setNoiseOffset) {
      elecNoiseOffsetLayerHistoName = m_elecNoiseOffsetHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << elecNoiseOffsetLayerHistoName << endmsg;
      data.m_histoElecNoiseOffset.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(elecNoiseOffsetLayerHistoName.c_str())));
      if (data.m_histoElecNoiseOffset.at(i).GetNbinsX() < 1) {
        error() << "Histogram  " << elecNoiseOffsetLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
    }
    if (m_addPileup) {
      pileupLayerHistoName = m_pileupHistoName + std::to_string(i + 1);
      debug() << "Getting histogram with a name " << pileupLayerHistoName << endmsg;
      data.m_histoPileupNoiseRMS.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(pileupLayerHistoName.c_str())));
      if (data.m_histoPileupNoiseRMS.at(i).GetNbinsX() < 1) {
        error() << "Histogram  " << pileupLayerHistoName
                << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
        return StatusCode::FAILURE;
      }
      if (m_setNoiseOffset == true) {
        pileupOffsetLayerHistoName = m_pileupOffsetHistoName + std::to_string(i + 1);
        debug() << "Getting histogram with a name " << pileupOffsetLayerHistoName << endmsg;
        data.m_histoPileupOffset.push_back(*dynamic_cast<TH1F*>(noiseFile->Get(pileupOffsetLayerHistoName.c_str())));
        if (data.m_histoElecNoiseOffset.at(i).GetNbinsX() < 1) {
          error() << "Histogram  " << pileupOffsetLayerHistoName
                  << " has 0 bins! check the file with noise and the name of the histogram!" << endmsg;
          return StatusCode::FAILURE;
        }
      }
    }
  }
  // Check if we have same number of histograms (all layers) for pileup and electronics noise
  if (data.m_histoElecNoiseRMS.size() == 0) {
    error() << "No histograms with noise found!!!!" << endmsg;
    return StatusCode::FAILURE;
  }
  if (m_addPileup) {
    if (data.m_histoElecNoiseRMS.size() != data.m_histoPileupNoiseRMS.size()) {
      error() << "Missing histograms! Different number of histograms for electronics noise and pileup!!!!" << endmsg;
      return StatusCode::FAILURE;
    }
  }

  K4RECCALORIMETER_CHECK( initBinning (data, *m_indexer));

  return StatusCode::SUCCESS;
}


StatusCode  ReadNoiseFromFileTool::initBinning (NoiseData& data,
                                                const k4::recCalo::ICaloIndexer& indexer) const
{
  /// PhiEta segmentation
  const auto* segmentation = 
    dynamic_cast<const dd4hep::DDSegmentation::FCCSWGridPhiEta_k4geo*>(
      m_geoSvc->getDetector()->readout(m_readoutName).segmentation().segmentation());
  if (segmentation == nullptr) {
    error() << "There is no phi-eta segmentation!!!!" << endmsg;
    return StatusCode::FAILURE;
  }

  data.m_bins.resize (indexer.cellIDs().size());

  auto deltaEta = [] (TH1F* h) -> std::pair<int, double>
  {
    if (!h) return std::make_pair (0, 0.);
    int Nbins = h->GetNbinsX();
    double delta = (h->GetBinLowEdge(Nbins) + h->GetBinWidth(Nbins) - h->GetBinLowEdge(1)) /
        Nbins;
    return std::make_pair (Nbins, delta);
  };
  auto [NbinsRMS, deltaEtaRMS] = deltaEta (&data.m_histoElecNoiseRMS.at(0));
  auto [NbinsOffset, deltaEtaOffset] = deltaEta (m_setNoiseOffset ? &data.m_histoElecNoiseOffset.at(0) : nullptr);

  for (uint64_t id : indexer.cellIDs()) {
    unsigned ndx = indexer.index (id);
    double cellEta = segmentation->eta(id);

    {
      int ibin = floor(fabs(cellEta) / deltaEtaRMS) + 1;
      if (ibin > NbinsRMS) {
        error() << "eta outside range of the RMS histograms! Cell eta: " << cellEta << " Nbins in histogram: " << NbinsRMS
                << endmsg;
        ibin = NbinsRMS;
      }
      data.m_bins.at(ndx).first = ibin;
    }
    {
      int ibin = floor(fabs(cellEta) / deltaEtaOffset) + 1;
      if (ibin > NbinsOffset) {
        error() << "eta outside range of the offset histograms! Cell eta: " << cellEta << " Nbins in histogram: " << NbinsOffset
                << endmsg;
        ibin = NbinsOffset;
      }
      data.m_bins.at(ndx).second = ibin;
    }

  }

  return StatusCode::SUCCESS;
}


double ReadNoiseFromFileTool::getNoiseRMSPerCell(uint64_t aCellId) const {
  // Get cell coordinates: bin and radial layer
  unsigned ndx = m_indexer->index (aCellId);
  int ibin = m_data->m_bins.at(ndx).first;

  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);
  return getNoiseRMSPerCell(ibin, cellLayer);
}

double ReadNoiseFromFileTool::getNoiseRMSPerCell(int ibin, unsigned cellLayer) const {
  double elecNoiseRMS = 0.;
  double pileupNoiseRMS = 0.;

  // All histograms have same binning, all bins with same size
  if (m_data->m_histoElecNoiseRMS.size() != 0) {
    // Check that there are not more layers than the constants are provided for
    if (cellLayer < m_data->m_histoElecNoiseRMS.size()) {
      elecNoiseRMS = m_data->m_histoElecNoiseRMS.at(cellLayer).GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoiseRMS = m_data->m_histoPileupNoiseRMS.at(cellLayer).GetBinContent(ibin);
      }
    } else {
      error()
          << "More radial layers than we have noise for!!!! Using the last layer for all histograms outside the range."
          << endmsg;
    }
  } else {
    error() << "No histograms with noise constants!!!!! " << endmsg;
  }

  // Total noise: electronics noise + pileup
  double totalNoiseRMS = sqrt(elecNoiseRMS * elecNoiseRMS + pileupNoiseRMS * pileupNoiseRMS) * m_scaleFactor;

  if (totalNoiseRMS < 1e-6) {
    warning() << "Zero noise: cell bin " << ibin << " layer " << cellLayer << " noise " << totalNoiseRMS << endmsg;
  }

  return totalNoiseRMS;
}

double ReadNoiseFromFileTool::getNoiseOffsetPerCell(uint64_t aCellId) const {

  if (!m_setNoiseOffset)
    return 0.;

  // Get cell coordinates: bin and radial layer
  unsigned ndx = m_indexer->index (aCellId);
  int ibin = m_data->m_bins.at(ndx).second;

  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);
  return getNoiseOffsetPerCell(ibin, cellLayer);
}

double ReadNoiseFromFileTool::getNoiseOffsetPerCell(int ibin, unsigned cellLayer) const {
  double elecNoiseOffset = 0.;
  double pileupNoiseOffset = 0.;

  // All histograms have same binning, all bins with same size
  if (m_data->m_histoElecNoiseOffset.size() != 0) {
    // Check that there are not more layers than the constants are provided for
    if (cellLayer < m_data->m_histoElecNoiseOffset.size()) {
      elecNoiseOffset = m_data->m_histoElecNoiseOffset.at(cellLayer).GetBinContent(ibin);
      if (m_addPileup) {
        pileupNoiseOffset = m_data->m_histoPileupOffset.at(cellLayer).GetBinContent(ibin);
      }
    } else {
      error()
          << "More radial layers than we have noise for!!!! Using the last layer for all histograms outside the range."
          << endmsg;
    }
  } else {
    error() << "No histograms with noise offset!!!!! " << endmsg;
  }

  // Total noise: electronics noise + pileup
  double totalNoiseOffset = sqrt(elecNoiseOffset * elecNoiseOffset + pileupNoiseOffset * pileupNoiseOffset) *
                            m_scaleFactor; // shouldnt the offset be summed linearly?
  // No warning is printed if offset is zero because that is the usual scenario

  return totalNoiseOffset;
}


std::pair<double, double>
ReadNoiseFromFileTool::getNoisePerCell(uint64_t aCellId) const
{
  // Get cell coordinates: bin and radial layer
  unsigned ndx = m_indexer->index (aCellId);
  const std::pair<unsigned, unsigned>& bins = m_data->m_bins.at(ndx);

  unsigned cellLayer = m_decoder->get(aCellId, m_index_activeField);

  double rms = getNoiseRMSPerCell(bins.first, cellLayer);
  double offset = m_setNoiseOffset ? getNoiseOffsetPerCell(bins.second, cellLayer) : 0;
  return std::make_pair (rms, offset);
}

