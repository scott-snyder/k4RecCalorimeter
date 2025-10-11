#include "ReadCaloCrosstalkMap.h"
#include "k4FWCore/k4_check.h"
#include "k4Interface/IGeoSvc.h"
#include "DD4hep/Detector.h"

#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

DECLARE_COMPONENT(ReadCaloCrosstalkMap)

StatusCode ReadCaloCrosstalkMap::initialize() {
  // prevent to initialize the tool if not intended (input file path empty)
  // otherwise things will crash if m_fileName is not available
  // not a perfect solution but tools seems to not be meant to be optional
  if (m_fileName == "") {
    debug() << "Empty 'fileName' provided, it means cross-talk map is not needed, exiting ReadCaloCrosstalkMap "
               "initialization"
            << endmsg;
    return StatusCode::SUCCESS;
  }

  info() << "Loading crosstalk map..." << endmsg;

  K4_GAUDI_CHECK( AlgTool::initialize() );
  K4_GAUDI_CHECK( m_constantsSvc.retrieve() );
  K4_GUADI_CHECK( m_indexerSvc.retrieve() );

  int detID = m_detID;
  if (detID < 0) {
    // If defaulted, get the ECAL_Barrel ID from the geometry service.
    ServiceHandle<IGeoSvc> geoSvc ("GeoSvc", name());
    K4RECCALORIMETER_CHECK( geoSvc.retrieve() );
    detID = geoSvc->getDetector()->constantAsDouble("DetID_ECAL_Barrel");
  }
  m_indexer = m_indexerSvc->indexer (detID);
  K4RECCALORIMETER_CHECK( m_indexer != nullptr );

  m_data = m_constantsSvc->getObj<CrosstalkData> (m_fileName);
  if (!m_data) {
    // Check if crosstalk file exists
    if (gSystem->AccessPathName(m_fileName.value().c_str())) {
      error() << "Provided file with the crosstalk map not found!" << endmsg;
      error() << "File path: " << m_fileName.value() << endmsg;
      return StatusCode::FAILURE;
    }
    std::unique_ptr<TFile> xtalkFile (TFile::Open(m_fileName.value().c_str(), "READ"));
    if (xtalkFile->IsZombie()) {
      error() << "Unable to read the file with the crosstalk map!" << endmsg;
      error() << "File path: " << m_fileName.value() << endmsg;
      return StatusCode::FAILURE;
    } else {
      info() << "Using the following file with the crosstalk map: " << m_fileName.value() << endmsg;
    }

    CrosstalkData data = readData (*xtalkFile);
    K4RECCALORIMETER_CHECK( m_constantsSvc->putObj (m_fileName, std::move (data)) );
    m_data = m_constantsSvc->getObj<CrosstalkData> (m_fileName);
    K4RECCALORIMETER_CHECK( m_data != nullptr );
  }

  return StatusCode::SUCCESS;
}


auto ReadCaloCrosstalkMap::readData (TFile& xtalkFile) const -> CrosstalkData
{
  CrosstalkData data;

  TTree* tree = nullptr;
  xtalkFile.GetObject("crosstalk_neighbours", tree);
  ULong64_t read_cellId;
  std::vector<uint64_t>* read_neighbours = nullptr;
  std::vector<double>* read_crosstalks = nullptr;

  tree->SetBranchAddress("cellId", &read_cellId);
  tree->SetBranchAddress("list_crosstalk_neighbours", &read_neighbours);
  tree->SetBranchAddress("list_crosstalks", &read_crosstalks);

  size_t ncells = m_indexer->cellIDs().size();
  data.m_neighbourIndices.resize (ncells);
  data.m_crosstalkIndices.resize (ncells);

  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);

    unsigned ndx = m_indexer->index (read_cellId);
    {
      size_t oldsz = data.m_neighbours.size();
      data.m_neighbours.insert (data.m_neighbours.end(), read_neighbours->begin(), read_neighbours->end());
      data.m_neighbourIndices[ndx] = std::make_pair (oldsz, data.m_neighbours.size()-oldsz);
    }
    {
      size_t oldsz = data.m_crosstalks.size();
      data.m_crosstalks.insert (data.m_crosstalks.end(), read_crosstalks->begin(), read_crosstalks->end());
      data.m_crosstalkIndices[ndx] = std::make_pair (oldsz, data.m_crosstalks.size()-oldsz);
    }
  }

  info() << "Crosstalk input: " << xtalkFile.GetName() << endmsg;
  info() << "Total number of cells = " << tree->GetEntries() << endmsg;
  //<< ", Size of crosstalk neighbours = " << data.m_mapNeighbours.size()
  //<< ", Size of coefficients = " << data.m_mapCrosstalks.size() << endmsg;
  delete tree;
  delete read_neighbours;
  delete read_crosstalks;
  xtalkFile.Close();

  return data;
}


std::span<const uint64_t>
ReadCaloCrosstalkMap::getNeighbours(uint64_t aCellId) const {
  return m_data->getNeighbours (m_indexer->index (aCellId));
}

std::span<const double>
ReadCaloCrosstalkMap::getCrosstalks(uint64_t aCellId) const {
  return m_data->getCrosstalks (m_indexer->index (aCellId));
}
