#include "TopoCaloNoisyCells.h"
#include "RecCaloCommon/k4RecCalorimeter_check.h"

#include "TBranch.h"
#include "TFile.h"
#include "TSystem.h"
#include "TTree.h"

DECLARE_COMPONENT(TopoCaloNoisyCells)

StatusCode TopoCaloNoisyCells::initialize() {
  K4RECCALORIMETER_CHECK( AlgTool::initialize() );
  K4RECCALORIMETER_CHECK( m_constantsSvc.retrieve() );
  K4RECCALORIMETER_CHECK( m_indexerSvc.retrieve() );

  m_indexer = m_indexerSvc->indexer (m_detID);
  K4RECCALORIMETER_CHECK( m_indexer != nullptr );

  m_data = m_constantsSvc->getObj<NoiseData> (m_fileName);
  if (!m_data) {
    // Check if file exists
    if (m_fileName.empty()) {
      error() << "Name of the file with the noisy cells not provided!" << endmsg;
      return StatusCode::FAILURE;
    }
    if (gSystem->AccessPathName(m_fileName.value().c_str())) {
      error() << "Provided file with the noisy cells not found!" << endmsg;
      error() << "File path: " << m_fileName.value() << endmsg;
      return StatusCode::FAILURE;
    }
    std::unique_ptr<TFile> inFile(TFile::Open(m_fileName.value().c_str(), "READ"));
    if (inFile->IsZombie()) {
      error() << "Unable to open the file with the noisy cells!" << endmsg;
      error() << "File path: " << m_fileName.value() << endmsg;
      return StatusCode::FAILURE;
    } else {
      info() << "Using the following file with the noisy cells: " << m_fileName.value() << endmsg;
    }

    NoiseData data = readData (*inFile);
    K4RECCALORIMETER_CHECK( m_constantsSvc->putObj (m_fileName, std::move (data)) );
    m_data = m_constantsSvc->getObj<NoiseData> (m_fileName);
    K4RECCALORIMETER_CHECK( m_data != nullptr );
  }

  return StatusCode::SUCCESS;
}


auto TopoCaloNoisyCells::readData (TFile& inFile) const -> NoiseData
{
  NoiseData data;

  TTree* tree = nullptr;
  inFile.GetObject("noisyCells", tree);
  ULong64_t readCellId;
  double readNoisyCells;
  double readNoisyCellsOffset;
  tree->SetBranchAddress("cellId", &readCellId);
  tree->SetBranchAddress("noiseLevel",
                         &readNoisyCells); // would be better to call branch noiseRMS rather than noiseLevel
  tree->SetBranchAddress("noiseOffset", &readNoisyCellsOffset);

  data.resize (m_indexer->cellIDs().size());

  for (uint i = 0; i < tree->GetEntries(); i++) {
    tree->GetEntry(i);
    unsigned ndx = m_indexer->index (readCellId);
    data.at(ndx) = std::make_pair (readNoisyCells, readNoisyCellsOffset);
  }
  delete tree;
  inFile.Close();

  return data;
}


double TopoCaloNoisyCells::getNoiseRMSPerCell(uint64_t aCellId) const
{
  unsigned ndx = m_indexer->index (aCellId);
  return m_data->at(ndx).first;
}


double TopoCaloNoisyCells::getNoiseOffsetPerCell(uint64_t aCellId) const
{
  unsigned ndx = m_indexer->index (aCellId);
  return m_data->at(ndx).second;
}


std::pair<double, double>
TopoCaloNoisyCells::getNoisePerCell(uint64_t aCellId) const
{
  unsigned ndx = m_indexer->index (aCellId);
  return m_data->at(ndx);
}

