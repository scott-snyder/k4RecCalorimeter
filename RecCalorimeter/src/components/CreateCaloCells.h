#ifndef RECCALORIMETER_CREATECALOCELLS_H
#define RECCALORIMETER_CREATECALOCELLS_H

// k4FWCore
#include "k4FWCore/DataHandle.h"
#include "k4FWCore/MetaDataHandle.h"
#include "k4Interface/ICalibrateCaloHitsTool.h"
#include "k4Interface/ICaloReadCrosstalkMap.h"
#include "k4Interface/ICalorimeterTool.h"
#include "k4Interface/INoiseCaloCellsTool.h"
#include "k4Interface/ICellPositionsTool.h"

// Gaudi
#include "Gaudi/Algorithm.h"
#include "GaudiKernel/ToolHandle.h"

// edm4hep
#include "edm4hep/CalorimeterHitCollection.h"
#include "edm4hep/CaloHitSimCaloHitLinkCollection.h"
#include "edm4hep/Constants.h"
#include "edm4hep/SimCalorimeterHitCollection.h"

// DD4hep
#include "DD4hep/Detector.h"
#include "DD4hep/Volumes.h"
#include "TGeoManager.h"

#include <variant>

class IGeoSvc;

/** @class CreateCaloCells
 *
 *  Algorithm for creating calorimeter cells from Geant4 hits.
 *  Tube geometry with PhiEta segmentation expected.
 *
 *  Flow of the program:
 *  1/ Merge Geant4 energy deposits with same cellID
 *  2/ Emulate cross-talk (if switched on)
 *  3/ Calibrate to electromagnetic scale (if calibration switched on)
 *  4/ Add random noise to each cell (if noise switched on)
 *  5/ Filter cells and remove those with energy below threshold (if noise +
 * filtering switched on)
 *
 *  Tools called:
 *    - CalibrateCaloHitsTool
 *    - NoiseCaloCellsTool
 *    - CaloReadCrosstalkMap
 *    - CalorimeterTool (for geometry)
 *
 *  @author Jana Faltova
 *  @author Anna Zaborowska
 *  @date   2016-09
 *
 */

class CreateCaloCells : public Gaudi::Algorithm {

public:
  CreateCaloCells(const std::string& name, ISvcLocator* svcLoc);

  virtual StatusCode initialize() override;

  virtual StatusCode execute(const EventContext&) const override;
  
  virtual StatusCode finalize() override;

private:
  static constexpr size_t INVALID = static_cast<size_t> (-1);
  using CellsIndexMap_t = std::unordered_map<uint64_t, size_t>;
  using CellsIndexPair_t = std::pair<size_t, size_t>;  // cell index, hit index

  struct CellsFullIndex
  {
    CellsFullIndex (const CellsIndexMap_t& cellsIndexMap)
      : m_cellsIndexMap (cellsIndexMap),
        m_indices (cellsIndexMap.size(), {INVALID, INVALID})
    {
    }

    CellsIndexPair_t& index (uint64_t cellid)
    {
      auto it = m_cellsIndexMap.find (cellid);
      if (it == m_cellsIndexMap.end()) {
        throw std::out_of_range ("bad cellid");
      }
      return m_indices.at (it->second);
    }

    const CellsIndexMap_t& m_cellsIndexMap;
    std::vector<CellsIndexPair_t> m_indices;
  };


  struct CellsSparseIndex
  {
    CellsIndexPair_t& index (uint64_t cellid)
    {
      return m_indices.try_emplace (cellid, CellsIndexPair_t{INVALID,INVALID}).first->second;
    }

    using CellsIndexPairMap_t = std::unordered_map<uint64_t, CellsIndexPair_t>;
    CellsIndexPairMap_t m_indices;
  };

  struct CellsIndex
  {
    CellsIndex (const CellsIndexMap_t& cellsIndexMap)
    {
      if (!cellsIndexMap.empty()) {
        m_indices.emplace<1> (cellsIndexMap);
      }
      else {
        m_indices.emplace<2>();
      }
    }

    CellsIndexPair_t& pair (uint64_t cellid)
    {
      if (m_indices.index() == 1) {
        return std::get<1> (m_indices).index (cellid);
      }
      return std::get<2> (m_indices).index (cellid);
    }
    size_t& index (uint64_t cellid, size_t ihit)
    {
      CellsIndexPair_t& p = pair (cellid);
      if (p.second == INVALID) p.second = ihit;
      return p.first;
    }
    size_t& index (uint64_t cellid)
    {
      return pair (cellid).first;
    }
    size_t& ihit (uint64_t cellid)
    {
      return pair (cellid).second;
    }

    std::variant<int, CellsFullIndex, CellsSparseIndex> m_indices;
  };

  struct CellsInfo
  {
    CellsInfo (size_t capacity)
    {
      m_cells.reserve (capacity);
    }

    size_t size() const
    {
      return m_cells.size();
    }

    size_t add (uint64_t cellID, double energy)
    {
      m_cells.emplace_back (cellID, energy);
      return m_cells.size() - 1;
    }

    uint64_t cellID (size_t icell) const
    {
      return m_cells.at(icell).first;
    }
      
    double& energy (size_t icell)
    {
      return m_cells.at(icell).second;
    }

    std::vector<std::pair<uint64_t, double> > m_cells;
  };

  /// Handle for the calorimeter cells crosstalk tool
  ToolHandle<ICaloReadCrosstalkMap> m_crosstalksTool
  {this, "crosstalksTool", "ReadCaloCrosstalkMap", "Handle for the cell crosstalk tool"};

  /// Handle for tool to calibrate Geant4 energy to EM scale tool
  mutable ToolHandle<ICalibrateCaloHitsTool> m_calibTool{"CalibrateCaloHitsTool", this};
  /// Handle for the calorimeter cells noise tool
  mutable ToolHandle<INoiseCaloCellsTool> m_noiseTool{"NoiseCaloCellsFlatTool", this};
  /// Handle for the geometry tool
  ToolHandle<ICalorimeterTool> m_geoTool{"TubeLayerPhiEtaCaloTool", this};
  ToolHandle<ICellPositionsTool> m_cellPos
    { this, "CellPositionsTool", "", "Cell positions tool.  If defaulted, position based on volume only." };

  /// Add crosstalk to cells?
  Gaudi::Property<bool> m_addCrosstalk{this, "addCrosstalk", false, "Add crosstalk effect?"};
  /// Calibrate to EM scale?
  Gaudi::Property<bool> m_doCellCalibration{this, "doCellCalibration", true, "Calibrate to EM scale?"};
  /// Add noise to cells?
  Gaudi::Property<bool> m_addCellNoise{this, "addCellNoise", true, "Add noise to cells?"};
  /// Save only cells with energy above threshold?
  Gaudi::Property<bool> m_filterCellNoise{this, "filterCellNoise", false,
                                          "Save only cells with energy above threshold?"};
  // Add position information to the cells? (based on Volumes, not cells, could be improved)
  Gaudi::Property<bool> m_addPosition{this, "addPosition", false, "Add position information to the cells?"};

  /// Handle for calo hits (input collection)
  mutable DataHandle<edm4hep::SimCalorimeterHitCollection> m_hits{"hits", Gaudi::DataHandle::Reader, this};
  /// Handle for the cellID encoding string of the input collection
  MetaDataHandle<std::string> m_hitsCellIDEncoding{m_hits, edm4hep::labels::CellIDEncoding, Gaudi::DataHandle::Reader};
  /// Handle for calo cells (output collection)
  mutable DataHandle<edm4hep::CalorimeterHitCollection> m_cells{"cells", Gaudi::DataHandle::Writer, this};
  MetaDataHandle<std::string> m_cellsCellIDEncoding{m_cells, edm4hep::labels::CellIDEncoding,
                                                    Gaudi::DataHandle::Writer};
  /// Handle for hit<->cell link (output collection)
  mutable DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection> m_links{"", Gaudi::DataHandle::Writer, this};
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalBarrelPhiEta", "Name of the detector readout"};
  /// Name of active volumes
  Gaudi::Property<std::string> m_activeVolumeName{this, "activeVolumeName", "_sensitive", "Name of the active volumes"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer",
                                                 "Name of active layers for sampling calorimeter"};
  /// Name of the bit-fields (in the readout) describing the volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{
      this, "fieldNames", {}, "Name of the bit-fields (in the readout) describing the volume"};
  /// Values of the fields that identify the volume to change segmentation (e.g.
  /// ID of the ECal)
  Gaudi::Property<std::vector<int>> m_fieldValues{this,
                                                  "fieldValues",
                                                  {},
                                                  "Value of the field that identifies the volume "
                                                  "to to change segmentation (e.g. ID of the "
                                                  "ECal)"};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;
  dd4hep::VolumeManager m_volman;
  /// Map of cell IDs to cell indices.
  /// This assigns to each cell a dense index in the range 0..ncells-1.
  CellsIndexMap_t m_cellsIndexMap;

  /// Maps of cell IDs (corresponding to DD4hep IDs) on final energies to be used for clustering
  mutable std::unordered_map<uint64_t, double> m_cellsMap;
  /// Maps of cell IDs (corresponding to DD4hep IDs) on transfer of signals due to crosstalk
  mutable std::unordered_map<uint64_t, double> m_CrosstalkCellsMap;
  /// Maps of cell IDs with zero energy, for all cells in calo (needed if addCellNoise and filterCellNoise are both set)
  mutable std::unordered_map<uint64_t, double> m_emptyCellsMap;
};

#endif /* RECCALORIMETER_CREATECALOCELLS_H */
