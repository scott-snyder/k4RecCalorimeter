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
#include "RecCaloCommon/ICaloCellIndexerSvc.h"

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


private:
  /// Index of a cell in the container.
  using index_t = k4::recCalo::ICaloIndexer::index_t;

  // Forward declaration.
  class CaloCells;


  /**
   * @brief Helper for accessing cell data.
   *
   * This provides a class-like interface for accessing data for a single cell.
   * The members are bound to the proper members of the vectors in CaloCells.
   * It is expected that members that are not actually used will be entirely
   * optimized away.
   */
  class CaloCell
  {
  public:
    /// Construct a cell interface object.
    CaloCell (uint64_t& the_cellid, double& the_energy, size_t& the_ihit)
      : cellid (the_cellid), energy (the_energy), ihit (the_ihit)
    {
    }

    uint64_t& cellid;
    double&   energy;
    size_t&   ihit;
  };


  /**
   * @brief Holder for cell data.
   *
   * This helper class holds information on cells during reconstruction.
   * For each cell, we store the cell id, the energy, and the index of the
   * original hit for this cell (if any).  A constraint is that numerous
   * interfaces want to operate on the cellid and energy data as a vector
   * of pairs, so we need to store those two values in that way (@c m_cells).
   *
   * Given that, we implement two different ways of storing the data,
   * full and sparse, depending on how many cells we expect to be processing.
   * These are identified by the @c m_mode variable.
   *
   * In full mode, all cells are expected to be present.  The @c m_cells vector
   * is initialized by writing in all possible cell ids.  To look up a cell
   * by cell id, we find the index within the cell id list and then index
   * into @c m_cells.  The hit indices are stored in @c m_ihits, which
   * is originally in 1-1 correspondence with @c m_cells.  No sorting
   * is needed in this case.  A complication, however, is that the cells
   * may be filtered, which breaks the correspondence between the cellid
   * index and the index in @c m_cells.  In that case, this is noted
   * by changing @c m_mode to @c FILTERED.  Then, looking up cells
   * by cell id is an error, and we find the ihit for a cell in @c m_cells
   * by looking up the cell's cell id index.
   *
   * In sparse mode, we expect only a small fraction of cells to be present.
   * In that case, @c m_cells is originally empty, and we add new cells
   * as needed.  An @c unordered_map (@c m_indices) maps from cell id to
   * both the cell's index in @c m_cells and the hit index for the cell.
   * At the end, the cells need to be sorted.  After that, we reset the
   * cell indices in @c m_indices (this also fixes up the effects of
   * any filtering).
   */
  class CaloCells
  {
  public:
    /// Invalid cell index.
    constexpr static index_t INVALID_ICELL = k4::recCalo::ICaloIndexer::INVALID;

    /// Invalid hit index.
    constexpr static size_t INVALID_IHIT = static_cast<size_t> (-1);


    /**
     * @brief Constructor.
     * @param full If true, we initialize for all cells.
     * @param cellids List of all possible cell ids.  Only used in full mode.
     * @param indexer Helper to convert from cell id to the indix within
     *                @c cellids.  Only used in full mode.
     */
    CaloCells (bool full,
               std::span<const uint64_t> cellids,
               const k4::recCalo::ICaloIndexer* indexer);


    /**
     * @brief Return the number of cells.
     */
    size_t size() const { return m_cells.size(); }


    /**
     * @brief Return the index within @c m_cells for a given cell id.
     * @param cellid The cell id to find.
     *
     * May return INVALID_ICELL if the cell id doesn't exist in the container.
     * Will throw if we have a filtered full container.
     */
    index_t indexByID (uint64_t cellid);


    /**
     * @brief Return a @c CaloCell wrapper for a given cell id, adding
     *        the cell if needed.
     * @param cellid The cell id to find.
     *
     * May add a new cell.
     * Will throw if we have a filtered full container.
     */
    CaloCell cellByID (uint64_t cellid);


    /**
     * @brief Return a @c CaloCell wrapper for a given either the
     *        cell id or the index within @c m_cells.
     * @param icell Index within @c m_cells, or @c INVALID_CELL.
     * @param cellid The cell id to find.
     *
     * Uses @c icell for a @c FULL container, otherwise @c cellid.
     */
    CaloCell cellByIndexOrID (index_t icell, uint64_t cellid);


    /**
     * @brief Return hit index for a given cell in @c m_cells.
     * @param icell Index of the cell in @c m_cells.
     *
     * Returns the corresponding hit index, or @c INVALID_IHIT.
     */
    size_t ihitByIndex (index_t icell) const;


    /**
     * @brief Sort the sells in order of increasing cell id.
     */
    void sort();


    /**
     * @brief Called after filtering.
     * @param orig_size The container size before filtering.
     */
    void afterFilter (size_t orig_size);


    /// How are we representing the data?
    enum {
      // All cells are present.
      FULL = 0,
      // FULL, but some cells have been removed.
      FILTERED = 1,
      // Sparse mode.
      SPARSE = 2
    } m_mode;   


    /// Helper to map from cell id to index within cell id list.
    // Only used in FULL/FILTERED modes.
    const k4::recCalo::ICaloIndexer* m_indexer = nullptr;

    /// Vector of cellid, energy pairs per cell.
    // In FULL mode, all cellids are present in sorted order.
    // If FILTERED, some may have been remove.
    // In SPARSE, only cells actually used are present.
    using CellData_t = std::vector<std::pair<uint64_t, double> >;
    CellData_t m_cells;

    /// Vector of original hit indices, indexed by cell id index.
    // Only used in FULL/FILTERED modes.
    std::vector<size_t> m_ihits;

    /// Map from cell id to icell, ihit pairs.
    // Only used in SPARSE mode.
    std::unordered_map<uint64_t, std::pair<index_t, size_t> > m_indices;
  };


  /// Build m_caloTypes, giving calorimeter type per system ID.
  void findCaloTypes();

  /// Do crosstalk processing.
  void addCrosstalk (CaloCells& cells) const;


  /// Handle for the calorimeter cells crosstalk tool
  ToolHandle<ICaloReadCrosstalkMap> m_crosstalksTool
  {this, "crosstalkTool", "ReadCaloCrosstalkMap", "Handle for the cell crosstalk tool"};

  /// Handle for tool to calibrate Geant4 energy to EM scale tool
  ToolHandle<ICalibrateCaloHitsTool> m_calibTool{"CalibrateCaloHitsTool", this};
  /// Handle for the calorimeter cells noise tool
  ToolHandle<INoiseCaloCellsTool> m_noiseTool{"NoiseCaloCellsFlatTool", this};
  /// Handle for the geometry tool
  ToolHandle<ICalorimeterTool> m_geoTool{"TubeLayerPhiEtaCaloTool", this};
  ToolHandle<ICellPositionsTool> m_cellPos
    { this, "positionsTool", "", "Cell positions tool.  If defaulted, position based on volume only." };

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
  mutable k4FWCore::DataHandle<edm4hep::SimCalorimeterHitCollection> m_hits{"hits", Gaudi::DataHandle::Reader, this};
  /// Handle for the cellID encoding string of the input collection
  k4FWCore::MetaDataHandle<std::string> m_hitsCellIDEncoding{m_hits, edm4hep::labels::CellIDEncoding,
                                                             Gaudi::DataHandle::Reader};
  /// Handle for calo cells (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CalorimeterHitCollection> m_cells{"cells", Gaudi::DataHandle::Writer, this};
  k4FWCore::MetaDataHandle<std::string> m_cellsCellIDEncoding{m_cells, edm4hep::labels::CellIDEncoding,
                                                              Gaudi::DataHandle::Writer};
  /// Handle for hit<->cell link (output collection)
  mutable k4FWCore::DataHandle<edm4hep::CaloHitSimCaloHitLinkCollection> m_links{"", Gaudi::DataHandle::Writer, this};
  /// Name of active volumes

  Gaudi::Property<float> m_discritMin {this, "discritMin", 0};
  Gaudi::Property<float> m_discritMax {this, "discritMax", 100};
  Gaudi::Property<int> m_discritN     {this, "discritN", -1};

  /// Pointer to the geometry service
  ServiceHandle<IGeoSvc> m_geoSvc;

  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc
  { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };

  /// Volume manager for our subdetector.
  dd4hep::VolumeManager m_volman;

  /// List of all cell ids for our subdetector.
  std::span<const uint64_t> m_cellIDs;

  /// Mapper from cell id to index within m_cellIDs.
  const k4::recCalo::ICaloIndexer* m_indexer = nullptr;

  /// Indexed by system ID, giving the calorimeter type word.
  /// Non-calorimeter system IDs are set to 0.
  /// We record this for all system IDs since we build this during
  /// initialization, before we know which system ID we're handling.
  std::vector<int> m_caloTypes;

  /// Cell ID decoder.
  dd4hep::DDSegmentation::BitFieldCoder m_decoder;

  /// Field indices for detector ID and layer.
  unsigned m_systemIndex = -1;
  unsigned m_layerIndex = -1;
};

#endif /* RECCALORIMETER_CREATECALOCELLS_H */
