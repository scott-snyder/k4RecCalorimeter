#ifndef RECCALORIMETER_READCALOXTALKMAP_H
#define RECCALORIMETER_READCALOXTALKMAP_H

// from Gaudi
#include "GaudiKernel/AlgTool.h"

// k4FWCore
#include "k4Interface/ICaloReadCrosstalkMap.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "RecCaloCommon/ICaloCellIndexerSvc.h"
#include <span>

class IGeoSvc;
class TFile;

/** @class ReadCaloCrosstalkMap Reconstruction/RecCalorimeter/src/components/ReadCaloCrosstalkMap.h
 *TopoCaloNeighbours.h
 *
 *  Tool that reads a ROOT file containing the TTree with branches "cellId", "list_crosstalk_neighbours" and
 *"list_crosstalks". This tools reads the tree, creates two maps, and allows a lookup of all crosstalk neighbours as
 *well as the corresponding crosstalk coefficients for a given cell.
 *
 *  @author Zhibo Wu
 */

class ReadCaloCrosstalkMap : public extends<AlgTool, ICaloReadCrosstalkMap> {
public:
  using base_class::base_class;
  virtual ~ReadCaloCrosstalkMap() = default;

  virtual StatusCode initialize() override final;

  /** Function to be called for the crosstalk neighbours of a cell.
   *   @param[in] aCellId, cellid of the cell of interest.
   *   @return vector of cellIDs, corresponding to the crosstalk neighbours.
   */
  virtual std::span<const uint64_t> getNeighbours(uint64_t aCellId) const override final;

  /** Function to be called for the crosstalk coefficients between the input cell and its neighbouring cells.
   *   @param[in] aCellId, cellid of the cell of interest.
   *   @return vector of crosstalk coefficients.
   */
  virtual std::span<const double> getCrosstalks(uint64_t aCellId) const override final;

private:
  /// Name of input root file that contains the TTree with cellID->vec<list_crosstalk_neighboursCellID> and
  /// cellId->vec<list_crosstalksCellID>
  Gaudi::Property<std::string> m_fileName{this, "fileName", "",
                                          "Name of the file that contains the crosstalk map. Leave the default empty "
                                          "to avoid crashes when cross-talk is not needed."};

  Gaudi::Property<int> m_detID
    { this, "detID", -1, "Subsystem ID for this detector." };
  ServiceHandle<k4::recCalo::ICaloCellConstantsSvc> m_constantsSvc
  { this, "CaloCellConstantsSvc", "k4::recCalo::CaloCellConstantsSvc", "" };
  ServiceHandle<k4::recCalo::ICaloCellIndexerSvc> m_indexerSvc
  { this, "CaloCellIndexerSvc", "k4::recCalo::CaloCellIndexerSvc", "" };

  /// Output maps to be used for the fast lookup in the creating calo-cells algorithm
  struct CrosstalkData {
    std::vector<std::pair<size_t, size_t> > m_neighbourIndices;
    std::vector<uint64_t> m_neighbours;
    std::vector<std::pair<size_t, size_t> > m_crosstalkIndices;
    std::vector<double> m_crosstalks;

    std::span<const uint64_t> getNeighbours (size_t cellNdx) const
    {
      const auto& indices = m_neighbourIndices[cellNdx];
      return std::span<const uint64_t> (m_neighbours.data()+indices.first,
                                        indices.second);
    }

    std::span<const double> getCrosstalks (size_t cellNdx) const
    {
      const auto& indices = m_crosstalkIndices[cellNdx];
      return std::span<const double> (m_crosstalks.data()+indices.first,
                                      indices.second);
    }
  };

  const CrosstalkData* m_data = nullptr;
  const k4::recCalo::ICaloIndexer* m_indexer = nullptr;

  CrosstalkData readData (TFile& xtalkFile) const;
};

#endif /* RECCALORIMETER_READCALOXTALKMAP_H */
