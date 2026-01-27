#ifndef RECFCCEECALORIMETER_TUBELAYERMODULETHETACALOTOOL_H
#define RECFCCEECALORIMETER_TUBELAYERMODULETHETACALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"
#include "RecCaloCommon/IDMap.h"
#include <optional>

/** @class TubeLayerModuleThetaCaloTool
 * k4RecCalorimeter/RecFCCeeCalorimeter/src/components/TubeLayerModuleThetaCaloTool.h TubeLayerModuleThetaCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  It assumes cylindrical geometry (layers) and phi-theta segmentation.
 *
 *  @author Anna Zaborowska
 *  @author Zhibo Wu
 */

class TubeLayerModuleThetaCaloTool : public CalorimeterToolBase {
public:
  TubeLayerModuleThetaCaloTool (const std::string& type,
                                const std::string& name,
                                const IInterface* parent);

  //using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~TubeLayerModuleThetaCaloTool() = default;

  virtual StatusCode initialize() override;

  virtual unsigned index(uint64_t cellID) const;


protected:
  /** Fill vector with all existing cells for this geometry.
   */
  virtual StatusCode collectCells(std::vector<uint64_t>& cells) const final override;


private:
  /// Name of active volumes
  Gaudi::Property<std::string> m_activeVolumeName{this, "activeVolumeName", "LAr_sensitive"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer"};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{this, "fieldNames"};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValues{this, "fieldValues"};
  /// Number of layers
  Gaudi::Property<unsigned int> m_activeVolumesNumber{this, "activeVolumesNumber", 0};

  using IDMap_t = k4::recCalo::IDMapN<unsigned, 3>;
  std::optional<IDMap_t> m_idmap;
};

#endif /* RECFCCEECALORIMETER_TUBELAYERMODULETHETACALOTOOL_H */
