#ifndef RECCALORIMETER_TUBELAYERPHIETACALOTOOL_H
#define RECCALORIMETER_TUBELAYERPHIETACALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class TubeLayerPhiEtaCaloTool Reconstruction/RecCalorimeter/src/components/TubeLayerPhiEtaCaloTool.h
 *TubeLayerPhiEtaCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  It assumes cylindrical geometry (layers) and phi-eta segmentation.
 * 
 *  Prepare a map of all existing cells in current geometry.
 *   Active layers (cylindrical tubes) are looked in the geometry manager by name ('\b activeVolumeName').
 *   Corresponding bitfield name is given in '\b activeFieldName'.
 *   If users wants to limit the number of active layers, it is possible by setting '\b activeVolumesNumber'.
 *   The total number of cells N = n_layer * n_eta * n_phi, where
 *   n_layer is number of layers (taken from geometry if activeVolumesNumber not set),
 *   n_eta is number of eta bins in that layer,
 *   n_phi is number of phi bins (the same for each layer).
 *   For more explanation please [see reconstruction documentation](@ref md_reconstruction_doc_reccalorimeter).
 *
 *  @author Anna Zaborowska
 */

class TubeLayerPhiEtaCaloTool : public CalorimeterToolBase
{
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~TubeLayerPhiEtaCaloTool() = default;

  virtual StatusCode initialize() override final;


protected:
  virtual StatusCode collectCells(std::function<void(uint64_t)> cellFunc) const override final;

private:
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", ""};
  /// Name of active volumes
  Gaudi::Property<std::string> m_activeVolumeName{this, "activeVolumeName", "LAr_sensitive"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::string> m_activeFieldName{this, "activeFieldName", "active_layer"};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{this, "fieldNames"};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValues{this, "fieldValues"};
  /// Temporary: for use with MergeLayer tool
  Gaudi::Property<unsigned int> m_activeVolumesNumber{this, "activeVolumesNumber", 0};
};

#endif /* RECCALORIMETER_TUBELAYERPHIETACALOTOOL_H */
