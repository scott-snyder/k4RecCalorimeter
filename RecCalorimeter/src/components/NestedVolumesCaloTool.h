#ifndef RECCALORIMETER_NESTEDVOLUMESCALOTOOL_H
#define RECCALORIMETER_NESTEDVOLUMESCALOTOOL_H

#include "RecCaloCommon/CalorimeterToolBase.h"

/** @class NestedVolumesCaloTool Reconstruction/RecCalorimeter/src/components/NestedVolumesCaloTool.h
 *NestedVolumesCaloTool.h
 *
 *  Tool for geometry-dependent settings of the digitisation.
 *  It assumes no segmentation is used. It may be used for nested volumes.
 *
 *  Prepare a collection of all existing cells in current geometry.
 *   Active volumes are looked in the geometry manager by name ('\b activeVolumeName').
 *   Corresponding bitfield name is given in '\b activeFieldName'.
 *   If more than one name is given, it is assumed that volumes are nested.
 *   For more explanation please [see reconstruction documentation](@ref md_reconstruction_doc_reccalorimeter).
 *
 *  @author Anna Zaborowska
 */
class NestedVolumesCaloTool : public CalorimeterToolBase
{
public:
  using CalorimeterToolBase::CalorimeterToolBase;
  virtual ~NestedVolumesCaloTool() = default;

  virtual StatusCode initialize() override final;


protected:
  virtual StatusCode collectCells(std::function<void(uint64_t)> cellFunc) const override final;

private:
  /// Name of the detector readout
  Gaudi::Property<std::string> m_readoutName{this, "readoutName", "ECalHitsPhiEta", "Name of the detector readout"};
  /// Name of active volumes (if different than all)
  Gaudi::Property<std::vector<std::string>> m_activeVolumeName{
      this, "activeVolumeName", {"LAr_sensitive"}, "Name of active volumes (if different than all)"};
  /// Name of active layers for sampling calorimeter
  Gaudi::Property<std::vector<std::string>> m_activeFieldName{
      this, "activeFieldName", {"active_layer"}, "Name of active layers for sampling calorimeter"};
  /// Name of the fields describing the segmented volume
  Gaudi::Property<std::vector<std::string>> m_fieldNames{
      this, "fieldNames", {}, "Name of the fields describing the segmented volume"};
  /// Values of the fields describing the segmented volume
  Gaudi::Property<std::vector<int>> m_fieldValues{
      this, "fieldValues", {}, "Values of the fields describing the segmented volume"};
};

#endif /* RECCALORIMETER_NESTEDVOLUMESCALOTOOL_H */
