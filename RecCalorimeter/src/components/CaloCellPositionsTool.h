// This file's extension implies that it's C, but it's really -*- C++ -*-.
/**
 * @file k4RecCalorimeter/RecCalorimeter/src/components/CaloCellPositionsTool.h
 * @author scott snyder <snyder@bnl.gov>
 * @date Apr, 2025
 * @brief Generic tool to find positions of calorimeter cells.
 */


#ifndef RECCALORIMETER_CALOCELLPOSITIONSTOOL_H
#define RECCALORIMETER_CALOCELLPOSITIONSTOOL_H


#include "GaudiKernel/AlgTool.h"
#include "k4Interface/ICellPositionsTool.h"
#include "DD4hep/Segmentations.h"
#include "DD4hep/Volumes.h"


class CaloCellPositionsTool
  : public extends<AlgTool, ICellPositionsTool>
{
public:
  using base_class::base_class;

  virtual StatusCode initialize() override;
  virtual void getPositions(const edm4hep::CalorimeterHitCollection& aCells,
                            edm4hep::CalorimeterHitCollection&       outputColl) const override;

  virtual dd4hep::Position xyzPosition(const uint64_t& aCellId) const override;
  virtual int              layerId(const uint64_t& aCellId)     const override;


private:
  Gaudi::Property<std::string> m_readoutName
  { this, "readoutName", "", "Name of the readout for this detector" };
  Gaudi::Property<std::string> m_layerFieldName
  { this, "layerFieldName", "layer", "Name of the decoder field for layer" };

  dd4hep::VolumeManager m_volman;
  dd4hep::Segmentation m_segmentation;
  const dd4hep::BitFieldCoder* m_decoder = nullptr;
  size_t m_layerFieldIdx = 0;
};


#endif // not RECCALORIMETER_CALOCELLPOSITIONSTOOL_H
