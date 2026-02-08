/**
 * @file RecCaloCommon/tests/MultiIndexer_test.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Feb, 2026
 * @brief Unit test for MultiIndexer.
 */

#undef NDEBUG
#include "RecCaloCommon/MultiIndexer.h"
#include "RecCaloCommon/IDMapIndexer.h"
#include "RecCaloCommon/ICaloCellConstantsSvc.h"
#include "DD4hep/IDDescriptor.h"
#include <vector>
#include <algorithm>


using mapkey_t = uint64_t; // libc defines key_t...
using mapkey_span = std::span<const mapkey_t>;


// Helpers to make lists of IDS.
#include "make_ecal_ids.icc"
#include "make_hcal_ids.icc"


class DummyCellConstantsSvc : public k4::recCalo::ICaloCellConstantsSvc
{
public:
  virtual const std::any* getAnyObj (const std::string& key) const override;
  virtual bool putAnyObj (const std::string& key, std::any&& obj) override;

  virtual std::vector<std::string> getInterfaceNames() const override { return std::vector<std::string>(); }
  virtual unsigned long addRef() const override { return 0; }
  virtual unsigned long release() const override { return 0; }
  virtual unsigned long refCount() const override { return 0; }
  virtual void const* i_cast( const InterfaceID& ) const override { return nullptr; }
  virtual unsigned long decRef() const override { return 0; }
  virtual StatusCode queryInterface(const InterfaceID&, void**) override { return StatusCode::FAILURE; }


private:
  std::map<std::string, std::any> m_objs;
};


const std::any* DummyCellConstantsSvc::getAnyObj (const std::string& key) const
{
  auto it = m_objs.find (key);
  if (it != m_objs.end()) return &it->second;
  return nullptr;
}


bool DummyCellConstantsSvc::putAnyObj (const std::string& key, std::any&& obj)
{
  return m_objs.try_emplace (key, std::move(obj)).second;
}


//************************************************************************


void test1 (mapkey_span ecal_ids, mapkey_span hcal_ids)
{
  using Indexer_t = k4::recCalo::IDMapIndexer<3>;
  using IDMap_t = Indexer_t::IDMap_t;
  using FieldDesc_t = Indexer_t::FieldDesc_t;

  dd4hep::IDDescriptor ecal_desc ("desc", ecal_descstr);
  std::vector<FieldDesc_t> ecal_fielddescs
    {
      IDMap_t::makeDesc (*ecal_desc.field ("layer")),
      IDMap_t::makeDesc (*ecal_desc.field ("theta")),
      IDMap_t::makeDesc (*ecal_desc.field ("module")),
    };

  dd4hep::IDDescriptor hcal_desc ("desc", hcal_descstr);
  std::vector<FieldDesc_t> hcal_fielddescs
    {
      IDMap_t::makeDesc (*hcal_desc.field ("layer")),
      IDMap_t::makeDesc (*hcal_desc.field ("theta")),
      IDMap_t::makeDesc (*hcal_desc.field ("phi")),
    };

  Indexer_t ecal_map (4, 4, ecal_fielddescs, ecal_ids);
  Indexer_t hcal_map (8, 4, hcal_fielddescs, hcal_ids);

  DummyCellConstantsSvc constsSvc;
  std::vector<const k4::recCalo::ICaloIndexer*> indexers { &ecal_map, &hcal_map };
  k4::recCalo::MultiIndexer map (4, indexers, constsSvc);
  assert (map.detIDs().size() == 2);
  assert (map.detIDs()[0] == 4);
  assert (map.detIDs()[1] == 8);
  assert (map.detIDBits() == 4);

  std::span<const uint64_t> cell_ids = map.cellIDs();
  assert (cell_ids.size() == ecal_ids.size() + hcal_ids.size());
  assert (std::equal (ecal_ids.begin(), ecal_ids.end(), cell_ids.begin()));
  assert (std::equal (hcal_ids.begin(), hcal_ids.end(), cell_ids.begin() + ecal_ids.size()));

  size_t ncell = cell_ids.size();
  for (size_t i = 0; i < ncell; i++) {
    assert (map.index(cell_ids[i]) == i);
  }
  assert (map.index (0) == Indexer_t::INVALID);


  //***************************************************
  // Testing for errors during construction.
  using MultiIndexerException = k4::recCalo::MultiIndexer::MultiIndexerException;

#define EXPECT_EXCEPTION(EXP, CODE) do {        \
  bool caught = false;                          \
  try {                                         \
    CODE;                                       \
  }                                             \
  catch (const MultiIndexerException& exc) {    \
    if (std::string (exc.what()) == EXP) {      \
      caught = true;                            \
    }                                           \
    else {                                      \
      std::cerr << exc.what() << "\n";          \
    }                                           \
  }                                             \
  assert (caught);                              \
} while(0)

  EXPECT_EXCEPTION( "MultiIndexerException: detIDBits 40 out of range; should be between 1 and 8",
                    k4::recCalo::MultiIndexer map2 (40, indexers, constsSvc) );

  std::vector<const k4::recCalo::ICaloIndexer*> indexers3 { &map };
  EXPECT_EXCEPTION( "MultiIndexerException: bad indexer returns detIDs of size 2",
                    k4::recCalo::MultiIndexer map3 (4, indexers3, constsSvc) );

  Indexer_t ecal_map4 (40, 10, ecal_fielddescs, ecal_ids);
  std::vector<const k4::recCalo::ICaloIndexer*> indexers4 { &ecal_map4 };
  EXPECT_EXCEPTION( "MultiIndexerException: bad indexer returns out-of-range detID 40",
                    k4::recCalo::MultiIndexer map4 (4, indexers4, constsSvc) );
}


int main()
{
  std::vector<mapkey_t> ecal_ids = make_ecal_ids();
  std::vector<mapkey_t> hcal_ids = make_hcal_ids();
  test1 (ecal_ids, hcal_ids);
  return 0;
}
