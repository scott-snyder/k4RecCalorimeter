/**
 * @file RecCaloCommon/tests/IDMapIndexer_test.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Unit test for IDMapIndexer.
 */

#undef NDEBUG
#include "RecCaloCommon/IDMapIndexer.h"
#include "DD4hep/IDDescriptor.h"
#include <vector>
#include <algorithm>
#include <span>
#include <cassert>


using mapkey_t = uint64_t; // libc defines key_t...
using mapkey_span = std::span<const mapkey_t>;


// Helper to make list of IDS.
#include "make_ecal_ids.icc"


//************************************************************************


void test1 (mapkey_span ids)
{
  using Indexer_t = k4::recCalo::IDMapIndexer<3>;
  using IDMap_t = Indexer_t::IDMap_t;
  using FieldDesc_t = Indexer_t::FieldDesc_t;
  dd4hep::IDDescriptor desc ("desc", ecal_descstr);

  std::vector<FieldDesc_t> fielddescs
    {
      IDMap_t::makeDesc (*desc.field ("layer")),
      IDMap_t::makeDesc (*desc.field ("theta")),
      IDMap_t::makeDesc (*desc.field ("module")),
    };

  Indexer_t map (4, 6, fielddescs, ids);
  assert (map.detIDs().size() == 1);
  assert (map.detIDs()[0] == 4);
  assert (map.detIDBits() == 6);
  assert (map.byteSize() > 100000);

  size_t ncell = ids.size();
  assert (map.cellIDs().size() == ncell);
  for (size_t i = 0; i < ncell; i++) {
    assert (map.cellIDs()[i] == ids[i]);
    assert (map.index(ids[i]) == i);
  }
  assert (map.index (0) == Indexer_t::INVALID);
}


int main()
{
  std::vector<mapkey_t> ids = make_ecal_ids();
  test1 (ids);
  return 0;
}
