/**
 * @file IDMapIndexer_test.cpp
 * @author scott snyder <snyder@bnl.gov>
 * @date Jan, 2026
 * @brief Unit test for IDMapIndexer.
 */

#undef NDEBUG
#include "RecCaloCommon/IDMapIndexer.h"
#include "DD4hep/IDDescriptor.h"
#include <vector>
#include <span>
#include <cassert>


using mapkey_t = uint64_t; // libc defines key_t...
using mapkey_span = std::span<const mapkey_t>;


//************************************************************************
// Helper to make a list of IDs for testing.
// Hardcoded to match the set of Allegro ECal barrel IDs as of this writing.
//

const unsigned int numLayers = 11;
const char* const descstr = "system:4,cryo:1,type:3,subtype:3,layer:8,module:11,theta:10";


std::vector<mapkey_t> make_ids()
{
  std::vector<mapkey_t> ids;
  ids.reserve (2041344);
  dd4hep::IDDescriptor desc ("desc", descstr);
  auto decoder = desc.decoder();

  static const unsigned module_high = 1534;
  static const unsigned module_step = 2;
  static const unsigned theta_low[11] = { 8, 12, 12, 16, 20, 20, 24, 28, 32, 36, 40 };
  static const unsigned theta_high[11] = { 788, 787, 784, 780, 776, 776, 772, 768, 764, 760, 756 };
  static const unsigned theta_step[11] = { 4, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4 };

  size_t module_index = decoder->index ("module");
  size_t theta_index  = decoder->index ("theta");

  for (unsigned int ilayer = 0; ilayer < numLayers; ilayer++) {
    mapkey_t id = 0;
    decoder->set (id, "system", 4);
    decoder->set (id, "cryo", 0);
    decoder->set (id, "type", 0);
    decoder->set (id, "subtype", 0);
    decoder->set (id, "layer", ilayer);
    decoder->set (id, theta_index, 0);
    decoder->set (id, module_index, 0);

    for (unsigned int mod = 0; mod <= module_high; mod += module_step) {
      decoder->set (id, module_index, mod);
      for (unsigned int theta = theta_low[ilayer]; theta <= theta_high[ilayer]; theta += theta_step[ilayer]) {
        decoder->set (id, theta_index, theta);
        ids.push_back (id);
      }
    }
  }
  std::ranges::sort (ids);
  return ids;
}



//************************************************************************


void test1 (mapkey_span ids)
{
  using Indexer_t = k4::recCalo::IDMapIndexer<3>;
  using IDMap_t = Indexer_t::IDMap_t;
  using FieldDesc_t = Indexer_t::FieldDesc_t;
  dd4hep::IDDescriptor desc ("desc", descstr);

  std::vector<FieldDesc_t> fielddescs
    {
      IDMap_t::makeDesc (*desc.field ("layer")),
      IDMap_t::makeDesc (*desc.field ("theta")),
      IDMap_t::makeDesc (*desc.field ("module")),
    };

  Indexer_t map (4, fielddescs, ids);
  assert (map.detIDs().size() == 1);
  assert (map.detIDs()[0] == 4);

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
  std::vector<mapkey_t> ids = make_ids();
  test1 (ids);
  return 0;
}
