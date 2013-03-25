#include <string>
#include <iostream>
#define BOOST_TEST_MODULE test_summarize_matched_meat
#define BOOST_TEST_DYN_LINK
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include "summarize_matched_meat.hh"

#define CHECK_CLOSE(x, y) \
  do { \
    double w = x; \
    double v = y; \
    double z = fabs(w - v); \
    if (z > 1e-5) { \
      std::string sx = boost::lexical_cast<std::string>(w); \
      std::string sy = boost::lexical_cast<std::string>(v); \
      std::string sz = boost::lexical_cast<std::string>(z); \
      BOOST_ERROR(sx + " != " + sy + " [fabs(x - y) = " + sz + " > 1e-5]"); \
    } \
  } while (0)

typedef std::vector<alignment_segment> Segs;
typedef std::vector<size_t>            Lens;
typedef std::vector<double>            Taus;

// one perfect alignment from a0 -> b0, where these are the only seqs present
BOOST_AUTO_TEST_CASE(sanity)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 100, 0, 100, {}, {}} } };
  size_t A_card = 1, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1}, A_card, B_card, Lens{100}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0);
  CHECK_CLOSE(recall.nucl, 1.0);
  CHECK_CLOSE(recall.tran, 1.0);
}
