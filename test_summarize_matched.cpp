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

size_t choose_2(size_t n) { return n*(n-1)/2; }
typedef std::vector<alignment_segment> Segs;
typedef std::vector<size_t>            Lens;
typedef std::vector<double>            Taus;

// one perfect alignment from a0 -> b0, where these are the only seqs present
BOOST_AUTO_TEST_CASE(sanity)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 99, 0, 99, {}, {}} }, nan("") };
  size_t A_card = 1, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1}, A_card, B_card, Lens{100}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0);
  CHECK_CLOSE(recall.nucl, 1.0);
  CHECK_CLOSE(recall.tran, 1.0);
}

// no alignments from a0 -> b0, where these are the only seqs present
BOOST_AUTO_TEST_CASE(sanity_2)
{
  size_t A_card = 1, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({}, A_card, B_card, Lens{100}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 0.0);
  CHECK_CLOSE(recall.nucl, 0.0);
  CHECK_CLOSE(recall.tran, 0.0);
}

// a0 and a1 together cover b0 perfectly
//  a: 00000
//          11111
//  b: ----------
BOOST_AUTO_TEST_CASE(perfect_two_to_one)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 499,   0, 499, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 499, 500, 999, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*(choose_2(500+1) + choose_2(500+1))/choose_2(1000+1));
  CHECK_CLOSE(recall.nucl, 1.0);
  CHECK_CLOSE(recall.tran, 1.0);
}

// a0 covers b0 perfectly by itself, and a1 covers part of b0
//  a: 0000000000
//          11111
//  b: ----------
BOOST_AUTO_TEST_CASE(clipped_two_to_one)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 999,   0, 999, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 499, 500, 999, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;

  stats_tuple recall;
  
  recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0);
  CHECK_CLOSE(recall.nucl, 1.0);
  CHECK_CLOSE(recall.tran, 1.0);

  recall = do_it_all_wrapper({al2, al1}, A_card, B_card, Lens{1000}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0);
  CHECK_CLOSE(recall.nucl, 1.0);
  CHECK_CLOSE(recall.tran, 1.0);
}

// a0 and a1 together cover b0 and partially each other
//  a: 00000
//         111111
//  b: ----------
BOOST_AUTO_TEST_CASE(overlapping_two_to_one)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 599,   0, 599, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 599, 400, 999, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*(choose_2(600+1) + choose_2(400+1))/choose_2(1000+1));
  CHECK_CLOSE(recall.nucl, 1.0);
  CHECK_CLOSE(recall.tran, 1.0);
}

// a0 and a1 together cover 970/1000 of b0
// the alignment from a0->b0 covers 490/500 of a0
// the alignment from a1->b0 covers 480/500 of a1
//  a: 00000
//          11111
//  b: ----------
BOOST_AUTO_TEST_CASE(two_to_one_partial_coverage)
{
  tagged_alignment al1{ 0, 0, Segs{ {1, 490,   1, 490, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {1, 480, 501, 980, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*(choose_2(490+1) + choose_2(480+1))/choose_2(1000+1));
  CHECK_CLOSE(recall.nucl, 1.0*(490 + 480)/1000);
  CHECK_CLOSE(recall.tran, 1.0);
}

#if 0
//  a0: ---x---
//  b0: -------
BOOST_AUTO_TEST_CASE(one_mismatch)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 99, 1000, 1099, {50}, {1050}} }, nan("") };
  std::vector<tagged_alignment> als{al1};
  size_t A_card = 1, B_card = 1;
  stats_tuple recall = do_it_all_wrapper(als, A_card, B_card, Lens{100}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*choose_2(99+1)/choose_2(100+1));
  CHECK_CLOSE(recall.nucl, 1.0*99/100);
  CHECK_CLOSE(recall.tran, 1.0);
}
#endif
