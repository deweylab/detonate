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

size_t choose_2(size_t n) { return n*(n-1)/2; }

double pair_denom(Lens lens, Taus taus)
{
  double denom = 0.0;
  for (size_t i = 0; i < lens.size(); ++i)
    denom += taus[i]*choose_2(lens[i]+1);
  return denom;
}

double nucl_denom(Lens lens, Taus taus)
{
  double denom = 0.0;
  for (size_t i = 0; i < lens.size(); ++i)
    denom += taus[i]*lens[i];
  return denom;
}

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
BOOST_AUTO_TEST_CASE(covered_two_to_one)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 999,   0, 999, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 499, 500, 999, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0);
  CHECK_CLOSE(recall.nucl, 1.0);
  CHECK_CLOSE(recall.tran, 1.0);
}

// a0 and a1 cover b0 completely together, but after subtracting a0->b0, the
// a1->b1 is too short so it is thrown away.
//  a: 0000000000
//          11111
//  b: ----------
BOOST_AUTO_TEST_CASE(mostly_covered_two_to_one)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 989,   0, 989, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 499, 500, 999, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*choose_2(990+1)/choose_2(1000+1));
  CHECK_CLOSE(recall.nucl, 1.0*990/1000);
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

// a1->b0 is bigger than a0->b0 overall, so it is dealt with first. The small
// fragment of a1->a0 is contained within a0->0, so a0->b0 is split.
//  a: 0000000000  
//        111       111111111111111
//  b: ----------------------------
BOOST_AUTO_TEST_CASE(two_to_one_split_alignment)
{
  tagged_alignment al1{ 0, 0, Segs{ {  0, 299,    0,  299, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {  0,  49,  100,  149, {}, {}},
                                    {100, 999, 1000, 1899, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{2000}, Taus{1.0});
  //                            a0->b0 part 1     a1->b0 part 1    a0->b0 part 2     a1->b0 part 2
  //                            0-99              100-149          150-299           100-999
  CHECK_CLOSE(recall.pair, 1.0*(choose_2(100+1) + choose_2(50+1) + choose_2(150+1) + choose_2(900+1))/choose_2(2000+1));
  CHECK_CLOSE(recall.nucl, 1.0*(         100    +          50    +          150    +          900   )/         2000   );
  CHECK_CLOSE(recall.tran, 0.0);
}

// b0 is long but has low weight, b1 is short but has high weight, so the
// alignment a0->b1 is considered before a0->b0. As a result, the alignment
// a0->b0 is split.
//  a: 0000000000
//  b: 0000000000  with very low expression
//        111      with very high expression
BOOST_AUTO_TEST_CASE(one_to_two_split_alignment)
{
  tagged_alignment al1{ 0, 0, Segs{ {  0, 499, 0, 499, {}, {}} }, nan("") };
  tagged_alignment al2{ 0, 1, Segs{ {200, 299, 0,  99, {}, {}} }, nan("") };
  size_t A_card = 1, B_card = 2;
  Lens lens{500, 100};
  Taus taus{0.01, 0.99};
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, lens, taus);
  CHECK_CLOSE(recall.pair, 1.0*(taus[0]*choose_2(200+1) + taus[1]*choose_2(100+1) + taus[0]*choose_2(200+1))/pair_denom(lens, taus));
  CHECK_CLOSE(recall.nucl, 1.0*(taus[0]*         200    + taus[1]*         100    + taus[0]*         200   )/nucl_denom(lens, taus));
  CHECK_CLOSE(recall.tran, taus[1]);
}

//  a0: ---x---
//  b0: -------
BOOST_AUTO_TEST_CASE(one_mismatch)
{
  tagged_alignment al1{ 0, 0, Segs{ {1000, 1099, 0, 99, {1050}, {50}} }, nan("") };
  size_t A_card = 1, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1}, A_card, B_card, Lens{100}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*choose_2(99+1)/choose_2(100+1));
  CHECK_CLOSE(recall.nucl, 1.0*99/100);
  CHECK_CLOSE(recall.tran, 1.0);
}

//  a0: -x-x---
//  b0: -------
BOOST_AUTO_TEST_CASE(two_mismatches)
{
  tagged_alignment al1{ 0, 0, Segs{ {1000, 1099, 0, 99, {1048, 1050}, {48, 50}} }, nan("") };
  size_t A_card = 1, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1}, A_card, B_card, Lens{100}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*choose_2(98+1)/choose_2(100+1));
  CHECK_CLOSE(recall.nucl, 1.0*98/100);
  CHECK_CLOSE(recall.tran, 1.0);
}

//  a0: --xx---   discarded
//  a1: -----x-
//  b0: -------
BOOST_AUTO_TEST_CASE(disjoint_mismatches)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 99, 0, 99, {48, 50}, {48, 50}} }, nan("") };
  tagged_alignment al2{ 0, 0, Segs{ {0, 99, 0, 99, {    52}, {    52}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{100}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*choose_2(99+1)/choose_2(100+1));
  CHECK_CLOSE(recall.nucl, 1.0*99/100);
  CHECK_CLOSE(recall.tran, 1.0);
}

//  a0: ----x-x-x-      used in full
//  a1:       --x----   clipped
//  b0: -------------
BOOST_AUTO_TEST_CASE(more_complicated_mismatches)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 99,   0,  99, {50, 80, 81}, {50, 80, 81}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 74,  75, 149, {81-75     }, {81        }} }, nan("") };
  size_t A_card = 2, B_card = 1;
  stats_tuple recall = do_it_all_wrapper({al1, al2}, A_card, B_card, Lens{150}, Taus{1.0});
  CHECK_CLOSE(recall.pair, 1.0*(choose_2(97+1) + choose_2(50+1))/choose_2(150+1));
  CHECK_CLOSE(recall.nucl, 1.0*(97+50)/150);
  CHECK_CLOSE(recall.tran, 1.0);
}

// A: 00000000   2222  11111111
// B: 000000000  1111  2222222222
BOOST_AUTO_TEST_CASE(several_B_elements)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 99, 0, 99, {}, {}} }, nan("") };
  tagged_alignment al2{ 2, 1, Segs{ {0, 49, 0, 49, {}, {}} }, nan("") };
  tagged_alignment al3{ 1, 2, Segs{ {0, 99, 0, 99, {}, {}} }, nan("") };
  size_t A_card = 3, B_card = 3;
  Lens lens{110, 50, 120};
  Taus taus{0.1, 0.3, 0.6};
  stats_tuple recall = do_it_all_wrapper({al1, al2, al3}, A_card, B_card, lens, taus);
  CHECK_CLOSE(recall.pair, (taus[0]*choose_2(100+1) + taus[1]*choose_2(50+1) + taus[2]*choose_2(100+1)) / pair_denom(lens, taus));
  CHECK_CLOSE(recall.nucl, (taus[0]*         100    + taus[1]*         50    + taus[2]*         100   ) / nucl_denom(lens, taus));
  CHECK_CLOSE(recall.tran, taus[1]);
}

//    discarded    discarded    kept
// A: 0000000000   0000000000   0000000000
// B: 000000000    1111         2222222222
BOOST_AUTO_TEST_CASE(several_B_elements_from_one_A_element)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 89, 0, 89, {}, {}} }, nan("") };
  tagged_alignment al2{ 0, 1, Segs{ {0, 49, 0, 49, {}, {}} }, nan("") };
  tagged_alignment al3{ 0, 2, Segs{ {0, 99, 0, 99, {}, {}} }, nan("") };
  size_t A_card = 1, B_card = 3;
  Lens lens{ 90,  50, 100};
  Taus taus{0.1, 0.3, 0.6};
  stats_tuple recall = do_it_all_wrapper({al1, al2, al3}, A_card, B_card, lens, taus);
  CHECK_CLOSE(recall.pair, (taus[2]*choose_2(100+1)) / pair_denom(lens, taus));
  CHECK_CLOSE(recall.nucl, (taus[2]*         100   ) / nucl_denom(lens, taus));
  CHECK_CLOSE(recall.tran, taus[2]);
}

//    discarded    discarded    kept
// A: 0000x0x00    1x11         2222222222
// B: 0000000000   0000000000   0000000000
BOOST_AUTO_TEST_CASE(several_A_elements_from_one_B_element)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 89, 0, 89, {60, 62}, {60, 62}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 49, 0, 49, {30    }, {30    }} }, nan("") };
  tagged_alignment al3{ 2, 0, Segs{ {0, 99, 0, 99, {      }, {      }} }, nan("") };
  size_t A_card = 3, B_card = 1;
  Lens lens{100};
  Taus taus{1.0};
  stats_tuple recall = do_it_all_wrapper({al1, al2, al3}, A_card, B_card, lens, taus);
  CHECK_CLOSE(recall.pair, (taus[0]*choose_2(100+1)) / pair_denom(lens, taus));
  CHECK_CLOSE(recall.nucl, (taus[0]*         100   ) / nucl_denom(lens, taus));
  CHECK_CLOSE(recall.tran, taus[0]);
}
