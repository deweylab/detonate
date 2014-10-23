// Copyright (c) 2013
// Nathanael Fillmore (University of Wisconsin-Madison)
// nathanae@cs.wisc.edu
//
// This file is part of REF-EVAL.
//
// REF-EVAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// REF-EVAL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with REF-EVAL.  If not, see <http://www.gnu.org/licenses/>.

#include <string>
#include <iostream>
#include <random>       // std::default_random_engine
#include <chrono>       // std::chrono::system_clock
#define BOOST_TEST_MODULE test_summarize_matched_meat
#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/foreach.hpp>
#include "re_matched.hh"

using namespace re::matched;

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
  double pair_recall = compute_recall<pair_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0);
  CHECK_CLOSE(nucl_recall, 1.0);
  CHECK_CLOSE(tran_recall, 1.0);
}

// no alignments from a0 -> b0, where these are the only seqs present
BOOST_AUTO_TEST_CASE(sanity_2)
{
  size_t A_card = 1, B_card = 1;
  double pair_recall = compute_recall<pair_helper>({}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 0.0);
  CHECK_CLOSE(nucl_recall, 0.0);
  CHECK_CLOSE(tran_recall, 0.0);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0*(choose_2(500+1) + choose_2(500+1))/choose_2(1000+1));
  CHECK_CLOSE(nucl_recall, 1.0);
  CHECK_CLOSE(tran_recall, 1.0);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0);
  CHECK_CLOSE(nucl_recall, 1.0);
  CHECK_CLOSE(tran_recall, 1.0);
}

// a0 and a1 cover b0 completely together, but after subtracting a0->b0, the
// a1->b1 is too short so it is thrown away.
//  a: 000000000
//          11111
//  b: ----------
BOOST_AUTO_TEST_CASE(mostly_covered_two_to_one)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 989,   0, 989, {}, {}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 499, 500, 999, {}, {}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 20);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 20);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 20);
  CHECK_CLOSE(pair_recall, 1.0*choose_2(990+1)/choose_2(1000+1));
  CHECK_CLOSE(nucl_recall, 1.0*990/1000);
  CHECK_CLOSE(tran_recall, 1.0);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0*(choose_2(600+1) + choose_2(400+1))/choose_2(1000+1));
  CHECK_CLOSE(nucl_recall, 1.0);
  CHECK_CLOSE(tran_recall, 1.0);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{1000}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0*(choose_2(490+1) + choose_2(480+1))/choose_2(1000+1));
  CHECK_CLOSE(nucl_recall, 1.0*(490 + 480)/1000);
  CHECK_CLOSE(tran_recall, 1.0);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{2000}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{2000}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{2000}, Taus{1.0}, 0);
  //                            a0->b0 part 1     a1->b0 part 1    a0->b0 part 2     a1->b0 part 2
  //                            0-99              100-149          150-299           100-999
  CHECK_CLOSE(pair_recall, 1.0*(choose_2(100+1) + choose_2(50+1) + choose_2(150+1) + choose_2(900+1))/choose_2(2000+1));
  CHECK_CLOSE(nucl_recall, 1.0*(         100    +          50    +          150    +          900   )/         2000   );
  CHECK_CLOSE(tran_recall, 0.0);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, lens, taus, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, lens, taus, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, lens, taus, 0);
  CHECK_CLOSE(pair_recall, 1.0*(taus[0]*choose_2(200+1) + taus[1]*choose_2(100+1) + taus[0]*choose_2(200+1))/pair_denom(lens, taus));
  CHECK_CLOSE(nucl_recall, 1.0*(taus[0]*         200    + taus[1]*         100    + taus[0]*         200   )/nucl_denom(lens, taus));
  CHECK_CLOSE(tran_recall, taus[1]);
}

//  a0: ---x---
//  b0: -------
BOOST_AUTO_TEST_CASE(one_mismatch)
{
  tagged_alignment al1{ 0, 0, Segs{ {1000, 1099, 0, 99, {1050}, {50}} }, nan("") };
  size_t A_card = 1, B_card = 1;
  double pair_recall = compute_recall<pair_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0*choose_2(99+1)/choose_2(100+1));
  CHECK_CLOSE(nucl_recall, 1.0*99/100);
  CHECK_CLOSE(tran_recall, 1.0);
}

//  a0: -x-x---
//  b0: -------
BOOST_AUTO_TEST_CASE(two_mismatches)
{
  tagged_alignment al1{ 0, 0, Segs{ {1000, 1099, 0, 99, {1048, 1050}, {48, 50}} }, nan("") };
  size_t A_card = 1, B_card = 1;
  double pair_recall = compute_recall<pair_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0*choose_2(98+1)/choose_2(100+1));
  CHECK_CLOSE(nucl_recall, 1.0*98/100);
  CHECK_CLOSE(tran_recall, 1.0);
}

//  a0: --xx---   discarded
//  a1: -----x-
//  b0: -------
BOOST_AUTO_TEST_CASE(disjoint_mismatches)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 99, 0, 99, {48, 50}, {48, 50}} }, nan("") };
  tagged_alignment al2{ 0, 0, Segs{ {0, 99, 0, 99, {    52}, {    52}} }, nan("") };
  size_t A_card = 2, B_card = 1;
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{100}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0*choose_2(99+1)/choose_2(100+1));
  CHECK_CLOSE(nucl_recall, 1.0*99/100);
  CHECK_CLOSE(tran_recall, 1.0);
}

//  a0: ----x-x-x-      used in full
//  a1:       --x----   clipped
//  b0: -------------
BOOST_AUTO_TEST_CASE(more_complicated_mismatches)
{
  tagged_alignment al1{ 0, 0, Segs{ {0, 99,   0,  99, {50, 80, 81}, {50, 80, 81}} }, nan("") };
  tagged_alignment al2{ 1, 0, Segs{ {0, 74,  75, 149, {81-75     }, {81        }} }, nan("") };
  size_t A_card = 2, B_card = 1;
  double pair_recall = compute_recall<pair_helper>({al1, al2}, A_card, B_card, Lens{150}, Taus{1.0}, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2}, A_card, B_card, Lens{150}, Taus{1.0}, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2}, A_card, B_card, Lens{150}, Taus{1.0}, 0);
  CHECK_CLOSE(pair_recall, 1.0*(choose_2(97+1) + choose_2(50+1))/choose_2(150+1));
  CHECK_CLOSE(nucl_recall, 1.0*(97+50)/150);
  CHECK_CLOSE(tran_recall, 1.0);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  CHECK_CLOSE(pair_recall, (taus[0]*choose_2(100+1) + taus[1]*choose_2(50+1) + taus[2]*choose_2(100+1)) / pair_denom(lens, taus));
  CHECK_CLOSE(nucl_recall, (taus[0]*         100    + taus[1]*         50    + taus[2]*         100   ) / nucl_denom(lens, taus));
  CHECK_CLOSE(tran_recall, taus[1]);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  CHECK_CLOSE(pair_recall, (taus[2]*choose_2(100+1)) / pair_denom(lens, taus));
  CHECK_CLOSE(nucl_recall, (taus[2]*         100   ) / nucl_denom(lens, taus));
  CHECK_CLOSE(tran_recall, taus[2]);
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
  double pair_recall = compute_recall<pair_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  double nucl_recall = compute_recall<nucl_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  double tran_recall = compute_recall<tran_helper>({al1, al2, al3}, A_card, B_card, lens, taus, 0);
  CHECK_CLOSE(pair_recall, (taus[0]*choose_2(100+1)) / pair_denom(lens, taus));
  CHECK_CLOSE(nucl_recall, (taus[0]*         100   ) / nucl_denom(lens, taus));
  CHECK_CLOSE(tran_recall, taus[0]);
}

//    b0                    b1
// B  -----------------     -------------
//         /    \   /  \   /  | / |
//        /      \ /    \ /   |/  |
//       /  x     /  y   \  z / w |
//      /        / \    / \  /|   |
// A    ---------   --------- ---------
//      a0          a1        a2
//
// Assume:
// - x > y > z > w
// - y - x < z
// - y - x is contained in z, wrt A
//
// Step 1: Process alignment x, resulting in
//
//    b0                    b1
// B  ------------------   --------------
//         /        /\ \   /  | / |
//        /        /  \ \ /   |/  |
//       /  x     /    \*\  z / w |        * = y - x
//      /        /      / \  /|   |
// A    ---------   --------- ---------
//      a0          a1        a2
//
// Step 2: Process alignment z, resulting in
//
//    b0                    b1
// B  -----------------    --------------
//         /        /      /    /||
//        /        /      /    / ||
//       /  x     /      /  z /  *|        * = w - z
//      /        /      /    /   ||
// A    ---------   --------- ---------
//      a0          a1        a2
//
// Now we have a 1-1 mathing. The intervals of B used to compute recall are:
//
//    b0                    b1
// B  -----[--------]--    [----][]------
//         /        /      /    /||
//        ...
BOOST_AUTO_TEST_CASE(complicated_ordering_1)
{
  tagged_alignment x{ 0, 0, Segs{ {  0, 499, 200, 699, {}, {}} }, nan("") }; // #id = 500
  tagged_alignment y{ 1, 0, Segs{ {150, 549, 500, 899, {}, {}} }, nan("") }; // #id = 400
  tagged_alignment z{ 1, 1, Segs{ {300, 599,   0, 299, {}, {}} }, nan("") }; // #id = 300
  tagged_alignment w{ 2, 1, Segs{ {  0, 249, 150, 399, {}, {}} }, nan("") }; // #id = 250
  // Note:
  // * y - x is [irrelevant] -> [700, 899], #id = 200
  // * y - x is [350, 549] -> [700, 899] which is contained within z, wrt A, when the difference is taken wrt B.
  // * w - z is [irrelevant] -> [300, 399], #id = 100
  size_t A_card = 3, B_card = 2;
  Lens lens{1000, 1000};
  Taus taus{0.5, 0.5};
  double pair_recall = compute_recall<pair_helper>({x, y, z, w}, A_card, B_card, lens, taus, 0);
  double nucl_recall = compute_recall<nucl_helper>({x, y, z, w}, A_card, B_card, lens, taus, 0);
  double tran_recall = compute_recall<tran_helper>({x, y, z, w}, A_card, B_card, lens, taus, 0);
  CHECK_CLOSE(pair_recall, (taus[0]*choose_2(500+1) + taus[1]*choose_2(300+1) + taus[1]*choose_2(100+1)) / pair_denom(lens, taus));
  CHECK_CLOSE(nucl_recall, (taus[0]*         500    + taus[1]*         300    + taus[1]*         100   ) / nucl_denom(lens, taus));
  CHECK_CLOSE(tran_recall, 0.0);
}

// Same initial picture:
//
//    b0                    b1
// B  -----------------     -------------
//         /    \   /  \   /  | / |
//        /      \ /    \ /   |/  |
//       /  x     /  y   \  z / w |
//      /        / \    / \  /|   |
// A    ---------   --------- ---------
//      a0          a1        a2
//
// Assume:
// - x > y > w > z  (w and z are interchanged)
// - y - x < w
// - y - x > z - w
// - y - x is contained in z, wrt A
//
// Step 1: Process alignment x, resulting in
//
//    b0                    b1
// B  ------------------   --------------
//         /        /\ \   /  | / |
//        /        /  \ \ /   |/  |
//       /  x     /    \*\  z / w |        * = y - x
//      /        /      / \  /|   |
// A    ---------   --------- ---------
//      a0          a1        a2
//
// Step 2: Process alignment w, resulting in
//
//    b0                    b1
// B  ------------------   --------------
//         /        /\ \   /  |   |
//        /        /  \ \ /  /|   |
//       /  x     /    \*\ +/ | w |        * = y - x
//      /        /      / \/  |   |        + = z - w
// A    ---------   --------- ---------
//      a0          a1        a2
//
// Step 3: Process alignment y - x, resulting in
//
//    b0                    b1
// B  ------------------   --------------
//         /        /\ \  +/  |   |
//        /        /  \*\//   |   |
//       /  x     /    \ \    | w |        * = y - x
//      /        /     // \   |   |        + = (z - w) - (y - x)
// A    ---------   --------- ---------
//      a0          a1        a2
//
// Now we have a 1-1 mathing. The intervals of B used to compute recall are:
//
//    b0                    b1
// B  -----[--------][-]   []-[---]------
//         /        /\ \  //  |   |        * = y - x
//          x         *   +     w          + = (z - w) - (y - x)
//        ...                              
BOOST_AUTO_TEST_CASE(complicated_ordering_2)
{
  tagged_alignment x{ 0, 0, Segs{ {  0, 499, 200, 699, {}, {}} }, nan("") }; // #id = 500
  tagged_alignment y{ 1, 0, Segs{ {150, 549, 500, 899, {}, {}} }, nan("") }; // #id = 400
  tagged_alignment z{ 1, 1, Segs{ {300, 549,   0, 249, {}, {}} }, nan("") }; // #id = 250
  tagged_alignment w{ 2, 1, Segs{ {  0, 299, 100, 399, {}, {}} }, nan("") }; // #id = 300
  // Note:
  // * y - x is [350, 549] -> [700, 899], #id = 200, which is contained within z, wrt A, when the difference is taken wrt B.
  // * z - w is [300, 399] -> [0, 99], #id = 100
  // * (z - w) - (y - x) is [300, 349] -> [0, 49], #id = 50
  size_t A_card = 3, B_card = 2;
  Lens lens{1000, 1000};
  Taus taus{0.5, 0.5};
  double pair_recall = compute_recall<pair_helper>({x, y, z, w}, A_card, B_card, lens, taus, 0);
  double nucl_recall = compute_recall<nucl_helper>({x, y, z, w}, A_card, B_card, lens, taus, 0);
  double tran_recall = compute_recall<tran_helper>({x, y, z, w}, A_card, B_card, lens, taus, 0);
  CHECK_CLOSE(pair_recall, (taus[0]*(choose_2(500+1)+choose_2(200+1)) + taus[1]*(choose_2(50+1)+choose_2(300+1))) / pair_denom(lens, taus));
  CHECK_CLOSE(nucl_recall, (taus[0]*(         500   +         200   ) + taus[1]*(         50   +         300   )) / nucl_denom(lens, taus));
  CHECK_CLOSE(tran_recall, 0.0);
}

BOOST_AUTO_TEST_CASE(allpairs)
{
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine rng(seed);

  for (size_t trial = 0; trial < 50; ++trial) {

    size_t n = 6;       // number of elements of A and of B
    size_t len = n*n*2; // length of each element of A and of B

    // E[i][j] is going to be the number of mismatches in the alignment from a[i] -> b[j].
    //
    // Simulate E randomly such that all nonzeros are distinct integers between
    // 1 and n*n.
    std::vector<std::vector<size_t>> E(n, std::vector<size_t>(n, 0)); // all zeros
    {
      std::vector<size_t> tmp;
      for (size_t k = 0; k < n*n; ++k)
        tmp.push_back(k+1);
      std::shuffle(tmp.begin(), tmp.end(), rng);
      std::vector<size_t>::iterator it = tmp.begin();
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          E[i][j] = *it;
          ++it;
        }
      }
    }

    // A[i][j] is going to be 1 if there is an alignment a[i] -> b[j], and 0 otherwise.
    std::vector<std::vector<size_t>> A(n, std::vector<size_t>(n, 1)); // all ones

    // D[i][j] is going to be 1 if there the alignment a[i] -> b[j] has already been selected.
    std::vector<std::vector<size_t>> D(n, std::vector<size_t>(n, 0)); // all zeros

    // Figure out the matching that should be selected.
    for (;;) {

      // Choose i,j = argmin E[i][j] subject to constraints.
      size_t best_i = 0, best_j = 0, best_E = std::numeric_limits<size_t>::max();
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          if (D[i][j] == 0 &&     // is not already done
              A[i][j] == 1 &&     // is still aligned
              E[i][j] < best_E) { // has fewer mismatches than current candidate
            best_i = i;
            best_j = j;
            best_E = E[i][j];
          }
        }
      }
      D[best_i][best_j] = 1;

      // Set all other entries of A in (i,j)'s row and column to zero.
      for (size_t i = 0; i < n; ++i) if (i != best_i) A[i][best_j] = 0;
      for (size_t j = 0; j < n; ++j) if (j != best_j) A[best_i][j] = 0;

      // If nnz(A) == n, we are done.
      size_t nnz = 0;
      for (size_t i = 0; i < n; ++i)
        for (size_t j = 0; j < n; ++j)
          if (A[i][j] != 0)
            ++nnz;
      if (nnz == n)
        break;
      BOOST_CHECK(nnz >= n);

    }

    // Compute predicted recall.
    double pred_pair_recall, pred_nucl_recall;
    {
      double pair_numer = 0.0, nucl_numer = 0.0;
      for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
          if (A[i][j]) {
            size_t num_id = len - E[i][j];
            double tau    = 1.0/n;
            pair_numer += tau*choose_2(num_id+1);
            nucl_numer += tau*         num_id   ;
          }
        }
      }
      double pair_denom = n*(1.0/n)*choose_2(len+1);
      double nucl_denom = n*(1.0/n)*         len   ;
      pred_pair_recall = pair_numer / pair_denom;
      pred_nucl_recall = nucl_numer / nucl_denom;
    }

    // Make all pairs of alignments, in which the alignment from a[i] -> b[j]
    // has E[i][j] mismatches.
    std::vector<tagged_alignment> alignments;
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        std::vector<size_t> mis;
        for (size_t k = 0; k < E[i][j]; ++k)
          mis.push_back(k);
        alignments.push_back(tagged_alignment { i, j, Segs{ {  0, len-1, 0, len-1, mis, mis} }, nan("") });
      }
    }
    size_t A_card = n, B_card = n;
    Lens lens(n, len);
    Taus taus(n, 1.0/n);
    
    // Compute hopefully true recall.
    double pair_recall = compute_recall<pair_helper>(alignments, A_card, B_card, lens, taus, 0);
    double nucl_recall = compute_recall<nucl_helper>(alignments, A_card, B_card, lens, taus, 0);

    // Compare the recalls.
    CHECK_CLOSE(pair_recall, pred_pair_recall);
    CHECK_CLOSE(nucl_recall, pred_nucl_recall);

  } // loop over trials
}
