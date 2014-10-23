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

#define BOOST_TEST_MODULE test_alignment_segment
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "alignment_segment.hh"

typedef std::vector<alignment_segment> VAS;

std::ostream& operator<<(std::ostream& out, const std::vector<alignment_segment>& v)
{
  out << "<vector>\n";
  BOOST_FOREACH(const alignment_segment& seg, v)
    out << "  " << seg << "\n";
  out << "</vector>\n";
  return out;
}

// test wrt b, then interchange a and b and test wrt a
void check_empty_both_ways(size_t seg1_start_a, size_t seg1_end_a, size_t seg1_start_b, size_t seg1_end_b,
                           size_t seg2_start_a, size_t seg2_end_a, size_t seg2_start_b, size_t seg2_end_b)
{
  {
    VAS segs1 { {seg1_start_a, seg1_end_a, seg1_start_b, seg1_end_b, {}, {}} };
    VAS segs2 { {seg2_start_a, seg2_end_a, seg2_start_b, seg2_end_b, {}, {}} };
    BOOST_CHECK(intersects<segment_ops_wrt_b>(segs1, segs2));

    VAS pred = subtract<segment_ops_wrt_b>(segs2, segs1);
    VAS trut {};
    BOOST_CHECK_EQUAL(pred, trut);
  }

  {
    VAS segs1 { {seg1_start_b, seg1_end_b, seg1_start_a, seg1_end_a, {}, {}} };
    VAS segs2 { {seg2_start_b, seg2_end_b, seg2_start_a, seg2_end_a, {}, {}} };
    BOOST_CHECK(intersects<segment_ops_wrt_a>(segs1, segs2));

    VAS pred = subtract<segment_ops_wrt_a>(segs2, segs1);
    VAS trut {};
    BOOST_CHECK_EQUAL(pred, trut);
  }
}

// test wrt b, then interchange a and b and test wrt a
void check_equal_both_ways(size_t seg1_start_a, size_t seg1_end_a, size_t seg1_start_b, size_t seg1_end_b,
                           size_t seg2_start_a, size_t seg2_end_a, size_t seg2_start_b, size_t seg2_end_b,
                           bool does_intersect,
                           size_t trut_start_a, size_t trut_end_a, size_t trut_start_b, size_t trut_end_b)
{
  {
    VAS segs1 { {seg1_start_a, seg1_end_a, seg1_start_b, seg1_end_b, {}, {}} };
    VAS segs2 { {seg2_start_a, seg2_end_a, seg2_start_b, seg2_end_b, {}, {}} };
    BOOST_CHECK_EQUAL(intersects<segment_ops_wrt_b>(segs1, segs2), does_intersect);

    VAS pred = subtract<segment_ops_wrt_b>(segs2, segs1);
    VAS trut { {trut_start_a, trut_end_a, trut_start_b, trut_end_b, {}, {}} };
    BOOST_CHECK_EQUAL(pred, trut);
  }

  {
    VAS segs1 { {seg1_start_b, seg1_end_b, seg1_start_a, seg1_end_a, {}, {}} };
    VAS segs2 { {seg2_start_b, seg2_end_b, seg2_start_a, seg2_end_a, {}, {}} };
    BOOST_CHECK_EQUAL(intersects<segment_ops_wrt_a>(segs1, segs2), does_intersect);

    VAS pred = subtract<segment_ops_wrt_a>(segs2, segs1);
    VAS trut { {trut_start_b, trut_end_b, trut_start_a, trut_end_a, {}, {}} };
    BOOST_CHECK_EQUAL(pred, trut);
  }
}

// test wrt b, then interchange a and b and test wrt a
void check_equal_both_ways(size_t seg1_start_a, size_t seg1_end_a, size_t seg1_start_b, size_t seg1_end_b,
                           size_t seg2_start_a, size_t seg2_end_a, size_t seg2_start_b, size_t seg2_end_b,
                           bool does_intersect,
                           size_t tru1_start_a, size_t tru1_end_a, size_t tru1_start_b, size_t tru1_end_b,
                           size_t tru2_start_a, size_t tru2_end_a, size_t tru2_start_b, size_t tru2_end_b)
{
  {
    VAS segs1 { {seg1_start_a, seg1_end_a, seg1_start_b, seg1_end_b, {}, {}} };
    VAS segs2 { {seg2_start_a, seg2_end_a, seg2_start_b, seg2_end_b, {}, {}} };
    BOOST_CHECK_EQUAL(intersects<segment_ops_wrt_b>(segs1, segs2), does_intersect);

    VAS pred = subtract<segment_ops_wrt_b>(segs2, segs1);
    VAS trut { {tru1_start_a, tru1_end_a, tru1_start_b, tru1_end_b, {}, {}},
               {tru2_start_a, tru2_end_a, tru2_start_b, tru2_end_b, {}, {}} };
    BOOST_CHECK_EQUAL(pred, trut);
  }

  {
    VAS segs1 { {seg1_start_b, seg1_end_b, seg1_start_a, seg1_end_a, {}, {}} };
    VAS segs2 { {seg2_start_b, seg2_end_b, seg2_start_a, seg2_end_a, {}, {}} };
    BOOST_CHECK_EQUAL(intersects<segment_ops_wrt_a>(segs1, segs2), does_intersect);

    VAS pred = subtract<segment_ops_wrt_a>(segs2, segs1);
    VAS trut { {tru1_start_b, tru1_end_b, tru1_start_a, tru1_end_a, {}, {}},
               {tru2_start_b, tru2_end_b, tru2_start_a, tru2_end_a, {}, {}} };
    BOOST_CHECK_EQUAL(pred, trut);
  }
}

BOOST_AUTO_TEST_CASE(forward_forward_exact)
{
  check_empty_both_ways(2, 4, 22, 24,
                        2, 4, 22, 24);
}

BOOST_AUTO_TEST_CASE(forward_forward_contained_1)
{
  check_empty_both_ways(3, 7, 22, 26,
                        2, 4, 23, 25);
}

BOOST_AUTO_TEST_CASE(forward_forward_contained_2)
{
  check_empty_both_ways(0, 0, 0, 0,
                        0, 0, 0, 0);
}

BOOST_AUTO_TEST_CASE(forward_forward_left_1)
{
  check_equal_both_ways(2, 4, 22, 24,
                        3, 5, 23, 25,
                        true,
                        5, 5, 25, 25);
}

BOOST_AUTO_TEST_CASE(forward_forward_left_2)
{
  check_equal_both_ways(2, 4, 22, 24,
                        4, 6, 24, 26,
                        true,
                        5, 6, 25, 26);
}

BOOST_AUTO_TEST_CASE(forward_forward_left_3)
{
  check_equal_both_ways(2, 4, 22, 24,
                        5, 7, 25, 27,
                        false,
                        5, 7, 25, 27);
}

BOOST_AUTO_TEST_CASE(forward_forward_right_1)
{
  check_equal_both_ways(2, 4, 22, 24,
                        1, 3, 21, 23,
                        true,
                        1, 1, 21, 21);
}

BOOST_AUTO_TEST_CASE(forward_forward_right_2)
{
  check_equal_both_ways(2, 4, 22, 24,
                        0, 2, 20, 22,
                        true,
                        0, 1, 20, 21);
}

BOOST_AUTO_TEST_CASE(forward_forward_right_3)
{
  check_equal_both_ways(3, 5, 23, 25,
                        0, 2, 20, 22,
                        false,
                        0, 2, 20, 22);
}

BOOST_AUTO_TEST_CASE(forward_forward_right_4_at_zero)
{
  check_equal_both_ways(3, 5, 1, 3,
                        0, 2, 0, 2,
                        true,
                        0, 0, 0, 0);
}

BOOST_AUTO_TEST_CASE(backward_forward_left_1)
{
  check_equal_both_ways(2, 4, 22, 24,
                        5, 3, 23, 25,
                        true,
                        3, 3, 25, 25);
}

BOOST_AUTO_TEST_CASE(backward_forward_left_2)
{
  check_equal_both_ways(2, 4, 22, 24,
                        6, 3, 23, 26,
                        true,
                        4, 3, 25, 26);
}

BOOST_AUTO_TEST_CASE(backward_forward_left_3)
{
  check_equal_both_ways(2, 4, 22, 24,
                        6, 3, 25, 28,
                        false,
                        6, 3, 25, 28);
}

BOOST_AUTO_TEST_CASE(backward_backward_left_1)
{
  check_equal_both_ways(2, 4, 22, 24,
                        5, 3, 25, 23,
                        true,
                        5, 5, 25, 25);
}

BOOST_AUTO_TEST_CASE(forward_backward_left_1)
{
  check_equal_both_ways(2, 4, 22, 24,
                        3, 5, 25, 23,
                        true,
                        3, 3, 25, 25);
}

BOOST_AUTO_TEST_CASE(forward_forward_split_1)
{
  check_equal_both_ways(100, 110,  1,  9,
                         50,  60,  0, 10,
                        true,
                         50,  50,  0,  0,
                         60,  60, 10, 10);
}

BOOST_AUTO_TEST_CASE(mismatches)
{
  VAS segs1 { {4, 6, 104, 106, {4, 5}, {104, 105}} };
  VAS segs2 { {0, 5, 100, 105, {3, 5}, {103, 105}} };
  BOOST_CHECK(intersects<segment_ops_wrt_b>(segs1, segs2));
  BOOST_CHECK(intersects<segment_ops_wrt_a>(segs1, segs2));

  VAS predb = subtract<segment_ops_wrt_b>(segs2, segs1);
  VAS trutb { {0, 3, 100, 103, {3},    {103}     } };
  BOOST_CHECK_EQUAL(predb, trutb);

  VAS preda = subtract<segment_ops_wrt_a>(segs2, segs1);
  VAS truta { {0, 3, 100, 103, {3},    {103}     } };
  BOOST_CHECK_EQUAL(preda, truta);
}

BOOST_AUTO_TEST_CASE(multiple_segs)
{
  VAS segs1 { { 2,  3, 102, 103, {},     {}        },    // 1a
              { 9, 10, 109, 110, {},     {}        },    // 1b
              {11, 20, 111, 120, {},     {}        },    // 1c
              {25, 26, 125, 126, {},     {}        } };  // 1d

  VAS segs2 { { 0,  2, 100, 102, {0},    {100}     },    // 2a: clipped by 1a
              { 3,  5, 103, 105, {3, 4}, {103, 104}},    // 2b: clipped by 1a
              {10, 11, 110, 111, {},     {}        },    // 2c: removed by 1b and 1c together
              {20, 30, 120, 130, {},     {}        } };  // 2d: clipped by 1c, then split by 1d

  BOOST_CHECK(intersects<segment_ops_wrt_b>(segs1, segs2));

  VAS trut  { { 0,  1, 100, 101, {0},    {100}     },    // 2a
              { 4,  5, 104, 105, {4},    {104}     },    // 2b
              /* empty */                                // 2c
              {21, 24, 121, 124, {},     {}        },    // 2d, part 1
              {27, 30, 127, 130, {},     {}        } };  // 2d, part 2

  VAS predb = subtract<segment_ops_wrt_b>(segs2, segs1);
  BOOST_CHECK_EQUAL(predb, trut);

  VAS preda = subtract<segment_ops_wrt_a>(segs2, segs1);
  BOOST_CHECK_EQUAL(preda, trut);
}
