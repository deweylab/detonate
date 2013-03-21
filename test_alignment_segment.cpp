#define BOOST_TEST_MODULE test_alignment_segment
#define BOOST_TEST_DYN_LINK
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

