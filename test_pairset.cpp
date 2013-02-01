#define BOOST_TEST_MODULE test_pairset
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "pairset.hh"

typedef boost::mpl::list<brute_force_pairset, big_matrix_pairset> test_types;

size_t choose_2(size_t n) { return n*(n+1)/2; }

BOOST_AUTO_TEST_CASE_TEMPLATE(simple, PairsetType, test_types)
{
  PairsetType ps(5);
  BOOST_CHECK_EQUAL(ps.size(), 0ul);
  ps.add_square(1, 1);
  BOOST_CHECK_EQUAL(ps.size(), 1ul);
  ps.remove_pair(1, 1);
  BOOST_CHECK_EQUAL(ps.size(), 0ul);
  ps.add_square(4, 4);
  BOOST_CHECK_EQUAL(ps.size(), 1ul);
  ps.add_square(0, 4);
  BOOST_CHECK_EQUAL(ps.size(), choose_2(5));
}

BOOST_AUTO_TEST_CASE_TEMPLATE(crossovers, PairsetType, test_types)
{
  PairsetType ps(100);
  BOOST_CHECK_EQUAL(ps.size(), 0ul);
  ps.add_square(0, 50-1);
  BOOST_CHECK_EQUAL(ps.size(), choose_2(50));
  ps.add_square(40, 60-1);
  BOOST_CHECK_EQUAL(ps.size(), choose_2(50) + choose_2(20) - choose_2(10));
  ps.remove_pair(35, 35);
  BOOST_CHECK_EQUAL(ps.size(), choose_2(50) + choose_2(20) - choose_2(10) - 1);
  ps.remove_pair(45, 45);
  BOOST_CHECK_EQUAL(ps.size(), choose_2(50) + choose_2(20) - choose_2(10) - 2);
  std::vector<size_t> x;
  x.push_back(35);
  x.push_back(46);
  ps.remove_all_pairs(x); // (35,35) [already gone], (35,46), (45,46)
  BOOST_CHECK_EQUAL(ps.size(), choose_2(50) + choose_2(20) - choose_2(10) - 2 - 2);
  ps.add_square(40, 60-1); // restore (45,46) and (46,46)
  BOOST_CHECK_EQUAL(ps.size(), choose_2(50) + choose_2(20) - choose_2(10) - 1 - 1);
}
