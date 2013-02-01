#define BOOST_TEST_MODULE test_pairset
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "mask.hh"

BOOST_AUTO_TEST_CASE(simple)
{
  mask m(5);
  BOOST_CHECK_EQUAL(m.num_ones(), 0ul);
  m.add_interval(1, 1);
  BOOST_CHECK_EQUAL(m.num_ones(), 1ul);
  std::vector<size_t> v;
  v.push_back(1);
  m.remove(v);
  BOOST_CHECK_EQUAL(m.num_ones(), 0ul);
  m.add_interval(4, 4);
  BOOST_CHECK_EQUAL(m.num_ones(), 1ul);
  m.add_interval(0, 4);
  BOOST_CHECK_EQUAL(m.num_ones(), 5ul);
}

BOOST_AUTO_TEST_CASE(crossover)
{
  mask m(100);
  BOOST_CHECK_EQUAL(m.num_ones(), 0ul);
  m.add_interval(0, 50-1);
  BOOST_CHECK_EQUAL(m.num_ones(), 50ul);
  m.add_interval(40, 60-1);
  BOOST_CHECK_EQUAL(m.num_ones(), 60ul);
  std::vector<size_t> v;
  v.push_back(35);
  v.push_back(45);
  m.remove(v);
  BOOST_CHECK_EQUAL(m.num_ones(), 58ul);
  m.add_interval(40, 60-1);
  BOOST_CHECK_EQUAL(m.num_ones(), 59ul);
}
