#define BOOST_TEST_MODULE test_mask
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include "mask.hh"

static size_t no_exceptions[0] = {};

BOOST_AUTO_TEST_CASE(simple)
{
  mask m(5);
  BOOST_CHECK_EQUAL(m.num_ones(), 0ul);

  m.add_interval_with_exceptions(1, 1, no_exceptions, no_exceptions);
  BOOST_CHECK_EQUAL(m.num_ones(), 1ul);

  size_t one_exception[1] = {3};
  m.add_interval_with_exceptions(3, 3, one_exception, one_exception+1);
  BOOST_CHECK_EQUAL(m.num_ones(), 1ul);

  m.add_interval_with_exceptions(2, 4, one_exception, one_exception+1);
  BOOST_CHECK_EQUAL(m.num_ones(), 3ul);

  m.add_interval_with_exceptions(0, 4, no_exceptions, no_exceptions);
  BOOST_CHECK_EQUAL(m.num_ones(), 5ul);
}

BOOST_AUTO_TEST_CASE(crossover)
{
  mask m(100);
  BOOST_CHECK_EQUAL(m.num_ones(), 0ul);

  m.add_interval_with_exceptions(0, 50-1, no_exceptions, no_exceptions);
  BOOST_CHECK_EQUAL(m.num_ones(), 50ul);

  m.add_interval_with_exceptions(70, 80-1, no_exceptions, no_exceptions);
  BOOST_CHECK_EQUAL(m.num_ones(), 60ul);

  m.add_interval_with_exceptions(40, 75-1, no_exceptions, no_exceptions);
  BOOST_CHECK_EQUAL(m.num_ones(), 80ul);
}
