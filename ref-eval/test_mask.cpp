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

#define BOOST_TEST_MODULE test_mask
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

BOOST_AUTO_TEST_CASE(reverse)
{
  mask f(100), r(100);
  BOOST_CHECK_EQUAL(f.num_ones(), 0ul);
  BOOST_CHECK_EQUAL(r.num_ones(), 0ul);

  f.add_interval_with_exceptions(0,   50-1, no_exceptions, no_exceptions);
  r.add_interval_with_exceptions(50-1,   0, no_exceptions, no_exceptions);
  BOOST_CHECK_EQUAL(f.num_ones(), 50ul);
  BOOST_CHECK_EQUAL(r.num_ones(), 50ul);

  size_t two_exceptions[2] = {78, 80};
  f.add_interval_with_exceptions(75,   90-1, two_exceptions, two_exceptions+2);
  r.add_interval_with_exceptions(90-1,   75, two_exceptions, two_exceptions+2);
  BOOST_CHECK_EQUAL(f.num_ones(), 50ul + 13);
  BOOST_CHECK_EQUAL(r.num_ones(), 50ul + 13);
}
