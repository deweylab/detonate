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
#define BOOST_TEST_MODULE test_lazycsv
#include <boost/test/unit_test.hpp>
#include "lazycsv.hh"

using namespace std;

BOOST_AUTO_TEST_CASE(basic)
{
  string s = "one\ttwo\tthree";
  lazycsv<3,'\t'> l(s);
  BOOST_CHECK_EQUAL(l.at<string>(0), string("one"));
  BOOST_CHECK_EQUAL(l.at<string>(1), string("two"));
  BOOST_CHECK_EQUAL(l.at<string>(2), string("three"));
}

BOOST_AUTO_TEST_CASE(basic_nl)
{
  string s = "one\ttwo\tthree\n";
  lazycsv<3,'\t'> l(s);
  BOOST_CHECK_EQUAL(l.at<string>(0), string("one"));
  BOOST_CHECK_EQUAL(l.at<string>(1), string("two"));
  BOOST_CHECK_EQUAL(l.at<string>(2), string("three\n"));
}

BOOST_AUTO_TEST_CASE(basic_numbers)
{
  string s = "1\t2.5\t-3";
  lazycsv<3,'\t'> l(s);
  BOOST_CHECK_EQUAL(l.at<size_t>(0), size_t(1));
  BOOST_CHECK_EQUAL(l.at<double>(1), 2.5);
  BOOST_CHECK_EQUAL(l.at<int   >(2), -3);
}

BOOST_AUTO_TEST_CASE(basic_empty_first)
{
  string s = "\ttwo\tthree";
  lazycsv<3,'\t'> l(s);
  BOOST_CHECK_EQUAL(l.at<string>(0), string(""));
  BOOST_CHECK_EQUAL(l.at<string>(1), string("two"));
  BOOST_CHECK_EQUAL(l.at<string>(2), string("three"));
}

BOOST_AUTO_TEST_CASE(basic_empty_second)
{
  string s = "one\t\tthree";
  lazycsv<3,'\t'> l(s);
  BOOST_CHECK_EQUAL(l.at<string>(0), string("one"));
  BOOST_CHECK_EQUAL(l.at<string>(1), string(""));
  BOOST_CHECK_EQUAL(l.at<string>(2), string("three"));
}

BOOST_AUTO_TEST_CASE(basic_empty_last)
{
  string s = "one\ttwo\t";
  lazycsv<3,'\t'> l(s);
  BOOST_CHECK_EQUAL(l.at<string>(0), string("one"));
  BOOST_CHECK_EQUAL(l.at<string>(1), string("two"));
  BOOST_CHECK_EQUAL(l.at<string>(2), string(""));
}

BOOST_AUTO_TEST_CASE(basic_empty_all)
{
  string s = "\t\t\t";
  lazycsv<4,'\t'> l(s);
  BOOST_CHECK_EQUAL(l.at<string>(0), string(""));
  BOOST_CHECK_EQUAL(l.at<string>(1), string(""));
  BOOST_CHECK_EQUAL(l.at<string>(2), string(""));
  BOOST_CHECK_EQUAL(l.at<string>(3), string(""));
}

BOOST_AUTO_TEST_CASE(basic_comma)
{
  string s = "one,two,three,four";
  lazycsv<4> l(s);
  BOOST_CHECK_EQUAL(l.at<string>(0), string("one"));
  BOOST_CHECK_EQUAL(l.at<string>(1), string("two"));
  BOOST_CHECK_EQUAL(l.at<string>(2), string("three"));
  BOOST_CHECK_EQUAL(l.at<string>(3), string("four"));
}
