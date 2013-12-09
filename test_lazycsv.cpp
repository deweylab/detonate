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
