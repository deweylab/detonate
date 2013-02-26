#include <iostream>
#include <fstream>
#include <string>
#define BOOST_TEST_MODULE test_line_stream
#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include "line_stream.hh"

using namespace std;

BOOST_AUTO_TEST_CASE(realistic_usage)
{
  string s = "line 1\nline 2\nline 3";
  const char *results[] = {"line 1", "line 2", "line 3"};
  boost::shared_ptr<istringstream> iss(new istringstream(s));
  line_stream ls(iss);
  string line;
  size_t i = 0;
  while (ls >> line)
    BOOST_CHECK_EQUAL(line, results[i++]);
  BOOST_CHECK_EQUAL(i, 3UL);
}

BOOST_AUTO_TEST_CASE(basic_no_trailing_eol)
{
  string s = "line 1\nline 2\nline 3";
  boost::shared_ptr<istringstream> iss(new istringstream(s));
  line_stream ls(iss);
  string line;
  ls >> line;
  BOOST_CHECK_EQUAL(line, "line 1");
  ls >> line;
  BOOST_CHECK_EQUAL(line, "line 2");
  ls >> line;
  BOOST_CHECK_EQUAL(line, "line 3");
  ls >> line;
  BOOST_CHECK_EQUAL(static_cast<bool>(ls), false);
  BOOST_CHECK_EQUAL(!ls, true);
}

BOOST_AUTO_TEST_CASE(basic_with_trailing_eol)
{
  string s = "line 1\nline 2\nline 3\n";
  boost::shared_ptr<istringstream> iss(new istringstream(s));
  line_stream ls(iss);
  string line;
  ls >> line;
  BOOST_CHECK_EQUAL(line, "line 1");
  ls >> line;
  BOOST_CHECK_EQUAL(line, "line 2");
  ls >> line;
  BOOST_CHECK_EQUAL(line, "line 3");
  ls >> line;
  BOOST_CHECK_EQUAL(static_cast<bool>(ls), false);
  BOOST_CHECK_EQUAL(!ls, true);
}

BOOST_AUTO_TEST_CASE(no_lines)
{
  string s = "";
  boost::shared_ptr<istringstream> iss(new istringstream(s));
  line_stream ls(iss);
  string line;
  ls >> line;
  BOOST_CHECK_EQUAL(static_cast<bool>(ls), false);
  BOOST_CHECK_EQUAL(!ls, true);
}

BOOST_AUTO_TEST_CASE(just_newline)
{
  string s = "\n";
  boost::shared_ptr<istringstream> iss(new istringstream(s));
  line_stream ls(iss);
  string line;
  ls >> line;
  BOOST_CHECK_EQUAL(line, "");
  ls >> line;
  BOOST_CHECK_EQUAL(static_cast<bool>(ls), false);
  BOOST_CHECK_EQUAL(!ls, true);
}
