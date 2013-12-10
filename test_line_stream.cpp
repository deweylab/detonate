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

#include <iostream>
#include <fstream>
#include <string>
#define BOOST_TEST_MODULE test_line_stream
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
