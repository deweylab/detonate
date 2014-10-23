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

#pragma once
#include <boost/lexical_cast.hpp>
#include <boost/array.hpp>

template<size_t num_fields, char sep = ','>
class lazycsv
{
public:
  lazycsv() {}
  lazycsv(const std::string& line) : line(line) { update_starts(); }

  // Get the value of the t'th field, converted to type T.
  template<typename T>
  inline T at(size_t t) const
  {
    size_t len = (starts[t+1] - 1) - starts[t]; // "- 1" so as not to include sep char
    return boost::lexical_cast<T>(line.substr(starts[t], len));
  }

  void parse_line(const std::string& line_)
  {
    line = line_;
    update_starts();
  }

  bool operator==(const lazycsv& other) const { return line == other.line; }
  bool operator!=(const lazycsv& other) const { return !(*this == other); }

private:
  std::string line;
  boost::array<size_t, num_fields + 1> starts;

  // Figure out where fields start (and end).
  //
  // Note: we record field start positions, not iterators pointing to those
  // positions, because it makes copying this object easier and more efficient.
  void update_starts()
  {
    #if 0
    // much slower for some reason
    starts[0] = 0;
    size_t f = 1;
    for (size_t i = 0; i < line.size(); ++i) {
      if (line[i] == sep) {
        starts[f] = i + 1;
        ++f;
      }
    }
    #else
    starts[0] = 0;
    size_t f = 1;
    std::string::size_type i = 0;
    for (i = line.find(sep, 0);
         i != std::string::npos;
         i = line.find(sep, i + 1)) {
      starts[f] = i+1;
      ++f;
    }
    #endif

    starts[num_fields] = line.size() + 1; // "+ 1" to imitate sep char
    if (f != num_fields)
      throw std::runtime_error("Invalid number of fields in line: '" + line + "'");
  }
};
