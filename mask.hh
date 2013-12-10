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
#include <vector>
#include <set>
#include <boost/foreach.hpp>

class mask
{
public:
  mask(size_t length) : data(length, false), num_1s(0) {}

  // Adds the elements [start, end] - exceptions to the mask
  template<typename ConstIterator>
  void add_interval_with_exceptions(size_t start, size_t end, ConstIterator exceptions_begin, ConstIterator exceptions_end)
  {
    if (start > end)
      std::swap(start, end);
    std::set<size_t> exceptions(exceptions_begin, exceptions_end);
    for (size_t i = start; i <= end; ++i) {
      if (exceptions.count(i) == 0) {
        if (!data[i]) // excluded -> included => increment num ones
          ++num_1s;
        data[i] = true;
      }
    }
  }

  size_t num_ones() { return num_1s; }

private:
  std::vector<bool> data;
  size_t num_1s;
};
