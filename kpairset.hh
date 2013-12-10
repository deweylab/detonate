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
#include <set>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

class kpairset
{
public:
  kpairset(size_t k, size_t top) : k(k), top(top), sz(0), data(top+1-k, false) {}

  size_t size() { return sz; }

  template<typename ConstIterator>
  void add_kpairs_with_exceptions(size_t lo, size_t hi, ConstIterator exceptions_begin, ConstIterator exceptions_end)
  {
    if (lo > hi)
      std::swap(lo, hi);

    // If this segment is shorter than the kmer length k, then we should skip
    // adding it.
    // 
    // E.g., k=3, lo=4, hi=6: 456  -> ok,     6-4+1 = 3 >= 3
    //       k=3, lo=4, hi=5: 45   -> not ok, 5-4+1 = 2 < 3
    //       k=4, lo=4, hi=6: 456  -> not ok, 6-4+1 = 3 < 4
    //       k=3, lo=4, hi=7: 4567 -> ok,     7-4+1 = 4 >= 3
    // So the requirement is that hi-lo+1 >= k, i.e., hi+1-lo >= k.
    if (hi+1-lo < k)
      return;

    std::set<size_t> exceptions;
    for (ConstIterator x = exceptions_begin; x != exceptions_end; ++x)
      exceptions.insert(*x);

    // E.g., k=3, lo=10, hi=15, then we can add (10,12), (11,13), (12,14), (13,15).
    // Note 13 = 15-3+1, i.e., hi-k+1, or (for size_t) hi+1-k.
    // Note 12 = 10+3-1, i.e, i+k-1.
    for (size_t i = lo; i <= hi+1-k; ++i) {
      if (i >= top+1-k)
        throw std::runtime_error("Out of bounds.");
      if (exceptions.count(i) == 0 && exceptions.count(i+k-1) == 0) {
        if (data[i] == false) {
          ++sz;
          data[i] = true;
        }
      }
    }
  }

private:
  size_t k, top, sz;
  std::vector<bool> data;
};
