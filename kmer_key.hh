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
#include "city.h"

typedef const char * kmer_key;

struct kmer_key_hash
{
  size_t kmerlen;

  kmer_key_hash(size_t kmerlen)
  : kmerlen(kmerlen)
  {}

  size_t operator()(const char *k) const
  {
    return CityHash64(k, kmerlen);
  }
};

struct kmer_key_equal_to
{
  size_t kmerlen;

  kmer_key_equal_to(size_t kmerlen)
  : kmerlen(kmerlen)
  {}

  bool operator()(const char *lhs, const char *rhs) const
  {
    for (size_t i = 0; i < kmerlen; ++i, ++lhs, ++rhs)
      if (*lhs != *rhs)
        return false;
    return true;
  }
};
