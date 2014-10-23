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

inline bool is_N(char c)
{
  return c == 'N' || c == 'n';
}

// Precondition: If at_beginning is true, then there are no preconditions. If
// at_beginning is false, then we assume that start[0:(k-2)] are all non-N,
// where k is the kmer length. In other words, in the initial kmer starting at
// start, the first (k-1) bases are all non-N, and only the kth base might be
// N.
//
// Other input: end points at one past the last valid start position of a kmer
// within the underlying string. kmerlen is the kmer length, k.
//
// Postcondition: The returned pointer, ret, points at the first kmer found >=
// start such that start[0:(k-1)] are all non-N, or ret is end if there is no
// such kmer.
const char *skip_Ns(const char *start, const char *end, size_t kmerlen,
                    bool at_beginning)
{
  const char *ret;

  // Step 1: Fast check (just check the last base of the kmer).
  //
  // If we are NOT at the beginning of the contig, then we can assume that all
  // but the last base of this kmer have been checked and found to contain only
  // non-N characters. (This would have been done by previous calls to
  // skip_Ns.) Hence, we can start with a fast check of just the last base.
  if (!at_beginning) {
    // If the last base of this kmer is non-N, then we don't need to do anything
    // special.
    if (!is_N(start[kmerlen - 1]))
      return start;

    // Otherwise, skip past this kmer, which has N in its last base, at position
    // k - 1.
    ret = start + kmerlen;
  }
  // If we ARE at the beginning of the contig, then we can't make any
  // assumptions whether or not there are Ns inside the initial kmer, since we
  // haven't checked these bases yet. So we don't do anything in this fast
  // step.
  else {
    ret = start;
  }

  // Step 2: Slow check (check each base of the kmer, moving the kmer forward
  // until each base is non-N).
  check_that_the_kmer_starting_at_ret_is_valid:

  // Make sure ret is still a valid place to start a kmer.
  if (ret >= end)
    return end;

  // We cannot assume that the kmer starting at ret contains non-Ns at any of
  // its positions. So, check whether each position contains an N. If so, jump
  // past the N and try again.
  for (size_t l = 0; l < kmerlen; ++l) {
    if (is_N(ret[l])) {
      ret = ret + l + 1;
      goto check_that_the_kmer_starting_at_ret_is_valid;
    }
  }

  return ret;
}

