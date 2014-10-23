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

#define BOOST_TEST_MODULE test_re_kc
#include <boost/test/unit_test.hpp>
#include "re_kc.hh"

using namespace re::kc;

BOOST_AUTO_TEST_CASE(example_1)
{
  //                 012345678
  const char *seq = "ACGNGCAT";
  size_t k = 3;
  const char *end = seq + 5;
  BOOST_CHECK_EQUAL(skip_Ns(seq + 0, end, k, false), seq + 0); // ACG -> ok
  BOOST_CHECK_EQUAL(skip_Ns(seq + 1, end, k, false), seq + 4); // CGN -> skip to GCA
  BOOST_CHECK_EQUAL(skip_Ns(seq + 4, end, k, false), seq + 4); // CGA -> ok
  BOOST_CHECK_EQUAL(skip_Ns(seq + 5, end, k, false), seq + 5); // CGA -> ok
}

BOOST_AUTO_TEST_CASE(example_2)
{
  //                 012345678
  const char *seq = "ACGNGCNT";
  size_t k = 3;
  const char *end = seq + 5;
  BOOST_CHECK_EQUAL(skip_Ns(seq + 0, end, k, false), seq + 0); // ACG -> ok
  BOOST_CHECK_EQUAL(skip_Ns(seq + 1, end, k, false), end);     // CGN -> skip to GCN, which contains N, hence to end
}

BOOST_AUTO_TEST_CASE(example_3)
{
  //                 012345678
  const char *seq = "ACGNGNAT";
  size_t k = 3;
  const char *end = seq + 5;
  BOOST_CHECK_EQUAL(skip_Ns(seq + 0, end, k, false), seq + 0); // ACG -> ok
  BOOST_CHECK_EQUAL(skip_Ns(seq + 1, end, k, false), end);     // CGN -> skip to GNA, which contains N, hence to end
}

BOOST_AUTO_TEST_CASE(example_4)
{
  //                 012345678
  const char *seq = "ACGNGNATT";
  size_t k = 3;
  const char *end = seq + 6;
  BOOST_CHECK_EQUAL(skip_Ns(seq + 0, end, k, false), seq + 0); // ACG -> ok
  BOOST_CHECK_EQUAL(skip_Ns(seq + 1, end, k, false), seq + 6); // CGN -> skip to GCN, which contains N, hence to ATT
  BOOST_CHECK_EQUAL(skip_Ns(seq + 6, end, k, false), seq + 6); // ATT -> ok
}

// There is only one N in this sequence, and it occurs within 76 bases of the
// end, so either skip_Ns should return the current position, or it should skip
// past that N, which will result in us being at the end of the line.
BOOST_AUTO_TEST_CASE(example_from_oases_mouse_chr1_cc_0)
{
  const char *seq = "ATAACAGTTACTTGTCTGGGTGTGTGTTAGGGGTCTCCGGTATGTCCTCGTATACCTAGGTTATGGGGGGGGGAGGGGCAAAGCAATGAGAAACATAACTAAATTTATGTGTAAGTTGGGGGTGTCANGAAAAAGCAACTTAG";
  size_t k = 76;
  size_t len = strlen(seq);
  const char *end = seq + len + 1 - k;
  for (size_t i = 0; i + k - 1 < len; ++i) {
    if (seq[i + k - 1] != 'N')
      BOOST_CHECK_EQUAL(skip_Ns(seq + i, end, k, i == 0), seq + i);
    else
      BOOST_CHECK_EQUAL(skip_Ns(seq + i, end, k, i == 0), end);
  }
}

// There is only one N in this sequence, and it occurs within the first 76
// bases of the sequence. So (assuming we're at the beginning of the sequence),
// skip_Ns should skip past that N, and from then on it should return the
// current position. (If we are not at the beginning of the sequence, then that
// N could not exist inside the first 76 bases of the sequence, so skip_Ns will
// ignore it.)
BOOST_AUTO_TEST_CASE(example_from_oases_mouse_chr1_cc_0_reverse_complemented)
{
  const char *seq = "CTAAGTTGCTTTTTCNTGACACCCCCAACTTACACATAAATTTAGTTATGTTTCTCATTGCTTTGCCCCTCCCCCCCCCATAACCTAGGTATACGAGGACATACCGGAGACCCCTAACACACACCCAGACAAGTAACTGTTAT";
  size_t k = 76;
  size_t len = strlen(seq);
  const char *end = seq + len + 1 - k;
  size_t after_N = 0;
  for (size_t i = 0; i < len; ++i) {
    if (seq[i] == 'N') {
      after_N = i + 1;
      break;
    }
  }
  // Check, assuming that we're at the beginning (each time).
  for (size_t i = 0; i + k - 1 < len; ++i) {
    if (i <= after_N)
      BOOST_CHECK_EQUAL(skip_Ns(seq + i, end, k, true), seq + after_N);
    else
      BOOST_CHECK_EQUAL(skip_Ns(seq + i, end, k, true), seq + i);
  }
  // Check, assuming that we are not at the beginning (each time).
  for (size_t i = 0; i + k - 1 < len; ++i)
    BOOST_CHECK_EQUAL(skip_Ns(seq + i, end, k, false), seq + i);
}
