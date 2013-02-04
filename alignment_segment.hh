#pragma once
#include <vector>
#include <iostream>
#include <boost/foreach.hpp>

struct alignment_segment
{
  size_t a_start, a_end, b_start, b_end; // the segment is [a_start, a_end] -> [b_start, b_end]
  std::vector<size_t> a_mismatches, b_mismatches;
};

std::ostream& operator<<(std::ostream& out, const alignment_segment& seg)
{
  out << "(a_start=" << seg.a_start
      << ", a_end=" << seg.a_end
      << ", b_start=" << seg.b_start
      << ", b_end=" << seg.b_end
      << ", a_mismatches=( ";
  BOOST_FOREACH(size_t x, seg.a_mismatches)
    out << x << " ";
  out << "), b_mismatches=( ";
  BOOST_FOREACH(size_t x, seg.b_mismatches)
    out << x << " ";
  out << ")";
  return out;
}
