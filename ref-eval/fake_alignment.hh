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
#include <string>
#include "line_stream.hh"
#include "alignment_segment.hh"

namespace detail
{
  class fake_alignment_input_stream;
}

////////////////////////////////////////////////////////////////////////////
// fake_alignment represents a single fake alignment, used for testing. It
// realizes the Alignment concept.
////////////////////////////////////////////////////////////////////////////

class fake_alignment
{
public:
  // Realization of Alignment concept
  void parse_line(const std::string& /*line*/) { throw std::runtime_error("Not implemented."); }
  std::string a_name()        const { return a_name_; }
  std::string b_name()        const { return b_name_; }
  double      frac_identity() const { return frac_identity_; }
  double      frac_indel()    const { return frac_indel_; }
  typedef detail::fake_alignment_input_stream   input_stream_type;
  typedef const std::vector<alignment_segment>& segments_type;
  segments_type segments(const std::string& /*a*/, const std::string& /*b*/) const { return alignment_segments; }

  // Data
  std::string a_name_, b_name_;
  double frac_identity_, frac_indel_;
  std::vector<alignment_segment> alignment_segments;
};

std::ostream& operator<<(std::ostream& out, const fake_alignment& al)
{
  out << "Alignment " << al.a_name() << " -> " << al.b_name()
      << ", frac_identity=" << al.frac_identity() << ", frac_indel=" << al.frac_indel() 
      << ", alignment_segments={ ";
  BOOST_FOREACH(const alignment_segment& seg, al.alignment_segments)
    out << "seg=" << seg << " ";
  out << "}";
  return out;
}

namespace detail
{
  class fake_alignment_input_stream
  {
  public:
    fake_alignment_input_stream(const std::vector<fake_alignment>& fake_alignments) : fake_alignments(fake_alignments), i(0), done(false) {}

    fake_alignment_input_stream& operator>>(fake_alignment& al)
    {
      if (i < fake_alignments.size())
        al = fake_alignments[i];
      else
        done = true;
      ++i;
      return *this;
    }

    inline operator bool() { return !done; }
    inline bool operator!() { return done; }

    std::vector<fake_alignment> fake_alignments;
    size_t i;
    bool done;
  };
} // namespace detail
