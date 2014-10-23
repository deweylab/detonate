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
#include "util.hh"
#include "deweylab/bio/formats/fasta.hh"

struct fasta
{
  size_t card;
  std::vector<std::string> seqs;
  std::vector<std::string> names;
  std::map<std::string, size_t> names_to_idxs;
  std::vector<size_t> lengths;
};

void read_fasta(fasta& fa, const std::string& filename)
{
  try {
    boost::shared_ptr<std::ifstream> ifs = open_or_throw(filename);
    deweylab::bio::formats::fasta::InputStream is(*ifs);
    deweylab::bio::formats::fasta::Record rec;
    size_t idx = 0;
    while (is >> rec) {
      if (fa.names_to_idxs.count(rec.id) != 0)
        throw std::runtime_error("Found duplicate sequence id in " + filename + ".");
      fa.names_to_idxs[rec.id] = idx;
      fa.names  .push_back(rec.id);
      fa.seqs   .push_back(rec.sequence);
      fa.lengths.push_back(rec.sequence.size());
      assert(fa.names[fa.names_to_idxs[rec.id]] == rec.id);
      ++idx;
    }
    fa.card = fa.seqs.size();
  } catch (const std::runtime_error& x) {
    throw std::runtime_error("Can't parse " + filename + ": " + x.what());
  }
}
