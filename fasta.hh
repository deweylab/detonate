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
