#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "deweylab/bio/formats/fasta.hh"

// need to declare "clock_t start" before tic
#define tic start = clock();
#define toc std::cerr << "Done in " << 1.0*(clock() - start)/CLOCKS_PER_SEC << " seconds." << std::endl;

void open_or_throw(std::ifstream& ifs, const std::string& filename)
{
  ifs.clear();
  ifs.open(filename.c_str());
  if (!ifs)
    throw std::runtime_error("Could not open file '" + filename + "'.");
}

void read_fasta_names_and_seqs(const std::string& filename,
                               std::vector<std::string>& seqs,
                               std::vector<std::string>& names,
                               std::map<std::string, size_t>& names_to_idxs)
{
  std::ifstream ifs;
  open_or_throw(ifs, filename);
  deweylab::bio::formats::fasta::InputStream is(ifs);
  deweylab::bio::formats::fasta::Record rec;
  size_t idx = 0;
  while (is >> rec) {
    assert(names_to_idxs.count(rec.id) == 0);
    names_to_idxs[rec.id] = idx;
    names.push_back(rec.id);
    seqs.push_back(rec.sequence);
    assert(names[names_to_idxs[rec.id]] == rec.id);
    ++idx;
  }
}

void read_transcript_expression(const std::string& filename,
                                std::vector<double>& expr,
                                const std::map<std::string, size_t>& names_to_idxs)
{
  std::ifstream ifs;
  open_or_throw(ifs, filename);
  std::string line;
  std::vector<std::string> col(8);
  { // header
    getline(ifs, line);
    std::stringstream ss(line);
    for (size_t i = 0; i < 8; ++i)
      getline(ss, col[i], '\t');
    assert(col[0] == "transcript_id");
    assert(col[5] == "TPM");
  }
  while (getline(ifs, line)) {
    std::stringstream ss(line);
    for (size_t i = 0; i < 8; ++i)
      getline(ss, col[i], '\t');
    assert(names_to_idxs.count(col[0]) != 0);
    size_t idx = names_to_idxs.find(col[0])->second;
    expr[idx] = boost::lexical_cast<double>(col[5]) / 1000000.0;
  }
}

void compute_nucl_expression(const std::vector<std::string>& A,
                             const std::vector<double>& tau_A,
                             std::vector<double>& nu_A)
{
  size_t n = tau_A.size();
  assert(tau_A.size() == n);
  assert(nu_A.size() == n);

  for (size_t i = 0; i < n; ++i)
    nu_A[i] = tau_A[i] * A[i].size();

  double denom = 0;
  for (size_t i = 0; i < n; ++i)
    denom += nu_A[i];

  for (size_t i = 0; i < n; ++i)
    nu_A[i] /= denom;
}

std::vector<size_t> make_random_permutation(size_t n)
{
  std::vector<size_t> x(n);
  for (size_t i = 0; i < n; ++i)
    x[i] = i;
  boost::random::random_device rng;
  boost::random::uniform_int_distribution<> index_dist(0, n - 1);
  for (size_t k = 0; k < n*10; ++k) {
    size_t i = index_dist(rng);
    size_t j = index_dist(rng);
    std::swap(x[i], x[j]);
  }
  return x;
}
