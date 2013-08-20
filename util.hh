#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "deweylab/bio/formats/fasta.hh"

// need to declare "clock_t start" before tic
#define tic start = clock();
#define toc std::cerr << "Done in " << 1.0*(clock() - start)/CLOCKS_PER_SEC << " seconds." << std::endl;

boost::shared_ptr<std::ifstream> open_or_throw(const std::string& filename)
{
  boost::shared_ptr<std::ifstream> ifs(new std::ifstream(filename.c_str()));
  if (!*ifs)
    throw std::runtime_error("Could not open file '" + filename + "'.");
  return ifs;
}

void read_fasta_names_and_seqs(const std::string& filename,
                               std::vector<std::string>& seqs,
                               std::vector<std::string>& names,
                               std::map<std::string, size_t>& names_to_idxs)
{
  boost::shared_ptr<std::ifstream> ifs = open_or_throw(filename);
  deweylab::bio::formats::fasta::InputStream is(*ifs);
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
  boost::shared_ptr<std::ifstream> ifs = open_or_throw(filename);
  std::string line;
  std::vector<std::string> col(8);
  { // header
    getline(*ifs, line);
    std::stringstream ss(line);
    for (size_t i = 0; i < 8; ++i)
      getline(ss, col[i], '\t');
    assert(col[0] == "transcript_id");
    assert(col[5] == "TPM");
  }
  while (getline(*ifs, line)) {
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
  //boost::random::random_device rng;
  boost::random::mt19937 rng;
  boost::random::uniform_int_distribution<> index_dist(0, n - 1);
  for (size_t k = 0; k < n*10; ++k) {
    size_t i = index_dist(rng);
    size_t j = index_dist(rng);
    std::swap(x[i], x[j]);
  }
  return x;
}

char complement(char c)
{
  switch (c)
  {
    case 'A': return 'T';
    case 'T': return 'A';
    case 'C': return 'G';
    case 'G': return 'C';
    case 'N': return 'N';

    case 'a': return 't';
    case 't': return 'a';
    case 'c': return 'g';
    case 'g': return 'c';
    case 'n': return 'n';

    case 'S': return 'S'; // compl(S) = compl({C,G}) = {G,C} = S
    case 'W': return 'W'; // compl(W) = compl({A,T}) = {T,A} = W
    case 'M': return 'K'; // compl(M) = compl({A,C}) = {T,G} = K
    case 'K': return 'M'; // compl(K) = compl({G,T}) = {C,A} = M
    case 'Y': return 'R'; // compl(Y) = compl({C,T}) = {G,A} = R
    case 'R': return 'Y'; // compl(R) = compl({A,G}) = {T,C} = Y
    case 'B': return 'V'; // compl(B) = compl({C,G,T}) = {G,C,A} = V
    case 'V': return 'B'; // compl(V) = compl({A,C,G}) = {T,G,C} = B
    case 'D': return 'H'; // compl(D) = compl({A,G,T}) = {T,C,A} = H
    case 'H': return 'D'; // compl(H) = compl({A,C,T}) = {T,G,A} = D

    case 's': return 's';
    case 'w': return 'w';
    case 'm': return 'k';
    case 'k': return 'm';
    case 'y': return 'r';
    case 'r': return 'y';
    case 'b': return 'v';
    case 'v': return 'b';
    case 'd': return 'h';
    case 'h': return 'd';

    default:
      throw std::runtime_error("Cannot complement invalid nucleotide.");
  }
}

std::string reverse_complement(const std::string& x)
{
  std::string y(x);
  transform(x.begin(), x.end(), y.begin(), complement);
  reverse(y.begin(), y.end());
  return y;
}

