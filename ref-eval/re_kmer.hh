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
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <boost/foreach.hpp>
#include <sparsehash/sparse_hash_map>
#include <sparsehash/dense_hash_map>
#include "skip_Ns.hh"
#include "util.hh"
#include "kmer_key.hh"

namespace re {
namespace kmer {

struct kmer_info
{
  double weights[2]; // sometimes normalized, sometimes not
  kmer_info()
  {
    weights[0] = 0;
    weights[1] = 0;
  }
  double& weight_in_A() { return weights[0]; }
  double& weight_in_B() { return weights[1]; }
  const double& weight_in_A() const { return weights[0]; }
  const double& weight_in_B() const { return weights[1]; }
};

typedef google::sparse_hash_map<kmer_key, kmer_info, kmer_key_hash, kmer_key_equal_to> sparse_kmer_map;
typedef google::dense_hash_map<kmer_key, kmer_info, kmer_key_hash, kmer_key_equal_to> dense_kmer_map;

template<typename Ht>
struct empty_key_initializer
{
  empty_key_initializer(Ht&, size_t)
  {}
};

template<>
struct empty_key_initializer<dense_kmer_map>
{
  std::string empty_string;
  const char *empty_key;

  empty_key_initializer(dense_kmer_map& ht, size_t kmerlen)
  : empty_string(kmerlen, ' '),
    empty_key(empty_string.c_str())
  {
    ht.set_empty_key(empty_key);
  }
};

// A_or_B is 0 if we're counting A's kmers, 1 if B's
template<typename Ht, size_t A_or_B>
void count_kmers(
    Ht& ht,
    const std::vector<std::string>& A,
    const std::vector<std::string>& A_rc,
    const std::vector<double>& tau_A,
    size_t kmerlen,
    bool strand_specific)
{
  // For each contig a in A:
  //   For each kmer r in a or reverse_complement(a):
  //     Add (1/2)*[1/(length(a) - k + 1)]*tau_A(a) to count_A(r)
  //   For each kmer r in reverse_complement(a):
  //     Add (1/2)*[1/(length(a) - k + 1)]*tau_A(a) to count_A(r)
  size_t num_strands = strand_specific ? 1 : 2;
  for (size_t i = 0; i < A.size(); ++i) {
    for (size_t which = 0; which < num_strands; ++which) {
      const std::string& a = which == 0 ? A[i] : A_rc[i];
      if (a.size() >= kmerlen) {
        double c = tau_A[i];
        const char *beg = a.c_str();
        const char *a_end = a.c_str() + a.size() + 1 - kmerlen;
        beg = skip_Ns(beg, a_end, kmerlen, true);
        for (; beg != a_end; ++beg) {
          beg = skip_Ns(beg, a_end, kmerlen, false);
          if (beg == a_end)
            break;
          ht[beg].weights[A_or_B] += c; // relies on default init to 0
        }
      }
    }
  }
}

template<typename Ht>
void normalize_kmer_distributions(Ht& ht)
{
  typedef std::pair<kmer_key const, kmer_info> X;
  double denom_A = 0.0;
  double denom_B = 0.0;

  BOOST_FOREACH(const X& x, ht) {
    denom_A += x.second.weight_in_A();
    denom_B += x.second.weight_in_B();
  }

  BOOST_FOREACH(X& x, ht) {
    x.second.weight_in_A() /= denom_A;
    x.second.weight_in_B() /= denom_B;
  }
}

size_t estimate_hashtable_size(
    const std::vector<std::string>& A,
    const std::vector<std::string>& B,
    size_t kmerlen,
    double hash_table_fudge_factor)
{
  size_t max_entries = 0;
  BOOST_FOREACH(const std::string& a, A)
    if (a.size() >= kmerlen)
      max_entries += static_cast<size_t>(0.5 +
        2 * (a.size() + 1 - kmerlen) / hash_table_fudge_factor);
  BOOST_FOREACH(const std::string& b, B)
    if (b.size() >= kmerlen)
      max_entries += static_cast<size_t>(0.5 +
        2 * (b.size() + 1 - kmerlen) / hash_table_fudge_factor);
  return max_entries;
}

template<typename Ht>
void compute_stats(
    const Ht& ht,
    const std::string& prefix)
{
  typedef std::pair<kmer_key, kmer_info> X;
  double KL_A_to_M = 0.0;
  double KL_B_to_M = 0.0;
  double hellinger = 0.0;
  double total_var = 0.0;

  BOOST_FOREACH(const X& x, ht) {
    double w_A = x.second.weight_in_A();
    double w_B = x.second.weight_in_B();

    double mean_prob = 0.5*(w_A + w_B);
    KL_A_to_M += w_A == 0 ? 0 : w_A * (log2(w_A) - log2(mean_prob));
    KL_B_to_M += w_B == 0 ? 0 : w_B * (log2(w_B) - log2(mean_prob));

    double y = sqrt(w_A) - sqrt(w_B);
    hellinger += y * y;
    
    total_var += fabs(w_A - w_B);
  }
  double JS = 0.5*KL_A_to_M + 0.5*KL_B_to_M;
  hellinger = sqrt(hellinger)/sqrt(2.0);
  total_var *= 0.5;

  std::cout << prefix << "_kmer_KL_A_to_M\t"       << KL_A_to_M << std::endl;
  std::cout << prefix << "_kmer_KL_B_to_M\t"       << KL_B_to_M << std::endl;
  std::cout << prefix << "_kmer_jensen_shannon\t"  << JS << std::endl;
  std::cout << prefix << "_kmer_hellinger\t"       << hellinger << std::endl;
  std::cout << prefix << "_kmer_total_variation\t" << total_var << std::endl;
}

template<typename Ht>
void main_2(
    const opts& o,
    const fasta& A,
    const fasta& B,
    const std::vector<std::string>& A_rc,
    const std::vector<std::string>& B_rc,
    const expr& tau_A,
    const expr& tau_B,
    const std::string& prefix)
{
  size_t max_entries = estimate_hashtable_size(A.seqs, B.seqs, o.kmerlen, o.hash_table_fudge_factor);
  std::cerr << "Initializing the hash table with space for " << max_entries << " entries..." << std::endl;
  Ht ht(max_entries, kmer_key_hash(o.kmerlen), kmer_key_equal_to(o.kmerlen));
  empty_key_initializer<Ht> eki(ht, o.kmerlen);

  std::cerr << "Populating the hash table..." << std::flush;
  count_kmers<Ht, 0>(ht, A.seqs, A_rc, tau_A, o.kmerlen, o.strand_specific);
  count_kmers<Ht, 1>(ht, B.seqs, B_rc, tau_B, o.kmerlen, o.strand_specific);
  std::cerr << "done; hash table contains " << ht.size() << " entries." << std::endl;

  std::cerr << "Normalizing the induced distributions..." << std::endl;
  normalize_kmer_distributions(ht);

  std::cerr << "Computing kmer Jensen-Shannon, Hellinger, and total variation scores..." << std::endl;
  compute_stats(ht, prefix);
}

template<typename Ht>
void main_1(
    const opts& o,
    const fasta& A,
    const fasta& B,
    const expr& tau_A,
    const expr& tau_B,
    const expr& unif_A,
    const expr& unif_B)
{
  std::cerr << "Reverse complementing the sequences..." << std::flush;
  std::vector<std::string> A_rc, B_rc;
  transform(A.seqs.begin(), A.seqs.end(), back_inserter(A_rc), reverse_complement);
  transform(B.seqs.begin(), B.seqs.end(), back_inserter(B_rc), reverse_complement);
  std::cerr << "done." << std::endl;

  if (o.weighted)
    main_2<Ht>(o, A, B, A_rc, B_rc, tau_A, tau_B, "weighted");

  if (o.unweighted)
    main_2<Ht>(o, A, B, A_rc, B_rc, unif_A, unif_B, "unweighted");
}

void main(
    const opts& o,
    const fasta& A,
    const fasta& B,
    const expr& tau_A,
    const expr& tau_B,
    const expr& unif_A,
    const expr& unif_B)
{
  if (o.kmer) {

    if (o.hash_table_type == "sparse")
      main_1<sparse_kmer_map>(o, A, B, tau_A, tau_B, unif_A, unif_B);
    else if (o.hash_table_type == "dense")
      main_1<dense_kmer_map>(o, A, B, tau_A, tau_B, unif_A, unif_B);
    else
      throw std::runtime_error("Unknown hash map type.");

  }
}

} // namespace kmer
} // namespace re
