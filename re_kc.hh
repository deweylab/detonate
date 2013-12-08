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
#include "util.hh"
#include "kmer_key.hh"

namespace re {
namespace kc {

struct kmer_info
{
  bool is_present_in_A;
  double weight_in_B;
  kmer_info()
  : is_present_in_A(false),
    weight_in_B(0.0)
  {}
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

template<typename Ht>
void count_kmers_in_A(
    Ht& ht,
    const vector<string>& A,
    const vector<string>& A_rc,
    size_t kmerlen,
    bool strand_specific)
{
  // For each contig a in A:
  //   For each kmer r in a or reverse_complement(a):
  //     Mark r as being present in A.
  size_t num_strands = strand_specific ? 1 : 2;
  for (size_t i = 0; i < A.size(); ++i) {
    for (size_t which = 0; which < num_strands; ++which) {
      const string& a = which == 0 ? A[i] : A_rc[i];
      if (a.size() >= kmerlen) {
        const char *beg = a.c_str();
        const char *a_end = a.c_str() + a.size() + 1 - kmerlen;
        for (; beg != a_end; ++beg)
          ht[beg].is_present_in_A = true;
      }
    }
  }
}

template<typename Ht>
void count_kmers_in_B(
    Ht& ht,
    const vector<string>& B,
    const vector<string>& B_rc,
    const vector<double>& tau_B,
    size_t kmerlen,
    bool strand_specific)
{
  // For each contig b in B:
  //   For each kmer r in b or reverse_complement(b):
  //     Add weight(b) to weight_in_B(r).
  size_t num_strands = strand_specific ? 1 : 2;
  for (size_t i = 0; i < B.size(); ++i) {
    for (size_t which = 0; which < num_strands; ++which) {
      const string& b = which == 0 ? B[i] : B_rc[i];
      if (b.size() >= kmerlen) {
        double c = tau_B[i];
        const char *beg = b.c_str();
        const char *b_end = b.c_str() + b.size() + 1 - kmerlen;
        for (; beg != b_end; ++beg)
          ht[beg].weight_in_B += c; // relies on default init to 0
      }
    }
  }
}

size_t estimate_hashtable_size(
    const vector<string>& A,
    const vector<string>& B,
    size_t kmerlen,
    double hash_table_fudge_factor)
{
  size_t max_entries = 0;
  BOOST_FOREACH(const string& a, A)
    if (a.size() >= kmerlen)
      max_entries += static_cast<size_t>(0.5 +
        2 * (a.size() + 1 - kmerlen) / hash_table_fudge_factor);
  BOOST_FOREACH(const string& b, B)
    if (b.size() >= kmerlen)
      max_entries += static_cast<size_t>(0.5 +
        2 * (b.size() + 1 - kmerlen) / hash_table_fudge_factor);
  return max_entries;
}

template<typename Ht>
double compute_kmer_recall(const Ht& ht)
{
  typedef std::pair<const kmer_key, kmer_info> X;
  double numer = 0.0, denom = 0.0;
  BOOST_FOREACH(const X& x, ht) {
    const kmer_info& i = x.second;
    if (i.is_present_in_A > 0)
      numer += i.weight_in_B;
    denom += i.weight_in_B;
  }
  return numer / denom;
}

double compute_inverse_compression_rate(
    const opts& o,
    const fasta& A)
{
  size_t num_bases_in_A = std::accumulate(A.lengths.begin(), A.lengths.end(), 0);
  return 1.0 * num_bases_in_A / (o.num_reads * o.readlen);
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

  size_t max_entries = estimate_hashtable_size(A.seqs, B.seqs, o.readlen, o.hash_table_fudge_factor);
  std::cerr << "Initializing the hash table with space for " << max_entries << " entries..." << std::flush;
  Ht ht(max_entries, kmer_key_hash(o.readlen), kmer_key_equal_to(o.readlen));
  empty_key_initializer<Ht> eki(ht, o.readlen);
  std::cerr << "done." << std::endl;

  std::cerr << "Populating the hash table..." << std::flush;
  count_kmers_in_A(ht, A.seqs, A_rc,        o.readlen, o.strand_specific);
  count_kmers_in_B(ht, B.seqs, B_rc, tau_B, o.readlen, o.strand_specific);
  std::cerr << "done; hash table contains " << ht.size() << " entries." << std::endl;

  double wkr = compute_kmer_recall(ht);
  double icr = compute_inverse_compression_rate(o, A);

  std::cout << "weighted_kmer_recall\t" << wkr << std::endl;
  std::cout << "inverse_compression_rate\t" << icr << std::endl;
  std::cout << "kmer_compression_score\t" << wkr - icr << std::endl;
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
  if (o.kc) {

    if (o.hash_table_type == "sparse")
      main_1<sparse_kmer_map>(o, A, B, tau_A, tau_B, unif_A, unif_B);
    else if (o.hash_table_type == "dense")
      main_1<dense_kmer_map>(o, A, B, tau_A, tau_B, unif_A, unif_B);
    else
      throw std::runtime_error("Unknown hash map type.");

  }
}

} // namespace kc
} // namespace re
