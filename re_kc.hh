#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <list>
#include <boost/foreach.hpp>
#include <sparsehash/sparse_hash_map>
//#include <sparsehash/dense_hash_map>
#include "util.hh"
#include "kmer_key.hh"
using namespace std;

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

typedef google::sparse_hash_map<kmer_key, kmer_info, kmer_key_hash, kmer_key_equal_to> kmer_hash_map;
//typedef google::dense_hash_map<kmer_key, kmer_info, kmer_key_hash, kmer_key_equal_to> kmer_hash_map;

void count_kmers_in_A(
    kmer_hash_map& ht,
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
        string::const_iterator beg = a.begin();
        string::const_iterator end = a.begin() + kmerlen;
        for (; end != a.end(); ++beg, ++end)
          ht[kmer_key(beg, end)].is_present_in_A = true;
        ht[kmer_key(beg, end)].is_present_in_A = true;
      }
    }
  }
}

void count_kmers_in_B(
    kmer_hash_map& ht,
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
        string::const_iterator beg = b.begin();
        string::const_iterator end = b.begin() + kmerlen;
        for (; end != b.end(); ++beg, ++end)
          ht[kmer_key(beg, end)].weight_in_B += c; // relies on default init to 0
        ht[kmer_key(beg, end)].weight_in_B += c; // relies on default init to 0
      }
    }
  }
}

void normalize_kmer_distributions(kmer_hash_map& ht)
{
  typedef std::pair<const kmer_key, kmer_info> X;

  double denom = 0;
  BOOST_FOREACH(const X& x, ht) {
    denom += x.second.weight_in_B;
  }

  BOOST_FOREACH(X& x, ht) {
    x.second.weight_in_B /= denom;
  }
}

size_t estimate_hashtable_size(
    const vector<string>& A,
    const vector<string>& B,
    size_t kmerlen)
{
  // The fudge factor is based partly on empirical observation and partly on
  // the assumption that most kmers will be shared between the oracleset and
  // the assembly.
  size_t fudge_factor = 2;
  size_t max_entries = 0;
  BOOST_FOREACH(const string& a, A)
    if (a.size() >= kmerlen)
      max_entries += 2 * (a.size() + 1 - kmerlen) / fudge_factor;
  BOOST_FOREACH(const string& a, B)
    if (a.size() >= kmerlen)
      max_entries += 2 * (a.size() + 1 - kmerlen) / fudge_factor;
  return max_entries;
}

kmer_hash_map build_hash_table(
    const vector<string>& A,
    const vector<string>& A_rc, 
    const vector<string>& B,
    const vector<string>& B_rc,
    const vector<double>& tau_B,
    size_t kmerlen,
    bool strand_specific)
{
  // Figure out roughly how many entries to expect in the hash table.
  size_t max_entries = estimate_hashtable_size(A, B, kmerlen);
  std::cerr << "...estimated number of hash table entries: " << max_entries << std::endl;

  // Create hash table.
  kmer_hash_map ht(max_entries, kmer_key_hash(), kmer_key_equal_to());

  // Add kmers to the hash table.
  count_kmers_in_A(ht, A, A_rc,        kmerlen, strand_specific);
  count_kmers_in_B(ht, B, B_rc, tau_B, kmerlen, strand_specific);
  std::cerr << "...actual number of hash table entries: " << ht.size() << std::endl;

  return ht;
}

double compute_kmer_recall(const kmer_hash_map& ht)
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
    if (false) {
    //if (vm.count("estimate-hashtable-size")) {
      std::cout << estimate_hashtable_size(A.seqs, B.seqs, o.readlen) * sizeof(kmer_info) << std::endl;
      return;
    }

    std::cerr << "Reverse complementing the sequences..." << std::flush;
    std::vector<std::string> A_rc, B_rc;
    transform(A.seqs.begin(), A.seqs.end(), back_inserter(A_rc), reverse_complement);
    transform(B.seqs.begin(), B.seqs.end(), back_inserter(B_rc), reverse_complement);
    std::cerr << "done." << std::endl;

    std::cerr << "Building the hash table..." << std::endl;
    kmer_hash_map ht = build_hash_table(A.seqs, A_rc, B.seqs, B_rc, tau_B, o.readlen, o.strand_specific);
    std::cerr << "...done." << std::endl;

    double wkr = compute_kmer_recall(ht);
    double icr = compute_inverse_compression_rate(o, A);

    std::cout << "weighted_kmer_recall\t" << wkr << std::endl;
    std::cout << "inverse_compression_rate\t" << icr << std::endl;
    std::cout << "kmer_compression_score\t" << wkr - icr << std::endl;

    // compute_kmer_recall(A.seqs, A_rc, B.seqs, B_rc, unif_B, "unweighted_kmer", "at_one",    o.readlen,   o.strand_specific);

    // expr dc_unif_B(B.card, 1.0);
    // compute_kmer_recall(A.seqs, A_rc, B.seqs, B_rc, dc_unif_B, "dc_unweighted_kmer", "at_one",    o.readlen,   o.strand_specific);

    // // compute_kmer_recall(A.seqs, A_rc, B.seqs, B_rc, tau_B, "weighted_kmer", "at_double", o.readlen*2, o.strand_specific);
    // // compute_kmer_recall(A.seqs, A_rc, B.seqs, B_rc, tau_B, "weighted_kmer", "at_half",   o.readlen/2, o.strand_specific);
  }
}

} // namespace kc
} // namespace re
