#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <list>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <sparsehash/sparse_hash_map>
//#include <sparsehash/dense_hash_map>
#include "util.hh"
using namespace std;

struct KmerInfo
{
  bool is_present_in_A;
  double weight_in_B;
  KmerInfo()
  : is_present_in_A(false),
    weight_in_B(0.0)
  {}
};

struct KmerKey
{
  string::const_iterator begin;
  string::const_iterator end;
  KmerKey() {}
  KmerKey(const string::const_iterator& begin, const string::const_iterator& end) : begin(begin), end(end) {}
};

namespace std
{
  namespace tr1
  {
    template<>
    struct hash<KmerKey>
    {
      size_t operator()(const KmerKey& k) const
      {
        return hash<string>()(string(k.begin, k.end));
      }
    };
  }

  template<>
  struct equal_to<KmerKey>
  {
    bool operator()(const KmerKey& lhs, const KmerKey& rhs) const
    {
      string::const_iterator i = lhs.begin;
      string::const_iterator j = rhs.begin;
      for (;;) {
        // If i == end and j == end, return true.
        // If i != end and j != end and *i == *j, continue;
        // In all other cases, return false.
        if (i == lhs.end) {
          if (j == rhs.end)
            return true;    // i == end and j == end              -> true
          else
            return false;   // i == end and j != end              -> false
        } else {
          if (j == rhs.end)
            return false;   // i != end and j == end              -> false
          else {
            if (*i != *j)
              return false; // i != end and j != end and *i != *j -> false
            else {
            }               // i != end and j != end and *i == *j -> continue
          }
        }
        ++i, ++j;
      }
    }
  };
}

struct KmerStats
{
  typedef KmerKey key_type;
  typedef KmerInfo value_type;
  typedef google::sparse_hash_map<key_type, value_type> container_type;
  //typedef google::dense_hash_map<key_type, value_type> container_type;
};

void count_kmers_in_A(
    KmerStats::container_type& stats,
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
          stats[KmerKey(beg, end)].is_present_in_A = true;
        stats[KmerKey(beg, end)].is_present_in_A = true;
      }
    }
  }
}

void count_kmers_in_B(
    KmerStats::container_type& stats,
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
          stats[KmerKey(beg, end)].weight_in_B += c; // relies on default init to 0
        stats[KmerKey(beg, end)].weight_in_B += c; // relies on default init to 0
      }
    }
  }
}

void normalize_kmer_distributions(KmerStats::container_type& stats)
{
  typedef std::pair<const KmerStats::key_type, KmerStats::value_type> X;

  double denom = 0;
  BOOST_FOREACH(const X& x, stats) {
    denom += x.second.weight_in_B;
  }

  BOOST_FOREACH(X& x, stats) {
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

KmerStats::container_type
compute_kmer_stats(
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
  std::cerr << "Max entries: " << max_entries << std::endl;

  // Create hash table.
  KmerStats::container_type stats(max_entries);

  // Actually compute the desired stats.
  count_kmers_in_A(stats, A, A_rc,        kmerlen, strand_specific);
  count_kmers_in_B(stats, B, B_rc, tau_B, kmerlen, strand_specific);
  std::cerr << "... hashtable size is " << stats.size() << std::endl;
  //normalize_kmer_distributions(stats);

  return stats;
}

void print_kmer_stats(
    KmerStats::container_type& stats,
    const string& prefix,
    const string& suffix)
{
  typedef std::pair<const KmerStats::key_type, KmerStats::value_type> X;

  // First, compute and print stats based on the UNnormalized distribution.

  {
    double simple_recall_numer = 0.0, simple_recall_denom = 0.0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      if (i.is_present_in_A > 0)
        simple_recall_numer += i.weight_in_B;
      simple_recall_denom += i.weight_in_B;
    }
    double simple_recall = simple_recall_numer / simple_recall_denom;
    cout << prefix << "_simple_recall_numer_" << suffix << "\t" << simple_recall_numer << endl;
    cout << prefix << "_simple_recall_denom_" << suffix << "\t" << simple_recall_denom << endl;
    cout << prefix << "_simple_recall_"       << suffix << "\t" << simple_recall       << endl;
  }

  // Second, compute and print stats based on the normalized distribution.
  normalize_kmer_distributions(stats);

  {
    double simple_recall = 0.0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      if (i.is_present_in_A > 0)
        simple_recall += i.weight_in_B;
    }
    cout << prefix << "_simple_recall_doublecheck_" << suffix << "\t" << simple_recall << endl;
  }
}

void run(
    const vector<string>& A,
    const vector<string>& A_rc, 
    const vector<string>& B,
    const vector<string>& B_rc,
    const vector<double>& tau_B,
    const string& prefix,
    const string& suffix,
    size_t kmerlen,
    bool strand_specific)
{
  KmerStats::container_type kmer_stats = compute_kmer_stats(A, A_rc, B, B_rc, tau_B, kmerlen, strand_specific);
  print_kmer_stats(kmer_stats, prefix, suffix);
}

void parse_options(boost::program_options::variables_map& vm, int argc, const char **argv)
{
  namespace po = boost::program_options;

  po::options_description desc("Options");
  desc.add_options()
    ("help,?", "Display this information.")
    ("A-seqs", po::value<std::string>()->required(), "The assembly sequences, in FASTA format.")
    ("B-seqs", po::value<std::string>()->required(), "The oracleset sequences, in FASTA format.")
    ("B-expr", po::value<std::string>()->required(), "The oracleset expression, as produced by RSEM in a file called *.isoforms.results.")
    ("readlen", po::value<size_t>()->required(),     "The read length.")
    ("estimate-hashtable-size",                      "Estimate hashtable size, in bytes.")
    ("strand-specific",                              "Ignore alignments that are to the reverse strand.")
  ;

  try {

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cerr << desc << std::endl;
      exit(1);
    }

    po::notify(vm);

  } catch (std::exception& x) {
    std::cerr << "Error: " << x.what() << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}

void main_1(const boost::program_options::variables_map& vm)
{
  bool strand_specific = vm.count("strand-specific");

  std::cerr << "Reading the sequences" << std::endl;
  std::vector<std::string> A, B;
  std::vector<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  read_fasta_names_and_seqs(vm["A-seqs"].as<std::string>(), A, A_names, A_names_to_idxs);
  read_fasta_names_and_seqs(vm["B-seqs"].as<std::string>(), B, B_names, B_names_to_idxs);

  size_t readlen = vm["readlen"].as<size_t>();

  if (vm.count("estimate-hashtable-size")) {
    std::cout << estimate_hashtable_size(A, B, readlen) * sizeof(KmerInfo) << std::endl;
    return;
  }

  std::cerr << "Reverse complementing the sequences" << std::endl;
  std::vector<std::string> A_rc, B_rc;
  transform(A.begin(), A.end(), back_inserter(A_rc), reverse_complement);
  transform(B.begin(), B.end(), back_inserter(B_rc), reverse_complement);

  std::cerr << "Reading transcript-level expression" << std::endl;
  std::vector<double> tau_B(B.size());
  std::string B_expr_fname = vm["B-expr"].as<std::string>();
  read_transcript_expression(B_expr_fname, tau_B, B_names_to_idxs);

  std::cout << "summarize_wkr_version_1\t0" << std::endl;
  run(A, A_rc, B, B_rc, tau_B, "weighted_kmer", "at_one",    readlen,   strand_specific);
  run(A, A_rc, B, B_rc, tau_B, "weighted_kmer", "at_double", readlen*2, strand_specific);
  run(A, A_rc, B, B_rc, tau_B, "weighted_kmer", "at_half",   readlen/2, strand_specific);
  std::cerr << "Done" << std::endl;
}
