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
  double probs[2]; // sometimes normalized, sometimes not
  KmerInfo()
  {
    probs[0] = 0;
    probs[1] = 0;
  }
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

void count_kmers(
    KmerStats::container_type& stats,
    size_t A_or_B, // 0 if we're counting A's mers, 1 if B's
    const vector<string>& A,
    const vector<string>& A_rc,
    const vector<double>& tau_A,
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
      const string& a = which == 0 ? A[i] : A_rc[i];
      if (a.size() >= kmerlen) {
        double c = tau_A[i];
        string::const_iterator beg = a.begin();
        string::const_iterator end = a.begin() + kmerlen;
        for (; end != a.end(); ++beg, ++end)
          stats[KmerKey(beg, end)].probs[A_or_B] += c; // relies on default init to 0
        stats[KmerKey(beg, end)].probs[A_or_B] += c; // relies on default init to 0
      }
    }
  }
}

void normalize_kmer_distributions(KmerStats::container_type& stats)
{
  typedef std::pair<const KmerStats::key_type, KmerStats::value_type> X;

  double denoms[2] = {0, 0};
  BOOST_FOREACH(const X& x, stats) {
    denoms[0] += x.second.probs[0];
    denoms[1] += x.second.probs[1];
  }

  BOOST_FOREACH(X& x, stats) {
    x.second.probs[0] /= denoms[0];
    x.second.probs[1] /= denoms[1];
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
    const vector<double>& tau_A,
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
  count_kmers(stats, 0, A, A_rc, tau_A, kmerlen, strand_specific);
  count_kmers(stats, 1, B, B_rc, tau_B, kmerlen, strand_specific);
  std::cerr << "... hashtable size is " << stats.size() << std::endl;
  //normalize_kmer_distributions(stats);

  return stats;
}

double compute_F1(double precis, double recall)
{
  if (precis == 0.0 && recall == 0.0)
    return 0.0;
  else
    return 2*precis*recall/(precis+recall);
}

void print_kmer_stats(
    KmerStats::container_type& stats,
    const string& prefix,
    const string& suffix)
{
  typedef std::pair<const KmerStats::key_type, KmerStats::value_type> X;

  // First, compute and print stats based on the UNnormalized distribution.

  {
    size_t num_nonzero_A = 0, num_nonzero_B = 0, num_shared = 0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerKey& k = x.first;
      const KmerInfo& i = x.second;
      if (i.probs[0] != 0)
        ++num_nonzero_A;
      if (i.probs[1] != 0)
        ++num_nonzero_B;
      if (i.probs[0] != 0 && i.probs[1] != 0)
        ++num_shared;
    }
    cout << prefix << "_num_nonzero_in_A_" << suffix << "\t" << num_nonzero_A << endl;
    cout << prefix << "_num_nonzero_in_B_" << suffix << "\t" << num_nonzero_B << endl;
    cout << prefix << "_num_shared_" << suffix << "\t" << num_shared << endl;
  }

  {
    double numer = 0.0, denom_precis = 0.0, denom_recall = 0.0, frac_present_numer = 0.0, frac_present_denom = 0.0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      numer += std::min(i.probs[0], i.probs[1]); // frac "retrieved" and "relevant"
      denom_precis += i.probs[0]; // frac "retrieved"
      denom_recall += i.probs[1]; // frac "relevant"
      if (i.probs[0] > 0)
        frac_present_numer += i.probs[1];
      frac_present_denom += i.probs[1];
    }
    double precis = numer / denom_precis;
    double recall = numer / denom_recall;
    double F1 = compute_F1(precis, recall);
    double frac_present = frac_present_numer / frac_present_denom;
    cout << prefix << "_precision_" << suffix << "\t" << precis << endl;
    cout << prefix << "_recall_"    << suffix << "\t" << recall << endl;
    cout << prefix << "_F1_"        << suffix << "\t" << F1     << endl;
    cout << prefix << "_frac_present_numer_" << suffix << "\t" << frac_present_numer << endl;
    cout << prefix << "_frac_present_denom_" << suffix << "\t" << frac_present_denom << endl;
    cout << prefix << "_frac_present_"       << suffix << "\t" << frac_present       << endl;
  }

  // Second, compute and print stats based on the normalized distribution.
  normalize_kmer_distributions(stats);

  {
    double KL_A_to_M = 0, KL_B_to_M = 0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      double mean_prob = 0.5*(i.probs[0] + i.probs[1]);
      KL_A_to_M += i.probs[0] == 0 ? 0 : i.probs[0] * (log2(i.probs[0]) - log2(mean_prob));
      KL_B_to_M += i.probs[1] == 0 ? 0 : i.probs[1] * (log2(i.probs[1]) - log2(mean_prob));
    }
    double JS = 0.5*KL_A_to_M + 0.5*KL_B_to_M;
    cout << prefix << "_probs_KL_A_to_M_" << suffix << "\t" << KL_A_to_M << endl;
    cout << prefix << "_probs_KL_B_to_M_" << suffix << "\t" << KL_B_to_M << endl;
    cout << prefix << "_probs_JS_" << suffix << "\t" << JS << endl;
  }

  {
    double probs_hellinger = 0.0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      double y = sqrt(i.probs[0]) - sqrt(i.probs[1]);
      probs_hellinger += y * y;
    }
    probs_hellinger = sqrt(probs_hellinger)/sqrt(2.0);
    cout << prefix << "_probs_hellinger_" << suffix << "\t" << probs_hellinger << endl;
  }

}

void run(
    const vector<string>& A,
    const vector<string>& A_rc, 
    const vector<double>& tau_A,
    const vector<string>& B,
    const vector<string>& B_rc,
    const vector<double>& tau_B,
    const string& prefix,
    const string& suffix,
    size_t kmerlen,
    bool strand_specific)
{
  KmerStats::container_type kmer_stats = compute_kmer_stats(A, A_rc, tau_A, B, B_rc, tau_B, kmerlen, strand_specific);
  print_kmer_stats(kmer_stats, prefix, suffix);
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

void parse_options(boost::program_options::variables_map& vm, int argc, const char **argv)
{
  namespace po = boost::program_options;

  po::options_description desc("Options");
  desc.add_options()
    ("help,?", "Display this information.")
    ("A-seqs", po::value<std::string>()->required(), "The assembly sequences, in FASTA format.")
    ("B-seqs", po::value<std::string>()->required(), "The oracleset sequences, in FASTA format.")
    ("A-expr", po::value<std::string>(),             "The assembly expression, as produced by RSEM in a file called *.isoforms.results.")
    ("B-expr", po::value<std::string>(),             "The oracleset expression, as produced by RSEM in a file called *.isoforms.results.")
    ("no-expr",                                      "Do not use expression at all. No weighted scores will be produced.")
    ("readlen", po::value<size_t>()->required(),     "The read length.")
    ("estimate-hashtable-size",                      "Estimate hashtable size, in bytes.")
    ("strand-specific",                              "Ignore alignments that are to the reverse strand.")
    ("at-half-and-double",                           "Also run with k = half and double the read length.")
  ;

  try {

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cerr << desc << std::endl;
      exit(1);
    }

    po::notify(vm);

    if (vm.count("no-expr")) {
      if (vm.count("A-expr") != 0 || vm.count("B-expr") != 0)
        throw po::error("If --no-expr is given, then --A-expr and --B-expr cannot be given.");
    } else {
      if (vm.count("A-expr") == 0 || vm.count("B-expr") == 0)
        throw po::error("If --no-expr is not given, then --A-expr and --B-expr must be given.");
    }

  } catch (std::exception& x) {
    std::cerr << "Error: " << x.what() << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}

void main_1(const boost::program_options::variables_map& vm)
{
  bool strand_specific = vm.count("strand-specific");
  bool at_half_and_double = vm.count("at-half-and-double");
  bool no_expr = vm.count("no-expr");

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
  std::vector<double> tau_A(A.size()), tau_B(B.size());
  if (!no_expr) {
    std::string A_expr_fname = vm["A-expr"].as<std::string>();
    std::string B_expr_fname = vm["B-expr"].as<std::string>();
    read_transcript_expression(A_expr_fname, tau_A, A_names_to_idxs);
    read_transcript_expression(B_expr_fname, tau_B, B_names_to_idxs);
  }

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_A(A.size(), 1.0/A.size());
  std::vector<double> unif_tau_B(B.size(), 1.0/B.size());

  std::cout << "summarize_kmer_version_6\t0" << std::endl;
  if (!no_expr) {
    run(  A, A_rc, tau_A, B, B_rc, tau_B, "weighted_kmer", "at_one",    readlen,   strand_specific);
    if (at_half_and_double) {
      run(A, A_rc, tau_A, B, B_rc, tau_B, "weighted_kmer", "at_half",   readlen/2, strand_specific);
      run(A, A_rc, tau_A, B, B_rc, tau_B, "weighted_kmer", "at_double", readlen*2, strand_specific);
    }
  }

  run(  A, A_rc, unif_tau_A, B, B_rc, unif_tau_B, "unweighted_kmer", "at_one",    readlen,   strand_specific);
  if (at_half_and_double) {
    run(A, A_rc, unif_tau_A, B, B_rc, unif_tau_B, "unweighted_kmer", "at_half",   readlen/2, strand_specific);
    run(A, A_rc, unif_tau_A, B, B_rc, unif_tau_B, "unweighted_kmer", "at_double", readlen*2, strand_specific);
  }
  std::cerr << "Done" << std::endl;
}
