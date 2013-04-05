#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <string>
#include <vector>
#include <list>
#include <omp.h>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/array.hpp>
#include <sparsehash/sparse_hash_map>
//#include <sparsehash/dense_hash_map>
#include "util.hh"
using namespace std;

struct Probs
{
  double for_A, for_B;
  Probs() : for_A(0.0), for_B(0.0) {}
};

//typedef boost::array<google::sparse_hash_map<size_t, Probs>, 16> KmerStats;

struct KmerStats
{
  boost::array<std::vector<Probs>, 16> data;

  double& get_or_init(size_t pair, size_t dist, size_t A_or_B)
  {
    if (dist >= data[pair].size())
      data[pair].resize(std::max(dist+1, data[pair].size()*2));
    if (A_or_B == 0)
      return data[pair][dist].for_A;
    else
      return data[pair][dist].for_B;
  }
};

inline size_t choose_2(size_t n) { return n*(n-1)/2; }

inline size_t encode(char c)
{
  if (c == 'C')
    return 0;
  else if (c == 'G')
    return 1;
  else if (c == 'A')
    return 2;
  else if (c == 'T')
    return 3;
  else
    return 4;
}

void count_kmers(
    KmerStats& stats,
    size_t A_or_B, // 0 if we're counting A's mers, 1 if B's
    const vector<string>& A,
    const vector<string>& A_rc,
    const vector<double>& tau_A)
{
  // For each X in {A, A_rc}:
  //   For each contig a in X:
  //     For each l in 1:length(a):
  //       For each r in l:length(a):
  //         Add (1/2)*[1/choose(length(a)+1,2)]*tau_A(a) to count_A(a[l], a[r], l-r+1)
  for (size_t which = 0; which < 2; ++which) {
    const vector<string>& X = (which == 0) ? A : A_rc;
    for (size_t i = 0; i < X.size(); ++i) {
      const string& a = X[i];
      double c = 0.5*(1.0/choose_2(a.size()+1))*tau_A[i];
      string::const_iterator a_l, a_r, a_end = a.end();
      for (a_l = a.begin(); a_l != a_end; ++a_l) {
        size_t dist = 1;
        for (a_r = a_l; a_r != a_end; ++a_r) {
          size_t i = encode(*a_l);
          size_t j = encode(*a_r);
          if (i < 4 && j < 4)
            stats.get_or_init(i+4*j, dist, A_or_B) += c; // relies on default init of Probs to {0.0, 0.0}
          ++dist;
        }
      }
    }
  }
}

void normalize_kmer_distributions(KmerStats& stats)
{
  double denom_for_A = 0.0, denom_for_B = 0.0;
  for (size_t pair = 0; pair < 16; ++pair) {
    BOOST_FOREACH(const Probs& probs, stats.data[pair]) {
      denom_for_A += probs.for_A;
      denom_for_B += probs.for_B;
    }
  }

  for (size_t pair = 0; pair < 16; ++pair) {
    BOOST_FOREACH(Probs& probs, stats.data[pair]) {
      probs.for_A /= denom_for_A;
      probs.for_B /= denom_for_B;
    }
  }
}

void compute_kmer_stats(
    KmerStats& stats,
    const vector<string>& A,
    const vector<string>& A_rc, 
    const vector<double>& tau_A,
    const vector<string>& B,
    const vector<string>& B_rc,
    const vector<double>& tau_B)
{
  count_kmers(stats, 0, A, A_rc, tau_A);
  count_kmers(stats, 1, B, B_rc, tau_B);
  normalize_kmer_distributions(stats);
}

void print_kmer_stats(
    const KmerStats& stats,
    const string& prefix)
{
  double KL_A_to_M = 0.0, KL_B_to_M = 0.0;
  for (size_t pair = 0; pair < 16; ++pair) {
    BOOST_FOREACH(const Probs& probs, stats.data[pair]) {
      double mean_prob = 0.5*(probs.for_A + probs.for_B);
      KL_A_to_M += probs.for_A == 0 ? 0 : probs.for_A * (log2(probs.for_A) - log2(mean_prob));
      KL_B_to_M += probs.for_B == 0 ? 0 : probs.for_B * (log2(probs.for_B) - log2(mean_prob));
    }
  }
  double JS = 0.5*KL_A_to_M + 0.5*KL_B_to_M;
  cout << prefix << "_KL_A_to_M" << "\t" << KL_A_to_M << endl;
  cout << prefix << "_KL_B_to_M" << "\t" << KL_B_to_M << endl;
  cout << prefix << "_JS"        << "\t" << JS        << endl;
}

void compute_and_print_kmer_stats(
    const vector<string>& A,
    const vector<string>& A_rc, 
    const vector<double>& tau_A,
    const vector<string>& B,
    const vector<string>& B_rc,
    const vector<double>& tau_B,
    const string& prefix)
{
  KmerStats kmer_stats;
  compute_kmer_stats(kmer_stats, A, A_rc, tau_A, B, B_rc, tau_B);
  print_kmer_stats(kmer_stats, prefix);
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
    ("A-expr", po::value<std::string>()->required(), "The assembly expression, as produced by RSEM in a file called *.isoforms.results.")
    ("B-expr", po::value<std::string>()->required(), "The oracleset expression, as produced by RSEM in a file called *.isoforms.results.")
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
  //clock_t start;

  std::cerr << "Reading the sequences" << std::endl;
  std::vector<std::string> A, B;
  std::vector<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  read_fasta_names_and_seqs(vm["A-seqs"].as<std::string>(), A, A_names, A_names_to_idxs);
  read_fasta_names_and_seqs(vm["B-seqs"].as<std::string>(), B, B_names, B_names_to_idxs);

  std::cerr << "Reverse complementing the sequences" << std::endl;
  std::vector<std::string> A_rc, B_rc;
  transform(A.begin(), A.end(), back_inserter(A_rc), reverse_complement);
  transform(B.begin(), B.end(), back_inserter(B_rc), reverse_complement);

  std::cerr << "Reading transcript-level expression" << std::endl;
  std::vector<double> tau_A(A.size()), tau_B(B.size());
  std::string A_expr_fname = vm["A-expr"].as<std::string>();
  std::string B_expr_fname = vm["B-expr"].as<std::string>();
  read_transcript_expression(A_expr_fname, tau_A, A_names_to_idxs);
  read_transcript_expression(B_expr_fname, tau_B, B_names_to_idxs);

  std::cerr << "Computing nucleotide-level expression" << std::endl;
  std::vector<double> nu_A(A.size()), nu_B(B.size());
  compute_nucl_expression(A, tau_A, nu_A);
  compute_nucl_expression(B, tau_B, nu_B);

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_A(A.size(), 1.0/A.size());
  std::vector<double> unif_tau_B(B.size(), 1.0/B.size());

  std::cerr << "Computing uniform nucleotide-level expression" << std::endl;
  std::vector<double> unif_nu_A(A.size()), unif_nu_B(B.size());
  compute_nucl_expression(A, unif_tau_A, unif_nu_A);
  compute_nucl_expression(B, unif_tau_B, unif_nu_B);

  std::cout << "summarize_kmerpair_version\t1" << std::endl;

  compute_and_print_kmer_stats(A, A_rc, tau_A, B, B_rc, tau_B, "weighted_kmerpair");
  compute_and_print_kmer_stats(A, A_rc, unif_tau_A, B, B_rc, unif_tau_B, "unweighted_kmerpair");
}
