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

struct for_A {};
struct for_B {};
struct wei {};
struct unw {};

struct Probs
{
  double data[2][2];

  Probs()
  {
    data[0][0] = data[0][1] = data[1][0] = data[1][1] = 0.0;
  }

  template<typename for_A_or_B, typename wei_or_unw> inline       double& at();
  template<typename for_A_or_B, typename wei_or_unw> inline const double& at() const;
};

template<> inline       double& Probs::at<for_A, wei>()       { return data[0][0]; }
template<> inline       double& Probs::at<for_A, unw>()       { return data[0][1]; }
template<> inline       double& Probs::at<for_B, wei>()       { return data[1][0]; }
template<> inline       double& Probs::at<for_B, unw>()       { return data[1][1]; }

template<> inline const double& Probs::at<for_A, wei>() const { return data[0][0]; }
template<> inline const double& Probs::at<for_A, unw>() const { return data[0][1]; }
template<> inline const double& Probs::at<for_B, wei>() const { return data[1][0]; }
template<> inline const double& Probs::at<for_B, unw>() const { return data[1][1]; }

struct KmerStats
{
  boost::array<std::vector<Probs>, 16> data;

  KmerStats(size_t max)
  {
    for (size_t pair = 0; pair < 16; ++pair)
      data[pair].resize(max);
  }
};

inline size_t choose_2(size_t n) { return n*(n-1)/2; }

size_t encode(char c)
{
  switch (c)
  {
    case 'C': return 0;
    case 'G': return 1;
    case 'A': return 2;
    case 'T': return 3;
    case 'c': return 0;
    case 'g': return 1;
    case 'a': return 2;
    case 't': return 3;
    default:  return 4;
  }
}

// A_or_B is 0 if we're counting A's mers, 1 if B's
template<typename for_A_or_B>
void count_kmers(
    KmerStats& stats,
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
      double c_weighted   = 0.5*(1.0/choose_2(a.size()+1))*tau_A[i];
      double c_unweighted = 0.5*(1.0/choose_2(a.size()+1));
      string::const_iterator a_l, a_r, a_end = a.end();
      for (a_l = a.begin(); a_l != a_end; ++a_l) {
        size_t dist = 0;
        for (a_r = a_l; a_r != a_end; ++a_r) {
          size_t i = encode(*a_l);
          size_t j = encode(*a_r);
          if (i < 4 && j < 4) {
            Probs& probs = stats.data[i+4*j][dist];
            probs.at<for_A_or_B, wei>() += c_weighted;   // relies on default init of Probs to all 0's
            probs.at<for_A_or_B, unw>() += c_unweighted; // relies on default init of Probs to all 0's
          }
          ++dist;
        }
      }
    }
  }
}

template<typename for_A_or_B, typename wei_or_unw>
void normalize_kmer_distributions(KmerStats& stats)
{
  double denom = 0.0;

  for (size_t pair = 0; pair < 16; ++pair)
    BOOST_FOREACH(const Probs& probs, stats.data[pair])
      denom += probs.at<for_A_or_B, wei_or_unw>();

  for (size_t pair = 0; pair < 16; ++pair)
    BOOST_FOREACH(Probs& probs, stats.data[pair])
      probs.at<for_A_or_B, wei_or_unw>() /= denom;
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
  count_kmers<for_A>(stats, A, A_rc, tau_A);
  count_kmers<for_B>(stats, B, B_rc, tau_B);
  normalize_kmer_distributions<for_A, wei>(stats);
  normalize_kmer_distributions<for_A, unw>(stats);
  normalize_kmer_distributions<for_B, wei>(stats);
  normalize_kmer_distributions<for_B, unw>(stats);
}

template<typename wei_or_unw>
void print_kmer_stats(
    const KmerStats& stats,
    const string& prefix)
{
  double KL_A_to_M = 0.0, KL_B_to_M = 0.0;
  for (size_t pair = 0; pair < 16; ++pair) {
    BOOST_FOREACH(const Probs& probs, stats.data[pair]) {
      double prob_for_A = probs.at<for_A, wei_or_unw>();
      double prob_for_B = probs.at<for_B, wei_or_unw>();
      double mean_prob = 0.5*(prob_for_A + prob_for_B);
      KL_A_to_M += prob_for_A == 0 ? 0 : prob_for_A * (log2(prob_for_A) - log2(mean_prob));
      KL_B_to_M += prob_for_B == 0 ? 0 : prob_for_B * (log2(prob_for_B) - log2(mean_prob));
    }
  }
  double JS = 0.5*KL_A_to_M + 0.5*KL_B_to_M;
  cout << prefix << "_KL_A_to_M" << "\t" << KL_A_to_M << endl;
  cout << prefix << "_KL_B_to_M" << "\t" << KL_B_to_M << endl;
  cout << prefix << "_JS"        << "\t" << JS        << endl;
}

size_t max_length(const vector<string>& A)
{
  size_t mx = 0;
  BOOST_FOREACH(const string& a, A)
    if (a.size() > mx)
      mx = a.size();
  return mx;
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
  KmerStats kmer_stats(std::max(max_length(A), max_length(B)));
  compute_kmer_stats(kmer_stats, A, A_rc, tau_A, B, B_rc, tau_B);
  print_kmer_stats<wei>(kmer_stats, "weighted_" + prefix);
  print_kmer_stats<unw>(kmer_stats, "unweighted_" + prefix);
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

  std::cout << "summarize_kmerpair_version\t1" << std::endl;

  compute_and_print_kmer_stats(A, A_rc, tau_A, B, B_rc, tau_B, "kmerpair");
}
