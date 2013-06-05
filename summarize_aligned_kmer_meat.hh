#pragma once
#include <iostream>
#include <fstream>
#include <sstream>
#include <iterator>
#include <vector>
#include <list>
#include <omp.h>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include "blast.hh"
#include "psl.hh"
#include "pairset.hh"
#include "mask.hh"
#include "util.hh"
#include <sparsehash/sparse_hash_map>

#define MIN_FRAC_ID (MIN_PCT_ID/100.0)

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
        double c = 0.5 * tau_A[i];
        string::const_iterator beg = a.begin();
        string::const_iterator end = a.begin() + kmerlen;
        for (; end != a.end(); ++beg, ++end) {
          KmerKey kmer_key(beg, end);
          stats[kmer_key].probs[A_or_B] += c; // relies on default init to 0
        }
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

double compute_F1(double precis, double recall)
{
  if (precis == 0.0 && recall == 0.0)
    return 0.0;
  else
    return 2*precis*recall/(precis+recall);
}

void print_kmer_stats(
    KmerStats::container_type& stats,
    const string& prefix)
{
  typedef std::pair<const KmerStats::key_type, KmerStats::value_type> X;

  // First, compute and print stats based on the UNnormalized distribution.

  {
    size_t num_shared = 0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      if (i.probs[0] != 0 && i.probs[1] != 0)
        ++num_shared;
    }
    cout << prefix << "_num_shared" << "\t" << num_shared << endl;
  }

  {
    double numer = 0.0, denom_precis = 0.0, denom_recall = 0.0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      numer += std::min(i.probs[0], i.probs[1]); // frac "retrieved" and "relevant"
      denom_precis += i.probs[0]; // frac "retrieved"
      denom_recall += i.probs[1]; // frac "relevant"
    }
    double precis = numer / denom_precis;
    double recall = numer / denom_recall;
    double F1 = compute_F1(precis, recall);
    cout << prefix << "_precision" << "\t" << precis << endl;
    cout << prefix << "_recall"    << "\t" << recall << endl;
    cout << prefix << "_F1"        << "\t" << F1     << endl;
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
    cout << prefix << "_probs_KL_A_to_M" << "\t" << KL_A_to_M << endl;
    cout << prefix << "_probs_KL_B_to_M" << "\t" << KL_B_to_M << endl;
    cout << prefix << "_probs_JS" << "\t" << JS << endl;
  }

  {
    double probs_hellinger = 0.0;
    BOOST_FOREACH(const X& x, stats) {
      const KmerInfo& i = x.second;
      double y = sqrt(i.probs[0]) - sqrt(i.probs[1]);
      probs_hellinger += y * y;
    }
    probs_hellinger = sqrt(probs_hellinger)/sqrt(2.0);
    cout << prefix << "_probs_hellinger" << "\t" << probs_hellinger << endl;
  }

}

template<typename AlignmentType>
struct BestTuple
{
  double frac_identity_wrt_a;
  double frac_identity_wrt_b;
  double frac_indel_wrt_a;
  double frac_indel_wrt_b;
  AlignmentType al;
  bool is_empty;
  BestTuple() : is_empty(true) {}
};

template<typename T>
inline bool is_good_enough_helper(const T& al)
{
  #if (GOOD_POLICY == 1)
  return al.frac_identity_wrt_a() >= MIN_FRAC_ID && al.frac_indel_wrt_a() <= 0.0;
  #elif (GOOD_POLICY == 2)
  return (al.frac_identity_wrt_a() >= MIN_FRAC_ID && al.frac_indel_wrt_a() <= 0.0) ||
         (al.frac_identity_wrt_b() >= MIN_FRAC_ID && al.frac_indel_wrt_b() <= 0.0);
  #elif (GOOD_POLICY == 3)
  return al.frac_identity_wrt_a() >= MIN_FRAC_ID || al.frac_identity_wrt_b() >= MIN_FRAC_ID;
  #elif (GOOD_POLICY == 4)
  return true;
  #elif (GOOD_POLICY == 5)
  return al.frac_identity_wrt_b() >= MIN_FRAC_ID;
  #elif (GOOD_POLICY == 6)
  return al.frac_identity_wrt_a() >= MIN_FRAC_ID && al.frac_identity_wrt_b() >= MIN_FRAC_ID;
  #else
  #error "need to define GOOD_POLICY"
  #endif
}

inline bool is_good_enough(const psl_alignment& al) { return is_good_enough_helper(al); }
inline bool is_good_enough(const blast_alignment& al) { return al.evalue() <= 1e-5 && is_good_enough_helper(al); }

// is al better than ref?
template<typename T, typename S>
inline bool is_better(const T& al, const S& ref)
{
  #if (BETTER_POLICY == 1)
  return (al.frac_identity_wrt_a() > ref.frac_identity_wrt_a
          || (al.frac_identity_wrt_a() == ref.frac_identity_wrt_a && al.frac_indel_wrt_a() == ref.frac_indel_wrt_a));
  #elif (BETTER_POLICY == 2)
  double d1 = al.frac_identity_wrt_a() - ref.frac_identity_wrt_a;
  double d2 = -(al.frac_indel_wrt_a()  - ref.frac_indel_wrt_a);
  double d3 = al.frac_identity_wrt_b() - ref.frac_identity_wrt_b;
  double d4 = -(al.frac_indel_wrt_b()  - ref.frac_indel_wrt_b);
  return d1 > 0 || (d1 == 0 && d2 > 0)
                || (d1 == 0 && d2 == 0 && d3 > 0)
                || (d1 == 0 && d2 == 0 && d3 == 0 && d4 > 0);
  #elif (BETTER_POLICY == 3)
  double d1 = al.frac_identity_wrt_a() - ref.frac_identity_wrt_a;
  double d2 = al.frac_identity_wrt_b() - ref.frac_identity_wrt_b;
  return d1 > 0 || (d1 == 0 && d2 > 0);
  #elif (BETTER_POLICY == 4)
  double d1 = 1.0*al.num_identity() - 1.0*ref.al.num_identity();
  double d2 = -(1.0*al.num_indel()  - 1.0*ref.al.num_indel());
  return d1 > 0 || (d1 == 0 && d2 > 0);
  #else
  #error "need to define BETTER_POLICY"
  #endif
}

inline bool is_valid(const psl_alignment& al, bool strand_specific)
{
  if (!strand_specific)
    return true;
  else
    return al.strand() == "+";
}

inline bool is_valid(const blast_alignment& /*al*/, bool strand_specific)
{
  if (!strand_specific)
    return true;
  else
    throw std::runtime_error("Strand-specificity has not been implemented yet for blast alignments.");
}

// For each a in A, figure out which alignment from a -> b (for some b in B) is
// best.
//
// Preconditions:
// - best_from_A needs to be default-initialized of size A.size()
template<typename AlignmentType>
void read_alignments_and_filter_by_best_from_A(std::vector<BestTuple<AlignmentType> >&    best_from_A, 
                                               typename AlignmentType::input_stream_type& input_stream,
                                               const std::map<std::string, size_t>&       A_names_to_idxs,
                                               bool                                       strand_specific)
{
  AlignmentType al;
  while (input_stream >> al) {
    size_t a_idx = A_names_to_idxs.find(al.a_name())->second;
    BestTuple<AlignmentType>& bt = best_from_A[a_idx];
    if (is_valid(al, strand_specific) && is_good_enough(al)) {
      if (bt.is_empty || is_better(al, bt)) {
        bt.frac_identity_wrt_a = al.frac_identity_wrt_a();
        bt.frac_identity_wrt_b = al.frac_identity_wrt_b();
        bt.frac_indel_wrt_a    = al.frac_indel_wrt_a();
        bt.frac_indel_wrt_b    = al.frac_indel_wrt_b();
        bt.al                  = al;
        bt.is_empty            = false;
      }
    }
  }
}

// For each b in B, collect all the best alignments from a's that have b as
// target. The purpose here is to facilitate (i) make parallel processing and
// (ii) reduced memory usage later on.
//
// Preconditions:
// - best_to_B needs to be default-initialized of size B.size()
// - best_from_A needs to be the result of read_alignments_and_filter_by_best_from_A
template<typename AlignmentType>
void cluster_best_alignments_to_B(std::vector<std::vector<const AlignmentType *> >& best_to_B,
                                  const std::vector<BestTuple<AlignmentType> >&     best_from_A,
                                  const std::map<std::string, size_t>&              B_names_to_idxs)
{
  BOOST_FOREACH(const BestTuple<AlignmentType>& bt, best_from_A) {
    if (!bt.is_empty) {
      size_t b_idx = B_names_to_idxs.find(bt.al.b_name())->second;
      best_to_B[b_idx].push_back(&(bt.al));
    }
  }
}

template<typename AlignmentType>
void filter_by_best_alignment_to_B(std::vector<std::vector<const AlignmentType *> >& best_to_B,
                                   const std::vector<BestTuple<AlignmentType> >&     best_from_A,
                                   const std::map<std::string, size_t>&              B_names_to_idxs)
{
  BOOST_FOREACH(const BestTuple<AlignmentType>& bt, best_from_A) {
    if (!bt.is_empty) {
      size_t b_idx = B_names_to_idxs.find(bt.al.b_name())->second;
      if (best_to_B[b_idx].size() == 0)
        best_to_B[b_idx].push_back(&(bt.al));
      else {
        assert(best_to_B[b_idx].size() == 1);
        const AlignmentType *old = best_to_B[b_idx][0];
        if (!is_better(*old, bt)) // note: not exactly the same as is_better(bt, old) would be
          best_to_B[b_idx][0] = &(bt.al);
      }
    }
  }
}

template<typename PairsetType, typename AlignmentType>
void compute_alignment_stats(KmerStats::container_type&                              stats,
                             size_t                                                  A_or_B,
                             const std::vector<std::vector<const AlignmentType *> >& best_to_B,
                             const std::vector<std::string>&                         A,
                             const std::vector<std::string>&                         B,
                             const std::vector<std::string>&                         B_rc,
                             const std::vector<double>&                              tau_B,
                             const std::map<std::string, size_t>&                    A_names_to_idxs,
                             size_t                                                  kmerlen)
{
  for (int b_idx = 0; b_idx < static_cast<int>(B.size()); ++b_idx) {
    BOOST_FOREACH(const AlignmentType *al, best_to_B[b_idx]) {
      size_t a_idx = A_names_to_idxs.find(al->a_name())->second;
      BOOST_FOREACH(const alignment_segment& seg, al->segments(A[a_idx], B[b_idx])) {

        size_t start = seg.b_start;
        size_t end   = seg.b_end;

        const std::string *b;
        if (start > end) {
          // E.g.,
          // - 76543210
          // + 01234567
          //      e  s
          // -->
          // s = (8 - 1) - 6 = 1
          // e = (8 - 1) - 3 = 4
          //assert(al->strand() == "-");
          start = (B[b_idx].size() - 1) - start;
          end   = (B[b_idx].size() - 1) - end;
          b = &(B_rc[b_idx]);
        } else {
          //assert(al->strand() == "+");
          b = &(B[b_idx]);
        }

        if (end + 1 - start >= kmerlen) {
          std::string::const_iterator start_it = b->begin() + start;
          std::string::const_iterator end_it   = b->begin() + end;
          std::string::const_iterator x = start_it;
          std::string::const_iterator y = start_it + kmerlen - 1;
          for (; y != end_it; ++x, ++y) {
            KmerKey kmer_key(x, y);
            stats[kmer_key].probs[A_or_B] += tau_B[b_idx]; // relies on default init to 0
          }
        }

      }
    }
  }
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
    ("A-to-B", po::value<std::string>()->required(), "The alignments of A to B.")
    ("B-to-A", po::value<std::string>()->required(), "The alignments of B to A.")
    ("alignment-type", po::value<std::string>()->required(), "The type of alignments used, either 'blast' or 'psl'.")
    ("strand-specific",                              "Ignore alignments that are to the reverse strand.")
    ("readlen", po::value<size_t>()->required(),     "The read length.")
  ;

  try {

    //po::variables_map vm;
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

    if (vm["alignment-type"].as<std::string>() != "blast" &&
        vm["alignment-type"].as<std::string>() != "psl")
      throw po::error("Option --alignment-type needs to be 'blast' or 'psl'.");

  } catch (std::exception& x) {
    std::cerr << "Error: " << x.what() << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}

template<typename AlignmentType>
void main_1(const boost::program_options::variables_map& vm)
{
  bool strand_specific = vm.count("strand-specific");
  size_t readlen = vm["readlen"].as<size_t>();

  std::cerr << "Reading the sequences" << std::endl;
  std::vector<std::string> A, B;
  std::vector<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  read_fasta_names_and_seqs(vm["A-seqs"].as<std::string>(), A, A_names, A_names_to_idxs);
  read_fasta_names_and_seqs(vm["B-seqs"].as<std::string>(), B, B_names, B_names_to_idxs);
  size_t A_card = A.size();
  size_t B_card = B.size();

  std::cerr << "Reverse complementing the sequences" << std::endl;
  std::vector<std::string> A_rc, B_rc;
  transform(A.begin(), A.end(), back_inserter(A_rc), reverse_complement);
  transform(B.begin(), B.end(), back_inserter(B_rc), reverse_complement);

  std::cerr << "Reading alignments and filtering them by A" << std::endl;
  std::vector<BestTuple<AlignmentType> > best_from_A(A_card), best_from_B(B_card);
  typename AlignmentType::input_stream_type A_to_B_is(open_or_throw(vm["A-to-B"].as<std::string>()));
  typename AlignmentType::input_stream_type B_to_A_is(open_or_throw(vm["B-to-A"].as<std::string>()));
  read_alignments_and_filter_by_best_from_A(best_from_A, A_to_B_is, A_names_to_idxs, strand_specific);
  read_alignments_and_filter_by_best_from_A(best_from_B, B_to_A_is, B_names_to_idxs, strand_specific);

  std::cerr << "Clustering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > clustered_best_to_B(B_card), clustered_best_to_A(A_card);
  cluster_best_alignments_to_B(clustered_best_to_B, best_from_A, B_names_to_idxs);
  cluster_best_alignments_to_B(clustered_best_to_A, best_from_B, A_names_to_idxs);

  std::cerr << "Filtering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > filtered_best_to_B(B_card), filtered_best_to_A(A_card);
  filter_by_best_alignment_to_B(filtered_best_to_B, best_from_A, B_names_to_idxs);
  filter_by_best_alignment_to_B(filtered_best_to_A, best_from_B, A_names_to_idxs);

  std::vector<double> real_tau_A(A_card), real_tau_B(B_card);
  if (!vm.count("no-expr")) {
    std::cerr << "Reading transcript-level expression for A and B" << std::endl;
    std::string A_expr_fname = vm["A-expr"].as<std::string>();
    std::string B_expr_fname = vm["B-expr"].as<std::string>();
    read_transcript_expression(A_expr_fname, real_tau_A, A_names_to_idxs);
    read_transcript_expression(B_expr_fname, real_tau_B, B_names_to_idxs);
  }

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_A(A_card, 1.0/A_card);
  std::vector<double> unif_tau_B(B_card, 1.0/B_card);

  std::cout << "summarize_aligned_kmer_version_1\t0" << std::endl;
  std::cout << "summarize_min_frac_identity_" << MIN_FRAC_ID << "\t0" << std::endl;

  if (!vm.count("no-expr")) {
    std::cerr << "Computing weighted clustered stats" << std::endl;
    KmerStats::container_type stats;
    compute_alignment_stats<smart_pairset>(stats, 1, clustered_best_to_B, A, B, B_rc, real_tau_B, A_names_to_idxs, readlen);
    compute_alignment_stats<smart_pairset>(stats, 0, clustered_best_to_A, B, A, A_rc, real_tau_A, B_names_to_idxs, readlen);
    print_kmer_stats(stats, "weighted_clustered_aligned_kmer");
  }

  {
    KmerStats::container_type stats;
    compute_alignment_stats<smart_pairset>(stats, 1, clustered_best_to_B, A, B, B_rc, unif_tau_B, A_names_to_idxs, readlen);
    compute_alignment_stats<smart_pairset>(stats, 0, clustered_best_to_A, B, A, A_rc, unif_tau_A, B_names_to_idxs, readlen);
    print_kmer_stats(stats, "unweighted_clustered_aligned_kmer");
  }

  if (!vm.count("no-expr")) {
    std::cerr << "Computing weighted filtered stats" << std::endl;
    KmerStats::container_type stats;
    compute_alignment_stats<smart_pairset>(stats, 1, filtered_best_to_B, A, B, B_rc, real_tau_B, A_names_to_idxs, readlen);
    compute_alignment_stats<smart_pairset>(stats, 0, filtered_best_to_A, B, A, A_rc, real_tau_A, B_names_to_idxs, readlen);
    print_kmer_stats(stats, "weighted_filtered_aligned_kmer");
  }

  {
    KmerStats::container_type stats;
    compute_alignment_stats<smart_pairset>(stats, 1, filtered_best_to_B, A, B, B_rc, unif_tau_B, A_names_to_idxs, readlen);
    compute_alignment_stats<smart_pairset>(stats, 0, filtered_best_to_A, B, A, A_rc, unif_tau_A, B_names_to_idxs, readlen);
    print_kmer_stats(stats, "unweighted_filtered_aligned_kmer");
  }
}
