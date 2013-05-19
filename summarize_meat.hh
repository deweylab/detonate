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

struct stats_tuple
{
  double pair, nucl, tran;
};

double compute_F1(double precis, double recall)
{
  if (precis == 0.0 && recall == 0.0)
    return 0.0;
  else
    return 2*precis*recall/(precis+recall);
}

void print_stats(const stats_tuple& precis, const stats_tuple& recall, const std::string& prefix)
{
  std::cout << prefix << "_pair_precision\t" << precis.pair                          << std::endl;
  std::cout << prefix << "_pair_recall\t"    << recall.pair                          << std::endl;
  std::cout << prefix << "_pair_F1\t"        << compute_F1(precis.pair, recall.pair) << std::endl;

  std::cout << prefix << "_nucl_precision\t" << precis.nucl                          << std::endl;
  std::cout << prefix << "_nucl_recall\t"    << recall.nucl                          << std::endl;
  std::cout << prefix << "_nucl_F1\t"        << compute_F1(precis.nucl, recall.nucl) << std::endl;

  std::cout << prefix << "_tran_precision\t" << precis.tran                          << std::endl;
  std::cout << prefix << "_tran_recall\t"    << recall.tran                          << std::endl;
  std::cout << prefix << "_tran_F1\t"        << compute_F1(precis.tran, recall.tran) << std::endl;
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

#define min_frac_id 0.95
//#define min_frac_id 0.8

template<typename T>
inline bool is_good_enough_helper(const T& al)
{
  #if (GOOD_POLICY == 1)
  return al.frac_identity_wrt_a() >= min_frac_id && al.frac_indel_wrt_a() <= 0.0;
  #elif (GOOD_POLICY == 2)
  return (al.frac_identity_wrt_a() >= min_frac_id && al.frac_indel_wrt_a() <= 0.0) ||
         (al.frac_identity_wrt_b() >= min_frac_id && al.frac_indel_wrt_b() <= 0.0);
  #elif (GOOD_POLICY == 3)
  return al.frac_identity_wrt_a() >= min_frac_id || al.frac_identity_wrt_b() >= min_frac_id;
  #elif (GOOD_POLICY == 4)
  return true;
  #elif (GOOD_POLICY == 5)
  return al.frac_identity_wrt_b() >= min_frac_id;
  #elif (GOOD_POLICY == 6)
  return al.frac_identity_wrt_a() >= min_frac_id && al.frac_identity_wrt_b() >= min_frac_id;
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

inline bool is_valid(const blast_alignment& /*al*/, bool /*strand_specific*/)
{
  return true;
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

// Preconditions:
// - best_from_A should be of size 0
template<typename AlignmentType>
void read_alignments_without_filtering(std::vector<BestTuple<AlignmentType> >&    best_from_A, 
                                       typename AlignmentType::input_stream_type& input_stream,
                                       bool                                       strand_specific)
{
  AlignmentType al;
  while (input_stream >> al) {
    if (is_valid(al, strand_specific) && is_good_enough(al)) {
      BestTuple<AlignmentType> bt;
      bt.frac_identity_wrt_a = al.frac_identity_wrt_a();
      bt.frac_identity_wrt_b = al.frac_identity_wrt_b();
      bt.frac_indel_wrt_a    = al.frac_indel_wrt_a();
      bt.frac_indel_wrt_b    = al.frac_indel_wrt_b();
      bt.al                  = al;
      bt.is_empty            = false;
      best_from_A.push_back(bt);
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
stats_tuple compute_alignment_stats(std::vector<double>& B_frac_ones,
                                    const std::vector<std::vector<const AlignmentType *> >& best_to_B,
                                    const std::vector<std::string>&                         A,
                                    const std::vector<std::string>&                         B,
                                    const std::vector<double>&                              tau_B,
                                    const std::map<std::string, size_t>&                    A_names_to_idxs)
{
  std::vector<size_t> perm = make_random_permutation(B.size()); // for better load balancing (in case e.g. seqs are ordered by length)

  double pair_recall_numer = 0.0, pair_recall_denom = 0.0,
         nucl_recall_numer = 0.0, nucl_recall_denom = 0.0,
         tran_recall = 0.0;

  #pragma omp parallel
  {
    double private_pair_recall_numer = 0.0, private_pair_recall_denom = 0.0,
           private_nucl_recall_numer = 0.0, private_nucl_recall_denom = 0.0,
           private_tran_recall       = 0.0;

    #pragma omp for
    for (int b_pre_idx = 0; b_pre_idx < static_cast<int>(B.size()); ++b_pre_idx) {

      size_t b_idx = perm[b_pre_idx];
      PairsetType b_pairset(B[b_idx].size());
      mask b_mask(B[b_idx].size());

      BOOST_FOREACH(const AlignmentType *al, best_to_B[b_idx]) {
        size_t a_idx = A_names_to_idxs.find(al->a_name())->second;
        BOOST_FOREACH(const alignment_segment& seg, al->segments(A[a_idx], B[b_idx])) {
          b_pairset.add_square_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
          b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
        }
      }

      private_pair_recall_numer += tau_B[b_idx] * b_pairset.size();
      private_pair_recall_denom += tau_B[b_idx] * B[b_idx].size() * (B[b_idx].size() + 1) / 2;

      private_nucl_recall_numer += tau_B[b_idx] * b_mask.num_ones();
      private_nucl_recall_denom += tau_B[b_idx] * B[b_idx].size();

      B_frac_ones[b_idx] = (1.0 * b_mask.num_ones()) / B[b_idx].size();
      if (B_frac_ones[b_idx] >= 0.95)
        private_tran_recall += tau_B[b_idx];

    }

    #pragma omp critical
    {
      pair_recall_numer += private_pair_recall_numer;
      pair_recall_denom += private_pair_recall_denom;
      nucl_recall_numer += private_nucl_recall_numer;
      nucl_recall_denom += private_nucl_recall_denom;
      tran_recall += private_tran_recall;
    }

  } // omp parallel

  stats_tuple recall;
  recall.pair = pair_recall_numer / pair_recall_denom;
  recall.nucl = nucl_recall_numer / nucl_recall_denom;
  recall.tran = tran_recall;
  return recall;
}

template<typename AlignmentType>
void induce_prot_expression(std::vector<double>&                                    tau_B,
                            const std::vector<std::vector<const AlignmentType *> >& best_to_B,
                            const std::vector<double>&                              tau_A,
                            const std::map<std::string, size_t>&                    A_names_to_idxs,
                            const std::vector<std::string>&                         A,
                            const std::vector<std::string>&                         B)

{
  // Let IAC(c,p) be the interval of the contig c that is aligned to protein c (an interval of nucleotides)
  // Let IAP(c,p) be the interval of the protein p that is aligned to contig c (an interval of amino acids)
  // Then, let
  // n(p) = sum_{c in C(p)} (expr(c) * |IAC(c,p)|) / (3 * |union_{c in C(p)} (IAP(c,p))|)
  // and normalize with z = sum_p n(p), expr(p) = n(p)/z, as before.
  #pragma omp parallel for
  for (int b_idx = 0; b_idx < static_cast<int>(tau_B.size()); ++b_idx) {
    // Compute numer = sum_{c in C(p)} (expr(c) * |IAC(c,p)|)
    // and     denom = (3 * |union_{c in C(p)} (IAP(c,p))|)
    double numer = 0.0;
    std::set<size_t> union_IAP;
    BOOST_FOREACH(const AlignmentType *al, best_to_B[b_idx]) {
      size_t a_idx = A_names_to_idxs.find(al->a_name())->second;
      size_t IAC_size = 0;
      BOOST_FOREACH(const alignment_segment& seg, al->segments(A[a_idx], B[b_idx])) {
        size_t a_start = seg.a_start, a_end = seg.a_end;
        size_t b_start = seg.b_start, b_end = seg.b_end;
        if (a_start > a_end) std::swap(a_start, a_end);
        if (b_start > b_end) std::swap(b_start, b_end);
        //assert(a_end - a_start + 1 == 3*(b_end - b_start + 1)); // i.e., a is nucl, b is prot
        IAC_size += a_end - a_start + 1;
        for (size_t i = b_start; i != b_end; ++i)
          union_IAP.insert(i);
      }
      numer += tau_A[a_idx] * IAC_size;
    }
    double denom = 3 * union_IAP.size();
    tau_B[b_idx] = (numer == 0.0) ? 0.0 : numer / denom;
  }
  // Normalize
  double z = 0.0;
  for (size_t b_idx = 0; b_idx < tau_B.size(); ++b_idx) z += tau_B[b_idx];
  for (size_t b_idx = 0; b_idx < tau_B.size(); ++b_idx) tau_B[b_idx] /= z;
}

void print_plot_output(std::vector<double>& B_frac_ones,
                       const std::string& plot_output_prefix,
                       const std::string& plot_output_suffix)
{
  std::string s = plot_output_prefix + "_" + plot_output_suffix;
  std::ofstream f(s.c_str());
  BOOST_FOREACH(double x, B_frac_ones)
    f << x << std::endl;
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
    ("B-expr", po::value<std::string>(),             "The oracleset expression, as produced by RSEM in a file called *.isoforms.results.")
    ("induce-B-expr",                                "Induce oracleset expression from assembly expression and alignments.")
    ("A-to-B", po::value<std::string>()->required(), "The alignments of A to B.")
    ("B-to-A", po::value<std::string>()->required(), "The alignments of B to A.")
    ("alignment-type", po::value<std::string>()->required(), "The type of alignments used, either 'blast' or 'psl'.")
    ("plot-output", po::value<std::string>()->required(),   "File where plot values will be written.")
    ("strand-specific",                              "Ignore alignments that are to the reverse strand.")
  ;

  try {

    //po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cerr << desc << std::endl;
      exit(1);
    }

    po::notify(vm);

    if (!vm.count("B-expr") && !vm.count("induce-B-expr"))
      throw po::error("Either --B-expr or --induce-B-expr is required.");

    if (vm["alignment-type"].as<std::string>() != "blast" &&
        vm["alignment-type"].as<std::string>() != "psl")
      throw po::error("Option --alignment-type needs to be 'blast' or 'psl'.");

    if (vm.count("induce-B-expr") && vm["alignment-type"].as<std::string>() != "blast")
      throw po::error("Option --induce-B-expr only makes sense if --alignment-type=blast, i.e., if A is of nucleotide type and B is of protein type.");

  } catch (std::exception& x) {
    std::cerr << "Error: " << x.what() << std::endl;
    std::cerr << desc << std::endl;
    exit(1);
  }
}

template<typename AlignmentType>
void main_1(const boost::program_options::variables_map& vm)
{
  std::string plot_output_prefix = vm["plot-output"].as<std::string>();
  bool strand_specific = vm.count("strand-specific");

  std::cerr << "Reading the sequences" << std::endl;
  std::vector<std::string> A, B;
  std::vector<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  read_fasta_names_and_seqs(vm["A-seqs"].as<std::string>(), A, A_names, A_names_to_idxs);
  read_fasta_names_and_seqs(vm["B-seqs"].as<std::string>(), B, B_names, B_names_to_idxs);
  size_t A_card = A.size();
  size_t B_card = B.size();

  std::cerr << "Reading alignments and filtering them by A" << std::endl;
  std::vector<BestTuple<AlignmentType> > best_from_A(A_card), best_from_B(B_card);
  typename AlignmentType::input_stream_type A_to_B_is(open_or_throw(vm["A-to-B"].as<std::string>()));
  typename AlignmentType::input_stream_type B_to_A_is(open_or_throw(vm["B-to-A"].as<std::string>()));
  read_alignments_and_filter_by_best_from_A(best_from_A, A_to_B_is, A_names_to_idxs, strand_specific);
  read_alignments_and_filter_by_best_from_A(best_from_B, B_to_A_is, B_names_to_idxs, strand_specific);

  // std::cerr << "Reading alignments without filtering them by A" << std::endl;
  // std::vector<BestTuple<AlignmentType> > all_from_A, all_from_B;
  // typename AlignmentType::input_stream_type A_to_B_is2(open_or_throw(vm["A-to-B"].as<std::string>()));
  // typename AlignmentType::input_stream_type B_to_A_is2(open_or_throw(vm["B-to-A"].as<std::string>()));
  // read_alignments_without_filtering(all_from_A, A_to_B_is2, strand_specific);
  // read_alignments_without_filtering(all_from_B, B_to_A_is2, strand_specific);

  std::cerr << "Clustering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > clustered_best_to_B(B_card), clustered_best_to_A(A_card);
  cluster_best_alignments_to_B(clustered_best_to_B, best_from_A, B_names_to_idxs);
  cluster_best_alignments_to_B(clustered_best_to_A, best_from_B, A_names_to_idxs);

  std::cerr << "Filtering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > filtered_best_to_B(B_card), filtered_best_to_A(A_card);
  filter_by_best_alignment_to_B(filtered_best_to_B, best_from_A, B_names_to_idxs);
  filter_by_best_alignment_to_B(filtered_best_to_A, best_from_B, A_names_to_idxs);

  // std::cerr << "Clustering 'all' alignments by B" << std::endl;
  // std::vector<std::vector<const AlignmentType *> > jumbled_to_B(B_card), jumbled_to_A(A_card);
  // cluster_best_alignments_to_B(jumbled_to_B, all_from_A, B_names_to_idxs);
  // cluster_best_alignments_to_B(jumbled_to_A, all_from_B, A_names_to_idxs);

  std::cerr << "Reading transcript-level expression for A" << std::endl;
  std::vector<double> real_tau_A(A_card), real_tau_B(B_card);
  std::string A_expr_fname = vm["A-expr"].as<std::string>();
  read_transcript_expression(A_expr_fname, real_tau_A, A_names_to_idxs);
  if (vm.count("induce-B-expr")) {
    std::cerr << "Inducing transcript-level expression for B" << std::endl;
    induce_prot_expression(real_tau_B, clustered_best_to_B, real_tau_A, A_names_to_idxs, A, B);
  } else {
    std::cerr << "Reading transcript-level expression for B" << std::endl;
    std::string B_expr_fname = vm["B-expr"].as<std::string>();
    read_transcript_expression(B_expr_fname, real_tau_B, B_names_to_idxs);
  }

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_A(A_card, 1.0/A_card);
  std::vector<double> unif_tau_B(B_card, 1.0/B_card);

  std::cout << "summarize_version_11\t0" << std::endl;
  std::cout << "summarize_min_frac_identity_" << min_frac_id << "\t0" << std::endl;

  stats_tuple recall, precis;
  std::vector<double> B_frac_ones(B_card), A_frac_ones(A_card);

  std::cerr << "Computing weighted clustered stats" << std::endl;
  recall = compute_alignment_stats<smart_pairset>(B_frac_ones, clustered_best_to_B, A, B, real_tau_B, A_names_to_idxs);
  precis = compute_alignment_stats<smart_pairset>(A_frac_ones, clustered_best_to_A, B, A, real_tau_A, B_names_to_idxs);
  print_stats(precis, recall, "weighted_clustered");

  std::cerr << "Computing unweighted clustered stats" << std::endl;
  recall = compute_alignment_stats<smart_pairset>(B_frac_ones, clustered_best_to_B, A, B, unif_tau_B, A_names_to_idxs);
  precis = compute_alignment_stats<smart_pairset>(A_frac_ones, clustered_best_to_A, B, A, unif_tau_A, B_names_to_idxs);
  print_stats(precis, recall, "unweighted_clustered");
  print_plot_output(B_frac_ones, plot_output_prefix, "clustered");

  // std::cerr << "Computing weighted jumbled stats" << std::endl;
  // recall = compute_alignment_stats<smart_pairset>(B_frac_ones, jumbled_to_B, A, B, real_tau_B, A_names_to_idxs);
  // precis = compute_alignment_stats<smart_pairset>(A_frac_ones, jumbled_to_A, B, A, real_tau_A, B_names_to_idxs);
  // print_stats(precis, recall, "weighted_jumbled");
  //
  // std::cerr << "Computing unweighted jumbled stats" << std::endl;
  // recall = compute_alignment_stats<smart_pairset>(B_frac_ones, jumbled_to_B, A, B, unif_tau_B, A_names_to_idxs);
  // precis = compute_alignment_stats<smart_pairset>(A_frac_ones, jumbled_to_A, B, A, unif_tau_A, B_names_to_idxs);
  // print_stats(precis, recall, "unweighted_jumbled");
  // print_plot_output(B_frac_ones, plot_output_prefix, "jumbled");

  if (vm.count("induce-B-expr")) {
    std::cerr << "Inducing transcript-level expression for B with filtered stats" << std::endl;
    induce_prot_expression(real_tau_B, filtered_best_to_B, real_tau_A, A_names_to_idxs, A, B);
  }

  std::cerr << "Computing weighted filtered stats" << std::endl;
  recall = compute_alignment_stats<smart_pairset>(B_frac_ones, filtered_best_to_B, A, B, real_tau_B, A_names_to_idxs);
  precis = compute_alignment_stats<smart_pairset>(A_frac_ones, filtered_best_to_A, B, A, real_tau_A, B_names_to_idxs);
  print_stats(precis, recall, "weighted_filtered");

  std::cerr << "Computing unweighted filtered stats" << std::endl;
  recall = compute_alignment_stats<smart_pairset>(B_frac_ones, filtered_best_to_B, A, B, unif_tau_B, A_names_to_idxs);
  precis = compute_alignment_stats<smart_pairset>(A_frac_ones, filtered_best_to_A, B, A, unif_tau_A, B_names_to_idxs);
  print_stats(precis, recall, "unweighted_filtered");
  print_plot_output(B_frac_ones, plot_output_prefix, "filtered"); 
}
