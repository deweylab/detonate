#pragma once
#define N_POLICY 2
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
#include "mask.hh"
#include "util.hh"

template<typename AlignmentType>
struct BestTuple
{
  int num_identity;
  AlignmentType al;
  bool is_empty;
  BestTuple() : is_empty(true) {}
};

inline bool is_good_enough(const psl_alignment& /*al*/) { throw std::runtime_error("is_good_enough not implemented for psl"); }
inline bool is_good_enough(const blast_alignment& al) { return al.evalue() <= 1e-5; }

// is al better than ref?
template<typename T, typename S>
inline bool is_better(const T& al, const S& ref)
{
  return al.num_identity() > ref.num_identity;
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
        bt.num_identity        = al.num_identity();
        bt.al                  = al;
        bt.is_empty            = false;
      }
    }
  }
}

// For each b in B, figure out which alignment from a -> b (for some a in A) is
// best.
//
// Preconditions:
// - best_to_B needs to be default-initialized of size B.size()
template<typename AlignmentType>
void read_alignments_and_filter_by_best_to_B(std::vector<BestTuple<AlignmentType> >&    best_to_B, 
                                             typename AlignmentType::input_stream_type& input_stream,
                                             const std::map<std::string, size_t>&       B_names_to_idxs,
                                             bool                                       strand_specific)
{
  AlignmentType al;
  while (input_stream >> al) {
    size_t b_idx = B_names_to_idxs.find(al.b_name())->second;
    BestTuple<AlignmentType>& bt = best_to_B[b_idx];
    if (is_valid(al, strand_specific) && is_good_enough(al)) {
      if (bt.is_empty || is_better(al, bt)) {
        bt.num_identity        = al.num_identity();
        bt.al                  = al;
        bt.is_empty            = false;
      }
    }
  }
}

// Preconditions:
// - reciprocal_alignments should be NULL-initialized, of size B.size().
template<typename AlignmentType>
void filter_only_reciprocal_alignments(std::vector<const AlignmentType *>&           reciprocal_alignments,
                                       const std::vector<BestTuple<AlignmentType> >& best_from_A,
                                       const std::vector<BestTuple<AlignmentType> >& best_to_B,
                                       const std::map<std::string, size_t>&          B_names_to_idxs)
{
  for (size_t a_idx = 0; a_idx < best_from_A.size(); ++a_idx) {
    const BestTuple<AlignmentType>& bt_a = best_from_A[a_idx];
    if (!bt_a.is_empty) {
      size_t b_idx = B_names_to_idxs.find(bt_a.al.b_name())->second;
      const BestTuple<AlignmentType>& bt_b = best_to_B[b_idx];
      if (!bt_b.is_empty) {
        if (bt_a.al == bt_b.al) {
          // If all the above are satisfied, then this is a reciprocal
          // alignment, so add it to the set of reciprocal alignments.
          assert(reciprocal_alignments[b_idx] == NULL);
          reciprocal_alignments[b_idx] = &(bt_a.al);
        }
      }
    }
  }
}

template<typename AlignmentType>
void compute_alignment_stats(std::vector<double>&                      B_frac_ones,
                             const std::vector<const AlignmentType *>& reciprocal_alignments,
                             const std::vector<std::string>&           A,
                             const std::vector<std::string>&           B,
                             const std::vector<double>&                tau_B,
                             const std::map<std::string, size_t>&      A_names_to_idxs)
{
  std::vector<size_t> empty_mismatches;
  std::vector<size_t> perm = make_random_permutation(B.size()); // for better load balancing (in case e.g. seqs are ordered by length)
  #pragma omp parallel
  {

    #pragma omp for
    for (int b_pre_idx = 0; b_pre_idx < static_cast<int>(B.size()); ++b_pre_idx) {

      size_t b_idx = perm[b_pre_idx];
      mask b_mask(B[b_idx].size());

      const AlignmentType *al = reciprocal_alignments[b_idx];
      if (al != NULL) {
        size_t a_idx = A_names_to_idxs.find(al->a_name())->second;
        BOOST_FOREACH(const alignment_segment& seg, al->segments(A[a_idx], B[b_idx])) {
          b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, empty_mismatches.begin(), empty_mismatches.end());
          //b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
        }
      }

      B_frac_ones[b_idx] = (1.0 * b_mask.num_ones()) / B[b_idx].size();

    }

  } // omp parallel
}

void print_plot_output(std::vector<double>& B_frac_ones,
                       const std::string& plot_output_fname)
{
  std::ofstream f(plot_output_fname.c_str());
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
    ("A-to-B", po::value<std::string>()->required(), "The alignments of A to B.")
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
  std::string plot_output_fname = vm["plot-output"].as<std::string>();
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
  std::vector<BestTuple<AlignmentType> > best_from_A(A_card);
  {
    typename AlignmentType::input_stream_type A_to_B_is(open_or_throw(vm["A-to-B"].as<std::string>()));
    read_alignments_and_filter_by_best_from_A(best_from_A, A_to_B_is, A_names_to_idxs, strand_specific);
  }

  std::cerr << "Reading alignments and filtering them by B" << std::endl;
  std::vector<BestTuple<AlignmentType> > best_to_B(B_card);
  {
    typename AlignmentType::input_stream_type A_to_B_is(open_or_throw(vm["A-to-B"].as<std::string>()));
    read_alignments_and_filter_by_best_to_B(best_to_B, A_to_B_is, B_names_to_idxs, strand_specific);
  }

  std::cerr << "Filtering reciprocal alignments" << std::endl;
  std::vector<const AlignmentType *> reciprocal_alignments(B_card, NULL);
  filter_only_reciprocal_alignments(reciprocal_alignments, best_from_A, best_to_B, B_names_to_idxs);

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_B(B_card, 1.0/B_card);

  std::cerr << "Computing unweighted filtered stats" << std::endl;
  std::vector<double> B_frac_ones(B_card);
  compute_alignment_stats(B_frac_ones, reciprocal_alignments, A, B, unif_tau_B, A_names_to_idxs);
  print_plot_output(B_frac_ones, plot_output_fname);
}
