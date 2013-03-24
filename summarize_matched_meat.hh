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
#define N_POLICY 2
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
struct tagged_alignment
{
  size_t a_idx, b_idx;
  std::vector<alignment_segment> segments;
  //typename AlignmentType::segments_type segments;
  double contribution;
  bool is_deleted;
};

template<typename AlignmentType>
struct compare_tagged_alignments
{
  bool operator()(const tagged_alignment<AlignmentType> *l1, const tagged_alignment<AlignmentType> *l2) const
  {
    return l1->contribution < l2->contribution;
  }
};

template<typename AlignmentType>
struct pair_helper
{
  const std::vector<std::string>& B;
  const std::vector<double>&      tau_B;
  double                          numer;

  pair_helper(const std::vector<std::string>& B,
              const std::vector<double>&      tau_B)
  : B(B), tau_B(tau_B), numer(0.0)
  {}

  double compute_contribution(const tagged_alignment<AlignmentType>& l) const
  {
    smart_pairset b_pairset(B[l.b_idx].size());
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_pairset.add_square_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_pairset.size();
  }

  void add_contribution_to_recall(const tagged_alignment<AlignmentType>& l) { numer += l.contribution; }

  double get_recall()
  {
    double denom = 0.0;
    for (size_t b_idx = 0; b_idx < B.size(); ++b_idx)
      denom += tau_B[b_idx] * B[b_idx].size() * (B[b_idx].size() + 1) / 2;
    return numer/denom;
  }
};

template<typename AlignmentType>
struct nucl_helper
{
  const std::vector<std::string>& B;
  const std::vector<double>&      tau_B;
  double                          numer;

  nucl_helper(const std::vector<std::string>& B,
              const std::vector<double>&      tau_B)
  : B(B), tau_B(tau_B), numer(0.0)
  {}

  double compute_contribution(const tagged_alignment<AlignmentType>& l) const
  {
    mask b_mask(B[l.b_idx].size());
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_mask.num_ones();
  }

  void add_contribution_to_recall(const tagged_alignment<AlignmentType>& l) { numer += l.contribution; }

  double get_recall()
  {
    double denom = 0.0;
    for (size_t b_idx = 0; b_idx < B.size(); ++b_idx)
      denom += tau_B[b_idx] * B[b_idx].size();
    return numer/denom;
  }
};

template<typename AlignmentType>
struct tran_helper
{
  const std::vector<std::string>& B;
  const std::vector<double>&      tau_B;
  std::vector<mask>               B_mask;
  std::vector<double>             B_frac_ones;

  tran_helper(const std::vector<std::string>& B,
              const std::vector<double>&      tau_B)
  : B(B), tau_B(tau_B), B_frac_ones(B.size())
  {
    BOOST_FOREACH(const std::string& b, B)
      B_mask.push_back(mask(b.size()));
  }

  double compute_contribution(const tagged_alignment<AlignmentType>& l) const
  {
    mask b_mask(B[l.b_idx].size());
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_mask.num_ones();
  }

  void add_contribution_to_recall(const tagged_alignment<AlignmentType>& l)
  {
    mask& m = B_mask[l.b_idx];
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      m.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
  }

  double get_recall()
  {
    double recall = 0.0;
    for (size_t b_idx = 0; b_idx < B.size(); ++b_idx) {
      B_frac_ones[b_idx] = (1.0 * B_mask[b_idx].num_ones()) / B[b_idx].size();
      if (B_frac_ones[b_idx] >= 0.95)
        recall += tau_B[b_idx];
    }
    return recall;
  }
};

#if 0
// Preconditions:
// - best_from_A should be of size 0
template<typename AlignmentType, typename HelperType>
void do_it_all_2(HelperType&                                helper,
                 typename AlignmentType::input_stream_type& input_stream,
                 const std::vector<std::string>&            A,
                 const std::vector<std::string>&            B,
                 const std::map<std::string, size_t>&       A_names_to_idxs,
                 const std::map<std::string, size_t>&       B_names_to_idxs)
{
  typedef tagged_alignment<AlignmentType> TAT;

  // Read alignments into a big vector L.
  AlignmentType al;
  TAT l;
  std::vector<TAT> L;
  while (input_stream >> al) {
    if (is_good_enough(al)) {
      l.a_idx = A_names_to_idxs.find(al.a_name())->second;
      l.b_idx = B_names_to_idxs.find(al.b_name())->second;
      typename AlignmentType::segments_type segs = al.segments(A[l.a_idx], B[l.b_idx]);
      l.segments.assign(segs.begin(), segs.end());
      l.contribution = helper.compute_contribution(l);
      l.is_deleted = false;
      L.push_back(l);
    }
  }

  // Cluster alignments by their targets a in A.
  size_t A_card = A_names_to_idxs.size();
  std::vector<std::vector<TAT *> > L_A(A_card);
  BOOST_FOREACH(TAT& l, L)
    L_A[l.a_idx].push_back(&l);

  // Cluster alignments by their targets b in B.
  size_t B_card = B_names_to_idxs.size();
  std::vector<std::vector<TAT *> > L_B(B_card);
  BOOST_FOREACH(TAT& l, L)
    L_B[l.b_idx].push_back(&l);

  // Make priority queue initially filled with all elements of L.
  compare_tagged_alignments<AlignmentType> comparator;
  std::priority_queue<
    TAT *,
    std::vector<TAT *>,
    compare_tagged_alignments<AlignmentType> > Q(L.begin(), L.end(), comparator);

  for (;;) {

    // Pop l1 from Q.
    TAT *l1 = Q.top();
    Q.pop();

    // Record l1's contribution.
    helper.add_contribution_to_recall(*l1);
    l1->is_deleted = true;

    // For each l2 in L_b whose coverage of b overlaps l's coverage of b:
    BOOST_FOREACH(TAT* l2, L_B[l1->b_idx]) {
      if (intersects<segment_ops_wrt_b>(l1->segments, l2->segments)) {

        // Mark l2 as deleted.
        l2->is_deleted = true;

        // Add a new alignment, l3, with segments "l2 - l1"
        L.push_back(TAT());
        TAT *l3 = &(L.back());
        l3->b_idx = l2->b_idx;
        l3->segments = subtract<segment_ops_wrt_b>(l2->segments, l1->segments);
        l3->contribution = helper.compute_contribution(*l3);
        l3->is_deleted = false;

        // If l3 is good enough, add it to Q. (Otherwise, get rid of it.)
        //if (is_good_enough(l3))
        if (true)
          Q.push(l3);
        else
          L.pop_back();
      }
    }

    // For each l2 in L_a whose coverage of a overlaps l's coverage of a:
    BOOST_FOREACH(TAT* l2, L_A[l1->a_idx]) {
      if (intersects<segment_ops_wrt_a>(l1->segments, l2->segments)) {

        // Mark l2 as deleted.
        l2->is_deleted = true;

        // Add a new alignment, l3, with segments "l2 - l1"
        L.push_back(TAT());
        TAT *l3 = &(L.back());
        l3->a_idx = l2->a_idx;
        l3->segments = subtract<segment_ops_wrt_a>(l2->segments, l1->segments);
        l3->contribution = helper.compute_contribution(*l3);
        l3->is_deleted = false;

        // If l3 is good enough, add it to Q. (Otherwise, get rid of it.)
        //if (is_good_enough(l3))
        if (true)
          Q.push(l3);
        else
          L.pop_back();
      }
    }

    // If Q is not empty, loop.
    if (Q.empty())
      break;
  }
}
#endif

inline bool is_good_enough_for_do_it_all(const std::vector<alignment_segment>& segs)
{
  size_t len = 0;
  BOOST_FOREACH(const alignment_segment& seg, segs) {
    if (seg.a_end >= seg.a_start)
      len += seg.a_end - seg.a_start + 1;
    else
      len += seg.a_start - seg.a_end + 1;
  }
  return len >= 30;
}

// Preconditions:
// - best_from_A should be of size 0
template<typename AlignmentType, typename HelperType>
void do_it_all(HelperType&                                helper,
               typename AlignmentType::input_stream_type& input_stream,
               const std::vector<std::string>&            A,
               const std::vector<std::string>&            B,
               const std::map<std::string, size_t>&       A_names_to_idxs,
               const std::map<std::string, size_t>&       B_names_to_idxs)
{
  typedef tagged_alignment<AlignmentType> TAT;

  // Read alignments, perfom initial filtering, and convert the alignments to
  // segments.
  AlignmentType al;
  TAT l;
  std::vector<TAT> L;
  while (input_stream >> al) {
    l.a_idx = A_names_to_idxs.find(al.a_name())->second;
    l.b_idx = B_names_to_idxs.find(al.b_name())->second;
    typename AlignmentType::segments_type segs = al.segments(A[l.a_idx], B[l.b_idx]);
    l.segments.assign(segs.begin(), segs.end());
    if (is_good_enough_for_do_it_all(l.segments))
      L.push_back(l);
  }

  // Compute contributions.
  #pragma omp parallel for
  for (int i = 0; i < static_cast<int>(L.size()); ++i) {
    TAT& l = L[i];
    l.contribution = helper.compute_contribution(l);
  }

  // Init vector of pointers to popped alignments.
  std::vector<std::vector<TAT *> > popped_by_A(A.size());
  std::vector<std::vector<TAT *> > popped_by_B(B.size());

  // Make priority queue initially filled with all elements of L.
  compare_tagged_alignments<AlignmentType> comparator;
  std::priority_queue<
    TAT *,
    std::vector<TAT *>,
    compare_tagged_alignments<AlignmentType> > Q(comparator);
  BOOST_FOREACH(TAT& l1, L)
    Q.push(&l1);

  for ( ; !Q.empty(); ) {

    // Pop l1 from Q.
    TAT *l1 = Q.top();
    Q.pop();

    // Subtract all previous alignments from l1.
    bool l1_has_changed = false;
    typename std::vector<TAT *>::const_iterator
        l2_a = popped_by_A[l1->a_idx].begin(), end_a = popped_by_A[l1->a_idx].end(),
        l2_b = popped_by_B[l1->b_idx].begin(), end_b = popped_by_B[l1->b_idx].end();
    for ( ; l2_a != end_a && l2_b != end_b; ) {

      bool l2_b_was_popped_first;
      if (l2_a != end_a && l2_b != end_b) {
        // If *l2_a < *l2_b, i.e., the current alignment in popped_by_A has
        // lower priority than the current alignment in popped_by_B, then *l2_b
        // was popped before *l2_a.
        l2_b_was_popped_first = comparator(*l2_a, *l2_b);
      } else {
        // Since we are still in the for loop (so one of the iterators is not
        // at the end), and we got to this else (so one of the iterators is at
        // the end), we just have to figure out which iterator is still not at
        // the end.
        l2_b_was_popped_first = (l2_b != end_b);
      }

      if (l2_b_was_popped_first) {
        if (intersects<segment_ops_wrt_b>(l1->segments, (*l2_b)->segments)) {
          subtract_in_place<segment_ops_wrt_b>(l1->segments, (*l2_b)->segments);
          l1_has_changed = true;
        }
        ++l2_b;
      } else {
        if (intersects<segment_ops_wrt_a>(l1->segments, (*l2_a)->segments)) {
          subtract_in_place<segment_ops_wrt_a>(l1->segments, (*l2_a)->segments);
          l1_has_changed = true;
        }
        ++l2_a;
      }

    }

    // If the alignment has changed, then put it back in the priority queue.
    if (l1_has_changed) {

      if (is_good_enough_for_do_it_all(l1->segments)) {
        l1->contribution = helper.compute_contribution(*l1);
        Q.push(l1);
      }

    // Otherwise, if the alignment has not changed, we actually process it.
    } else {

      // Record l1's contribution.
      helper.add_contribution_to_recall(*l1);

      // Add l1 to list of popped alignments.
      popped_by_A[l1->a_idx].push_back(l1);
      popped_by_B[l1->b_idx].push_back(l1);

    }

  }
}

template<typename AlignmentType>
stats_tuple do_it_all_wrapper(const std::string                    input_fname,
                              const std::vector<std::string>&      A,
                              const std::vector<std::string>&      B,
                              const std::vector<double>&           tau_B,
                              const std::map<std::string, size_t>& A_names_to_idxs,
                              const std::map<std::string, size_t>& B_names_to_idxs)
{
  stats_tuple st;

  {
    pair_helper<AlignmentType> ph(B, tau_B);
    typename AlignmentType::input_stream_type is(open_or_throw(input_fname));
    do_it_all<AlignmentType>(ph, is, A, B, A_names_to_idxs, B_names_to_idxs);
    st.pair = ph.get_recall();
  }

  {
    nucl_helper<AlignmentType> nh(B, tau_B);
    typename AlignmentType::input_stream_type is(open_or_throw(input_fname));
    do_it_all<AlignmentType>(nh, is, A, B, A_names_to_idxs, B_names_to_idxs);
    st.nucl = nh.get_recall();
  }

  {
    tran_helper<AlignmentType> th(B, tau_B);
    typename AlignmentType::input_stream_type is(open_or_throw(input_fname));
    do_it_all<AlignmentType>(th, is, A, B, A_names_to_idxs, B_names_to_idxs);
    st.tran = th.get_recall();
  }

  return st;
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
    ("A-to-B", po::value<std::string>()->required(), "The alignments of A to B.")
    ("B-to-A", po::value<std::string>()->required(), "The alignments of B to A.")
    ("alignment-type", po::value<std::string>()->required(), "The type of alignments used, either 'blast' or 'psl'.")
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
  //clock_t start;

  std::cerr << "Reading the sequences" << std::endl;
  std::vector<std::string> A, B;
  std::vector<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  read_fasta_names_and_seqs(vm["A-seqs"].as<std::string>(), A, A_names, A_names_to_idxs);
  read_fasta_names_and_seqs(vm["B-seqs"].as<std::string>(), B, B_names, B_names_to_idxs);
  size_t A_card = A.size();
  size_t B_card = B.size();

  std::cerr << "Reading transcript-level expression for A and B" << std::endl;
  std::vector<double> real_tau_A(A_card), real_tau_B(B_card);
  std::string A_expr_fname = vm["A-expr"].as<std::string>();
  std::string B_expr_fname = vm["B-expr"].as<std::string>();
  read_transcript_expression(A_expr_fname, real_tau_A, A_names_to_idxs);
  read_transcript_expression(B_expr_fname, real_tau_B, B_names_to_idxs);

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_A(A_card, 1.0/A_card);
  std::vector<double> unif_tau_B(B_card, 1.0/B_card);

  std::cout << "summarize_matched_version_8\t0" << std::endl;

  stats_tuple recall, precis;
  recall = do_it_all_wrapper<AlignmentType>(vm["A-to-B"].as<std::string>(), A, B, real_tau_B, A_names_to_idxs, B_names_to_idxs);
  precis = do_it_all_wrapper<AlignmentType>(vm["B-to-A"].as<std::string>(), B, A, real_tau_A, B_names_to_idxs, A_names_to_idxs);
  print_stats(precis, recall, "weighted_matched");

  recall = do_it_all_wrapper<AlignmentType>(vm["A-to-B"].as<std::string>(), A, B, unif_tau_B, A_names_to_idxs, B_names_to_idxs);
  precis = do_it_all_wrapper<AlignmentType>(vm["B-to-A"].as<std::string>(), B, A, unif_tau_A, B_names_to_idxs, A_names_to_idxs);
  print_stats(precis, recall, "unweighted_matched");

}
