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

struct Stats
{
  double precis, recall, F1;
  void update_F1()
  {
    if (precis == 0.0 && recall == 0.0)
      F1 = 0.0;
    else
      F1 = 2*precis*recall/(precis+recall);
  }
};

void print_stats(const Stats& stats, const std::string& prefix)
{
  std::cout << prefix << "_precision\t" << stats.precis << std::endl;
  std::cout << prefix << "_recall\t"    << stats.recall << std::endl;
  std::cout << prefix << "_F1\t"        << stats.F1     << std::endl;
}

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

template<typename T>
inline bool is_good_enough(const T& al)
{
  #if (GOOD_POLICY == 1)
  return al.frac_identity_wrt_a() >= 0.95 && al.frac_indel_wrt_a() <= 0.0;
  #elif (GOOD_POLICY == 2)
  return (al.frac_identity_wrt_a() >= 0.95 && al.frac_indel_wrt_a() <= 0.0) ||
         (al.frac_identity_wrt_b() >= 0.95 && al.frac_indel_wrt_b() <= 0.0);
  #elif (GOOD_POLICY == 3)
  return al.frac_identity_wrt_a() >= 0.95 || al.frac_identity_wrt_b() >= 0.95;
  #elif (GOOD_POLICY == 4)
  return true;
  #else
  #error "need to define GOOD_POLICY"
  #endif
}

template<>
inline bool is_good_enough(const blast_alignment& al) { return al.evalue() <= 1e-5; }

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

// For each a in A, figure out which alignment from a -> b (for some b in B) is
// best.
//
// Preconditions:
// - best_from_A needs to be default-initialized of size A.size()
template<typename AlignmentType>
void read_alignments_and_filter_by_best_from_A(std::vector<BestTuple<AlignmentType> >&    best_from_A, 
                                               typename AlignmentType::input_stream_type& input_stream,
                                               const std::map<std::string, size_t>&       A_names_to_idxs)
{
  AlignmentType al;
  while (input_stream >> al) {
    size_t a_idx = A_names_to_idxs.find(al.a_name())->second;
    BestTuple<AlignmentType>& bt = best_from_A[a_idx];
    if (is_good_enough(al)) {
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
                                       typename AlignmentType::input_stream_type& input_stream)
{
  AlignmentType al;
  while (input_stream >> al) {
    if (is_good_enough(al)) {
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
    if (is_good_enough(al)) {
      l.a_idx = A_names_to_idxs.find(al.a_name())->second;
      l.b_idx = B_names_to_idxs.find(al.b_name())->second;
      typename AlignmentType::segments_type segs = al.segments(A[l.a_idx], B[l.b_idx]);
      l.segments.assign(segs.begin(), segs.end());
      if (is_good_enough_for_do_it_all(l.segments))
        L.push_back(l);
    }
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
  //clock_t start;

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
  read_alignments_and_filter_by_best_from_A(best_from_A, A_to_B_is, A_names_to_idxs);
  read_alignments_and_filter_by_best_from_A(best_from_B, B_to_A_is, B_names_to_idxs);

  std::cerr << "Reading alignments without filtering them by A" << std::endl;
  std::vector<BestTuple<AlignmentType> > all_from_A, all_from_B;
  typename AlignmentType::input_stream_type A_to_B_is2(open_or_throw(vm["A-to-B"].as<std::string>()));
  typename AlignmentType::input_stream_type B_to_A_is2(open_or_throw(vm["B-to-A"].as<std::string>()));
  read_alignments_without_filtering(all_from_A, A_to_B_is2);
  read_alignments_without_filtering(all_from_B, B_to_A_is2);

  std::cerr << "Clustering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > clustered_best_to_B(B_card), clustered_best_to_A(A_card);
  cluster_best_alignments_to_B(clustered_best_to_B, best_from_A, B_names_to_idxs);
  cluster_best_alignments_to_B(clustered_best_to_A, best_from_B, A_names_to_idxs);

  std::cerr << "Filtering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > filtered_best_to_B(B_card), filtered_best_to_A(A_card);
  filter_by_best_alignment_to_B(filtered_best_to_B, best_from_A, B_names_to_idxs);
  filter_by_best_alignment_to_B(filtered_best_to_A, best_from_B, A_names_to_idxs);

  std::cerr << "Clustering 'all' alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > jumbled_to_B(B_card), jumbled_to_A(A_card);
  cluster_best_alignments_to_B(jumbled_to_B, all_from_A, B_names_to_idxs);
  cluster_best_alignments_to_B(jumbled_to_A, all_from_B, A_names_to_idxs);

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

  std::cout << "summarize_version_8\t0" << std::endl;

  stats_tuple recall, precis;
  recall = do_it_all_wrapper<AlignmentType>(vm["A-to-B"].as<std::string>(), A, B, real_tau_B, A_names_to_idxs, B_names_to_idxs);
  precis = do_it_all_wrapper<AlignmentType>(vm["B-to-A"].as<std::string>(), B, A, real_tau_A, B_names_to_idxs, A_names_to_idxs);
  print_stats(precis, recall, "weighted_matched");

  recall = do_it_all_wrapper<AlignmentType>(vm["A-to-B"].as<std::string>(), A, B, unif_tau_B, A_names_to_idxs, B_names_to_idxs);
  precis = do_it_all_wrapper<AlignmentType>(vm["B-to-A"].as<std::string>(), B, A, unif_tau_A, B_names_to_idxs, A_names_to_idxs);
  print_stats(precis, recall, "unweighted_matched");

  #if 0
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

  std::cerr << "Computing weighted jumbled stats" << std::endl;
  recall = compute_alignment_stats<smart_pairset>(B_frac_ones, jumbled_to_B, A, B, real_tau_B, A_names_to_idxs);
  precis = compute_alignment_stats<smart_pairset>(A_frac_ones, jumbled_to_A, B, A, real_tau_A, B_names_to_idxs);
  print_stats(precis, recall, "weighted_jumbled");

  std::cerr << "Computing unweighted jumbled stats" << std::endl;
  recall = compute_alignment_stats<smart_pairset>(B_frac_ones, jumbled_to_B, A, B, unif_tau_B, A_names_to_idxs);
  precis = compute_alignment_stats<smart_pairset>(A_frac_ones, jumbled_to_A, B, A, unif_tau_A, B_names_to_idxs);
  print_stats(precis, recall, "unweighted_jumbled");

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

  std::ofstream plot_out(vm["plot-output"].as<std::string>().c_str());
  BOOST_FOREACH(double x, B_frac_ones)
    plot_out << x << std::endl;
  #endif

}
