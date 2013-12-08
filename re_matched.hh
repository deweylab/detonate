#pragma once
#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include "blast.hh"
#include "psl.hh"
#include "pairset.hh"
#include "kpairset.hh"
#include "kmerset.hh"
#include "mask.hh"
#include "blast.hh"
#include "psl.hh"
#include "util.hh"

namespace re {
namespace matched {

struct tagged_alignment
{
  size_t a_idx, b_idx;
  std::vector<alignment_segment> segments;
  double contribution;
};

struct compare_tagged_alignments
{
  bool operator()(const tagged_alignment *l1, const tagged_alignment *l2) const
  {
    return l1->contribution < l2->contribution;
  }
};

struct pair_helper
{
  const std::vector<size_t>& B_lengths;
  const std::vector<double>& tau_B;
  double                     numer;

  pair_helper(const std::vector<size_t>& B_lengths,
              const std::vector<double>& tau_B,
              size_t /*k*/)
  : B_lengths(B_lengths), tau_B(tau_B), numer(0.0)
  {}

  double compute_contribution(const tagged_alignment& l) const
  {
    smart_pairset b_pairset(B_lengths[l.b_idx]);
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_pairset.add_square_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_pairset.size();
  }

  void add_contribution_to_recall(const tagged_alignment& l) { numer += l.contribution; }

  double get_recall()
  {
    double denom = 0.0;
    for (size_t b_idx = 0; b_idx < B_lengths.size(); ++b_idx)
      denom += tau_B[b_idx] * B_lengths[b_idx] * (B_lengths[b_idx] + 1) / 2;
    return numer/denom;
  }
};

struct kpair_helper
{
  const std::vector<size_t>& B_lengths;
  const std::vector<double>& tau_B;
  size_t                     k;
  double                     numer;

  kpair_helper(const std::vector<size_t>& B_lengths,
               const std::vector<double>& tau_B,
               size_t k)
  : B_lengths(B_lengths), tau_B(tau_B), k(k), numer(0.0)
  {}

  double compute_contribution(const tagged_alignment& l) const
  {
    kpairset b_kpairset(k, B_lengths[l.b_idx]);
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_kpairset.add_kpairs_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_kpairset.size();
  }

  void add_contribution_to_recall(const tagged_alignment& l) { numer += l.contribution; }

  double get_recall()
  {
    // How many k-pairs are there in a seq of length n?
    // There are n-k+2 such k-pairs.
    // 
    // E.g., k=3, n=5
    // 01234  -> (0,2), (1,3), (2,4)
    // Indeed, 5-3+1 = 3
    double denom = 0.0;
    for (size_t b_idx = 0; b_idx < B_lengths.size(); ++b_idx) {
      if (B_lengths[b_idx] >= k)
        denom += tau_B[b_idx] * (B_lengths[b_idx] - k + 1);
    }
    return numer/denom;
  }
};

struct kmer_helper
{
  const std::vector<size_t>& B_lengths;
  const std::vector<double>& tau_B;
  size_t                     k;
  double                     numer;

  kmer_helper(const std::vector<size_t>& B_lengths,
              const std::vector<double>& tau_B,
              size_t k)
  : B_lengths(B_lengths), tau_B(tau_B), k(k), numer(0.0)
  {}

  double compute_contribution(const tagged_alignment& l) const
  {
    kmerset b_kmerset(k, B_lengths[l.b_idx]);
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_kmerset.add_kmers_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_kmerset.size();
  }

  void add_contribution_to_recall(const tagged_alignment& l) { numer += l.contribution; }

  double get_recall()
  {
    // How many kmers are there in a seq of length n?
    // There are n-k+2 such kmers.
    // 
    // E.g., k=3, n=5
    // 01234  -> (0,1,2), (1,2,3), (2,3,4)
    // Indeed, 5-3+1 = 3
    double denom = 0.0;
    for (size_t b_idx = 0; b_idx < B_lengths.size(); ++b_idx) {
      if (B_lengths[b_idx] >= k)
        denom += tau_B[b_idx] * (B_lengths[b_idx] - k + 1);
    }
    return numer/denom;
  }
};

struct nucl_helper
{
  const std::vector<size_t>& B_lengths;
  const std::vector<double>& tau_B;
  double                     numer;

  nucl_helper(const std::vector<size_t>& B_lengths,
              const std::vector<double>& tau_B,
              size_t /*k*/)
  : B_lengths(B_lengths), tau_B(tau_B), numer(0.0)
  {}

  double compute_contribution(const tagged_alignment& l) const
  {
    mask b_mask(B_lengths[l.b_idx]);
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_mask.num_ones();
  }

  void add_contribution_to_recall(const tagged_alignment& l) { numer += l.contribution; }

  double get_recall()
  {
    double denom = 0.0;
    for (size_t b_idx = 0; b_idx < B_lengths.size(); ++b_idx)
      denom += tau_B[b_idx] * B_lengths[b_idx];
    return numer/denom;
  }
};

struct tran_helper
{
  const std::vector<size_t>& B_lengths;
  const std::vector<double>& tau_B;
  std::vector<mask>          B_mask;
  std::vector<double>        B_frac_ones;

  tran_helper(const std::vector<size_t>& B_lengths,
              const std::vector<double>& tau_B,
              size_t /*k*/)
  : B_lengths(B_lengths), tau_B(tau_B), B_frac_ones(B_lengths.size())
  {
    BOOST_FOREACH(size_t b_length, B_lengths)
      B_mask.push_back(mask(b_length));
  }

  double compute_contribution(const tagged_alignment& l) const
  {
    mask b_mask(B_lengths[l.b_idx]);
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
    return tau_B[l.b_idx] * b_mask.num_ones();
  }

  void add_contribution_to_recall(const tagged_alignment& l)
  {
    mask& m = B_mask[l.b_idx];
    BOOST_FOREACH(const alignment_segment& seg, l.segments)
      m.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
  }

  double get_recall()
  {
    double recall = 0.0;
    for (size_t b_idx = 0; b_idx < B_lengths.size(); ++b_idx) {
      B_frac_ones[b_idx] = (1.0 * B_mask[b_idx].num_ones()) / B_lengths[b_idx];
      if (B_frac_ones[b_idx] >= 0.95)
        recall += tau_B[b_idx];
    }
    return recall;
  }
};

inline bool is_good_enough(const std::vector<alignment_segment>& segs, size_t min_seg_len)
{
  size_t len = 0;
  BOOST_FOREACH(const alignment_segment& seg, segs) {
    if (seg.a_end >= seg.a_start)
      len += seg.a_end - seg.a_start + 1;
    else
      len += seg.a_start - seg.a_end + 1;
  }
  return len >= min_seg_len;
}

// Reads alignments, perfoms initial filtering, and converts the alignments to
// segments. However, does *not* compute the initial contributions.
template<typename Al>
void read_alignments(std::vector<tagged_alignment>& alignments,
                     const std::string& filename,
                     const std::vector<std::string>& A,
                     const std::vector<std::string>& B,
                     const std::map<std::string, size_t>& A_names_to_idxs,
                     const std::map<std::string, size_t>& B_names_to_idxs,
                     bool strand_specific,
                     size_t readlen)
{
  try {
    typename Al::input_stream_type input_stream(open_or_throw(filename));
    Al al;
    tagged_alignment l;
    std::map<std::string, size_t>::const_iterator it;
    while (input_stream >> al) {
      if (al.is_on_valid_strand(strand_specific)) {
        // Extract a_name and look up its idx.
        it = A_names_to_idxs.find(al.a_name());
        if (it == A_names_to_idxs.end())
          throw std::runtime_error("Sequence name " + al.a_name() + " is not found in the corresponding fasta file.");
        l.a_idx = it->second;
        // Extract b_name and look up its idx.
        it = B_names_to_idxs.find(al.b_name());
        if (it == A_names_to_idxs.end())
          throw std::runtime_error("Sequence name " + al.b_name() + " is not found in the corresponding fasta file.");
        l.b_idx = it->second;
        // Extract the alignment segments.
        typename Al::segments_type segs = al.segments(A[l.a_idx], B[l.b_idx]);
        l.segments.assign(segs.begin(), segs.end());
        if (is_good_enough(l.segments, readlen))
          alignments.push_back(l);
      }
    }
  } catch (const std::runtime_error& x) {
    throw std::runtime_error("Can't parse " + filename + ": " + x.what());
  }
}

// Preconditions:
// - best_from_A should be of size 0
template<typename HelperType>
void process_alignments(HelperType&                   helper,
                        std::vector<tagged_alignment> alignments, /* passed by value intentionally */
                        size_t                        A_card,
                        size_t                        B_card,
                        size_t                        readlen)
{
  // Compute contributions.
  #pragma omp parallel for
  for (int i = 0; i < static_cast<int>(alignments.size()); ++i) {
    tagged_alignment& l = alignments[i];
    l.contribution = helper.compute_contribution(l);
  }

  // Init vector of pointers to popped alignments.
  std::vector<std::vector<tagged_alignment *> > popped_by_A(A_card);
  std::vector<std::vector<tagged_alignment *> > popped_by_B(B_card);

  // Make priority queue initially filled with all alignments.
  compare_tagged_alignments comparator;
  std::priority_queue<
    tagged_alignment *,
    std::vector<tagged_alignment *>,
    compare_tagged_alignments> Q(comparator);
  BOOST_FOREACH(tagged_alignment& l1, alignments)
    Q.push(&l1);

  while (!Q.empty()) {

    // Pop l1 from Q.
    tagged_alignment *l1 = Q.top();
    Q.pop();

    // Subtract all previous alignments from l1.
    bool l1_has_changed = false;
    typename std::vector<tagged_alignment *>::const_iterator
        l2_a = popped_by_A[l1->a_idx].begin(), end_a = popped_by_A[l1->a_idx].end(),
        l2_b = popped_by_B[l1->b_idx].begin(), end_b = popped_by_B[l1->b_idx].end();
    while (l2_a != end_a || l2_b != end_b) {

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

      if (is_good_enough(l1->segments, readlen)) {
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

template<typename Helper>
double compute_recall(const opts& o,
                      const fasta& A,
                      const fasta& B,
                      const expr& tau_B,
                      const std::vector<tagged_alignment>& A_to_B)
{
  Helper h(B.lengths, tau_B, o.readlen);
  process_alignments(h, A_to_B, A.card, B.card, o.readlen);
  return h.get_recall();
}

template<typename Helper>
void compute(const opts& o,
             const fasta& A,
             const fasta& B,
             const expr& tau_A,
             const expr& tau_B,
             std::vector<tagged_alignment> A_to_B,
             std::vector<tagged_alignment> B_to_A,
             const std::string& prefix)
{
  double precis = compute_recall<Helper>(o, B, A, tau_A, B_to_A);
  std::cout << prefix << "precision\t" << precis << std::endl;

  double recall = compute_recall<Helper>(o, A, B, tau_B, A_to_B);
  std::cout << prefix << "recall\t" << recall << std::endl;

  double F1 = compute_F1(precis, recall);
  std::cout << prefix << "F1\t" << F1 << std::endl;
}

template<typename Helper>
void main_2(const opts& o,
            const fasta& A,
            const fasta& B,
            const expr& tau_A,
            const expr& tau_B,
            const expr& unif_A,
            const expr& unif_B,
            std::vector<tagged_alignment> A_to_B,
            std::vector<tagged_alignment> B_to_A,
            const std::string& prefix)
{
  if (o.weighted)   compute<Helper>(o, A, B, tau_A,  tau_B,  A_to_B, B_to_A, "weighted_" + prefix);
  if (o.unweighted) compute<Helper>(o, A, B, unif_A, unif_B, A_to_B, B_to_A, "unweighted_" + prefix);
}

template<typename Al>
void main_1(const opts& o,
            const fasta& A,
            const fasta& B,
            const expr& tau_A,
            const expr& tau_B,
            const expr& unif_A,
            const expr& unif_B)
{
  std::cerr << "Reading the alignments and extracting intervals..." << std::flush;
  std::vector<tagged_alignment> A_to_B, B_to_A;
  read_alignments<Al>(A_to_B, o.A_to_B, A.seqs, B.seqs, A.names_to_idxs, B.names_to_idxs, o.strand_specific, o.readlen);
  read_alignments<Al>(B_to_A, o.B_to_A, B.seqs, A.seqs, B.names_to_idxs, A.names_to_idxs, o.strand_specific, o.readlen);
  std::cerr << "done." << std::endl;

  if (o.nucl) main_2<nucl_helper>(o, A, B, tau_A, tau_B, unif_A, unif_B, A_to_B, B_to_A, "nucl_");
  if (o.pair) main_2<pair_helper>(o, A, B, tau_A, tau_B, unif_A, unif_B, A_to_B, B_to_A, "pair_");
  //if (o.kpair) main_2<kpair_helper>(o, A, B, tau_A, tau_B, unif_A, unif_B, A_to_B, B_to_A, "kpair_");
  //if (o.kmer) main_2<kmer_helper>(o, A, B, tau_A, tau_B, unif_A, unif_B, A_to_B, B_to_A, "kmer_");
  //if (o.tran) main_2<kmer_helper>(o, A, B, tau_A, tau_B, unif_A, unif_B, A_to_B, B_to_A, "tran_");
}

void main(const opts& o,
          const fasta& A,
          const fasta& B,
          const expr& tau_A,
          const expr& tau_B,
          const expr& unif_A,
          const expr& unif_B)
{
  if (o.nucl || o.pair) {
    if (o.alignment_type == "blast")
      main_1<blast_alignment>(o, A, B, tau_A, tau_B, unif_A, unif_B);
    else if (o.alignment_type == "psl")
      main_1<psl_alignment>  (o, A, B, tau_A, tau_B, unif_A, unif_B);
  }
}

} // namespace matched
} // namespace re
