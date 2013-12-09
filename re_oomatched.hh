#pragma once
#include <iostream>
#include <vector>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#include "blast.hh"
#include "psl.hh"
#include "util.hh"

namespace re {
namespace oomatched {

struct result
{
  double weighted;
  double unweighted;
};

template<typename V, typename S>
std::string join(V vec, S sep)
{
  std::ostringstream oss;
  typename V::iterator b = vec.begin(), e = vec.end();
  if (b != e)
    oss << *b;
  ++b;
  for (; b != e; ++b)
    oss << sep << *b;
  return oss.str();
}

template<typename Al>
result compute_recall(
    const opts& o,
    typename Al::input_stream_type& input_stream,
    const fasta& A,
    const fasta& B,
    const std::vector<double>& tau_B)
{
  // Create the graph, and add a node for each contig and oracleset element.
  // Also create a reverse mapping from nodes to indices.
  lemon::SmartGraph graph;
  std::vector<lemon::SmartGraph::Node> A_nodes(A.card), B_nodes(B.card);
  std::map<lemon::SmartGraph::Node, size_t> A_nodes_to_idxs, B_nodes_to_idxs;
  for (size_t i = 0; i < A.card; ++i) {
    A_nodes[i] = graph.addNode();
    A_nodes_to_idxs[A_nodes[i]] = i;
  }
  for (size_t i = 0; i < B.card; ++i) {
    B_nodes[i] = graph.addNode();
    B_nodes_to_idxs[B_nodes[i]] = i;
  }

  // Add edges to the graph (and detetermine their weights) based on the given
  // alignments.
  lemon::SmartGraph::EdgeMap<double> wei_map(graph);
  lemon::SmartGraph::EdgeMap<boost::shared_ptr<Al> > al_map(graph);
  Al al;
  while (input_stream >> al) {
    if (al.is_on_valid_strand(o.strand_specific) &&
        al.frac_identity_wrt_a() >= o.contig_min_frac_identity && 
        al.frac_identity_wrt_b() >= o.contig_min_frac_identity &&
        al.frac_indel_wrt_a() <= o.contig_max_frac_indel &&
        al.frac_indel_wrt_b() <= o.contig_max_frac_indel) {
      size_t a_idx = A.names_to_idxs.find(al.a_name())->second;
      size_t b_idx = B.names_to_idxs.find(al.b_name())->second;
      lemon::SmartGraph::Edge edge = graph.addEdge(A_nodes[a_idx], B_nodes[b_idx]);
      if (o.weighted)
        wei_map[edge] = tau_B[b_idx];
      if (al_map[edge] == NULL || al.num_identity() >= al_map[edge]->num_identity())
        al_map[edge] = boost::make_shared<Al>(al);
    }
  }

  // Run the matching procedure.
  lemon::MaxWeightedMatching<lemon::SmartGraph, lemon::SmartGraph::EdgeMap<double> > wei_mm(graph, wei_map);
  lemon::MaxMatching<lemon::SmartGraph> unw_mm(graph);
  if (o.weighted)              wei_mm.run();
  if (o.unweighted || o.paper) unw_mm.run();

  // Compute the recall.
  result recall;
  if (o.weighted)              recall.weighted   = wei_mm.matchingWeight();
  if (o.unweighted || o.paper) recall.unweighted = 1.0*unw_mm.matchingSize()/B.card;

  // // Output the weighted matching.
  // ofstream wei_fo(wei_mm_fname.c_str());
  // for (size_t b_idx = 0; b_idx < B.card; ++b_idx) {
  //   lemon::SmartGraph::Node a_node = wei_mm.mate(B_nodes[b_idx]);
  //   lemon::SmartGraph::Edge edge = wei_mm.matching(B_nodes[b_idx]);
  //   if (a_node == lemon::INVALID)
  //     wei_fo << B.names[b_idx] << "\tNA\tNA\tNA\tNA\tNA" << std::endl;
  //   else {
  //     size_t a_idx = A_nodes_to_idxs[a_node];
  //     wei_fo << B.names[b_idx] << "\t" << A.names[a_idx] << "\t" << wei_map[edge] << "\t"
  //            << join(al_map[edge]->block_sizes(), ",") << "\t"
  //            << join(al_map[edge]->q_starts(), ",") << "\t"
  //            << join(al_map[edge]->t_starts(), ",")
  //            << std::endl;
  //   }
  // }
  // 
  // // Output the unweighted matching.
  // ofstream unw_fo(unw_mm_fname.c_str());
  // for (size_t b_idx = 0; b_idx < B.card; ++b_idx) {
  //   lemon::SmartGraph::Node a_node = unw_mm.mate(B_nodes[b_idx]);
  //   lemon::SmartGraph::Edge edge = unw_mm.matching(B_nodes[b_idx]);
  //   if (a_node == lemon::INVALID)
  //     unw_fo << B.names[b_idx] << "\tNA\tNA\tNA\tNA\tNA" << std::endl;
  //   else {
  //     size_t a_idx = A_nodes_to_idxs[a_node];
  //     unw_fo << B.names[b_idx] << "\t" << A.names[a_idx] << "\t" << 1.0/B.card << "\t"
  //            << join(al_map[edge]->block_sizes(), ",") << "\t"
  //            << join(al_map[edge]->q_starts(), ",") << "\t"
  //            << join(al_map[edge]->t_starts(), ",")
  //            << std::endl;
  //   }
  // }

  return recall;
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
  typename Al::input_stream_type A_to_B_is(open_or_throw(o.A_to_B));
  typename Al::input_stream_type B_to_A_is(open_or_throw(o.B_to_A));

  result recall = compute_recall<Al>(o, A_to_B_is, A, B, tau_B);
  result precis = compute_recall<Al>(o, B_to_A_is, B, A, tau_A);

  if (o.weighted) {
    std::cout << "weighted_contig_recall\t"   << recall.weighted   << std::endl;
    std::cout << "weighted_contig_precision\t"   << precis.weighted   << std::endl;
    std::cout << "weighted_contig_F1\t"   << compute_F1(precis.weighted,   recall.weighted)   << std::endl;
  }

  if (o.unweighted || o.paper) {
    std::cout << "unweighted_contig_recall\t" << recall.unweighted << std::endl;
    std::cout << "unweighted_contig_precision\t" << precis.unweighted << std::endl;
    std::cout << "unweighted_contig_F1\t" << compute_F1(precis.unweighted, recall.unweighted) << std::endl;
  }
}

void main(const opts& o,
          const fasta& A,
          const fasta& B,
          const expr& tau_A,
          const expr& tau_B,
          const expr& unif_A,
          const expr& unif_B)
{
  if (o.contig || o.paper) {
    if (o.alignment_type == "blast")
      main_1<blast_alignment>(o, A, B, tau_A, tau_B, unif_A, unif_B);
      //throw std::runtime_error("tran is not implemented for blast alignments yet.");
    else if (o.alignment_type == "psl")
      main_1<psl_alignment>  (o, A, B, tau_A, tau_B, unif_A, unif_B);
  }
}

} // namespace oomatched
} // namespace re
