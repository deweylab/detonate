// Copyright (c) 2013
// Nathanael Fillmore (University of Wisconsin-Madison)
// nathanae@cs.wisc.edu
//
// This file is part of REF-EVAL.
//
// REF-EVAL is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// REF-EVAL is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with REF-EVAL.  If not, see <http://www.gnu.org/licenses/>.

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
  result()
  : weighted(-1),
    unweighted(-1)
  {}
};

template<typename V, typename S>
std::string join(V vec, S sep)
{
  std::ostringstream oss;
  typename V::iterator b = vec.begin(), e = vec.end();
  if (b != e) {
    oss << *b;
    ++b;
    for (; b != e; ++b)
      oss << sep << *b;
  }
  return oss.str();
}

template<typename MM>
void output_matching(const std::string& fname,
                     const std::string& pr_string,
                     const fasta& A,
                     const fasta& B,
                     const std::vector<lemon::SmartGraph::Node>& B_nodes,
                     std::map<lemon::SmartGraph::Node, size_t>& A_nodes_to_idxs,
                     const lemon::SmartGraph& graph,
                     const MM& mm,
                     const lemon::SmartGraph::EdgeMap<double> *wei_map)
{
  std::ofstream fo(fname.c_str());
  if (pr_string == "recall")
    fo << "b_name\ta_name\tedge_weight\tcandidate_a_names" << std::endl;
  else // for precision, A and B have been interchanged
    fo << "a_name\tb_name\tedge_weight\tcandidate_b_names" << std::endl;
  for (size_t b_idx = 0; b_idx < B.card; ++b_idx) {
    lemon::SmartGraph::Node b_node = B_nodes[b_idx];
    lemon::SmartGraph::Node a_node = mm.mate(b_node);
    lemon::SmartGraph::Edge edge = mm.matching(b_node);
    // Output the first three columns.
    if (a_node == lemon::INVALID) {
      fo << B.names[b_idx] << "\tNA\tNA\t";
    } else {
      size_t a_idx = A_nodes_to_idxs[a_node];
      fo << B.names[b_idx] << "\t"
         << A.names[a_idx] << "\t"
         << (wei_map ? (*wei_map)[edge] : 1.0/B.card) << "\t";
    }
    // Output the last two columns.
    std::vector<std::string> alt_names;
    for (lemon::SmartGraph::OutArcIt it(graph, b_node); it != lemon::INVALID; ++it) {
      assert(graph.source(it) == b_node);
      lemon::SmartGraph::Node alternate = graph.target(it);
      alt_names.push_back(A.names[A_nodes_to_idxs[alternate]]);
    }
    if (alt_names.size() == 0)
      fo << "NA";
    else
      fo << join(alt_names, ",");
    fo << std::endl;
  }
}

template<typename Al>
result compute_recall(const opts& o,
                      typename Al::input_stream_type& input_stream,
                      const fasta& A,
                      const fasta& B,
                      const std::vector<double>& tau_B,
                      const std::vector<size_t>& num_non_N_in_A,
                      const std::vector<size_t>& num_non_N_in_B,
                      const std::string pr_string)
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
    if (al.is_on_valid_strand(o.strand_specific)) {
        // al.frac_identity_wrt_a() >= o.min_frac_identity && 
        // al.frac_identity_wrt_b() >= o.min_frac_identity &&
        // al.frac_indel_wrt_a() <= o.max_frac_indel &&
        // al.frac_indel_wrt_b() <= o.max_frac_indel)
      size_t a_idx = A.names_to_idxs.find(al.a_name())->second;
      size_t b_idx = B.names_to_idxs.find(al.b_name())->second;
      if (1.0*al.num_identity_wrt_a()/num_non_N_in_A[a_idx] >= o.min_frac_identity && 
          1.0*al.num_identity_wrt_b()/num_non_N_in_B[b_idx] >= o.min_frac_identity &&
          1.0*al.frac_indel_wrt_a()/num_non_N_in_A[a_idx] <= o.max_frac_indel &&
          1.0*al.frac_indel_wrt_b()/num_non_N_in_B[b_idx] <= o.max_frac_indel) {
        lemon::SmartGraph::Edge edge = graph.addEdge(A_nodes[a_idx], B_nodes[b_idx]);
        if (o.weighted)
          wei_map[edge] = tau_B[b_idx];
        if (al_map[edge] == NULL || al.num_identity() >= al_map[edge]->num_identity())
          al_map[edge] = boost::make_shared<Al>(al);
      }
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

  // Output the weighted matching.
  if (o.trace != "" && o.weighted) {
    std::ostringstream fname;
    fname << o.trace << ".weighted_contig_" << pr_string << "_matching";
    output_matching(fname.str(), pr_string, A, B, B_nodes, A_nodes_to_idxs, graph, wei_mm, &wei_map);
  }

  // Output the unweighted matching.
  if (o.trace != "" && (o.unweighted || o.paper)) {
    std::ostringstream fname;
    fname << o.trace << ".unweighted_contig_" << pr_string << "_matching";
    output_matching(fname.str(), pr_string, A, B, B_nodes, A_nodes_to_idxs, graph, unw_mm, NULL);
  }
  return recall;
}

void compute_num_non_N(std::vector<size_t>& num_non_N, const fasta& A)
{
  for (size_t i = 0; i < A.card; ++i) {
    size_t n = 0;
    const std::string& a = A.seqs[i];
    for (std::string::const_iterator it = a.begin(); it != a.end(); ++it)
      if (*it != 'N' && *it != 'n')
        ++n;
    num_non_N[i] = n;
  }
}

template<typename Al>
void main_1(const opts& o,
            const fasta& A,
            const fasta& B,
            const expr& tau_A,
            const expr& tau_B)
{
  typename Al::input_stream_type A_to_B_is(open_or_throw(o.A_to_B));
  typename Al::input_stream_type B_to_A_is(open_or_throw(o.B_to_A));

  std::cerr << "Computing number of non-N bases..." << std::endl;
  std::vector<size_t> num_non_N_in_A(A.card);
  std::vector<size_t> num_non_N_in_B(B.card);
  compute_num_non_N(num_non_N_in_A, A);
  compute_num_non_N(num_non_N_in_B, B);

  std::cerr << "Computing contig precision, recall, and F1 scores..." << std::endl;
  result recall = compute_recall<Al>(o, A_to_B_is, A, B, tau_B, num_non_N_in_A, num_non_N_in_B, "recall");
  result precis = compute_recall<Al>(o, B_to_A_is, B, A, tau_A, num_non_N_in_B, num_non_N_in_A, "precision");

  if (o.weighted) {
    std::cout << "weighted_contig_recall\t" << recall.weighted << std::endl;
    std::cout << "weighted_contig_precision\t" << precis.weighted << std::endl;
    std::cout << "weighted_contig_F1\t" << compute_F1(precis.weighted, recall.weighted) << std::endl;
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
          const expr& tau_B)
{
  if (o.contig || o.paper) {
    if (o.alignment_type == "blast")
      main_1<blast_alignment>(o, A, B, tau_A, tau_B);
      //throw std::runtime_error("tran is not implemented for blast alignments yet.");
    else if (o.alignment_type == "psl")
      main_1<psl_alignment>  (o, A, B, tau_A, tau_B);
  }
}

} // namespace oomatched
} // namespace re
