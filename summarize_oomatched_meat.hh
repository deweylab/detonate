#pragma once
#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>
#define N_POLICY 2
#include "blast.hh"
#include "psl.hh"
#include "util.hh"

double compute_F1(double precis, double recall)
{
  if (precis == 0.0 && recall == 0.0)
    return 0.0;
  else
    return 2*precis*recall/(precis+recall);
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

template<typename AlignmentType>
std::pair<double, double> compute_recall(
    typename AlignmentType::input_stream_type& input_stream,
    size_t                                     A_card,
    size_t                                     B_card,
    const std::map<std::string, size_t>&       A_names_to_idxs,
    const std::map<std::string, size_t>&       B_names_to_idxs,
    const std::vector<double>&                 tau_B,
    bool                                       strand_specific,
    double                                     thresh)
{
  lemon::SmartGraph graph;
  std::vector<lemon::SmartGraph::Node> A_nodes(A_card), B_nodes(B_card);
  for (size_t i = 0; i < A_card; ++i)
    A_nodes[i] = graph.addNode();
  for (size_t i = 0; i < B_card; ++i)
    B_nodes[i] = graph.addNode();

  lemon::SmartGraph::EdgeMap<double> wei_map(graph);

  AlignmentType al;
  while (input_stream >> al) {
    if (is_valid(al, strand_specific) &&
        al.frac_identity_wrt_a() >= thresh && 
        al.frac_identity_wrt_b() >= thresh) {
      size_t a_idx = A_names_to_idxs.find(al.a_name())->second;
      size_t b_idx = B_names_to_idxs.find(al.b_name())->second;
      lemon::SmartGraph::Edge edge = graph.addEdge(A_nodes[a_idx], B_nodes[b_idx]);
      wei_map[edge] = tau_B[b_idx];
    }
  }

  lemon::MaxWeightedMatching<lemon::SmartGraph, lemon::SmartGraph::EdgeMap<double> > wei_mwm(graph, wei_map);
  lemon::MaxMatching<lemon::SmartGraph> unw_mwm(graph);
  wei_mwm.run();
  unw_mwm.run();

  double wei_recall = wei_mwm.matchingWeight();
  double unw_recall = 1.0*unw_mwm.matchingSize()/B_card;

  return std::make_pair(wei_recall, unw_recall);
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
    ("thresh", po::value<double>()->required(), "The threshold for frac_identity (wrt both a and b) below which the alignment is not counted.")
  ;

  try {

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

  std::cerr << "Reading the sequences" << std::endl;
  std::vector<std::string> A, B;
  std::vector<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  read_fasta_names_and_seqs(vm["A-seqs"].as<std::string>(), A, A_names, A_names_to_idxs);
  read_fasta_names_and_seqs(vm["B-seqs"].as<std::string>(), B, B_names, B_names_to_idxs);

  std::cerr << "Computing sequence statistics" << std::endl;
  size_t A_card = A.size(), B_card = B.size();
  std::vector<size_t> A_lengths, B_lengths;
  BOOST_FOREACH(const std::string& a, A) A_lengths.push_back(a.size());
  BOOST_FOREACH(const std::string& b, B) B_lengths.push_back(b.size());

  std::vector<double> tau_A(A_card), tau_B(B_card);
  if (!vm.count("no-expr")) {
    std::cerr << "Reading transcript-level expression for A and B" << std::endl;
    std::string A_expr_fname = vm["A-expr"].as<std::string>();
    std::string B_expr_fname = vm["B-expr"].as<std::string>();
    read_transcript_expression(A_expr_fname, tau_A, A_names_to_idxs);
    read_transcript_expression(B_expr_fname, tau_B, B_names_to_idxs);
  }

  std::cerr << "Reading the alignments" << std::endl;
  typename AlignmentType::input_stream_type A_to_B_is(open_or_throw(vm["A-to-B"].as<std::string>()));
  typename AlignmentType::input_stream_type B_to_A_is(open_or_throw(vm["B-to-A"].as<std::string>()));

  double thresh = vm["thresh"].as<double>();
  std::pair<double, double> recall = compute_recall<AlignmentType>(A_to_B_is, A_card, B_card, A_names_to_idxs, B_names_to_idxs, tau_B, strand_specific, thresh);
  std::pair<double, double> precis = compute_recall<AlignmentType>(B_to_A_is, B_card, A_card, B_names_to_idxs, A_names_to_idxs, tau_A, strand_specific, thresh);
  
  double wei_F1 = compute_F1(precis.first,  recall.first);
  double unw_F1 = compute_F1(precis.second, recall.second);

  std::cout << "summarize_oomatched_version_1\t0\n";

  if (!vm.count("no-expr")) {
    std::cout << "weighted_oomatched_tran_recall\t"    << recall.first << "\n"
              << "weighted_oomatched_tran_precision\t" << precis.first << "\n"
              << "weighted_oomatched_tran_F1\t"        << wei_F1       << "\n";
  }

  std::cout << "unweighted_oomatched_tran_recall\t"    << recall.second << "\n"
            << "unweighted_oomatched_tran_precision\t" << precis.second << "\n"
            << "unweighted_oomatched_tran_F1\t"        << unw_F1        << "\n";
            
  std::cout << std::flush;

  std::cerr << "Done!" << std::endl;
}
