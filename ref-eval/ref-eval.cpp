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

#include <iostream>
#include <sstream>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include "opts.hh"
#include "fasta.hh"
#include "expr.hh"
#include "util.hh"
#include "re_matched.hh"
#include "re_oomatched.hh"
#include "re_kmer.hh"
#include "re_kc.hh"
#include "re_help.hh"

boost::program_options::options_description describe_options()
{
  namespace po = boost::program_options;
  po::options_description desc;
  desc.add_options()
    ("help,?", "Display this information.\n")
    ("scores", po::value<std::string>())
    ("weighted", po::value<std::string>())
    ("paper", "Flag")
    ("A-seqs", po::value<std::string>())
    ("B-seqs", po::value<std::string>())
    ("A-expr", po::value<std::string>())
    ("B-expr", po::value<std::string>())
    ("A-to-B", po::value<std::string>())
    ("B-to-A", po::value<std::string>())
    ("alignment-type", po::value<std::string>())
    ("strand-specific", "Flag")
    ("readlen", po::value<size_t>())
    ("num-reads", po::value<size_t>())
    ("kmerlen", po::value<size_t>())
    ("min-frac-identity", po::value<double>())
    ("max-frac-indel", po::value<double>())
    ("min-segment-len", po::value<size_t>())
    ("hash-table-type", po::value<std::string>())
    ("hash-table-fudge-factor", po::value<double>())
    ("trace", po::value<std::string>())
  ;
  return desc;
}

void parse_options(opts& o, boost::program_options::variables_map& vm)
{
  namespace po = boost::program_options;

  // Check that either --scores and --weighted, or --paper, are given.
  if (vm.count("paper")) {
    if (vm.count("scores") || vm.count("weighted"))
      throw po::error("--scores and --weighted are incompatible with --paper.");
  }
  else if (vm.count("scores") && vm["scores"].as<std::string>() == "kc") {
    if (vm.count("weighted"))
      throw po::error("--weighted is not needed if only --scores=kc is given.");
  }
  else if (vm.count("scores") || vm.count("weighted")) {
    if (!vm.count("scores"))
      throw po::error("If --weighted is given, --scores is also required.");
    if (!vm.count("weighted"))
      throw po::error("If --scores is given, --weighted is also required.");
  }
  else {
    throw po::error("Either --scores and --weighted or --paper is required.");
  }

  // Parse scores.
  if (vm.count("scores")) {
    std::string scores = vm["scores"].as<std::string>();
    if (scores.size() == 0)
      throw po::error("Invalid empty value for --scores.");
    std::stringstream ss(scores);
    std::string buf;
    while (getline(ss, buf, ',')) {
      if      (buf == "nucl")   { o.nucl   = true; o.alignment_based = true; }
      else if (buf == "contig") { o.contig = true; o.alignment_based = true; }
      else if (buf == "pair")   { o.pair   = true; o.alignment_based = true; }
      else if (buf == "kmer")   { o.kmer   = true; o.alignment_free  = true; }
      else if (buf == "kc")     { o.kc     = true; o.alignment_free  = true; }
      else throw po::error("Invalid value for --scores: " + buf);
    }
  }

  // Parse weighted.
  if (vm.count("weighted")) {
    std::string weighted = vm["weighted"].as<std::string>();
    if      (weighted == "yes")  o.weighted   = true;
    else if (weighted == "no")   o.unweighted = true;
    else if (weighted == "both") o.weighted = o.unweighted = true;
    else throw po::error("Invalid value for --weighted: " + weighted);
  }

  // Parse paper.
  if (vm.count("paper")) {
    o.paper = true;
    o.alignment_based = true;
    o.alignment_free = true;
  }

  // Parse sequences.
  if (!vm.count("A-seqs")) throw po::error("--A-seqs is required.");
  if (!vm.count("B-seqs")) throw po::error("--B-seqs is required.");
  o.A_seqs = vm["A-seqs"].as<std::string>();
  o.B_seqs = vm["B-seqs"].as<std::string>();

  // Parse expression.
  if (o.weighted) {
    if (!vm.count("A-expr"))
      throw po::error("--A-expr is required for weighted variants of scores.");
    if (!vm.count("B-expr"))
      throw po::error("--B-expr is required for weighted variants of scores.");
    o.A_expr = vm["A-expr"].as<std::string>();
    o.B_expr = vm["B-expr"].as<std::string>();
  }
  else if (o.kc || o.paper) {
    if (vm.count("A-expr"))
      throw po::error("--A-expr is not needed except for weighted variants of scores.");
    if (!vm.count("B-expr"))
      throw po::error("--B-expr is required for for the kc score.");
    o.B_expr = vm["B-expr"].as<std::string>();
  }
  else {
    if (vm.count("A-expr"))
      throw po::error("--A-expr is not needed except for weighted variants of scores.");
    if (vm.count("B-expr"))
      throw po::error("--B-expr is not needed except for weighted variants of scores and the kc score.");
  }

  // Parse alignments.
  if (o.alignment_based) {
    if (!vm.count("A-to-B"))
      throw po::error("--A-to-B is required for alignment-based scores.");
    if (!vm.count("B-to-A"))
      throw po::error("--B-to-A is required for alignment-based scores.");
    o.A_to_B = vm["A-to-B"].as<std::string>();
    o.B_to_A = vm["B-to-A"].as<std::string>();
    if (!vm.count("alignment-type"))
      o.alignment_type = "psl";
    else {
      o.alignment_type = vm["alignment-type"].as<std::string>();
      if (o.alignment_type != "blast" && o.alignment_type != "psl")
        throw po::error("Invalid value for --alignment-type: " + o.alignment_type);
      if (o.alignment_type == "blast")
        std::cerr << "Warning: Support for --alignment-type=blast is experimental and has not been thoroughly tested yet." << std::endl;
    }
  } else {
    if (vm.count("A-to-B"))
      throw po::error("--A-to-B is not needed except for alignment-based scores.");
    if (vm.count("B-to-A"))
      throw po::error("--B-to-A is not needed except for alignment-based scores.");
    if (vm.count("alignment-type"))
      throw po::error("--alignment-type is not needed except for alignment-based scores.");
  }

  // Parse strand-specific.
  if (vm.count("strand-specific"))
    o.strand_specific = true;

  // Parse read length.
  if (o.kc || o.paper) {
    if (!vm.count("readlen"))
      throw po::error("--readlen is required for kc scores.");
    o.readlen = vm["readlen"].as<size_t>();
  } else {
    if (vm.count("readlen"))
      throw po::error("--readlen is not needed except for kc scores.");
  }

  // Parse num reads.
  if (o.kc || o.paper) {
    if (!vm.count("num-reads"))
      throw po::error("--num-reads is required for kc scores.");
    o.num_reads = vm["num-reads"].as<size_t>();
  } else {
    if (vm.count("num-reads"))
      throw po::error("--num-reads is not needed except for kc scores.");
  }

  // Parse kmer length.
  if (o.kc || o.kmer || o.paper) {
    if (!vm.count("kmerlen"))
      throw po::error("--kmerlen is required for kc and kmer scores.");
    o.kmerlen = vm["kmerlen"].as<size_t>();
  } else {
    if (vm.count("kmerlen"))
      throw po::error("--kmerlen is not needed except for kc and kmer scores.");
  }

  // Parse min-frac-identity.
  if (vm.count("min-frac-identity"))
    o.min_frac_identity = vm["min-frac-identity"].as<double>();
  else
    o.min_frac_identity = 0.99;

  // Parse max-frac-indel.
  if (vm.count("max-frac-indel"))
    o.max_frac_indel = vm["max-frac-indel"].as<double>();
  else
    o.max_frac_indel = 0.01;

  // Parse min segment length.
  if (vm.count("min-segment-len"))
    o.min_segment_len = vm["min-segment-len"].as<size_t>();
  else
    o.min_segment_len = 100;

  // Parse hash-table-type.
  if (o.kc || o.kmer || o.paper) {
    if (vm.count("hash-table-type")) {
      o.hash_table_type = vm["hash-table-type"].as<std::string>();
      if (o.hash_table_type != "sparse" && o.hash_table_type != "dense")
        throw po::error("Invalid value for --hash-table-type: " + o.hash_table_type);
    } else {
      o.hash_table_type = "sparse";
    }
  } else {
    if (vm.count("hash-table-type"))
      throw po::error("--hash-table-type is not needed except for kmer and kc scores.");
  }

  // Parse hash-table-fudge-factor.
  if (o.kc || o.kmer || o.paper) {
    if (vm.count("hash-table-fudge-factor")) {
      o.hash_table_fudge_factor = vm["hash-table-fudge-factor"].as<double>();
      if (o.hash_table_fudge_factor < 0)
        //throw po::error("Invalid value for --hash-table-fudge-factor: " + boost::lexical_cast<std::string>(o.hash_table_fudge_factor));
        throw po::error("Invalid value for --hash-table-fudge-factor: " + vm["hash-table-fudge-factor"].as<std::string>());
    } else {
      o.hash_table_fudge_factor = 2.0;
    }
  } else {
    if (vm.count("hash-table-fudge-factor"))
      throw po::error("--hash-table-fudge-factor is not needed except for kmer and kc scores.");
  }

  // Parse trace.
  if (vm.count("trace")) {
    o.trace = vm["trace"].as<std::string>();
  }
}

void print_help()
{
  std::cout << get_help_string() << std::endl;
}

int main(int argc, const char **argv)
{
  try {

    std::ios::sync_with_stdio(false);
    namespace po = boost::program_options;

    po::options_description desc = describe_options();
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (argc == 1 || vm.count("help")) {
      print_help();
      exit(0);
    }

    opts o;
    parse_options(o, vm);
    notify(vm);

    std::cerr << "Reading the sequences..." << std::endl;
    fasta A, B;
    read_fasta(A, o.A_seqs);
    read_fasta(B, o.B_seqs);

    expr tau_A, tau_B;
    if (o.weighted) {
      std::cerr << "Reading the expression..." << std::endl;
      tau_A.resize(A.card);
      tau_B.resize(B.card);
      read_rsem_expr(tau_A, o.A_expr, A);
      read_rsem_expr(tau_B, o.B_expr, B);
    } else if (o.kc || o.paper) {
      std::cerr << "Reading the expression..." << std::endl;
      tau_B.resize(B.card);
      read_rsem_expr(tau_B, o.B_expr, B);
    }

    expr unif_A, unif_B;
    if (o.unweighted || o.paper) {
      unif_A.assign(A.card, 1.0/A.card);
      unif_B.assign(B.card, 1.0/B.card);
    }

    re::matched  ::main(o, A, B, tau_A, tau_B, unif_A, unif_B);
    re::oomatched::main(o, A, B, tau_A, tau_B);
    re::kc       ::main(o, A, B,        tau_B);
    re::kmer     ::main(o, A, B, tau_A, tau_B, unif_A, unif_B);

    std::cerr << "Done computing all scores." << std::endl;

  } catch (const boost::program_options::error& x) {

    std::cerr << std::endl;
    std::cerr << argv[0] << ": Error: " << x.what() << std::endl;
    std::cerr << "Check " << argv[0] << " --help for more information." << std::endl;
    return 1;

  } catch (const std::exception& x) {

    std::cerr << std::endl;
    std::cerr << argv[0] << ": Error: " << x.what() << std::endl;
    return 1;

  }

  return 0;
}
