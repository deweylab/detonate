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
#include "re_kc.hh"

void add_options(boost::program_options::options_description& desc)
{
  namespace po = boost::program_options;
  desc.add_options()
    ("help,?", "Display this information.\n")
    ("scores", po::value<std::string>(),
        "The score groups to compute, separated by commas "
        "(e.g., --scores=nucl,contig,kc):\n"
        "\n"
        "  Alignment-based score groups:\n"
        "  - nucl:   \tnucleotide precision, recall, and F1.\n"
        "  - contig: \tcontig precision, recall, and F1.\n"
        "  - pair:   \tpair precision, recall, and F1.\n"
        "\n"
        "  Alignment-free score groups:\n"
        "  - kmer:   \tkmer Kullback-Leibler divergence, "
                      "Jensen-Shannon divergence, and "
                      "Hellinger distance.\n"
        "  - kc:     \tkmer recall, number of nucleotides, and "
                      "kmer compression score.\n"
        "\n"
        "Required unless --paper is given.\n")
    ("weighted", po::value<std::string>(),
        "A string indicating whether to compute weighted or "
        "unweighted variants of scores, or both "
        "(e.g., --weighted=yes):\n"
        "\n"
        "  - yes:  \tcompute weighted variants of scores.\n"
        "  - no:   \tcompute unweighted variants of scores.\n"
        "  - both: \tcompute both weighted and unweighted "
                    "variants of scores.\n"
        "\n"
        "The weighted variants require --A-expr and --B-expr "
        "to be specified.\n"
        "\n"
        "Required unless --paper is given.\n")
    ("paper", 
        "If this flag is present, the reference-based scores "
        "described in the main text of the paper corresponding "
        "to this software will be computed. These scores are as "
        "follows:\n"
        "\n"
        "  Alignment-based scores:\n"
        "  - unweighted nucleotide F1\n"
        "  - unweighted contig F1\n"
        "\n"
        "  Alignment-free score groups:\n"
        "  - weighted kmer compression score\n"
        "\n"
        "For obvious reasons, the --scores and --weighted options "
        "are incompatible with this flag.\n")
    ("A-seqs", po::value<std::string>(),
        "The assembly sequences, in FASTA format. Required.\n")
    ("B-seqs", po::value<std::string>(),
        "The oracleset sequences, in FASTA format. Required.\n")
    ("A-expr", po::value<std::string>(),
        "The assembly expression, for use in weighted scores, as "
        "produced by RSEM in a file called *.isoforms.results. "
        "Required for weighted variants of scores.\n")
    ("B-expr", po::value<std::string>(),
        "The oracleset expression, for use in weighted scores, as "
        "produced by RSEM in a file called *.isoforms.results. "
        "Required for weighted variants of scores.\n")
    ("A-to-B", po::value<std::string>(),
        "The alignments of the assembly to the oracleset. The file "
        "format is specified by --alignment-type. Required for "
        "alignment-based scores.\n")
    ("B-to-A", po::value<std::string>(),
        "The alignments of the oracleset to the assembly. The file "
        "format is specified by --alignment-type. Required for "
        "alignment-based scores.\n")
    ("alignment-type", po::value<std::string>(),
        "The type of alignments used, either \"blast\" or \"psl\". "
        "Required for alignment-based scores.\n")
    ("strand-specific",
        "If this flag is present, ignore alignments or kmer matches "
        "that are to the reverse strand.\n")
    ("readlen", po::value<size_t>(),
        //"The read length. Required for kmer and kc scores.\n"
        "The read length of the reads used to build the assembly. Required.\n") // XXX really shouldn't be
    ("num_reads", po::value<size_t>(),
        "The number of reads used to build the assembly. Required "
        "for kc scores.\n")
    ("contig-min-frac-identity", po::value<double>(),
        "This option only applies to contig scores. Alignments with "
        "fraction identity less than this threshold are ignored. "
        "The fraction identity of an alignment is min(x/y, x/z), "
        "where \n"
        "\n"
        "  - x \tis the number of bases that are identical in the "
                 "assembly sequence and the oracleset sequence, "
                 "according to the alignment,\n"
        "  - y \tis the number of bases in the assembly sequence, and\n"
        "  - z \tis the number of bases in the oracleset sequence.\n"
        "\n"
        "Default: 0.99.\n")
    ("contig-max-frac-indel", po::value<double>(),
        "This option only applies to contig scores. Alignments with "
        "fraction indel greater than this threshold are ignored. "
        "For psl alignments, the fraction indel of an alignment "
        "is max(w/y, x/z), where \n"
        "\n"
        "  - w \tis the number of bases that are inserted in the "
                 "assembly sequence, according to the alignment "
                 "(\"Q gap bases\"),\n"
        "  - x \tis the number of bases that are inserted in the "
                 "oracleset sequence, according to the alignment "
                 "(\"T gap bases\"),\n"
        "  - y \tis the number of bases in the assembly sequence, and\n"
        "  - z \tis the number of bases in the oracleset sequence.\n"
        "\n"
        "For blast alignments, the fraction indel of an alignment "
        "is max(x/y, x/z), where \n"
        "\n"
        "  - x \tis the number of gaps bases that are inserted in the "
                 "oracleset sequence, according to the alignment "
                 "(\"gaps\"),\n"
        "  - y \tis the number of bases in the assembly sequence, and\n"
        "  - z \tis the number of bases in the oracleset sequence.\n"
        "\n"
        "Default: 0.01.")
  ;
}

void parse_options(opts& o, boost::program_options::variables_map& vm)
{
  namespace po = boost::program_options;

  // Check that either --scores and --weighted, or --paper, are given.
  if (vm.count("paper")) {
    if (vm.count("scores") || vm.count("weighted"))
      throw po::error("--scores and --weighted are incompatible with --paper.");
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

  // Parse weighted.
  std::string weighted = vm["weighted"].as<std::string>();
  if      (weighted == "yes")  o.weighted   = true;
  else if (weighted == "no")   o.unweighted = true;
  else if (weighted == "both") o.weighted = o.unweighted = true;
  else throw po::error("Invalid value for --weighted: " + weighted);

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
  } else {
    if (vm.count("A-expr"))
      throw po::error("--A-expr is not needed except for weighted variants of scores.");
    if (vm.count("B-expr"))
      throw po::error("--B-expr is not needed except for weighted variants of scores.");
  }

  // Parse alignments.
  if (o.alignment_based) {
    if (!vm.count("A-to-B"))
      throw po::error("--A-to-B is required for alignment-based scores.");
    if (!vm.count("B-to-A"))
      throw po::error("--B-to-A is required for alignment-based scores.");
    if (!vm.count("alignment-type"))
      throw po::error("--alignment-type is required for alignment-based scores.");
    o.A_to_B = vm["A-to-B"].as<std::string>();
    o.B_to_A = vm["B-to-A"].as<std::string>();
    o.alignment_type = vm["alignment-type"].as<std::string>();
    if (o.alignment_type != "blast" && o.alignment_type != "psl")
      throw po::error("Invalid value for --alignment-type: " + o.alignment_type);
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
  // if (o.kmer || o.kc) {
  //   if (!vm.count("readlen"))
  //     throw po::error("--readlen is required for kmer and kc scores.");
  //   o.readlen = vm["readlen"].as<size_t>();
  // } else {
  //   if (vm.count("readlen"))
  //     throw po::error("--readlen is not needed except for kmer and kc scores.");
  // }
  if (!vm.count("readlen"))
    throw po::error("--readlen is required.");
  o.readlen = vm["readlen"].as<size_t>();

  // Parse num reads.
  if (o.kc) {
    if (!vm.count("num_reads"))
      throw po::error("--num_reads is required for kc scores.");
    o.num_reads = vm["num_reads"].as<size_t>();
  } else {
    if (vm.count("num_reads"))
      throw po::error("--num_reads is not needed except for kc scores.");
  }

  // Parse contig-min-frac-identity.
  if (vm.count("contig-min-frac-identity"))
    o.contig_min_frac_identity = vm["contig-min-frac-identity"].as<double>();
  else
    o.contig_min_frac_identity = 0.99;

  // Parse contig-max-frac-indel.
  if (vm.count("contig-max-frac-indel"))
    o.contig_max_frac_indel = vm["contig-max-frac-indel"].as<double>();
  else
    o.contig_max_frac_indel = 0.01;
}

int main(int argc, const char **argv)
{
  try {

    std::ios::sync_with_stdio(false);
    namespace po = boost::program_options;

    po::options_description desc("Arguments");
    add_options(desc);
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cerr << desc << std::endl;
      exit(0);
    }

    opts o;
    parse_options(o, vm);
    notify(vm);

    std::cerr << "Reading the sequences..." << std::flush;
    fasta A, B;
    read_fasta(A, o.A_seqs);
    read_fasta(B, o.B_seqs);
    std::cerr << "done." << std::endl;

    expr tau_A, tau_B;
    if (o.weighted) {
      std::cerr << "Reading the expression..." << std::flush;
      tau_A.resize(A.card);
      tau_B.resize(B.card);
      read_rsem_expr(tau_A, o.A_expr, A);
      read_rsem_expr(tau_B, o.B_expr, B);
      std::cerr << "done." << std::endl;
    }

    expr unif_A, unif_B;
    if (o.unweighted) {
      unif_A.assign(A.card, 1.0/A.card);
      unif_B.assign(B.card, 1.0/B.card);
    }

    re::matched  ::main(o, A, B, tau_A, tau_B, unif_A, unif_B);
    re::oomatched::main(o, A, B, tau_A, tau_B, unif_A, unif_B);
    re::kc       ::main(o, A, B, tau_A, tau_B, unif_A, unif_B);

    std::cerr << "Done!" << std::endl;

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
