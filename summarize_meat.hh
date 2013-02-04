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
#include <boost/random/random_device.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include "deweylab/bio/formats/fasta.hh"
#include "blast.hh"
#include "psl.hh"
#include "pairset.hh"
#include "mask.hh"

// need to define "clock_t start" before tic
#define tic start = clock();
#define toc std::cerr << "Done in " << 1.0*(clock() - start)/CLOCKS_PER_SEC << " seconds." << std::endl;

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

void open_or_throw(std::ifstream& ifs, const std::string& filename)
{
  ifs.clear();
  ifs.open(filename.c_str());
  if (!ifs)
    throw std::runtime_error("Could not open file '" + filename + "'.");
}

void read_fasta_names_and_seqs(const std::string& filename,
                               std::vector<std::string>& seqs,
                               std::vector<std::string>& names,
                               std::map<std::string, size_t>& names_to_idxs)
{
  std::ifstream ifs;
  open_or_throw(ifs, filename);
  deweylab::bio::formats::fasta::InputStream is(ifs);
  deweylab::bio::formats::fasta::Record rec;
  size_t idx = 0;
  while (is >> rec) {
    assert(names_to_idxs.count(rec.id) == 0);
    names_to_idxs[rec.id] = idx;
    names.push_back(rec.id);
    seqs.push_back(rec.sequence);
    assert(names[names_to_idxs[rec.id]] == rec.id);
    ++idx;
  }
}

void read_transcript_expression(const std::string& filename,
                                std::vector<double>& expr,
                                const std::map<std::string, size_t>& names_to_idxs)
{
  std::ifstream ifs;
  open_or_throw(ifs, filename);
  std::string line;
  std::vector<std::string> col(8);
  { // header
    getline(ifs, line);
    std::stringstream ss(line);
    for (size_t i = 0; i < 8; ++i)
      getline(ss, col[i], '\t');
    assert(col[0] == "transcript_id");
    assert(col[5] == "TPM");
  }
  while (getline(ifs, line)) {
    std::stringstream ss(line);
    for (size_t i = 0; i < 8; ++i)
      getline(ss, col[i], '\t');
    assert(names_to_idxs.count(col[0]) != 0);
    size_t idx = names_to_idxs.find(col[0])->second;
    expr[idx] = boost::lexical_cast<double>(col[5]) / 1000000.0;
  }
}

void compute_nucl_expression(const std::vector<std::string>& A,
                             const std::vector<double>& tau_A,
                             std::vector<double>& nu_A)
{
  size_t n = tau_A.size();
  assert(tau_A.size() == n);
  assert(nu_A.size() == n);

  for (size_t i = 0; i < n; ++i)
    nu_A[i] = tau_A[i] * A[i].size();

  double denom = 0;
  for (size_t i = 0; i < n; ++i)
    denom += nu_A[i];

  for (size_t i = 0; i < n; ++i)
    nu_A[i] /= denom;
}

template<typename AlignmentType>
struct BestTuple
{
  double frac_identity, frac_indel;
  AlignmentType al;
  BestTuple() : frac_identity(-1.0) {}
};

inline bool is_good_enough(const psl_alignment& al) { return al.frac_identity() > 0.95; }
inline bool is_good_enough(const blast_alignment& al) { return al.evalue() < 1e-5; }

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
      double frac_identity = al.frac_identity();
      double frac_indel    = al.frac_indel();
      // relies on default init of bt.frac_identity to -1.0
      if (frac_identity > bt.frac_identity
          || (frac_identity == bt.frac_identity && frac_indel < bt.frac_indel)) {
        bt.frac_identity = frac_identity;
        bt.frac_indel    = frac_indel;
        bt.al            = al;
      }
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
  typename std::vector<BestTuple<AlignmentType> >::const_iterator it;
  for (it = best_from_A.begin(); it != best_from_A.end(); ++it) {
    if (it->frac_identity >= 0) { // the purpose here is to exclude default-constructed BestTuples, i.e., non-alignments
      size_t b_idx = B_names_to_idxs.find(it->al.b_name())->second;
      best_to_B[b_idx].push_back(&(it->al));
    }
  }
}

template<typename AlignmentType>
void filter_by_best_alignment_to_B(std::vector<std::vector<const AlignmentType *> >& best_to_B,
                                   const std::vector<BestTuple<AlignmentType> >&     best_from_A,
                                   const std::map<std::string, size_t>&              B_names_to_idxs)
{
  typename std::vector<BestTuple<AlignmentType> >::const_iterator it;
  for (it = best_from_A.begin(); it != best_from_A.end(); ++it) {
    if (it->frac_identity >= 0) { // the purpose here is to exclude default-constructed BestTuples, i.e., non-alignments
      size_t b_idx = B_names_to_idxs.find(it->al.b_name())->second;
      if (best_to_B[b_idx].size() == 0)
        best_to_B[b_idx].push_back(&(it->al));
      else {
        assert(best_to_B[b_idx].size() == 1);
        if (best_to_B[b_idx][0]->frac_identity() < it->frac_identity)
          best_to_B[b_idx][0] = &(it->al);
      }
    }
  }
}

std::vector<size_t> make_random_permutation(size_t n)
{
  std::vector<size_t> x(n);
  for (size_t i = 0; i < n; ++i)
    x[i] = i;
  boost::random::random_device rng;
  boost::random::uniform_int_distribution<> index_dist(0, n - 1);
  for (size_t k = 0; k < n*10; ++k) {
    size_t i = index_dist(rng);
    size_t j = index_dist(rng);
    std::swap(x[i], x[j]);
  }
  return x;
}

template<typename PairsetType, typename AlignmentType>
void compute_alignment_stats(Stats& pair, Stats& nucl, Stats& tran,
                             std::vector<double>& b_frac_ones,
                             const std::vector<std::vector<const AlignmentType *> >& best_to_B,
                             const std::vector<std::string>&                         A,
                             const std::vector<double>&                              nu_A,
                             const std::vector<double>&                              tau_A,
                             const std::vector<std::string>&                         B,
                             const std::vector<double>&                              nu_B,
                             const std::vector<double>&                              tau_B,
                             const std::map<std::string, size_t>&                    A_names_to_idxs)
{
  std::vector<size_t> perm = make_random_permutation(nu_B.size()); // for better load balancing (in case e.g. seqs are ordered by length)

  pair.precis = 0.0, pair.recall = 0.0,
  nucl.precis = 0.0, nucl.recall = 0.0,
  tran.precis = 0.0, tran.recall = 0.0;

  #pragma omp parallel
  {
    double private_pair_precis = 0.0, private_pair_recall = 0.0,
           private_nucl_precis = 0.0, private_nucl_recall = 0.0,
           private_tran_precis = 0.0, private_tran_recall = 0.0;

    #pragma omp for
    for (int b_pre_idx = 0; b_pre_idx < static_cast<int>(nu_B.size()); ++b_pre_idx) {

      size_t b_idx = perm[b_pre_idx];
      PairsetType b_pairset(B[b_idx].size());
      mask b_mask(B[b_idx].size());

      BOOST_FOREACH(const AlignmentType *al, best_to_B[b_idx]) {

        size_t a_idx = A_names_to_idxs.find(al->a_name())->second;
        PairsetType a_pairset(A[a_idx].size());
        mask a_mask(A[a_idx].size());

        BOOST_FOREACH(const alignment_segment& seg, al->segments(A[a_idx], B[b_idx])) {
          a_pairset.add_square_with_exceptions(seg.a_start, seg.a_end, seg.a_mismatches.begin(), seg.a_mismatches.end());
          b_pairset.add_square_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
          a_mask.add_interval_with_exceptions(seg.a_start, seg.a_end, seg.a_mismatches.begin(), seg.a_mismatches.end());
          b_mask.add_interval_with_exceptions(seg.b_start, seg.b_end, seg.b_mismatches.begin(), seg.b_mismatches.end());
        }

        int a_num_total_bases = A[a_idx].size();
        int a_num_total_pairs = a_num_total_bases * (a_num_total_bases + 1) / 2;
        double a_frac_ones = 1.0 * a_mask.num_ones() / a_num_total_bases;
        private_pair_precis += (1.0 * a_pairset.size() / a_num_total_pairs) * nu_A[a_idx];
        private_nucl_precis += a_frac_ones * nu_A[a_idx];
        if (a_frac_ones >= 0.95) private_tran_precis += tau_A[a_idx];

      }

      int b_num_total_bases = B[b_idx].size();
      int b_num_total_pairs = b_num_total_bases * (b_num_total_bases + 1) / 2;
      b_frac_ones[b_idx] = (1.0 * b_mask.num_ones()) / b_num_total_bases;
      //std::cout << (1.0 * b_mask.num_ones()) / (1.0 * b_num_total_bases) << std::endl;
      private_pair_recall += (1.0 * b_pairset.size() / b_num_total_pairs) * nu_B[b_idx];
      private_nucl_recall += b_frac_ones[b_idx] * nu_B[b_idx];
      if (b_frac_ones[b_idx] >= 0.95) private_tran_recall += tau_B[b_idx];

    }

    #pragma omp critical
    {
      pair.precis += private_pair_precis; pair.recall += private_pair_recall;
      nucl.precis += private_nucl_precis; nucl.recall += private_nucl_recall;
      tran.precis += private_tran_precis; tran.recall += private_tran_recall;
    }

  } // omp parallel

  pair.update_F1();
  nucl.update_F1();
  tran.update_F1();
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
    ("induce-B-expr",                                "The oracleset expression, as produced by RSEM in a file called *.isoforms.results.")
    ("A-to-B", po::value<std::string>()->required(), "The alignments of A to B.")
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
  clock_t start;

  // Read the sequences and make sequence names-to-idxs maps
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
    std::string fname = vm["A-to-B"].as<std::string>();
    std::ifstream ifs;
    open_or_throw(ifs, fname);
    typename AlignmentType::input_stream_type input_stream(ifs);
    tic; read_alignments_and_filter_by_best_from_A(best_from_A, input_stream, A_names_to_idxs); toc;
  }

  std::cerr << "Clustering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > best_to_B(B_card);
  cluster_best_alignments_to_B(best_to_B, best_from_A, B_names_to_idxs);

  std::cerr << "Filtering alignments by B" << std::endl;
  std::vector<std::vector<const AlignmentType *> > single_best_to_B(B_card);
  filter_by_best_alignment_to_B(single_best_to_B, best_from_A, B_names_to_idxs);

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_A(A_card, 1.0/A_card);
  std::vector<double> unif_tau_B(B_card, 1.0/B_card);

  std::cerr << "Computing uniform nucleotide-level expression" << std::endl;
  std::vector<double> unif_nu_A(A_card), unif_nu_B(B_card);
  compute_nucl_expression(A, unif_tau_A, unif_nu_A);
  compute_nucl_expression(B, unif_tau_B, unif_nu_B);

  std::cerr << "Reading transcript-level expression for A" << std::endl;
  std::vector<double> tau_A(A_card), tau_B(B_card);
  std::string A_expr_fname = vm["A-expr"].as<std::string>();
  read_transcript_expression(A_expr_fname, tau_A, A_names_to_idxs);
  if (vm.count("induce-B-expr")) {
    std::cerr << "Inducing transcript-level expression for B" << std::endl;
    induce_prot_expression(tau_B, best_to_B, tau_A, A_names_to_idxs, A, B);
  } else {
    std::cerr << "Reading transcript-level expression for B" << std::endl;
    std::string B_expr_fname = vm["B-expr"].as<std::string>();
    read_transcript_expression(B_expr_fname, tau_B, B_names_to_idxs);
  }

  std::cerr << "Computing nucleotide-level expression" << std::endl;
  std::vector<double> nu_A(A_card), nu_B(B_card);
  compute_nucl_expression(A, tau_A, nu_A);
  compute_nucl_expression(B, tau_B, nu_B);

  std::cout << "summarize_version\t1" << std::endl;

  Stats pair, nucl, tran;
  std::vector<double> b_frac_ones(B_card);

  std::cerr << "Computing weighted stats" << std::endl;
  tic;
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones, best_to_B, A,      nu_A,      tau_A, B,      nu_B,      tau_B, A_names_to_idxs);
  print_stats(pair, "weighted_clustered_pair");
  print_stats(nucl, "weighted_clustered_nucl");
  print_stats(tran, "weighted_clustered_tran");
  toc;

  std::cerr << "Computing unweighted stats" << std::endl;
  tic;
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones, best_to_B, A, unif_nu_A, unif_tau_A, B, unif_nu_B, unif_tau_B, A_names_to_idxs);
  print_stats(pair, "unweighted_clustered_pair");
  print_stats(nucl, "unweighted_clustered_nucl");
  print_stats(tran, "unweighted_clustered_tran");
  toc;

  std::cerr << "Computing weighted stats" << std::endl;
  tic;
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones, single_best_to_B, A,      nu_A,      tau_A, B,      nu_B,      tau_B, A_names_to_idxs);
  print_stats(pair, "weighted_filtered_pair");
  print_stats(nucl, "weighted_filtered_nucl");
  print_stats(tran, "weighted_filtered_tran");
  toc;

  std::cerr << "Computing unweighted stats" << std::endl;
  tic;
  compute_alignment_stats<smart_pairset>(pair, nucl, tran, b_frac_ones, single_best_to_B, A, unif_nu_A, unif_tau_A, B, unif_nu_B, unif_tau_B, A_names_to_idxs);
  print_stats(pair, "unweighted_filtered_pair");
  print_stats(nucl, "unweighted_filtered_nucl");
  print_stats(tran, "unweighted_filtered_tran");
  toc;

  std::ofstream plot_out(vm["plot-output"].as<std::string>().c_str());
  for (size_t b_idx = 0; b_idx < b_frac_ones.size(); ++b_idx)
    plot_out << b_frac_ones[b_idx] << std::endl;
}
