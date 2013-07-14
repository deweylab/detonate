#pragma once
#include <iostream>
#include <vector>
#include <boost/program_options.hpp>
#include <boost/foreach.hpp>
#include <boost/make_shared.hpp>
#include <boost/algorithm/string/join.hpp>
#include <sparsehash/sparse_hash_map>
#define N_POLICY 2
#include "blast.hh"
#include "psl.hh"
#include "util.hh"

struct kmer_key
{
  string::const_iterator begin;
  string::const_iterator end;
  kmer_key() {}
  kmer_key(const string::const_iterator& begin, const string::const_iterator& end) : begin(begin), end(end) {}
};

struct kmer_val
{
  double numer, denom;
  kmer_val() : numer(0), denom(0) {}
};

namespace std
{
  namespace tr1
  {
    template<>
    struct hash<kmer_key>
    {
      size_t operator()(const kmer_key& k) const
      {
        return hash<string>()(string(k.begin, k.end));
      }
    };
  }

  template<>
  struct equal_to<kmer_key>
  {
    bool operator()(const kmer_key& lhs, const kmer_key& rhs) const
    {
      string::const_iterator i = lhs.begin;
      string::const_iterator j = rhs.begin;
      for (;;) {
        // If i == end and j == end, return true.
        // If i != end and j != end and *i == *j, continue;
        // In all other cases, return false.
        if (i == lhs.end) {
          if (j == rhs.end)
            return true;    // i == end and j == end              -> true
          else
            return false;   // i == end and j != end              -> false
        } else {
          if (j == rhs.end)
            return false;   // i != end and j == end              -> false
          else {
            if (*i != *j)
              return false; // i != end and j != end and *i != *j -> false
            else {
            }               // i != end and j != end and *i == *j -> continue
          }
        }
        ++i, ++j;
      }
    }
  };
}

size_t estimate_hashtable_size(
  const std::vector<std::string>& B,
  size_t kmerlen)
{
  // The fudge factor is based partly on empirical observation and partly on
  // the assumption that most kmers will be shared between the oracleset and
  // the assembly.
  size_t fudge_factor = 2;
  size_t max_entries = 0;
  BOOST_FOREACH(const string& a, B)
    max_entries += 2 * (a.size() + 1 - kmerlen) / fudge_factor;
  return max_entries;
}

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

template<typename AlignmentType>
void compute_and_print_stats(
    typename AlignmentType::input_stream_type& input_stream,
    const std::vector<std::string>&            A,
    const std::vector<std::string>&            B,
    const std::map<std::string, size_t>&       A_names_to_idxs,
    const std::map<std::string, size_t>&       B_names_to_idxs,
    const std::vector<double>&                 tau_B,
    bool                                       strand_specific,
    size_t                                     kmerlen,
    const std::string&                         prefix,
    const std::string&                         suffix)
{
  // Make masks. For example, if k = 3, and a sequence is of length 5, then the
  // mask should be of length seq_len-k+1 = 5-3+1 = 3, since 3 kmer positions
  // are possible.
  //   seq:   01234
  //   kmers: 012
  //           123
  //            234
  std::vector<std::vector<bool> > masks(B.size());
  for (size_t b_idx = 0; b_idx < B.size(); ++b_idx)
    if (B[b_idx].size() >= kmerlen)
      masks[b_idx].resize(B[b_idx].size() + 1 - kmerlen, false);

  // Find good kmers.
  AlignmentType al;
  while (input_stream >> al) {
    if (is_valid(al, strand_specific)) {
      size_t a_idx = A_names_to_idxs.find(al.a_name())->second;
      size_t b_idx = B_names_to_idxs.find(al.b_name())->second;
      std::vector<bool>& m = masks[b_idx];
      BOOST_FOREACH(const alignment_segment& seg, al.segments(A[a_idx], B[b_idx])) {

        if (seg.b_start > seg.b_end)
          throw std::runtime_error("Fatal error: b_start > b_end.");

        // If this segment is shorter than the kmer length k, then we should
        // skip adding it.
        // 
        // E.g., k=3, lo=4, hi=6: 456  -> ok,     6-4+1 = 3 >= 3
        //       k=3, lo=4, hi=5: 45   -> not ok, 5-4+1 = 2 < 3
        //       k=4, lo=4, hi=6: 456  -> not ok, 6-4+1 = 3 < 4
        //       k=3, lo=4, hi=7: 4567 -> ok,     7-4+1 = 4 >= 3
        // So the requirement is that hi-lo+1 >= k, i.e., hi+1-lo >= k.
        if (seg.b_end + 1 - seg.b_start < kmerlen)
          continue;

        // seq:   01234 - here b_start = 0, b_end = 4
        // kmers: 012
        //         123
        //          234
        // Thus, the last kmer starts at position 2 = 4-3+1 = b_end-k+1.
        for (size_t i = seg.b_start; i <= seg.b_end + 1 - kmerlen; ++i)
          m[i] = true;

      }
    }
  }

  // Count good kmers and possible kmers.
  size_t est_size = estimate_hashtable_size(B, kmerlen);
  google::sparse_hash_map<kmer_key, kmer_val> ht(est_size);
  for (size_t b_idx = 0; b_idx < B.size(); ++b_idx) {
    const string& b = B[b_idx];
    const std::vector<bool>& m = masks[b_idx];
    if (b.size() >= kmerlen) {
      double c = tau_B[b_idx];
      string::const_iterator beg = b.begin();
      string::const_iterator end = b.begin() + kmerlen;
      std::vector<bool>::const_iterator m_it = m.begin();
      for (; end != b.end(); ++beg, ++end, ++m_it) {
        kmer_key key(beg, end);
        kmer_val& val = ht[key];
        if (*m_it)
          val.numer += c;
        val.denom += c;
      }
    }
  }

  // Compute statistics.
  typedef std::pair<const kmer_key, kmer_val> X;
  double numer = 0.0, denom = 0.0;
  BOOST_FOREACH(const X& x, ht) {
    const kmer_val& val = x.second;
    numer += val.numer;
    denom += val.denom;
  }
  double ratio = numer / denom;

  // Output statistics.
  cout << prefix << "_aligned_kmer_frac_present_numer_" << suffix << "\t" << numer << endl;
  cout << prefix << "_aligned_kmer_frac_present_denom_" << suffix << "\t" << denom << endl;
  cout << prefix << "_aligned_kmer_frac_present_"       << suffix << "\t" << ratio << endl;
}
      
void parse_options(boost::program_options::variables_map& vm, int argc, const char **argv)
{
  namespace po = boost::program_options;

  po::options_description desc("Options");
  desc.add_options()
    ("help,?", "Display this information.")
    ("A-seqs", po::value<std::string>()->required(), "The assembly sequences, in FASTA format.")
    ("B-seqs", po::value<std::string>()->required(), "The oracleset sequences, in FASTA format.")
    ("B-expr", po::value<std::string>(),             "The oracleset expression, as produced by RSEM in a file called *.isoforms.results.")
    ("no-expr",                                      "Do not use expression at all. No weighted scores will be produced.")
    ("A-to-B", po::value<std::string>()->required(), "The alignments of A to B.")
    ("alignment-type", po::value<std::string>()->required(), "The type of alignments used, either 'blast' or 'psl'.")
    ("strand-specific",                              "Ignore alignments that are to the reverse strand.")
    ("readlen", po::value<size_t>()->required(),     "The read length.")
  ;

  try {

    po::store(po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help")) {
      std::cerr << desc << std::endl;
      exit(1);
    }

    po::notify(vm);

    if (vm.count("no-expr")) {
      if (vm.count("B-expr") != 0)
        throw po::error("If --no-expr is given, then --B-expr cannot be given.");
    } else {
      if (vm.count("B-expr") == 0)
        throw po::error("If --no-expr is not given, then --B-expr must be given.");
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
  size_t readlen = vm["readlen"].as<size_t>();

  std::cerr << "Reading the sequences" << std::endl;
  std::vector<std::string> A, B;
  std::vector<std::string> A_names, B_names;
  std::map<std::string, size_t> A_names_to_idxs, B_names_to_idxs;
  read_fasta_names_and_seqs(vm["A-seqs"].as<std::string>(), A, A_names, A_names_to_idxs);
  read_fasta_names_and_seqs(vm["B-seqs"].as<std::string>(), B, B_names, B_names_to_idxs);

  std::vector<double> tau_B(B.size());
  if (!vm.count("no-expr")) {
    std::cerr << "Reading transcript-level expression for A and B" << std::endl;
    std::string B_expr_fname = vm["B-expr"].as<std::string>();
    read_transcript_expression(B_expr_fname, tau_B, B_names_to_idxs);
  }

  std::cerr << "Computing uniform transcript-level expression" << std::endl;
  std::vector<double> unif_tau_B(B.size(), 1.0);

  std::cout << "summarize_aligned_kmer_version_1\t0\n";

  if (!vm.count("no-expr")) {
    typename AlignmentType::input_stream_type A_to_B_is(open_or_throw(vm["A-to-B"].as<std::string>()));
    compute_and_print_stats<AlignmentType>(A_to_B_is, A, B, A_names_to_idxs, B_names_to_idxs, tau_B, strand_specific, readlen, "weighted", "at_one");
  }

  {
    typename AlignmentType::input_stream_type A_to_B_is(open_or_throw(vm["A-to-B"].as<std::string>()));
    compute_and_print_stats<AlignmentType>(A_to_B_is, A, B, A_names_to_idxs, B_names_to_idxs, unif_tau_B, strand_specific, readlen, "unweighted", "at_one");
  }

  std::cerr << "Done!" << std::endl;
}
