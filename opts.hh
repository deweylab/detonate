#pragma once

struct opts
{
  // Scores
  bool nucl;
  bool contig;
  bool pair;
  bool kpair;
  bool kmer;
  bool kc;

  // Score meta-groups
  bool alignment_based;
  bool alignment_free;

  // Weighted
  bool weighted;
  bool unweighted;

  // Paper
  bool paper;

  // Sequences
  std::string A_seqs;
  std::string B_seqs;

  // Expression
  std::string A_expr;
  std::string B_expr;

  // Alignments
  std::string A_to_B;
  std::string B_to_A;
  std::string alignment_type;

  // Strand-specific
  bool strand_specific;

  // Read length and num reads
  size_t readlen;
  size_t num_reads;

  // Contig thresholds
  double contig_min_frac_identity;
  double contig_max_frac_indel;

  // Hash table
  std::string hash_table_type;
  double hash_table_fudge_factor;

  opts()
  :
    // Scores
    nucl(false),
    contig(false),
    pair(false),
    kpair(false),
    kmer(false),
    kc(false),

    // Score meta-groups
    alignment_based(false),
    alignment_free(false),

    // Weighted
    weighted(false),
    unweighted(false),

    // Paper
    paper(false),

    // Sequences
    A_seqs(""),
    B_seqs(""),

    // Expression
    A_expr(""),
    B_expr(""),

    // Alignments
    A_to_B(""),
    B_to_A(""),
    alignment_type(""),

    // Strand-specific
    strand_specific(false),

    // Read length
    readlen(-1),
    num_reads(0),

    // Contig thresholds
    contig_min_frac_identity(2.0),
    contig_max_frac_indel(-1.0),

    // Hash table
    hash_table_type(""),
    hash_table_fudge_factor(-1.0)
  {}
};

std::ostream& operator<<(std::ostream& os, const opts& o)
{
  if (o.nucl)   os << "nucl"   << "\n";
  if (o.contig) os << "contig" << "\n";
  if (o.pair)   os << "pair"   << "\n";
  if (o.kpair)  os << "kpair"   << "\n";
  if (o.kmer)   os << "kmer"   << "\n";
  if (o.kc)     os << "kc"     << "\n";

  if (o.alignment_based) os << "alignment_based" << "\n";
  if (o.alignment_free)  os << "alignment_free"  << "\n";

  if (o.weighted)   os << "weighted"   << "\n";
  if (o.unweighted) os << "unweighted" << "\n";

  if (o.paper) os << "paper" << "\n";

  if (o.A_seqs.size()) os << "A_seqs: " << o.A_seqs << "\n";
  if (o.B_seqs.size()) os << "B_seqs: " << o.B_seqs << "\n";

  if (o.A_expr.size()) os << "A_expr: " << o.A_expr << "\n";
  if (o.B_expr.size()) os << "B_expr: " << o.B_expr << "\n";

  if (o.A_to_B.size())         os << "A_to_B: "         << o.A_to_B << "\n";
  if (o.B_to_A.size())         os << "B_to_A: "         << o.B_to_A << "\n";
  if (o.alignment_type.size()) os << "alignment_type: " << o.alignment_type << "\n";

  if (o.strand_specific) os << "strand_specific"        << "\n";
  if (o.readlen)         os << "readlen: " << o.readlen << "\n";

  return os;
}
