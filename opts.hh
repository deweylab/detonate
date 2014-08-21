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

  // kmer length
  size_t kmerlen;

  // Contig thresholds
  double min_frac_identity;
  double max_frac_indel;

  // Minimum segment length
  size_t min_segment_len;

  // Hash table
  std::string hash_table_type;
  double hash_table_fudge_factor;

  // Trace output
  std::string trace;

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

    // Read length and num reads
    readlen(-1),
    num_reads(0),

    // kmer length
    kmerlen(-1),

    // Contig thresholds
    min_frac_identity(2.0),
    max_frac_indel(-1.0),

    // Minimum segment length
    min_segment_len(0),

    // Hash table
    hash_table_type(""),
    hash_table_fudge_factor(-1.0),

    // Trace output
    trace("")
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
