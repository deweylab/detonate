// Copyright (c) 2013
// Nathanael Fillmore and Bo Li (University of Wisconsin-Madison)
// nathanae@cs.wisc.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include <algorithm>
#include <vector>
#include <list>
#include <string>
#include <stdio.h>
#include <time.h>
#include <boost/program_options.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/discrete_distribution.hpp>
#include "sam/bam.h"
#include "fasta.hh"
#include "util.hh"
#include "re_eta_help.hh"

struct probs_and_lidxs
{
  std::vector<size_t> lidxs;
  std::vector<double> probs;
};

struct ReadStruct
{
  // The read id, the transcript id (of the transcript the read is aligned to),
  // the leftmost position of the read within its transcript, and the read
  // length.
  int rid, tid, pos, readlen;

  ReadStruct(int rid, int tid, int pos, int readlen)
  {
    this->rid = rid;
    this->tid = tid;
    this->pos = pos;
    this->readlen = readlen;
  }

  bool operator< (const ReadStruct& o) const
  {
    if (tid != o.tid) return tid < o.tid;
    return pos < o.pos;
  }
};

struct ContigStruct
{
  // The transcript id (of the transcript this contig is a subsequence of), the
  // leftmost position of the contig within its transcript, and the contig
  // length.
  int tid, pos, len;

  // The index of *some* leftward contig in the scaffold. The parent is not
  // necessarily, however, either immediately leftward (i.e., the left-adjacent
  // contig in the scaffold) or as far left as possible (i.e., the root, or the
  // first contig in the scaffold). The find_and_normalize() function
  // normalizes this relation.
  int parent;

  // The next (immediately rightward) contig in the scaffold.
  int next;

  ContigStruct(int tid, int pos, int len, int parent)
  {
    this->tid = tid;
    this->pos = pos;
    this->len = len;
    this->parent = parent;
    this->next = -1;
  }
};

// Determine the read length of the mate read, assuming that there are
// no insertions or deletions in the alignment of the read or mate.
// This is adapted from Bio-SamTools-1.39,
// http://search.cpan.org/~lds/Bio-SamTools/
int mate_len(const bam1_t *rec)
{
  // Extract the insert length and check that it is nonzero.
  assert(rec->core.isize != 0);
  int ins_len = rec->core.isize;

  // Extract the (non-mate) read's start position, and check that it
  // is valid and that the read is indeed mapped.
  assert(rec->core.pos >= 0);
  assert(!(BAM_FUNMAP & rec->core.flag));
  int start = rec->core.pos + 1;

  // Extract the the mate's start position, and check that it is valid
  // and that the mate is indeed mapped.
  assert(rec->core.mpos >= 0);
  assert(!(BAM_FMUNMAP & rec->core.flag));
  int mate_start = rec->core.mpos + 1;

  // Actually compute the mate's length.
  int mate_len;
  if (ins_len > 0) {
    mate_len = ins_len + (start - mate_start);
  } else {
    mate_len = mate_start - (start + ins_len);
  }
  assert(mate_len > 0);

  return mate_len;
}

bool is_ok(const bam1_t *rec, bool paired)
{
  // In paired-end mode, we exclude alignments that are focused on the second
  // (right) mate of each pair. The reason for doing this is that the alignment
  // focused on the first (left) mate also contains information about the
  // second mate's position, and for sampling alignments, it is simpler to have
  // one alignment per read pair than to have two alignments per read pair.
  if (paired) {
    if (!(rec->core.flag & BAM_FPROPER_PAIR)) {
      throw std::runtime_error("Found a read that is not in a proper pair, "
                               "even though we are running in paired-end "
                               "mode. Currently, this is not supported. If "
                               "this feature would be of interest to you, "
                               "please contact us.");
    }
    if (rec->core.flag & BAM_FREAD2) // this is read2
      return false;
  }

  // Exclude noise reads.
  if (rec->core.qual == 0) // this is a noise read
    return false;
  if (rec->core.tid < 0)
    throw std::runtime_error("Found an alignment with qual > 0 "
        "but tid < 0. Was your BAM file produced by RSEM?");

  return true;
}

// Calculate the posterior probability that this read follows this alignment,
// as calculated by RSEM.
double compute_prob(const bam1_t *rec)
{
  // According to rsem-calculate-expression -?, ``The MAPQ field of each
  // alignment is set to min(100, floor(-10 * log10(1.0 - w) + 0.5)), where w
  // is the posterior probability of that alignment being the true mapping of
  // a read.''
  //
  // To go backwards, we solve:
  //        min(100, floor(-10 * log10(1.0 - w) + 0.5) = q
  // i.e.,           round(-10 * log10(1.0 - w))       = q
  // i.e.,                 -10 * log10(1.0 - w)        = q
  // i.e.,                       log10(1.0 - w)        = -q/10
  // i.e.,                             1.0 - w         = 10^(-q/10)
  // i.e.,                                 - w         = 10^(-q/10) - 1.0
  // i.e.,                                   w         = -( 10^(-q/10) - 1.0 )
  //
  // As a doublecheck: In sam_rsem_cvt.h, there is the following function:
  //
  //   uint8_t getMAPQ(double val) {
  //       double err = 1.0 - val;
  //       if (err <= 1e-10) return 100;
  //       return (uint8_t)(-10 * log10(err) + .5); // round it
  //   }
  //
  // If we plug in val = -(10^(-q/10) - 1.0), we get
  //
  //    -10 * log10(err)
  // =  -10 * log10(1.0 - val)
  // =  -10 * log10(1.0 - (-(10^(-q/10) - 1.0)))
  // =  -10 * log10(1.0 + 10^(-q/10) - 1.0)
  // =  -10 * log10(      10^(-q/10)      )
  // =  -10 * (               -q/10       )
  // =  -     (               -q          )
  // =                         q           
  //
  // as expected if the above is correct.
  double prob = -(pow(10.0, -rec->core.qual/10.0) - 1.0);
  if (prob < 0.0 || prob > 1.0)
    throw std::runtime_error("Posterior probability of an alignment is < 0 or > 1.");
  return prob;
}

size_t extract_alignment_probs(std::map<std::string, probs_and_lidxs>& read_to_probs_and_lidxs,
                               const std::string& bam_fname,
                               double min_alignment_prob,
                               bool paired)
{
  // Open BAM file and read header.
  bamFile fi = bam_open(bam_fname.c_str(), "r");
  if (!fi)
    throw std::runtime_error("Cannot open BAM file: " + bam_fname);
  bam_header_t *header = bam_header_read(fi);

  // Iterate over lines of the bamfile.
  bam1_t *rec = bam_init1();
  size_t lidx;
  for (lidx = 0;; ++lidx) { // alignment idx
    
    // Read next record. If no more records, we're done.
    if (bam_read1(fi, rec) < 0)
      break;

    // Check that this isn't an alignment to the noise transcript, etc.
    if (!is_ok(rec, paired))
      continue;

    // Calculate the posterior probability that this read follows this
    // alignment, as calculated by RSEM.
    double prob = compute_prob(rec);

    // Check that the posterior probability is big enough.
    if (prob < min_alignment_prob)
      continue;

    // Actually add the alignment.
    const char *read_name = bam1_qname(rec);
    probs_and_lidxs& pal = read_to_probs_and_lidxs[read_name];
    pal.probs.push_back(prob);
    pal.lidxs.push_back(lidx);

  }

  // Clean up.
  bam_destroy1(rec);
  bam_header_destroy(header);

  // Return the number of alignments.
  return lidx;
}

void choose_alignments(std::vector<bool>& is_lidx_chosen,
                       boost::random::mt19937& rng,
                       const std::map<std::string, probs_and_lidxs>& read_to_probs_and_lidxs,
                       const std::string& alignment_policy)
{
  std::map<std::string, probs_and_lidxs>::const_iterator it;

  if (alignment_policy == "sample") {
    for (it = read_to_probs_and_lidxs.begin();
         it != read_to_probs_and_lidxs.end();
         ++it) {
      boost::random::discrete_distribution<> dd(it->second.probs);
      size_t i = dd(rng);
      size_t chosen_lidx = it->second.lidxs[i];
      is_lidx_chosen[chosen_lidx] = true;
    }
  }

  else if (alignment_policy == "best") {
    for (it = read_to_probs_and_lidxs.begin();
         it != read_to_probs_and_lidxs.end();
         ++it) {
      size_t num_alignments = it->second.probs.size();
      size_t best = 0;
      if (num_alignments > 1) {
        for (size_t i = 1; i < num_alignments; ++i)
          if (it->second.probs[i] > it->second.probs[best])
            best = i;
      }
      size_t chosen_lidx = it->second.lidxs[best];
      is_lidx_chosen[chosen_lidx] = true;
    }
  }

  else if (alignment_policy == "all") {
    for (it = read_to_probs_and_lidxs.begin();
         it != read_to_probs_and_lidxs.end();
         ++it) {
      size_t num_alignments = it->second.probs.size();
      for (size_t i = 0; i < num_alignments; ++i)
        is_lidx_chosen[it->second.lidxs[i]] = true;
    }
  }

  else {
    throw std::runtime_error("Unknown alignment policy: " + alignment_policy);
  }
}

void extract_reads(std::vector<ReadStruct>& reads,
                   const std::vector<bool>& is_lidx_chosen,
                   const std::string& bam_fname,
                   bool paired)
{
  // Open BAM file and read header.
  bamFile fi = bam_open(bam_fname.c_str(), "r");
  if (!fi)
    throw std::runtime_error("Cannot open BAM file: " + bam_fname);
  bam_header_t *header = bam_header_read(fi);

  // Iterate over lines of the bamfile.
  bam1_t *rec = bam_init1();
  for (size_t lidx = 0;; ++lidx) { // alignment idx

    // Read next record. If no more records, we're done.
    if (bam_read1(fi, rec) < 0)
      break;

    // Check that this alignment has been chosen.
    if (!is_lidx_chosen[lidx])
      continue;

    // Extract the transcript idx, the leftmost position of the read, and
    // the read length.
    int tid = rec->core.tid;
    int pos = rec->core.pos;
    int readlen = rec->core.l_qseq;

    // Record the above info about this read. This is the main info that will
    // be relevant for making the contigs later on.
    //
    // We use the alignment idx (lidx) as the read id because the purpose of
    // the read id is simply to associate the read and its mate. Note that if
    // all alignments are used (--alignment-policy=all), then a read might
    // occur more than once, but for the purposes of assembly, each alignment
    // of the read should be added with a different read id, since otherwise we
    // would try to scaffold between the various alignments of the same read,
    // which doesn't make sense.
    reads.push_back(ReadStruct(lidx, tid, pos, readlen));

    // In the paired-end case, extract info about the mate read, and add it to
    // the list of reads, too.
    if (paired) {
      // Extract the mate's transcript idx, the leftmost position of the mate,
      // and the mate's length.
      int mtid = rec->core.mtid;
      int mpos = rec->core.mpos;
      int mreadlen = mate_len(rec);

      // Check that the read and its mate align to the same transcript.
      if (mtid != tid) {
        std::ostringstream ss;
        ss << "Read " << bam1_qname(rec) << " is aligned to transcript " << tid
           << " but its mate is aligned to transcript " << mtid
           << ". Were your alignments produced by RSEM?";
        throw std::runtime_error(ss.str());
      }

      // Check that the read and its mate have the same length. This doesn't
      // necessarily need to be the case, but it should be the case in our
      // initial experiment.
      if (mreadlen != readlen) {
        std::cerr << "Warning: Read " << bam1_qname(rec)  << " has length " 
                  << readlen << ", but its mate has imputed length "
                  << mreadlen << ". You might want to check that this is "
                  << "correct." << std::endl;
      }

      reads.push_back(ReadStruct(lidx, mtid, mpos, mreadlen));
    }

  }

  // Clean up.
  bam_destroy1(rec);
  bam_header_destroy(header);
}

// Returns the root of the scaffold involving u.
int find(int u, const std::vector<ContigStruct>& contigs)
{
  while (contigs[u].parent != u)
    u = contigs[u].parent;
  return u;
}

// Returns the root of the scaffold involving u, and updates all the contigs
// between u and the root to point at the root as their parent.
int find_and_normalize(int u, std::vector<ContigStruct>& contigs)
{
  // Remember the given contig.
  int x = u;
  // Set u as the root, i.e., as the first contig in the scaffold.
  u = find(u, contigs);
  // Go from the given contig toward the root, and for each contig
  // encountered, set the root as its parent.
  while (contigs[x].parent != u) {
    int tmp = contigs[x].parent;
    contigs[x].parent = u;
    x = tmp;
  }
  // Return the root.
  return u;
}

void assemble(std::vector<ContigStruct>& contigs,
              std::vector<ReadStruct>& reads,
              int min_overlap,
              bool paired,
              std::vector<int> read2contig)
{
  int s = -1;
  size_t num_reads = reads.size();

  if (paired) {
    read2contig.clear();
    read2contig.resize(num_reads, -1);
  }

  // For each read (each mate separately, in the case of paired-end data), ...
  for (size_t i = 0; i < num_reads; i++) {
    // The reads should be sorted in the BAM file, by transcript and then
    // position within the transcript, so the last contig's start position
    // should be before this read's position. The last contig should also end
    // before this read ends.
    if (contigs.size() > 0 && contigs[s].tid == reads[i].tid) {
      assert(contigs[s].pos <= reads[i].pos);
      assert(contigs[s].pos + contigs[s].len <= reads[i].pos + reads[i].readlen);
    }
    // If (1) there are no contigs, or (2) this read comes from a new
    // transcript (different from that of the previous contig, and hence from
    // that of all other previous contigs, since the reads are sorted by
    // transcript), or (3) this read starts after too much of a gap from the
    // current contig, then start a new contig.
    if (contigs.size() == 0 ||
        contigs[s].tid != reads[i].tid ||
        contigs[s].pos + contigs[s].len - reads[i].pos < min_overlap) {
      ++s;
      contigs.push_back(ContigStruct(reads[i].tid, reads[i].pos, reads[i].readlen, s));
    }
    // Otherwise, extend the current contig with this read.
    else {
      contigs[s].len = std::max(contigs[s].len,
                                reads[i].pos + reads[i].readlen -
                                contigs[s].pos);
    }

    if (paired) {
      // If this read (or more properly, its mate) has not been already assigned
      // to a contig, then assign it to the current contig.
      if (read2contig[reads[i].rid] < 0) {
        read2contig[reads[i].rid] = s;
      }
      // If this read (or more properly, its mate) *has* been already assigned
      // to a contig, then build a scaffold between that contig and the contig
      // we built or extended above based on the current mate. Of course, don't
      // do this if the two contigs are the same (fx == fy).
      else {
        int fx = find_and_normalize(read2contig[reads[i].rid], contigs);
        int fy = find_and_normalize(s, contigs);

        if (fx < fy) { contigs[fy].parent = fx; }
        else if (fx > fy) { contigs[fx].parent = fy; }
      }
    }
  }

  ++s; // there are s contigs

  if (paired) {
    // For each contig i, starting from the last one and going backwards,
    // normalize the path (make all contigs have the scaffold's root as their
    // parent) and let fx be the scaffold's root. If the current contig i is not
    // the root, then point i's next node as the root's current next node, and
    // point the root's next node as the current node. In pictures, assuming
    // that so far we've processed contigs j, k, and l in this contig:
    //
    //    root ---------> j -> k -> l -> -1
    //
    // is changed to
    //
    //    root ----> i -> j -> k -> l -> -1
    for (int i = s - 1; i >= 0; i--) {
      int fx = find_and_normalize(i, contigs);
      if (fx != i) {
        contigs[i].next = contigs[fx].next;
        contigs[fx].next = i;
      }
    }
  }
}

void output(const std::string& output_fname,
            const std::vector<ContigStruct>& contigs,
            const fasta& ref_fa,
            const std::string& bam_fname,
            bool paired)
{
  // Open BAM file and read header.
  bamFile fi = bam_open(bam_fname.c_str(), "r");
  if (!fi)
    throw std::runtime_error("Cannot open BAM file: " + bam_fname);
  bam_header_t *header = bam_header_read(fi);

  // Open file to write output.
  std::ofstream out(output_fname.c_str());
  if (!out)
    throw std::runtime_error("Cannot open file to write assembly: " + output_fname);

  // Initialize the number of scaffolds, for paired-end mode.
  int ns = 0;

  // Assemble each contig.
  for (int i = 0; i < static_cast<int>(contigs.size()); i++) {
    // In paired-end mode, skip contigs that are not the first element of their
    // scaffold.
    if (paired)
      if (find(i, contigs) != i)
        continue;

    // We have recorded, in contigs[i].tid, the idx of the transcript that this
    // contig comes from, relative to the ordering of transcripts in the BAM
    // file. However, the transcript sequences are stored in the FASTA file,
    // and the order here could be different. Thus, here, we will figure out
    // the idx of the transcript relative to the order in the FASTA file.
    std::string tname(header->target_name[contigs[i].tid]);
    std::map<std::string, size_t>::const_iterator it;
    it = ref_fa.names_to_idxs.find(tname);
    if (it == ref_fa.names_to_idxs.end()) {
      throw std::runtime_error("The following transcript name does not "
          "correspond to a reference sequence: " + tname);
    }
    size_t tidx = it->second;
    if (static_cast<int>(tidx) != contigs[i].tid) {
      throw std::runtime_error("tidx and tid don't match");
    }

    // Get the transcript sequence and extract the contig. (In the paired-end
    // case, this is the first contig of the scaffold.)
    const std::string& refseq = ref_fa.seqs[contigs[i].tid];
    if (contigs[i].pos + contigs[i].len > static_cast<int>(refseq.length())) {
      std::ostringstream ss;
      ss << "Contig extends past end of the transcript sequence. "
         << "Info for debugging: "
         << "tname="  << tname << ", " 
         << "pos="    << contigs[i].pos << ", " 
         << "length=" << contigs[i].len << ", " 
         << "transcript_length=" << refseq.length();
      throw std::runtime_error(ss.str());
    }
    std::string seq = refseq.substr(contigs[i].pos, contigs[i].len);
    int tpos = contigs[i].pos + contigs[i].len - 1;

    // In single-end mode, simply output the contig.
    if (!paired) {
      // Output the contig in FASTA format.
      out << ">contig_" << i
          << " tname="  << tname
          << " pos="    << contigs[i].pos 
          << " length=" << contigs[i].len << "\n";
      out << seq << "\n";
    }

    // In paired-end mode, finish constructing the scaffold, and then output it.
    else {

      // Record the start and length of the first contig of the scaffold, for
      // inclusion in the scaffold title.
      std::ostringstream pos_ss;
      std::ostringstream len_ss;
      pos_ss << contigs[i].pos;
      len_ss << contigs[i].len;

      // For each remaining contig in the scaffold, ...
      int nc = 1;
      int p = contigs[i].next;
      while (p >= 0) {
        ++nc;

        // Figure out how many N's should be between the end of the previous
        // contig and the beginning of this contig. Then add them to the
        // scaffold. Note that without the "max(., 0)", numN could be negative,
        // if the min_overlap is greater than 0. For example, with min_overlap
        // = 2:
        //
        // transcript:   ---------------------------
        // reads:           -----
        //                   -----
        //                       -----
        // contigs:         ------
        //                       -----
        // 
        // Here, numN would be -1 without the "max(., 0)".
        int numN = std::max(contigs[p].pos - tpos - 1, 0);
        if (numN > 0)
          seq.append(numN, 'n');
        tpos += numN;

        // Recompute the contig length again in a convoluted way, and if this
        // length is nonzero, add the contig sequence to the scaffold. The
        // convolutions are important in examples like the one above, where two
        // contigs overlap. In the case of no overlap, the calculation
        // simplifies, since in this case
        //
        //   tpos == old_tpos + numN
        //        == old_tpos + contigs[p].pos - old_tpos - 1
        //        == contigs[p].pos - 1
        //
        // so
        //
        //   rlen == contigs[p].pos + contigs[p].len - tpos - 1
        //        == contigs[p].pos + contigs[p].len - (contigs[p].pos - 1) - 1
        //        == contigs[p].pos + contigs[p].len - contigs[p].pos + 1 - 1
        //        == contigs[p].len
        //
        // and therefore in this case (namely, where the current contig doesn't
        // overlap the previous one) the below is equivalent to:
        //
        //   seq += refseq.substr(contigs[p].pos, contigs[p].len);
        //   tpos += contigs[p].len;
        //
        // In the case of overlap, only the "new", non-overlapping part of the
        // current contig will be added to the scaffold sequence. For example,
        // continuing the example above:
        //
        // transcript:   ---------------------------
        // contigs:         ------
        //                       -----
        // scaffold:        ----------
        int rlen = contigs[p].pos + contigs[p].len  - tpos - 1;
        if (rlen > 0) {
          seq += refseq.substr(tpos + 1, rlen);
          tpos += rlen;
        }

        // Record the start and length of this contig, for inclusion in the
        // scaffold title.
        pos_ss << "," << contigs[p].pos;
        len_ss << "," << contigs[p].len;

        p = contigs[p].next;
      }

      // Output the scaffold in FASTA format.
      out << ">scaffold_" << ns
          << " tname="  << tname
          << " pos="    << pos_ss.str()
          << " length=" << len_ss.str() << "\n";
      out << seq << "\n";
      ++ns;

    }

  }

  // Clean up.
  bam_header_destroy(header);
}

void print_help()
{
  std::cout << get_help_string() << std::endl;
}

boost::program_options::options_description describe_options()
{
  namespace po = boost::program_options;
  po::options_description desc;
  desc.add_options()
    ("help,?", "Display this information.")
    ("reference", po::value<std::string>())
    ("expression", po::value<std::string>())
    ("assembly", po::value<std::string>())
    ("paired-end", "Flag")
    ("min-overlap", po::value<std::string>())
    ("min-alignment-prob", po::value<double>())
    ("alignment-policy", po::value<std::string>())
    ;
  return desc;
}

struct opts
{
  std::string reference;
  std::string expression;
  std::string assembly;
  bool paired;
  size_t min_overlap_begin;
  size_t min_overlap_end;
  double min_alignment_prob;
  std::string alignment_policy;
};

void parse_options(opts& o, const boost::program_options::variables_map& vm)
{
  namespace po = boost::program_options;

  if (!vm.count("reference"))
    throw po::error("--reference is required.");
  o.reference = vm["reference"].as<std::string>();

  if (!vm.count("expression"))
    throw po::error("--expression is required.");
  o.expression = vm["expression"].as<std::string>();

  if (!vm.count("assembly"))
    throw po::error("--assembly is required.");
  o.assembly = vm["assembly"].as<std::string>();

  o.paired = vm.count("paired-end");

  std::string mo;
  if (vm.count("min-overlap"))
    mo = vm["min-overlap"].as<std::string>();
  else
    mo = "0";

  size_t tpos = mo.find_first_of(',');
  if (tpos == std::string::npos) {
    o.min_overlap_begin = boost::lexical_cast<size_t>(mo);
    o.min_overlap_end   = o.min_overlap_begin;
  }
  else {
    o.min_overlap_begin = boost::lexical_cast<size_t>(mo.substr(0, tpos));
    o.min_overlap_end   = boost::lexical_cast<size_t>(mo.substr(tpos + 1));
  }

  if (vm.count("alignment-policy")) {
    o.alignment_policy = vm["alignment-policy"].as<std::string>();
    if (o.alignment_policy != "sample" &&
        o.alignment_policy != "best" &&
        o.alignment_policy != "all")
      throw po::error("Invalid value for --alignment-policy: " + o.alignment_policy);
  } else {
    o.alignment_policy = "sample";
  }

  if (vm.count("min-alignment-prob")) {
    o.min_alignment_prob = vm["min-alignment-prob"].as<double>();
    if (o.min_alignment_prob < 0.0 || o.min_alignment_prob > 1.0)
      throw po::error("Invalid value for --min-alignment-prob: " + vm["min-alignment-prob"].as<std::string>());
  } else {
    o.min_alignment_prob = 0.0;
  }

}

int main(int argc, const char **argv)
{
  namespace po = boost::program_options;

  try {

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

    // Read transcripts and extract a map from transcript names to transcript
    // idxs. The idxs start at 1, with 0 being the noise transcript.
    std::cerr << curtime() << "Reading the transcripts..." << std::flush;
    fasta ref_fa;
    std::string ref_fname = o.reference + ".transcripts.fa";
    read_fasta(ref_fa, ref_fname);
    std::cerr << "done." << std::endl;

    // Read alignments and extract the posterior probability of each alignment.
    std::cerr << curtime() << "Reading the alignments and extracting the posterior probability of each alignment..." << std::flush;
    std::map<std::string, probs_and_lidxs> read_to_probs_and_lidxs;
    std::string bam_fname = o.expression + ".transcript.sorted.bam";
    size_t num_alignments = extract_alignment_probs(read_to_probs_and_lidxs, bam_fname, o.min_alignment_prob, o.paired);
    std::cerr << "done." << std::endl;

    // Select a subset of alignments (currently: by sampling one alignment for
    // each read, according to the posterior probability of each alignment).
    std::cerr << curtime() << "Selecting a subset of alignments..." << std::flush;
    std::vector<bool> is_lidx_chosen(num_alignments);
    boost::random::mt19937 rng(time(0));
    choose_alignments(is_lidx_chosen, rng, read_to_probs_and_lidxs, o.alignment_policy);
    std::cerr << "done." << std::endl;

    // Read the alignments again and, for the reads chosen above, extract info
    // about the read locations relative to the reference transcripts.
    std::cerr << curtime() << "Reading the alignments and extracting info about read locations..." << std::flush;
    std::vector<ReadStruct> reads;
    extract_reads(reads, is_lidx_chosen, bam_fname, o.paired);
    std::cerr << "done." << std::endl;

    // Sort the read by the location they align to.
    std::cerr << curtime() << "Sorting the reads by their location..." << std::flush;
    std::sort(reads.begin(), reads.end());
    std::cerr << "done." << std::endl;

    for (size_t mo = o.min_overlap_begin; mo <= o.min_overlap_end; ++mo) {
      // Assemble the reads into contigs.
      std::cerr << curtime() << "Assembling the reads into contigs at minimum overlap " << mo << "..." << std::flush;
      std::vector<ContigStruct> contigs;
      std::vector<int> read2contig;
      assemble(contigs, reads, mo, o.paired, read2contig);

      // Output the contigs.
      std::ostringstream ss;
      ss << o.assembly << "_" << mo << ".fa";
      std::cerr << "writing " << contigs.size() << " contigs..." << std::flush;
      output(ss.str(), contigs, ref_fa, bam_fname, o.paired);
      std::cerr << "done." << std::endl;
    }

    std::cerr << curtime() << "Done with everything." << std::endl;
    return 0;

  } catch (const boost::program_options::error& x) {

    std::cerr << std::endl;
    std::cerr << argv[0] << ": Error: " << x.what() << std::endl;
    std::cerr << "Check " << argv[0] << " --help for more information." << std::endl;
    return 1;

  } /* catch (const std::exception& x) {

    std::cerr << std::endl;
    std::cerr << argv[0] << ": Error: " << x.what() << std::endl;
    return 1;

  } */
}
