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

struct probs_and_lidxs {
  std::vector<size_t> lidxs;
  std::vector<double> probs;
};

struct ReadStruct
{
  int rid, tid, pos, readlen;
  int cid, rpos; // id in contigs and right most end's contig position

  ReadStruct(int rid, int tid, int pos, int readlen) {
    this->rid = rid;
    this->tid = tid;
    this->pos = pos;
    this->readlen = readlen;
    cid = rpos = -1;
  }

  bool operator< (const ReadStruct& o) const {
    if (tid != o.tid) return tid < o.tid;
    return pos < o.pos;
  }
};

struct ContigStruct
{
  int tid, pos, len;

  ContigStruct(int tid, int pos, int len) {
    this->tid = tid;
    this->pos = pos;
    this->len = len;
  }
};

struct BreakPointStruct
{
  int cid, pos, overlap;

  BreakPointStruct(int cid, int pos, int overlap) {
    this->cid = cid;
    this->pos = pos;
    this->overlap = overlap;
  }
};

bool is_ok(const bam1_t *rec)
{
  // // Exclude ``right'' mates.
  // // we only need to use one of the pairs of reads; let's choose the
  // // left/upstream ones
  // // XXX should we really do this?
  // if (rec->core.flag & BAM_FREAD2) // this is read2
  //   return false;

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

size_t extract_alignment_probs(
    std::map<std::string, probs_and_lidxs>& read_to_probs_and_lidxs,
    const std::string& bam_fname,
    double min_alignment_prob)
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
    if (!is_ok(rec))
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

void choose_alignments(
    std::vector<bool>& is_lidx_chosen,
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
}

void extract_reads(
    std::vector<ReadStruct>& reads,
    const std::vector<bool>& is_lidx_chosen,
    const std::string& bam_fname)
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

    // Record the above info about this read. This is the main info that
    // will be relevant for making the contigs later on.
    reads.push_back(ReadStruct(-1, tid, pos, readlen));

  }

  // Clean up.
  bam_destroy1(rec);
  bam_header_destroy(header);
}

void assemble(
  std::vector<ContigStruct>& contigs,
  std::vector<ReadStruct>& reads,
  int min_overlap)
{
  int s = -1;
  size_t num_reads = reads.size();

  for (size_t i = 0; i < num_reads; i++) {
    if (contigs.size() > 0 && contigs[s].tid == reads[i].tid) {
      assert(contigs[s].pos <= reads[i].pos);
      assert(contigs[s].pos + contigs[s].len <= reads[i].pos + reads[i].readlen);
    }
    if (contigs.size() == 0 ||
        contigs[s].tid != reads[i].tid ||
        contigs[s].pos + contigs[s].len - reads[i].pos < min_overlap) {
      ++s;
      contigs.push_back(ContigStruct(reads[i].tid, reads[i].pos, reads[i].readlen));
    }
    else {
      contigs[s].len = reads[i].pos + reads[i].readlen - contigs[s].pos;
    }

    reads[i].cid = s;
    reads[i].rpos = contigs[s].len - 1;
  }

  ++s; // there are s contigs

  // if (breakpoints != NULL) {
  //   for (int i = 0; i < num_reads - 1; i++) {
  //     // read is the end of a contig
  //     if (reads[i].rpos + 1 == contigs[reads[i].cid].len)
  //       continue;
  //     // there is another read coincide with this read
  //     if (reads[i].rpos == reads[i + 1].rpos)
  //       continue;
  //     assert(reads[i].cid == reads[i + 1].cid);
  //     breakpoints->push_back(BreakPointStruct(reads[i].cid, reads[i].rpos, readLen - (reads[i + 1].rpos - reads[i].rpos)));
  //   }
  // }
}

void output(
  const std::string& output_fname,
  const std::vector<ContigStruct>& contigs,
  const fasta& ref_fa,
  const std::string& bam_fname)
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

  // Assemble each contig.
  for (size_t i = 0; i < contigs.size(); i++) {

    // We have recorded, in contigs[i].tid, the idx of the transcript that this
    // contig comes from, relative to the ordering of transcripts in the BAM
    // file. However, the transcript sequences are stored in the FASTA file,
    // and the order here could be different. Thus, here, we will figure out
    // the idx of the transcript relative to the order in the FASTA file.
    std::string tname(header->target_name[contigs[i].tid]);
    std::map<std::string, size_t>::const_iterator it;
    it = ref_fa.names_to_idxs.find(tname);
    if (it == ref_fa.names_to_idxs.end())
      throw std::runtime_error("The following transcript name does not "
          "correspond to a reference sequence: " + tname);
    size_t tidx = it->second;

    // Get the transcript sequence and extract the contig.
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

    // Output the contig in FASTA format.
    out << ">contig_" << i
        << " tname="  << tname
        << " pos="    << contigs[i].pos 
        << " length=" << contigs[i].len << "\n";
    out << seq << "\n";

  }

  // Clean up.
  bam_header_destroy(header);

  // if (outputBreakPoints) {
  //   sprintf(breakPointsF, "%s_%d.bp", breakPointsFN, windowSize);
  //   fo = fopen(breakPointsF, "w");
  //   if (!fo)
  //     throw std::runtime_error(std::string("Cannot open file to write breakpoints: ") + breakPointsF);
  //   int nbp = (int)breakPoints.size();
  //   fprintf(fo, "%d\n", nbp);
  //   for (int i = 0; i < nbp; i++) 
  //     fprintf(fo, "%d %d %d\n", breakPoints[i].cid, breakPoints[i].pos, breakPoints[i].overlap);
  //   fclose(fo);
  //   printf("All break points are written to %s!\n", breakPointsF);
  // }
}

void print_help(std::ostream& os)
{
  std::cout <<
  "Overview\n"
  "--------\n"
  "\n"
  "This program constructs an estimate of the \"true\" assembly of a set of\n"
  "reads, relative to a set of reference sequences, based on alignment\n"
  "information produced by RSEM.\n"
  "\n"
  "As defined by DETONATE [1], the \"true\" assembly of a set of reads, relative\n"
  "to the true transcript sequences that the reads were generated from, is the\n"
  "set of contiguous subsequences of the transcript sequences that are covered\n"
  "by reads, with the reads overlapping by at least --min-overlap bases, when\n"
  "each read is aligned to its true location of origin.\n"
  "\n"
  "In practice, we do not know the reads' true locations of origin; in fact, we\n"
  "do not even know (precisely) the true transcript sequences that the reads\n"
  "were generated from. Instead, we start with a set of reference sequences\n"
  "(e.g., an Ensembl reference), and we align each read to these reference\n"
  "sequences using RSEM/Bowtie. The current program chooses a subset of these\n"
  "alignments, according to a policy specified by --alignment-policy and\n"
  "--min-alignment-prob, and then builds our best guess as to the \"true\"\n"
  "assembly based on these alignments.\n"
  "\n"
  "This program was mostly written by Nathanael Fillmore, with key parts of the\n"
  "code written by Bo Li, and with additional input from Colin Dewey. Please\n"
  "send questions and bug reports to the DETONATE mailing list (TODO@wisc.edu).\n"
  "\n"
  "[1] Bo Li*, Nathanael Fillmore*, Yongsheng Bai, Mike Collins, James A.\n"
  "    Thompson, Ron Stewart, Colin N. Dewey. Evaluation of de novo \n"
  "    transcriptome assemblies from RNA-Seq data.\n"
  "\n"
  "\n"
  "Example usage\n"
  "-------------\n"
  "\n"
  "First, use a recent version of RSEM to quantify the expression of the\n"
  "full-length transcripts relative to the reads. For example:\n"
  "\n"
  "$ rsem-prepare-reference --gtf mm9.ensembl63.filtered.gtf mm9.fa ref\n"
  "$ rsem-calculate-expression --num-threads 24 reads.fq ref expr\n"
  "\n"
  "Second, use this program to build the \"true\" assembly:\n"
  "\n"
  "$ ./ref-eval-build-true-assembly --reference ref --expression expr --assembly cc\n"
  "\n"
  "This will output a file, cc_0.fa, that contains the \"true\" assembly.\n"
  "\n"
  "\n"
  "General options\n"
  "---------------\n"
  "\n"
  "  -? [ --help ] \n"
  "\n"
  "        Display this information.\n"
  "\n"
  "\n"
  "Options that specify input and output\n"
  "-------------------------------------\n"
  "\n"
  "  --reference arg\n"
  "\n"
  "        The prefix of the reference built by rsem-prepare-reference. Required.\n"
  "\n"
  "  --expression arg\n"
  "\n"
  "        The prefix of the expression built by rsem-calculate-expression. \n"
  "        Required.\n"
  "\n"
  "  --assembly arg\n"
  "\n"
  "        A prefix to write the \"true\" assembly or sequence of assemblies to. The\n"
  "        suffix \"_x.fa\" will be appended to this prefix, where x is the minimum \n"
  "        overlap size. Required.\n"
  "\n"
  "\n"
  "Options that change the output\n"
  "------------------------------\n"
  "\n"
  "  --min-overlap arg\n"
  "\n"
  "        Either:\n"
  "\n"
  "        (1) An integer that specifies how much overlap between two reads is \n"
  "            required to merge two reads. For example, if --min-overlap=3, then \n"
  "            only reads whose chosen alignments overlap by at least 3 bases will\n"
  "            be joined into contigs. If --min-overlap=0, then only reads whose \n"
  "            chosen alignments are contiguous (or overlap by a positive amount) \n"
  "            will be joined into contigs.\n"
  "\n"
  "        Or:\n"
  "\n"
  "        (2) A pair of integers, separated by commas, specifying a range of \n"
  "            overlap sizes, as described above. For example, if \n"
  "            --min-overlap=2,4 is given, then three assemblies will be produced,\n"
  "            corresponding to --min-overlap=2, --min-overlap=3, and \n"
  "            --min-overlap=4. You might use this option to compute ideal \n"
  "            assemblies at all overlap sizes, e.g., --min-overlap=0,76 for \n"
  "            76-length reads.\n"
  "\n"
  "        Default: 0.\n"
  "\n"
  "  --min-alignment-prob arg\n"
  "\n"
  "        A number between 0 and 1 (inclusive). Any alignment (of a read to a \n"
  "        reference transcript) with posterior probability, as calculated by \n"
  "        RSEM, strictly less than this value will be discarded. Noise reads, \n"
  "        with posterior probability exactly 0, are always discarded. \n"
  "        Default: 0.\n"
  "\n"
  "  --alignment-policy arg\n"
  "\n"
  "        The policy used to choose which alignment(s) of each read to use in \n"
  "        constructing the \"true\" assembly. Options:\n"
  "\n"
  "        - sample: For each read, sample a single alignment (to some reference \n"
  "                  transcript) according to the posterior probability that the \n"
  "                  read follows each alignment, as calculated by RSEM.\n"
  "        - best:   For each read, choose the alignment that maximizes the \n"
  "                  posterior probability mentioned above. Ties are broken \n"
  "                  arbitrarily but deterministically (the first alignment in the\n"
  "                  BAM file is used).\n"
  "        - all:    For each read, use all its alignments. Some reads might end \n"
  "                  up with more than one alignment. In that case, contigs will \n"
  "                  be made assuming that the read aligns to each place. (In \n"
  "                  other words, the read is effectively duplicated, with one \n"
  "                  copy per alignment.)\n"
  "\n"
  "        This policy is applied after the thresholding implied by \n"
  "        --min-alignment-prob. For example, if \"--min-alignment-prob=0.10 \n"
  "        --alignment-policy=sample\" is given, then (first) all alignments with \n"
  "        posterior probability less than 0.10 will be discarded, and (second), \n"
  "        for each read, an alignment will be sampled from among the remaining \n"
  "        alignments, with the posterior distribution renormalized as \n"
  "        appropriate. As another example, if \"--min-alignment-prob=0.90 \n"
  "        --alignment-policy=all\" is given, then all alignments with posterior \n"
  "        probability at least 0.90 will be used.\n"
  "\n"
  "        Default: sample.\n"
  << std::flush;
}

boost::program_options::options_description make_options_old()
{
  namespace po = boost::program_options;

  po::options_description desc(
      "General options\n"
      "---------------\n", 80, 72);
  desc.add_options()
    ("help,?", "Display this information.\n");

  po::options_description data(
      "Input/output specification\n"
      "--------------------------\n", 80, 72);
  data.add_options()
    ("reference", po::value<std::string>(),
        "\nThe prefix of the reference built by rsem-prepare-reference. "
        "Required.\n")
    ("expression", po::value<std::string>(),
        "\nThe prefix of the expression built by rsem-calculate-expression. "
        "Required.\n")
    ("assembly", po::value<std::string>(),
        "\nA prefix to write the \"true\" assembly or sequence of assemblies "
        "to. The suffix \"_x.fa\" will be appended to this prefix, where "
        "x is the overlap size. Required.\n")
        //"Either:\n"
        //"(1) \tThe file name to write the assembly to, if --min-overlap is "
        //      "a single integer.\n"
        //"Or:\n"
        //"(2) \tA prefix to write the sequence of assemblies to, if "
        //      "--min-overlap is a range of integers. For example, if "
        //      "--output=./truth_ --min-overlap=2,4 is given, "
        //      "then the assembly with overlap 2 will be written to "
        //      "./truth_2.fa, etc.\n")
    ;
  desc.add(data);

  po::options_description params(
      "Parameters that change the output\n"
      "---------------------------------\n", 80, 72);
  params.add_options()
    ("min-overlap", po::value<std::string>(),
        "\nEither:\n"
        "\n"
        "(1) \tAn integer that specifies how much overlap between two reads "
              "is required to merge two reads. For example, if "
              "--min-overlap=3, then only reads whose chosen alignments "
              "overlap by at least 3 bases will be joined into contigs. "
              "If --min-overlap=0, then only reads whose chosen alignments "
              "are contiguous (or overlap by a positive amount) will be "
              "joined into contigs.\n"
        "\n"
        "Or:\n"
        "\n"
        "(2) \tA pair of integers, separated by commas, specifying a range "
              "of overlap sizes, as described above. For example, if "
              "--min-overlap=2,4 is given, then three assemblies will be "
              "produced, corresponding to --min-overlap=2, --min-overlap=3, "
              "and --min-overlap=4. You might use this option to compute "
              "ideal assemblies at all overlap sizes, e.g., "
              "--min-overlap=0,76 for 76-length reads.\n"
        "\n"
        "Default: 0.\n")
    ("min-alignment-prob", po::value<double>(),
        "\nA number between 0 and 1 (inclusive). "
        "Any alignment (of a read to a reference transcript) with "
        "posterior probability, as calculated by RSEM, less than this "
        "value will be discarded. Noise reads, with posterior probability "
        "exactly 0, are always discarded. Default: 0.\n")
    ("alignment-policy", po::value<std::string>(),
        "\nThe policy used to choose which alignment(s) of each read to "
        "use in constructing the \"true\" assembly. Options:\n"
        "\n"
        "- sample: \tFor each read, sample a single alignment (to some "
                    "reference transcript) according to the posterior "
                    "probability that the read follows each alignment, "
                    "as calculated by RSEM.\n"
        "- best:   \tFor each read, choose the alignment that maximizes "
                    "the posterior probability mentioned above. Ties "
                    "are broken arbitrarily but deterministically (the "
                    "first alignment in the BAM file is used).\n"
        "- all:    \tFor each read, use all its alignments. Some reads "
                    "might end up with more than one alignment. In that "
                    "case, contigs will be made assuming that the read "
                    "aligns to each place. (In other words, the read is "
                    "effectively duplicated, with one copy per alignment.)\n"
        "\n"
        "This policy is applied after the thresholding implied by "
        "--min-alignment-prob. For example, if \"--min-alignment-prob=0.10 "
        "--alignment-policy=sample\" is given, then (first) all alignments "
        "with posterior probability less than 0.10 will be discarded, and "
        "(second), for each read, an alignment will be sampled from among "
        "the remaining alignments, with the posterior distribution "
        "renormalized as appropriate. As another example, if "
        "\"--min-alignment-prob=0.90 --alignment-policy=all\" is given, "
        "then all alignments with posterior probability at least 0.90 "
        "will be used.\n"
        "\n"
        "Default: sample.\n")
    // ("breakpoints", po::value<std::string>(),
    //     "A prefix to write the breakpoints to. Optional. If not given, "
    //     "the breakpoints will not be written.\n")
    ;
  desc.add(params);

  return desc;
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
    if (o.alignment_policy != "random" &&
        o.alignment_policy != "best" &&
        o.alignment_policy != "all")
      throw po::error("Invalid value for --alignment-policy: " + o.alignment_policy);
  } else {
    o.alignment_policy = "random";
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
      print_help(std::cout);
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
    size_t num_alignments = extract_alignment_probs(read_to_probs_and_lidxs, bam_fname, o.min_alignment_prob);
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
    extract_reads(reads, is_lidx_chosen, bam_fname);
    std::cerr << "done." << std::endl;

    for (size_t mo = o.min_overlap_begin; mo <= o.min_overlap_end; ++mo) {
      // Assemble the reads into contigs.
      std::cerr << curtime() << "Assembling the reads into contigs at minimum overlap " << mo << "..." << std::flush;
      std::vector<ContigStruct> contigs;
      assemble(contigs, reads, mo);

      // Output the contigs.
      std::ostringstream ss;
      ss << o.assembly << "_" << mo << ".fa";
      std::cerr << "writing " << contigs.size() << " contigs..." << std::flush;
      output(ss.str(), contigs, ref_fa, bam_fname);
      std::cerr << "done." << std::endl;
    }

    std::cerr << curtime() << "Done with everything." << std::endl;
    return 0;

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
}
